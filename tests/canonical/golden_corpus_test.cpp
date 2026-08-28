#include <doctest/doctest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "chem/canonical/canonical_smiles.hpp"
#include "chem/core/composition.hpp"
#include "chem/core/molecular_graph.hpp"
#include "chem/parsing/formula_parser.hpp"
#include "chem/parsing/smiles_parser.hpp"
#include "graph_iso.hpp"

namespace {

using chem::CompositionMap;
using chem::MolecularGraph;

struct CorpusEntry {
  std::string id;
  std::string input_format;
  std::string input;
  std::string expected;
  std::uint32_t spec_version;
};

std::vector<std::string> split(std::string_view text, char delimiter) {
  std::vector<std::string> parts;
  std::size_t start = 0;
  while (start <= text.size()) {
    const std::size_t comma = text.find(delimiter, start);
    if (comma == std::string_view::npos) {
      parts.emplace_back(text.substr(start));
      break;
    }
    parts.emplace_back(text.substr(start, comma - start));
    start = comma + 1;
  }
  return parts;
}

std::vector<CorpusEntry> loadCorpus() {
  std::vector<CorpusEntry> entries;
  std::ifstream file(CHEM_GOLDEN_CORPUS_CSV_PATH);
  REQUIRE(file.is_open());

  std::string line;
  bool first = true;
  while (std::getline(file, line)) {
    if (first) { // header
      first = false;
      continue;
    }
    if (line.empty()) {
      continue;
    }
    const std::vector<std::string> fields = split(line, ',');
    REQUIRE(fields.size() == 5);
    CorpusEntry entry;
    entry.id = fields[0];
    entry.input_format = fields[1];
    entry.input = fields[2];
    entry.expected = fields[3];
    entry.spec_version = static_cast<std::uint32_t>(std::stoul(fields[4]));
    entries.push_back(std::move(entry));
  }
  return entries;
}

MolecularGraph parseEntry(const CorpusEntry& entry) {
  if (entry.input_format == "smiles") {
    return chem::parseSmiles(entry.input);
  }
  return chem::parseFormula(entry.input);
}

// Sum of stored (implicit+explicit) hydrogens over the whole graph.
int storedHydrogenCount(const MolecularGraph& graph) {
  int total = 0;
  for (const auto& atom : graph.atoms()) {
    total += atom.implicit_h + atom.explicit_h;
  }
  return total;
}

// Splits a dot-joined canonical string into its component strings.
std::vector<std::string> splitComponents(std::string_view canonical) {
  return split(canonical, '.');
}

// Merges dot-separated canonical components into one disconnected graph so the
// whole molecule can be compared against parse(input). parseSmiles rejects
// dots, so each component is parsed independently and its atoms/bonds are
// concatenated; components are disjoint, so index offsets simply shift.
MolecularGraph mergeComponents(const std::vector<std::string>& components) {
  MolecularGraph merged;
  for (const std::string& component : components) {
    const MolecularGraph parsed = chem::parseSmiles(component);
    const std::ptrdiff_t base = static_cast<std::ptrdiff_t>(merged.atoms().size());
    for (const auto& atom : parsed.atoms()) {
      merged.addAtom(atom);
    }
    for (const auto& bond : parsed.bonds()) {
      merged.addBond(base + bond.a, base + bond.b, bond.order);
    }
  }
  return merged;
}

} // namespace

namespace graph_iso = chem::test::graph_iso;

// ===========================================================================
// FR-10 / FR-11b: corpus contents and canonical equality
// ===========================================================================

TEST_CASE("FR-10a: golden corpus has exactly 30 entries") {
  const std::vector<CorpusEntry> entries = loadCorpus();
  CHECK(entries.size() == 30);
}

TEST_CASE("FR-11b: every corpus entry canonicalizes to its frozen string") {
  const std::vector<CorpusEntry> entries = loadCorpus();
  REQUIRE(entries.size() == 30);

  for (const CorpusEntry& entry : entries) {
    const MolecularGraph graph = parseEntry(entry);
    const std::string canonical = chem::canonicalSmiles(graph);
    CHECK_MESSAGE(canonical == entry.expected, "entry ", entry.id, ": got ", canonical);
    CHECK_MESSAGE(entry.spec_version == chem::kCanonicalSpecVersion, "entry ", entry.id);
  }
}

// ===========================================================================
// FR-9 / FR-12b: round-trip invariant over the whole corpus
// ===========================================================================

TEST_CASE("FR-12b: round-trip — canonical is a fixed point and composition is preserved") {
  const std::vector<CorpusEntry> entries = loadCorpus();
  REQUIRE(entries.size() == 30);

  for (const CorpusEntry& entry : entries) {
    const MolecularGraph original = parseEntry(entry);

    if (entry.expected.find('.') == std::string::npos) {
      // Connected entry: parseSmiles can round-trip the full string.
      const MolecularGraph reparsed = chem::parseSmiles(entry.expected);
      const std::string recanonical = chem::canonicalSmiles(reparsed);
      CHECK_MESSAGE(recanonical == entry.expected, "fixed point failed for ", entry.id);

      // Composition equality holds whenever the stored graph carries its
      // hydrogens explicitly (the spec's §8 formula path, e.g. C -> "C", is
      // accepted as a string-only fixed point).
      if (storedHydrogenCount(original) > 0) {
        const CompositionMap original_comp = chem::composition(original);
        const CompositionMap reparsed_comp = chem::composition(reparsed);
        CHECK_MESSAGE(original_comp == reparsed_comp, "composition drift for ", entry.id);
      }
    } else {
      // Disconnected form: parseSmiles rejects dots, so round-trip each
      // component independently; the sorted join must reproduce the frozen
      // string (FR-8b).
      const std::vector<std::string> components = splitComponents(entry.expected);
      std::vector<std::string> recanonicalized;
      for (const std::string& component : components) {
        const MolecularGraph reparsed = chem::parseSmiles(component);
        recanonicalized.push_back(chem::canonicalSmiles(reparsed));
      }
      std::sort(recanonicalized.begin(), recanonicalized.end());
      std::string joined;
      for (std::size_t i = 0; i < recanonicalized.size(); ++i) {
        if (i != 0) {
          joined.push_back('.');
        }
        joined += recanonicalized[i];
      }
      CHECK_MESSAGE(joined == entry.expected, "component fixed point failed for ", entry.id);
    }
  }
}

// ===========================================================================
// FR-12c: determinism smoke — the whole corpus twice, byte-equal
// ===========================================================================

TEST_CASE("FR-12c: full corpus canonicalization is byte-deterministic on a second run") {
  const std::vector<CorpusEntry> entries = loadCorpus();
  REQUIRE(entries.size() == 30);

  std::vector<std::string> first_run;
  first_run.reserve(entries.size());
  for (const CorpusEntry& entry : entries) {
    first_run.push_back(chem::canonicalSmiles(parseEntry(entry)));
  }
  for (std::size_t i = 0; i < entries.size(); ++i) {
    const std::string second = chem::canonicalSmiles(parseEntry(entries[i]));
    CHECK_MESSAGE(second == first_run[i], "determinism drift for ", entries[i].id);
  }
}

// ===========================================================================
// FR-12e: explicit graph-isomorphism guard
// ===========================================================================

// FR-11b (canonical string equality) plus FR-12b (fixed point) tie both graphs
// to the same canonical class, but never compare the graphs themselves. This
// guard closes that gap: parse(input) and parse(expected) must be isomorphic,
// so a canonicalizer that collapses two distinct graphs onto one string (a
// completeness regression) fails here even though the string checks pass.
TEST_CASE("FR-12e: every corpus input graph is isomorphic to its canonical graph") {
  const std::vector<CorpusEntry> entries = loadCorpus();
  REQUIRE(entries.size() == 30);

  for (const CorpusEntry& entry : entries) {
    const MolecularGraph input_graph = parseEntry(entry);

    MolecularGraph expected_graph;
    if (entry.expected.find('.') == std::string::npos) {
      expected_graph = chem::parseSmiles(entry.expected);
    } else {
      expected_graph = mergeComponents(splitComponents(entry.expected));
    }

    const graph_iso::SimpleGraph input_view = graph_iso::fromGraph(input_graph);
    const graph_iso::SimpleGraph expected_view = graph_iso::fromGraph(expected_graph);
    CHECK_MESSAGE(graph_iso::isomorphic(input_view, expected_view),
                  "input not isomorphic to canonical for entry ", entry.id);
  }
}