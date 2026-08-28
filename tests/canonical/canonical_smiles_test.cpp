#include <doctest/doctest.h>

#include <array>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include "chem/canonical/canonical_smiles.hpp"
#include "chem/core/errors.hpp"
#include "chem/core/molecular_graph.hpp"
#include "chem/parsing/formula_parser.hpp"
#include "chem/parsing/smiles_parser.hpp"

namespace {

using chem::Atom;
using chem::BondOrder;
using chem::Element;
using chem::MolecularGraph;
using chem::ValidationError;

// Reinserts a graph's atoms in a different order and rewires the bonds to the
// same logical molecule (FR-12a order-independence helper).
// `perm` lists logical atom indices in insertion order.
MolecularGraph permuteGraph(const MolecularGraph& src, const std::vector<int>& perm) {
  std::vector<int> new_index(static_cast<std::size_t>(src.atoms().size()), -1);
  MolecularGraph result;
  for (const int logical : perm) {
    new_index[static_cast<std::size_t>(logical)] =
        static_cast<int>(result.addAtom(src.atoms()[static_cast<std::size_t>(logical)]));
  }
  for (const auto& bond : src.bonds()) {
    result.addBond(new_index[static_cast<std::size_t>(bond.a)],
                   new_index[static_cast<std::size_t>(bond.b)], bond.order);
  }
  return result;
}

std::string can(std::string_view smiles) {
  return chem::canonicalSmiles(chem::parseSmiles(smiles));
}

std::string canFormula(std::string_view formula) {
  return chem::canonicalSmiles(chem::parseFormula(formula));
}

} // namespace

// ===========================================================================
// FR-1: Entry point and version constant
// ===========================================================================

TEST_CASE("kCanonicalSpecVersion is 1") { CHECK(chem::kCanonicalSpecVersion == 1); }

TEST_CASE("FR-1b: empty graph throws ValidationError") {
  const MolecularGraph empty;
  CHECK_THROWS_AS(static_cast<void>(chem::canonicalSmiles(empty)), const ValidationError&);
}

TEST_CASE("FR-1a: canonicalization is deterministic (twice, byte-equal)") {
  const std::array<std::string_view, 6> inputs = {
      "C1=CC=CC=C1", "CC(=O)O", "C1CCC2CCCCC2C1", "[NH4+]", "CS(=O)(=O)O", "CCO",
  };
  for (const std::string_view input : inputs) {
    const MolecularGraph g = chem::parseSmiles(input);
    CHECK(chem::canonicalSmiles(g) == chem::canonicalSmiles(g));
  }
}

// ===========================================================================
// FR-6: Atom emission rules
// ===========================================================================

TEST_CASE("FR-6: bare organic atoms emit bare (C, O)") {
  CHECK(can("C") == "C");
  CHECK(can("O") == "O");
  CHECK(can("CCO") == "C(O)C");
}

TEST_CASE("FR-6: equivalent bracket spellings normalize to bare form") {
  // design.md §4.1 / FR-12d: [OH2] and O collapse to one string.
  CHECK(can("[OH2]") == "O");
  CHECK(can("[CH4]") == "C");
  CHECK(can("[C]") == "C"); // formula path single carbon has no stored H
  CHECK(can("[O]") == "O");
}

TEST_CASE("FR-6: charge forces bracketed emission ([NH4+], [OH-])") {
  CHECK(can("[NH4+]") == "[NH4+]");
  CHECK(can("[OH-]") == "[OH-]");
  CHECK(can("CC(=O)[O-]") == "C([O-])(=O)C");
}

TEST_CASE("FR-6: isotope forces bracketed emission ([13CH4], [2H]O[2H])") {
  CHECK(can("[13CH4]") == "[13CH4]");
  CHECK(can("[2H]O[2H]") == "O([2H])[2H]");
}

TEST_CASE("FR-6: over-valent-in-bare atoms stay bracketed") {
  // [OH3] (oxonium) cannot be reproduced as bare O (valence model gives 2 H).
  CHECK(can("[OH3+]") == "[OH3+]");
}

TEST_CASE("FR-6: non-organic elements always bracket ([Fe+3])") {
  CHECK(can("[Fe+3]") == "[Fe+3]");
}

// ===========================================================================
// FR-7: Bond emission
// ===========================================================================

TEST_CASE("FR-7: single bonds omitted, double/triple explicit") {
  CHECK(can("C=C") == "C=C");
  CHECK(can("C#N") == "C#N");
  CHECK(can("C=O") == "C=O");
  CHECK(can("CO") == "CO");
  // Butadiene: min-string root is the central carbon, so the canonical is
  // C(=C)C=C rather than C=CC=C; both round-trip, the former is lexicographically smaller.
  CHECK(can("C=C-C=C") == "C(=C)C=C");
}

TEST_CASE("FR-7: branch bond symbol inside parens, ring symbol at both ends") {
  CHECK(can("CC(=O)O") == "C(=O)(O)C");
  CHECK(can("C1=CC=CC=C1") == "C1=CC=CC=C1");  // single ring closures omit symbol
  CHECK(can("C=1C=CC=CC=1") == "C1=CC=CC=C1"); // kekule rotation collapses
}

// ===========================================================================
// FR-8: Disconnected graphs
// ===========================================================================

TEST_CASE("FR-8: formula H2O canonicalizes to O.[H].[H] (sorted components)") {
  CHECK(canFormula("H2O") == "O.[H].[H]");
}

TEST_CASE("FR-8: formula methane canonicalizes to C plus explicit H") {
  CHECK(canFormula("CH4") == "C.[H].[H].[H].[H]");
}

TEST_CASE("FR-8: components joined in ascending lexicographic order") {
  // Ethanol formula multiset: two C, one O, six H. Bare "C" and "O" sort before
  // bracketed "[H]" ('C' 0x43 < 'O' 0x4F < '[' 0x5B).
  CHECK(canFormula("C2H6O") == "C.C.O.[H].[H].[H].[H].[H].[H]");
}

// ===========================================================================
// FR-12a: Order independence (permuted insertion orders)
// ===========================================================================

TEST_CASE("FR-12a: benzene insertion-order permutations canonicalize identically") {
  const MolecularGraph benzene = chem::parseSmiles("C1=CC=CC=C1");
  const std::string canonical = chem::canonicalSmiles(benzene);
  const std::vector<std::vector<int>> permutations = {
      {0, 1, 2, 3, 4, 5},
      {5, 4, 3, 2, 1, 0},
      {2, 1, 0, 5, 4, 3},
      {3, 4, 5, 0, 1, 2},
  };
  for (const auto& perm : permutations) {
    const MolecularGraph rotated = permuteGraph(benzene, perm);
    CHECK(chem::canonicalSmiles(rotated) == canonical);
  }
}

TEST_CASE("FR-12a: acetic acid insertion-order permutations canonicalize identically") {
  const MolecularGraph acetic = chem::parseSmiles("CC(=O)O");
  const std::string canonical = chem::canonicalSmiles(acetic);
  const std::vector<std::vector<int>> permutations = {
      {0, 1, 2, 3},
      {1, 0, 2, 3},
      {1, 2, 3, 0},
      {3, 2, 1, 0},
  };
  for (const auto& perm : permutations) {
    const MolecularGraph permuted = permuteGraph(acetic, perm);
    CHECK(chem::canonicalSmiles(permuted) == canonical);
  }
}

TEST_CASE("FR-12a: neopentane methyl-order permutations canonicalize identically") {
  const MolecularGraph neopentane = chem::parseSmiles("CC(C)(C)C");
  const std::string canonical = chem::canonicalSmiles(neopentane);
  const std::vector<std::vector<int>> permutations = {
      {0, 1, 2, 3, 4},
      {2, 3, 4, 1, 0},
      {4, 3, 2, 1, 0},
  };
  for (const auto& perm : permutations) {
    const MolecularGraph permuted = permuteGraph(neopentane, perm);
    CHECK(chem::canonicalSmiles(permuted) == canonical);
  }
}

// ===========================================================================
// FR-12d: Equivalence classes
// ===========================================================================

TEST_CASE("FR-12d: O and [OH2] collapse to one string") {
  CHECK(can("O") == can("[OH2]"));
  CHECK(can("[OH2]") == "O");
}

TEST_CASE("FR-12d: benzene kekule rotations collapse to one string") {
  CHECK(can("C1=CC=CC=C1") == can("C=1C=CC=CC=1"));
  CHECK(can("C1=CC=CC=C1") == "C1=CC=CC=C1");
}

TEST_CASE("FR-12d: hand-built ethanol graph equals parsed CCO") {
  MolecularGraph handmade;
  const auto c_methyl = handmade.addAtom(Atom{.element = Element(6), .implicit_h = 3});
  const auto c_methylene = handmade.addAtom(Atom{.element = Element(6), .implicit_h = 2});
  const auto o_hydroxyl = handmade.addAtom(Atom{.element = Element(8), .implicit_h = 1});
  handmade.addBond(c_methyl, c_methylene, BondOrder::kSingle);
  handmade.addBond(c_methylene, o_hydroxyl, BondOrder::kSingle);

  CHECK(chem::canonicalSmiles(handmade) == chem::canonicalSmiles(chem::parseSmiles("CCO")));
  CHECK(chem::canonicalSmiles(handmade) == "C(O)C");
}

// ===========================================================================
// FR-5d: ring-closure digit limits
// ===========================================================================

// Builds a molecule whose root is the unique lowest-class atom (B) holding n
// ring-closure back-edges that all open simultaneously at the root.
MolecularGraph buildStarOfClosures(int n) {
  MolecularGraph g;
  const auto root = g.addAtom(Atom{.element = Element(5)}); // B: atomic 5 < C 6
  const auto hub = g.addAtom(Atom{.element = Element(6)});
  g.addBond(root, hub, BondOrder::kSingle);
  for (int i = 0; i < n; ++i) {
    const auto leaf = g.addAtom(Atom{.element = Element(6)});
    g.addBond(hub, leaf, BondOrder::kSingle);
    g.addBond(root, leaf, BondOrder::kSingle); // back-edge: ring closure at root
  }
  return g;
}

TEST_CASE("FR-5d: more than 99 simultaneous ring closures throw ValidationError") {
  const MolecularGraph star = buildStarOfClosures(100);
  CHECK_THROWS_AS(static_cast<void>(chem::canonicalSmiles(star)), const ValidationError&);
}

TEST_CASE("FR-5d: ten simultaneous closures emit %10 on the output side") {
  const std::string canonical = chem::canonicalSmiles(buildStarOfClosures(10));
  CHECK(canonical.find("%10") != std::string::npos);
  CHECK(canonical.find("%11") == std::string::npos);
}
