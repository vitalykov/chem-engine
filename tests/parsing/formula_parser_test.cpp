#include <doctest/doctest.h>

#include <string>
#include <string_view>

#include "chem/core/composition.hpp"
#include "chem/core/errors.hpp"
#include "chem/core/molecular_graph.hpp"
#include "chem/parsing/formula_parser.hpp"

namespace {

using chem::MolecularGraph;
using chem::ParseError;

struct ParsedFormula {
  MolecularGraph graph;
  chem::CompositionMap counts;
};

ParsedFormula parse_and_count(std::string_view input) {
  ParsedFormula result;
  result.graph = chem::parse_formula(input);
  result.counts = chem::composition(result.graph);
  return result;
}

bool same_atoms(const MolecularGraph& lhs, const MolecularGraph& rhs) {
  if (lhs.atoms().size() != rhs.atoms().size() || lhs.bonds().size() != rhs.bonds().size()) {
    return false;
  }
  for (std::size_t i = 0; i < lhs.atoms().size(); ++i) {
    const auto& a = lhs.atoms()[i];
    const auto& b = rhs.atoms()[i];
    if (a.element.atomic_number() != b.element.atomic_number() || a.charge != b.charge ||
        a.implicit_h != b.implicit_h || a.explicit_h != b.explicit_h) {
      return false;
    }
  }
  return true;
}

void expect_rejection(std::string_view input, std::string_view message_fragment) {
  try {
    static_cast<void>(chem::parse_formula(input));
    FAIL("expected ParseError for input: ", input);
  } catch (const ParseError& e) {
    const std::string message = e.what();
    CHECK_MESSAGE(message.find(message_fragment) != std::string::npos, "input '", input,
                  "' produced message: ", message);
  }
}

TEST_CASE("Single element parses to one atom without bonds") {
  const ParsedFormula parsed = parse_and_count("O");
  REQUIRE(parsed.graph.atoms().size() == 1);
  CHECK(parsed.graph.atoms()[0].element.atomic_number() == 8);
  CHECK(parsed.graph.bonds().empty());
  const chem::CompositionMap expected{{8, 1}};
  CHECK(parsed.counts == expected);
}

TEST_CASE("Repeated elements aggregate counts") {
  const ParsedFormula parsed = parse_and_count("H2O");
  CHECK(parsed.graph.atoms().size() == 3);
  const chem::CompositionMap expected{{1, 2}, {8, 1}};
  CHECK(parsed.counts == expected);
  CHECK(chem::molar_mass(parsed.graph) == doctest::Approx(18.015).epsilon(0.001 / 18.015));
}

TEST_CASE("Groups with multipliers") {
  const ParsedFormula calcium_hydroxide = parse_and_count("Ca(OH)2");
  CHECK(calcium_hydroxide.graph.atoms().size() == 5);
  const chem::CompositionMap hydroxide{{20, 1}, {8, 2}, {1, 2}};
  CHECK(calcium_hydroxide.counts == hydroxide);

  const ParsedFormula ferric_nitrate = parse_and_count("Fe(NO3)3");
  CHECK(ferric_nitrate.graph.atoms().size() == 13);
  const chem::CompositionMap nitrate{{26, 1}, {7, 3}, {8, 9}};
  CHECK(ferric_nitrate.counts == nitrate);
}

TEST_CASE("Nested groups multiply through all levels") {
  const ParsedFormula parsed = parse_and_count("Mg(NO2(OH))2");
  const chem::CompositionMap expected{{12, 1}, {7, 2}, {8, 6}, {1, 2}};
  CHECK(parsed.counts == expected);

  const ParsedFormula deep = parse_and_count("K4(ON(SO3)2)2");
  const chem::CompositionMap deep_expected{{19, 4}, {8, 14}, {7, 2}, {16, 4}};
  CHECK(deep.counts == deep_expected);
}

TEST_CASE("Group without multiplier defaults to one") {
  const ParsedFormula parsed = parse_and_count("((H))");
  const chem::CompositionMap expected{{1, 1}};
  CHECK(parsed.counts == expected);
}

TEST_CASE("Multi-digit counts") {
  const ParsedFormula glucose = parse_and_count("C6H12O6");
  const chem::CompositionMap expected{{6, 6}, {1, 12}, {8, 6}};
  CHECK(glucose.counts == expected);

  const ParsedFormula many = parse_and_count("W123");
  const chem::CompositionMap tungsten{{74, 123}};
  CHECK(many.counts == tungsten);
}

TEST_CASE("Parsing is deterministic") {
  const MolecularGraph first = chem::parse_formula("Ca(OH)2");
  const MolecularGraph second = chem::parse_formula("Ca(OH)2");
  CHECK(same_atoms(first, second));
}

TEST_CASE("Empty input is rejected") { expect_rejection("", "empty"); }

TEST_CASE("Unknown element symbols are rejected quoting the symbol") {
  expect_rejection("Xx", "Xx");
  expect_rejection("Zz", "Zz");
  expect_rejection("Hx", "Hx");
}

TEST_CASE("Characters outside the grammar are rejected at their position") {
  expect_rejection("H2O ", "position 4");
  expect_rejection("h2O", "'h'");
  expect_rejection("2H", "'2'");
  expect_rejection("H.O", "'.'");
  expect_rejection("[OH-]", "'['");
  expect_rejection("C+H", "'+'");
  expect_rejection("Cl-", "'-'");
}

TEST_CASE("Middle-dot hydrate notation is rejected") { expect_rejection("CuSO4·5H2O", "byte 0x"); }

TEST_CASE("Formulas exceeding the atom limit are rejected without allocating") {
  expect_rejection("H99999999", "atom limit");
  expect_rejection("(H)999999999", "atom limit");
  expect_rejection("((((H)99)99)99)99)", "atom limit");
}

TEST_CASE("Structural errors are rejected") {
  expect_rejection("H)", "unbalanced closing parenthesis");
  expect_rejection("(H", "unclosed group");
  expect_rejection("()", "empty group");
  expect_rejection("()2", "empty group");
  expect_rejection("H0", "zero count");
  expect_rejection("H007", "leading zeros");
  expect_rejection("(OH)0", "zero count");
  expect_rejection(" ", "'");
}

} // namespace
