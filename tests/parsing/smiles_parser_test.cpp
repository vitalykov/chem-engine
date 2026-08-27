#include <doctest/doctest.h>

#include <array>
#include <string>
#include <string_view>

#include "chem/core/composition.hpp"
#include "chem/core/element.hpp"
#include "chem/core/errors.hpp"
#include "chem/core/molecular_graph.hpp"
#include "chem/parsing/formula_parser.hpp"
#include "chem/parsing/smiles_parser.hpp"

namespace {

using chem::Atom;
using chem::BondOrder;
using chem::CompositionMap;
using chem::Element;
using chem::MolecularGraph;
using chem::ParseError;

bool sameGraph(const MolecularGraph& lhs, const MolecularGraph& rhs) {
  if (lhs.atoms().size() != rhs.atoms().size() || lhs.bonds().size() != rhs.bonds().size()) {
    return false;
  }
  for (std::size_t i = 0; i < lhs.atoms().size(); ++i) {
    const Atom& a = lhs.atoms()[i];
    const Atom& b = rhs.atoms()[i];
    if (a.element.atomicNumber() != b.element.atomicNumber() || a.charge != b.charge ||
        a.implicit_h != b.implicit_h || a.explicit_h != b.explicit_h || a.isotope != b.isotope) {
      return false;
    }
  }
  for (std::size_t i = 0; i < lhs.bonds().size(); ++i) {
    const auto& la = lhs.bonds()[i];
    const auto& rb = rhs.bonds()[i];
    if (la.a != rb.a || la.b != rb.b || la.order != rb.order) {
      return false;
    }
  }
  return true;
}

void expectRejection(std::string_view input, std::string_view message_fragment) {
  try {
    static_cast<void>(chem::parseSmiles(input));
    FAIL("expected ParseError for input: ", input);
  } catch (const ParseError& e) {
    const std::string message = e.what();
    CHECK_MESSAGE(message.find(message_fragment) != std::string::npos, "input '", input,
                  "' produced message: ", message);
  }
}

} // namespace

// ===========================================================================
// FR-12a: Acceptance
// ===========================================================================

TEST_CASE("Methane 'C' parses to one carbon with four implicit H") {
  const MolecularGraph g = chem::parseSmiles("C");
  REQUIRE(g.atoms().size() == 1);
  CHECK(g.atoms()[0].element.atomicNumber() == 6);
  CHECK(g.atoms()[0].implicit_h == 4);
  CHECK(g.bonds().empty());
}

TEST_CASE("Water 'O' parses to one oxygen with two implicit H") {
  const MolecularGraph g = chem::parseSmiles("O");
  REQUIRE(g.atoms().size() == 1);
  CHECK(g.atoms()[0].element.atomicNumber() == 8);
  CHECK(g.atoms()[0].implicit_h == 2);
}

TEST_CASE("Bracket water '[OH2]' has explicit H, no implicit H") {
  const MolecularGraph g = chem::parseSmiles("[OH2]");
  REQUIRE(g.atoms().size() == 1);
  CHECK(g.atoms()[0].element.atomicNumber() == 8);
  CHECK(g.atoms()[0].explicit_h == 2);
  CHECK(g.atoms()[0].implicit_h == 0);
  // Same composition as bare "O"
  const MolecularGraph bare = chem::parseSmiles("O");
  CHECK(chem::composition(g) == chem::composition(bare));
}

TEST_CASE("Acetic acid 'CC(=O)O' parses with correct bonds") {
  const MolecularGraph g = chem::parseSmiles("CC(=O)O");
  REQUIRE(g.atoms().size() == 4);
  CHECK(g.atoms()[0].element.atomicNumber() == 6); // C
  CHECK(g.atoms()[1].element.atomicNumber() == 6); // C
  CHECK(g.atoms()[2].element.atomicNumber() == 8); // O
  CHECK(g.atoms()[3].element.atomicNumber() == 8); // O

  REQUIRE(g.bonds().size() == 3);
  // C0-C1 single, C1=O2 double, C1-O3 single
  CHECK(g.bonds()[0].a == 0);
  CHECK(g.bonds()[0].b == 1);
  CHECK(g.bonds()[0].order == BondOrder::kSingle);
  CHECK(g.bonds()[1].a == 1);
  CHECK(g.bonds()[1].b == 2);
  CHECK(g.bonds()[1].order == BondOrder::kDouble);
  CHECK(g.bonds()[2].a == 1);
  CHECK(g.bonds()[2].b == 3);
  CHECK(g.bonds()[2].order == BondOrder::kSingle);

  // C0 has 3 implicit H, C1 has 0 (bond sum 4), O2 has 0, O3 has 1
  CHECK(g.atoms()[0].implicit_h == 3);
  CHECK(g.atoms()[1].implicit_h == 0);
  CHECK(g.atoms()[2].implicit_h == 0);
  CHECK(g.atoms()[3].implicit_h == 1);
}

TEST_CASE("Benzene kekule 'C1=CC=CC=C1' parses as six-membered ring") {
  const MolecularGraph g = chem::parseSmiles("C1=CC=CC=C1");
  REQUIRE(g.atoms().size() == 6);
  for (const Atom& a : g.atoms()) {
    CHECK(a.element.atomicNumber() == 6);
    CHECK(a.implicit_h == 1); // each C has bond sum 3, valence 4
  }
  REQUIRE(g.bonds().size() == 6); // 5 chain + 1 ring closure
  for (const auto& b : g.bonds()) {
    const bool ok = b.order == BondOrder::kSingle || b.order == BondOrder::kDouble;
    CHECK(ok);
  }
}

TEST_CASE("Fused rings 'C1CCC2CCCCC2C1' parse without error") {
  const MolecularGraph g = chem::parseSmiles("C1CCC2CCCCC2C1");
  REQUIRE(g.atoms().size() == 10);
  CHECK(g.bonds().size() == 11); // 9 chain bonds + 2 ring closures
  // All carbons. Fusion atoms (indices 3 and 8) have 3 heavy bonds -> implicit_h 1;
  // the remaining atoms have 2 heavy bonds -> implicit_h 2.
  for (std::size_t i = 0; i < g.atoms().size(); ++i) {
    CHECK(g.atoms()[i].element.atomicNumber() == 6);
    const int expected_h = (i == 3 || i == 8) ? 1 : 2;
    CHECK(g.atoms()[i].implicit_h == expected_h);
  }
}

TEST_CASE("Multi-digit ring closure '%10' works") {
  const MolecularGraph g = chem::parseSmiles("C%10CCCCC%10");
  REQUIRE(g.atoms().size() == 6);
  REQUIRE(g.bonds().size() == 6); // 5 chain + 1 closure
  // The closure bond connects atom 0 and atom 5
  bool found_closure = false;
  for (const auto& b : g.bonds()) {
    if ((b.a == 0 && b.b == 5) || (b.a == 5 && b.b == 0)) {
      found_closure = true;
      CHECK(b.order == BondOrder::kSingle);
    }
  }
  CHECK(found_closure);
}

TEST_CASE("Isotope label '[13CH4]' stores isotope 13") {
  const MolecularGraph g = chem::parseSmiles("[13CH4]");
  REQUIRE(g.atoms().size() == 1);
  CHECK(g.atoms()[0].element.atomicNumber() == 6);
  CHECK(g.atoms()[0].isotope == 13);
  CHECK(g.atoms()[0].explicit_h == 4);
  CHECK(g.atoms()[0].implicit_h == 0);
}

TEST_CASE("Ammonium '[NH4+]' has charge +1 and explicit H 4") {
  const MolecularGraph g = chem::parseSmiles("[NH4+]");
  REQUIRE(g.atoms().size() == 1);
  CHECK(g.atoms()[0].element.atomicNumber() == 7);
  CHECK(g.atoms()[0].charge == 1);
  CHECK(g.atoms()[0].explicit_h == 4);
  CHECK(g.atoms()[0].implicit_h == 0);
}

TEST_CASE("Iron '[Fe+3]' has charge +3") {
  const MolecularGraph g = chem::parseSmiles("[Fe+3]");
  REQUIRE(g.atoms().size() == 1);
  CHECK(g.atoms()[0].element.atomicNumber() == 26);
  CHECK(g.atoms()[0].charge == 3);
}

TEST_CASE("Sulfoxide 'CS(=O)(=O)O' — S selects valence 6") {
  const MolecularGraph g = chem::parseSmiles("CS(=O)(=O)O");
  REQUIRE(g.atoms().size() == 5);
  // C(0)-S(1), S(1)=O(2), S(1)=O(3), S(1)-O(4)
  // S bond sum = 1+2+2+1 = 6 -> valence 6, implicit_h = 0
  CHECK(g.atoms()[1].element.atomicNumber() == 16);
  CHECK(g.atoms()[1].implicit_h == 0);
  // O(2) and O(3) have bond sum 2 -> implicit_h 0
  CHECK(g.atoms()[2].implicit_h == 0);
  CHECK(g.atoms()[3].implicit_h == 0);
  // O(4) has bond sum 1 -> implicit_h 1
  CHECK(g.atoms()[4].implicit_h == 1);
  // C(0) has bond sum 1 -> implicit_h 3
  CHECK(g.atoms()[0].implicit_h == 3);
}

// ===========================================================================
// FR-12b: Graph correctness — bond orders, H totals, charges, isotopes
// ===========================================================================

TEST_CASE("Explicit bond tokens produce correct orders") {
  const MolecularGraph g = chem::parseSmiles("C=O");
  REQUIRE(g.bonds().size() == 1);
  CHECK(g.bonds()[0].order == BondOrder::kDouble);
  CHECK(g.atoms()[0].implicit_h == 2); // C bond sum 2 -> 4-2=2
  CHECK(g.atoms()[1].implicit_h == 0); // O bond sum 2 -> 2-2=0
}

TEST_CASE("Triple bond 'C#N' parses") {
  const MolecularGraph g = chem::parseSmiles("C#N");
  REQUIRE(g.bonds().size() == 1);
  CHECK(g.bonds()[0].order == BondOrder::kTriple);
  CHECK(g.atoms()[0].implicit_h == 1); // C bond sum 3 -> 4-3=1
  CHECK(g.atoms()[1].implicit_h == 0); // N bond sum 3 -> valence 3, 3-3=0
}

TEST_CASE("Branch nesting 'C(C(C)C)C' produces correct atom count") {
  const MolecularGraph g = chem::parseSmiles("C(C(C)C)C");
  REQUIRE(g.atoms().size() == 5);
  CHECK(g.bonds().size() == 4);
}

TEST_CASE("Chlorine 'Cl' parses to one atom with implicit H 1") {
  const MolecularGraph g = chem::parseSmiles("Cl");
  REQUIRE(g.atoms().size() == 1);
  CHECK(g.atoms()[0].element.atomicNumber() == 17);
  CHECK(g.atoms()[0].implicit_h == 1);
}

TEST_CASE("Bromine 'Br' parses to one atom with implicit H 1") {
  const MolecularGraph g = chem::parseSmiles("Br");
  REQUIRE(g.atoms().size() == 1);
  CHECK(g.atoms()[0].element.atomicNumber() == 35);
  CHECK(g.atoms()[0].implicit_h == 1);
}

TEST_CASE("'[Fe++]' is valid with charge +2") {
  const MolecularGraph g = chem::parseSmiles("[Fe++]");
  REQUIRE(g.atoms().size() == 1);
  CHECK(g.atoms()[0].charge == 2);
}

TEST_CASE("'[O-2]' has charge -2") {
  const MolecularGraph g = chem::parseSmiles("[O-2]");
  REQUIRE(g.atoms().size() == 1);
  CHECK(g.atoms()[0].charge == -2);
}

TEST_CASE("'[OH-]' has charge -1 and explicit H 1") {
  const MolecularGraph g = chem::parseSmiles("[OH-]");
  REQUIRE(g.atoms().size() == 1);
  CHECK(g.atoms()[0].charge == -1);
  CHECK(g.atoms()[0].explicit_h == 1);
}

TEST_CASE("Tetrafluoromethane 'C(F)(F)(F)F' — F implicit_h 0 each") {
  const MolecularGraph g = chem::parseSmiles("C(F)(F)(F)F");
  REQUIRE(g.atoms().size() == 5);
  CHECK(g.atoms()[0].implicit_h == 0); // C bond sum 4
  for (const Atom& f : g.atoms().subspan(1)) {
    CHECK(f.implicit_h == 0); // F bond sum 1, valence 1
  }
}

// ===========================================================================
// FR-12c: Valence model — each row + smallest-valence selection + over-valence
// ===========================================================================

TEST_CASE("Boron 'B' implicit H 3") {
  const MolecularGraph g = chem::parseSmiles("B");
  CHECK(g.atoms()[0].implicit_h == 3);
}

TEST_CASE("Phosphine 'P' implicit H 3 (smallest valence)") {
  const MolecularGraph g = chem::parseSmiles("P");
  CHECK(g.atoms()[0].implicit_h == 3);
}

TEST_CASE("Sulfur 'S' implicit H 2 (smallest valence, not 4 or 6)") {
  const MolecularGraph g = chem::parseSmiles("S");
  CHECK(g.atoms()[0].implicit_h == 2);
}

TEST_CASE("Nitrogen 'N' implicit H 3") {
  const MolecularGraph g = chem::parseSmiles("N");
  CHECK(g.atoms()[0].implicit_h == 3);
}

TEST_CASE("Fluorine 'F' implicit H 1") {
  const MolecularGraph g = chem::parseSmiles("F");
  CHECK(g.atoms()[0].implicit_h == 1);
}

TEST_CASE("Iodine 'I' implicit H 1") {
  const MolecularGraph g = chem::parseSmiles("I");
  CHECK(g.atoms()[0].implicit_h == 1);
}

TEST_CASE("S at bond sum 3 selects valence 4 (not 6)") {
  // CS(=O): S has C(single)+O(double) = 3 -> valence 4
  const MolecularGraph g = chem::parseSmiles("CS(=O)");
  REQUIRE(g.atoms().size() == 3);
  CHECK(g.atoms()[1].element.atomicNumber() == 16);
  CHECK(g.atoms()[1].implicit_h == 1); // valence 4 - bond sum 3 = 1
}

TEST_CASE("Over-valent 'CO(C)C' is rejected") { expectRejection("CO(C)C", "over-valent"); }

// ===========================================================================
// FR-12d: Rejection classes FR-10a-h
// ===========================================================================

TEST_CASE("FR-10a: Lowercase aromatic atoms rejected") {
  expectRejection("c1ccccc1", "aromatic");
  expectRejection("co", "aromatic");
}

TEST_CASE("FR-10b: Stereo markers and directional bonds rejected") {
  expectRejection("C/C=C/C", "directional");
  expectRejection("C\\C=C\\C", "directional");
  expectRejection("[C@H](N)C", "stereo");
  expectRejection("[C@@H](N)C", "stereo");
}

TEST_CASE("FR-10c: Disconnected components rejected") {
  expectRejection("CC.CC", "disconnected");
  expectRejection("[Na+].[Cl-]", "disconnected");
}

TEST_CASE("FR-10d: Reaction SMILES rejected") { expectRejection("CC>>CC", "reaction"); }

TEST_CASE("FR-10e: Wildcards and atom-class labels rejected") {
  expectRejection("C*C", "wildcard");
  expectRejection("[*]", "wildcard");
  expectRejection("[C:1]", "atom-class");
}

TEST_CASE("FR-10f: Unknown symbols and invalid bracket contents rejected") {
  expectRejection("[Hx]", "unknown element");
  expectRejection("[]", "element symbol");
  expectRejection("[123]", "element symbol");
  expectRejection("[OH-", "unterminated");
  expectRejection("Fe", "brackets");
}

TEST_CASE("FR-10g: Structural errors rejected") {
  expectRejection("", "empty");
  expectRejection("C1CC", "unclosed ring");
  expectRejection("C C", "whitespace");
  expectRejection("C	C", "whitespace");
}

TEST_CASE("FR-10h: No silent skipping — unrecognized bytes are errors") {
  expectRejection("C$", "unexpected");
}

TEST_CASE("Edge: '[Fe+++]' rejected — run form limited to ++/--") {
  expectRejection("[Fe+++]", "charge");
}

TEST_CASE("Edge: 'C11' rejected — same-atom immediate ring reuse") {
  expectRejection("C11", "ring");
}

TEST_CASE("Edge: 'C%1CC%1' rejected — % requires two digits") {
  expectRejection("C%1CC%1", "two digits");
}

TEST_CASE("Edge: 'C(C' rejected — unclosed branch") { expectRejection("C(C", "unclosed branch"); }

TEST_CASE("Edge: 'CC)' rejected — unbalanced closing paren") {
  expectRejection("CC)", "unbalanced");
}

TEST_CASE("Edge: '()' rejected — empty branch") { expectRejection("C()", "empty branch"); }

TEST_CASE("Edge: trailing bond 'C=' rejected") { expectRejection("C=", "trailing bond"); }

TEST_CASE("Edge: '[Fe+9]' valid (max charge magnitude)") {
  const MolecularGraph g = chem::parseSmiles("[Fe+9]");
  CHECK(g.atoms()[0].charge == 9);
}

TEST_CASE("Edge: '[Fe+10]' rejected (charge magnitude > 9)") {
  expectRejection("[Fe+10]", "charge magnitude");
}

TEST_CASE("Edge: '[1000C]' rejected — isotope > 999") { expectRejection("[1000C]", "isotope"); }

TEST_CASE("Edge: ring bond order mismatch 'C=1CC-1' rejected") {
  expectRejection("C=1CC-1", "mismatched");
}

TEST_CASE("Edge: single-sided ring bond order 'C1CC=1' is valid and sets double") {
  const MolecularGraph g = chem::parseSmiles("C1CC=1");
  REQUIRE(g.atoms().size() == 3);
  REQUIRE(g.bonds().size() == 3); // 2 chain + 1 closure
  bool found_closure = false;
  for (const auto& b : g.bonds()) {
    if ((b.a == 0 && b.b == 2) || (b.a == 2 && b.b == 0)) {
      found_closure = true;
      CHECK(b.order == BondOrder::kDouble);
    }
  }
  CHECK(found_closure);
}

// ===========================================================================
// FR-12e: Purity — parse twice, compare field-by-field
// ===========================================================================

TEST_CASE("Purity: parsing the same input twice yields identical graphs") {
  const std::array<std::string_view, 15> inputs = {
      "C",
      "O",
      "[OH2]",
      "CC(=O)O",
      "C1=CC=CC=C1",
      "C1CCC2CCCCC2C1",
      "C%10CCCCC%10",
      "[13CH4]",
      "[NH4+]",
      "[Fe+3]",
      "CS(=O)(=O)O",
      "C(F)(F)(F)F",
      "C#N",
      "Br",
      "[Fe++]",
  };
  for (std::string_view input : inputs) {
    const MolecularGraph first = chem::parseSmiles(input);
    const MolecularGraph second = chem::parseSmiles(input);
    CHECK_MESSAGE(sameGraph(first, second), "purity failed for: ", input);
  }
}

// ===========================================================================
// FR-12f: M1 regression — formula parser graphs still validate, isotope default 0
// ===========================================================================

TEST_CASE("M1 regression: formula-parser atoms have isotope 0") {
  const MolecularGraph g = chem::parseFormula("H2O");
  for (const Atom& a : g.atoms()) {
    CHECK(a.isotope == 0);
  }
}

TEST_CASE("M1 regression: composition and molar mass unchanged") {
  const MolecularGraph water = chem::parseFormula("H2O");
  const CompositionMap expected{{1, 2}, {8, 1}};
  CHECK(chem::composition(water) == expected);
  CHECK(chem::molarMass(water) == doctest::Approx(18.015).epsilon(0.001 / 18.015));
}
