#include <doctest/doctest.h>

#include <map>

#include "chem/core/composition.hpp"
#include "chem/core/element.hpp"
#include "chem/core/molecular_graph.hpp"

namespace {

using chem::Atom;
using chem::BondOrder;
using chem::Element;
using chem::MolecularGraph;

MolecularGraph waterGraph() {
  MolecularGraph graph;
  graph.addAtom(Atom{.element = Element("O")});
  graph.addAtom(Atom{.element = Element("H")});
  graph.addAtom(Atom{.element = Element("H")});
  return graph;
}

TEST_CASE("addAtom returns sequential indices and stores fields") {
  MolecularGraph graph;
  CHECK(graph.addAtom(Atom{.element = Element("O")}) == 0);
  CHECK(graph.addAtom(Atom{.element = Element("H")}) == 1);

  REQUIRE(graph.atoms().size() == 2);
  CHECK(graph.atoms()[0].element.atomicNumber() == 8);
  CHECK(graph.atoms()[1].element.atomicNumber() == 1);
}

TEST_CASE("Atom hydrogen fields default to zero") {
  const Atom atom{.element = Element("C")};
  CHECK(atom.charge == 0);
  CHECK(atom.implicit_h == 0);
  CHECK(atom.explicit_h == 0);
}

TEST_CASE("addBond stores endpoints and order") {
  MolecularGraph graph;
  graph.addAtom(Atom{.element = Element("O")});
  graph.addAtom(Atom{.element = Element("O")});

  graph.addBond(0, 1, BondOrder::kSingle);
  REQUIRE(graph.bonds().size() == 1);
  CHECK(graph.bonds()[0].a == 0);
  CHECK(graph.bonds()[0].b == 1);
  CHECK(graph.bonds()[0].order == BondOrder::kSingle);
}

TEST_CASE("Composition aggregates repeated atoms") {
  const MolecularGraph graph = waterGraph();
  const chem::CompositionMap counts = chem::composition(graph);
  const chem::CompositionMap expected{{1, 2}, {8, 1}};
  CHECK(counts == expected);
}

TEST_CASE("Composition includes implicit and explicit hydrogens") {
  MolecularGraph graph;
  graph.addAtom(Atom{.element = Element("O"), .charge = 0, .implicit_h = 1, .explicit_h = 1});
  graph.addAtom(Atom{.element = Element("Cl"), .implicit_h = 2, .explicit_h = 3});

  const chem::CompositionMap counts = chem::composition(graph);
  const chem::CompositionMap expected{{1, 7}, {8, 1}, {17, 1}};
  CHECK(counts == expected);
}

TEST_CASE("Empty graph has empty composition and zero mass") {
  const MolecularGraph graph;
  CHECK(chem::composition(graph).empty());
  CHECK(chem::molarMass(graph) == doctest::Approx(0.0));
}

TEST_CASE("Molar mass of water is approximately 18.015 g/mol") {
  CHECK(chem::molarMass(waterGraph()) == doctest::Approx(18.015).epsilon(0.001 / 18.015));
}

TEST_CASE("Molar mass matches independent recomputation over elements") {
  MolecularGraph glucose;
  for (int i = 0; i < 6; ++i) {
    glucose.addAtom(Atom{.element = Element("C")});
  }
  for (int i = 0; i < 12; ++i) {
    glucose.addAtom(Atom{.element = Element("H")});
  }
  for (int i = 0; i < 6; ++i) {
    glucose.addAtom(Atom{.element = Element("O")});
  }

  double expected = 6 * Element("C").standardWeight();
  expected += 12 * Element("H").standardWeight();
  expected += 6 * Element("O").standardWeight();

  const double mass = chem::molarMass(glucose);
  CHECK(mass == doctest::Approx(expected).epsilon(1e-12));
  CHECK(mass == doctest::Approx(180.156).epsilon(0.001 / 180.156));
}

} // namespace
