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

MolecularGraph water_graph() {
  MolecularGraph graph;
  graph.add_atom(Atom{.element = Element("O")});
  graph.add_atom(Atom{.element = Element("H")});
  graph.add_atom(Atom{.element = Element("H")});
  return graph;
}

TEST_CASE("add_atom returns sequential indices and stores fields") {
  MolecularGraph graph;
  CHECK(graph.add_atom(Atom{.element = Element("O")}) == 0);
  CHECK(graph.add_atom(Atom{.element = Element("H")}) == 1);

  REQUIRE(graph.atoms().size() == 2);
  CHECK(graph.atoms()[0].element.atomic_number() == 8);
  CHECK(graph.atoms()[1].element.atomic_number() == 1);
}

TEST_CASE("Atom hydrogen fields default to zero") {
  const Atom atom{.element = Element("C")};
  CHECK(atom.charge == 0);
  CHECK(atom.implicit_h == 0);
  CHECK(atom.explicit_h == 0);
}

TEST_CASE("add_bond stores endpoints and order") {
  MolecularGraph graph;
  graph.add_atom(Atom{.element = Element("O")});
  graph.add_atom(Atom{.element = Element("O")});

  graph.add_bond(0, 1, BondOrder::Single);
  REQUIRE(graph.bonds().size() == 1);
  CHECK(graph.bonds()[0].a == 0);
  CHECK(graph.bonds()[0].b == 1);
  CHECK(graph.bonds()[0].order == BondOrder::Single);
}

TEST_CASE("Composition aggregates repeated atoms") {
  const MolecularGraph graph = water_graph();
  const chem::CompositionMap counts = chem::composition(graph);
  const chem::CompositionMap expected{{1, 2}, {8, 1}};
  CHECK(counts == expected);
}

TEST_CASE("Composition includes implicit and explicit hydrogens") {
  MolecularGraph graph;
  graph.add_atom(Atom{.element = Element("O"), .charge = 0, .implicit_h = 1, .explicit_h = 1});
  graph.add_atom(Atom{.element = Element("Cl"), .implicit_h = 2, .explicit_h = 3});

  const chem::CompositionMap counts = chem::composition(graph);
  const chem::CompositionMap expected{{1, 7}, {8, 1}, {17, 1}};
  CHECK(counts == expected);
}

TEST_CASE("Empty graph has empty composition and zero mass") {
  const MolecularGraph graph;
  CHECK(chem::composition(graph).empty());
  CHECK(chem::molar_mass(graph) == doctest::Approx(0.0));
}

TEST_CASE("Molar mass of water is approximately 18.015 g/mol") {
  CHECK(chem::molar_mass(water_graph()) == doctest::Approx(18.015).epsilon(0.001 / 18.015));
}

TEST_CASE("Molar mass matches independent recomputation over elements") {
  MolecularGraph glucose;
  for (int i = 0; i < 6; ++i) {
    glucose.add_atom(Atom{.element = Element("C")});
  }
  for (int i = 0; i < 12; ++i) {
    glucose.add_atom(Atom{.element = Element("H")});
  }
  for (int i = 0; i < 6; ++i) {
    glucose.add_atom(Atom{.element = Element("O")});
  }

  double expected = 6 * Element("C").standard_weight();
  expected += 12 * Element("H").standard_weight();
  expected += 6 * Element("O").standard_weight();

  const double mass = chem::molar_mass(glucose);
  CHECK(mass == doctest::Approx(expected).epsilon(1e-12));
  CHECK(mass == doctest::Approx(180.156).epsilon(0.001 / 180.156));
}

} // namespace
