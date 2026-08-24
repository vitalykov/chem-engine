#include "chem/core/molecular_graph.hpp"

#include <cassert>
#include <iterator>

namespace chem {

std::ptrdiff_t MolecularGraph::add_atom(Atom atom) {
  const auto index = static_cast<std::ptrdiff_t>(atoms_.size());
  atoms_.push_back(atom);
  return index;
}

void MolecularGraph::add_bond(std::ptrdiff_t a, std::ptrdiff_t b, BondOrder order) {
  assert(a >= 0 && a < std::ssize(atoms_));
  assert(b >= 0 && b < std::ssize(atoms_));
  bonds_.push_back(Bond{a, b, order});
}

std::span<const Atom> MolecularGraph::atoms() const noexcept { return atoms_; }

std::span<const Bond> MolecularGraph::bonds() const noexcept { return bonds_; }

} // namespace chem
