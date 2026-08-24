#include "chem/core/composition.hpp"

namespace chem {

CompositionMap composition(const MolecularGraph& graph) {
  CompositionMap counts;
  for (const Atom& atom : graph.atoms()) {
    ++counts[atom.element.atomic_number()];
    const int attached_h = atom.implicit_h + atom.explicit_h;
    if (attached_h != 0) {
      counts[1] += attached_h;
    }
  }
  return counts;
}

double molar_mass(const MolecularGraph& graph) {
  double total = 0.0;
  for (const auto& [atomic_number, count] : composition(graph)) {
    total += Element(atomic_number).standard_weight() * static_cast<double>(count);
  }
  return total;
}

} // namespace chem
