#include "compound.hpp"

namespace chem {

Compound::Compound(std::initializer_list<Atom> atoms, const std::string& name) : atoms_{atoms}, name_{name} {
  for (const auto& atom : atoms) {
    mol_mass_ += atom.GetMass();
  }
}

}  // namespace chem