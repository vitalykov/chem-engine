#pragma once

#include <initializer_list>
#include <vector>
#include <functional>
#include <string>
#include <string_view>

#include "atom.hpp"

namespace chem {

class Compound {
public:
  Compound(std::initializer_list<Atom> atoms, const std::string& name);
  Compound(std::string_view formula);
  Compound(std::string_view formula, std::string_view name);
  inline std::vector<Atom> GetAtoms() const { return atoms_; }
  inline auto MolMass() const { return mol_mass_; }
  inline auto Name() const { return name_; }
  inline bool operator==(const Compound& other) const { return mol_mass_ == other.mol_mass_; }

private:
  std::vector<Atom> atoms_;
  double mol_mass_ {0.0};
  std::string name_;
};

}  // namespace chem

template<>
struct std::hash<chem::Compound> {
  std::size_t operator()(const chem::Compound& mol) const noexcept {
    return std::hash<double>{}(mol.MolMass());
  }
};
