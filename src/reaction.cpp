#include "reaction.hpp"

#include <algorithm>
#include <iostream>

namespace chem {

Reaction::Reaction(std::initializer_list<Compound> reagents,
                   std::initializer_list<Compound> products, double k)
    : k_{k} {
  for (const auto &mol : reagents) {
    ++reagents_[mol];
  }
  for (const auto &mol : products) {
    ++products_[mol];
  }
}

void Reaction::Print() const {
  std::cout << "Reagents" << '\n';
  for (const auto &[mol, count] : reagents_) {
    std::cout << count << ' ' << mol.MolMass() << '\n';
  }
  std::cout << "Products" << '\n';
  for (const auto &[mol, count] : products_) {
    std::cout << count << ' ' << mol.MolMass() << '\n';
  }
}

} // namespace chem