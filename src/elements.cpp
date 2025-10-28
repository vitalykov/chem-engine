#include "elements.hpp"

#include <string>
#include <fstream>

namespace chem {

Elements::Elements() {
  std::ifstream file {"elements"};
  std::string symbol;
  int n;
  double mass;
  while (file >> symbol >> n >> mass) {
    elements_[symbol] = ElementInfo{.index = n, .mass = mass};
  }
}

// Elements::ElementInfo Elements::GetElementInfo(const std::string& symbol) const {
//   return elements_.at(symbol);
// }

// int Elements::GetElementNumber(const std::string& symbol) const {
//   return elements_.at(symbol).index;
// }

// double Elements::GetElementMass(const std::string& symbol) const {
//   return elements_.at(symbol).mass;
// }
  
}  // namespace chem
