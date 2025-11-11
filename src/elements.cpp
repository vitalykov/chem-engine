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
  
}  // namespace chem
