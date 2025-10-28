#include "atom.hpp"

#include "elements.hpp"

namespace chem {

Atom::Atom(const std::string& symbol) {
  static auto elements {Elements()};
  index_ = elements.GetElementNumber(symbol);
  mass_ = elements.GetElementMass(symbol);
}

}  // namespace chem
