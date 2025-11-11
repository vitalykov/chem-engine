#include "atom.hpp"

#include "elements.hpp"

namespace chem {

Atom::Atom(std::string_view symbol) {
  static auto elements {Elements()};
  index_ = elements.Index(symbol);
  mass_ = elements.Mass(symbol);
}

}  // namespace chem
