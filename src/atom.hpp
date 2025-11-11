#pragma once

#include <string>

namespace chem {

class Atom {
public:
  Atom(std::string_view symbol);
  inline auto Mass() const { return mass_; }
  inline auto Index() const { return index_; }

private:
  int index_;
  double mass_;
};

}  // namespace chem
