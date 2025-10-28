#pragma once

#include <string>

namespace chem {

class Atom {
public:
  Atom(const std::string& symbol);
  inline auto GetMass() const { return mass_; }
  inline auto GetIndex() const { return index_; }

private:
  int index_;
  double mass_;
};

}  // namespace chem
