#pragma once

#include <initializer_list>
#include <unordered_map>
#include <string>
#include <string_view>

#include "compound.hpp"

namespace chem {

class Reaction {
public:
  Reaction(std::initializer_list<Compound> reagents,
           std::initializer_list<Compound> products, double k = 0.01);
  Reaction(std::string_view equation, double k = 0.01);

  inline const auto& Reagents() const noexcept { return reagents_; }
  inline const auto& Products() const noexcept { return products_; }
  inline const auto& RateConst() const noexcept { return k_; }

  void Print() const;

private:
  std::unordered_map<Compound, double> reagents_;
  std::unordered_map<Compound, double> products_;
  double k_;
};

} // namespace chem