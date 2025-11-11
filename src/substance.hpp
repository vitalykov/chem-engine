#pragma once

#include "compound.hpp"

namespace chem {

class Substance {
public:
  Substance(const Compound& compound, double moles) : compound_{compound}, moles_{moles} {}
  inline const Compound& GetCompound() const noexcept { return compound_; }
  inline double Moles() const noexcept { return moles_; }
  inline auto Name() const noexcept { return GetCompound().Name(); }
  inline bool operator==(const Substance& other) const { return GetCompound() == other.GetCompound(); }

private:
  const Compound compound_;
  double moles_;
};

}  // namespace chem_engine

template<>
struct std::hash<chem::Substance> {
  std::size_t operator()(const chem::Substance& substance) const noexcept {
    return std::hash<double>{}(substance.GetCompound().Mass());
  }
};
