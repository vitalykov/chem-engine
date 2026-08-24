#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "chem/core/element.hpp"

namespace chem {

enum class BondOrder : std::uint8_t { kSingle, kDouble, kTriple };

struct Atom {
  Element element;
  int charge = 0;
  int implicit_h = 0;
  int explicit_h = 0;
};

struct Bond {
  std::ptrdiff_t a;
  std::ptrdiff_t b;
  BondOrder order;
};

class MolecularGraph {
public:
  MolecularGraph() = default;

  std::ptrdiff_t addAtom(Atom atom);
  void addBond(std::ptrdiff_t a, std::ptrdiff_t b, BondOrder order);

  [[nodiscard]] std::span<const Atom> atoms() const noexcept;
  [[nodiscard]] std::span<const Bond> bonds() const noexcept;

private:
  std::vector<Atom> atoms_;
  std::vector<Bond> bonds_;
};

} // namespace chem
