#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "chem/core/element.hpp"

namespace chem {

enum class BondOrder : std::uint8_t { Single, Double, Triple };

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

  std::ptrdiff_t add_atom(Atom atom);
  void add_bond(std::ptrdiff_t a, std::ptrdiff_t b, BondOrder order);

  [[nodiscard]] std::span<const Atom> atoms() const noexcept;
  [[nodiscard]] std::span<const Bond> bonds() const noexcept;

private:
  std::vector<Atom> atoms_;
  std::vector<Bond> bonds_;
};

} // namespace chem
