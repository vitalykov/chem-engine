#include "flask.hpp"

#include <cmath>
#include <iostream>
#include <fstream>

namespace chem {

Flask::Flask(std::initializer_list<Substance> substances,
        std::initializer_list<Reaction> reactions, double volume)
        : reactions_{reactions}, volume_{volume} {
  for (const auto& substance : substances) {
    concentrations_[substance.GetCompound()] += substance.Moles() / volume;
  }
  for (const auto& r : reactions) {
    for (const auto& [compound, coeff] : r.Reagents()) {
      concentrations_[compound];
    }
    for (const auto& [compound, coeff] : r.Products()) {
      concentrations_[compound];
    }
  }
}

void Flask::RunReactions(double dt, double duration) {
  constexpr char delim {','};
  auto out_file {std::ofstream("data.csv")};
  out_file << "t" << delim;
  for (const auto& [compound, _] : concentrations_) {
    out_file << compound.Name() << delim;
  }
  out_file << '\n';
  for (double t = 0.0; t < duration; t += dt) {
    out_file << t << delim;
    for (const auto& reaction : reactions_) {
      for (const auto& [_, conc] : concentrations_) {
        out_file << conc << delim;
      }
      out_file << '\n';
      auto k {reaction.K()};
      double conc_product {1.0};
      for (const auto& [compound, coeff] : reaction.Reagents()) {
        conc_product *= std::pow(concentrations_[compound], coeff);
      }
      for (const auto& [compound, coeff] : reaction.Reagents()) {
        concentrations_[compound] -= coeff * k * conc_product;
      }
      for (const auto& [compound, coeff] : reaction.Products()) {
        concentrations_[compound] += coeff * k * conc_product;
      }
    }
  }
}

}  // namespace chem