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
}

void Flask::RunReactions(double dt, double duration) {
  for (const auto& r : reactions_) {
    for (const auto& [compound, coeff] : r.Reagents()) {
      concentrations_[compound];
    }
    for (const auto& [compound, coeff] : r.Products()) {
      concentrations_[compound];
    }
  }

  constexpr char delim {','};
  auto out_file {std::ofstream("data.csv")};
  out_file << "t" << delim;
  for (const auto& [compound, _] : concentrations_) {
    out_file << compound.Name() << delim;
  }
  out_file << '\n';

  auto& old_conc {concentrations_};
  auto new_conc {concentrations_};
  for (double t = 0.0; t < duration; t += dt) {
    out_file << t << delim;
    for (const auto& [_, conc] : old_conc) {
      out_file << conc << delim;
    }
    out_file << '\n';
    for (const auto& reaction : reactions_) {
      auto k {reaction.RateConst()};
      double conc_mult {1.0};
      for (const auto& [compound, coeff] : reaction.Reagents()) {
        conc_mult *= std::pow(old_conc[compound], coeff);
      }
      auto rate {k * conc_mult};
      for (const auto& [compound, coeff] : reaction.Reagents()) {
        new_conc[compound] -= coeff * rate;
      }
      for (const auto& [compound, coeff] : reaction.Products()) {
        new_conc[compound] += coeff * rate;
      }
    }

    old_conc.swap(new_conc);
    for (auto old_it {old_conc.begin()}, new_it {new_conc.begin()};
         new_it != new_conc.end();
         ++new_it, ++old_it) {
      new_it->second = old_it->second;
    }
  }
}

}  // namespace chem