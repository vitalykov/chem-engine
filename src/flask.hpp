#pragma once

#include <initializer_list>
#include <vector>
#include <unordered_map>
#include <string_view>

#include "substance.hpp"
#include "reaction.hpp"

namespace chem {

class Flask {
public:
  Flask(std::initializer_list<Substance> substances,
        std::initializer_list<Reaction> reactions, double volume = 1.0);
  Flask() : volume_{1.0} {}

  inline void AddSubstance(const Substance &substance) {
    concentrations_[substance.GetCompound()] += substance.Moles() / volume_;
  }
  inline void AddSubstance(std::string_view formula, double moles) {
    concentrations_[Compound(formula)] += moles / volume_;
  }

  inline void AddReaction(const Reaction &reaction) {
    reactions_.push_back(reaction);
  }
  inline void AddReaction(std::string_view reaction, double k = 0.01) {
    reactions_.push_back(Reaction(reaction, k));
  }

  void RunReactions(double dt = 0.1, double duration = 10.0);

private:
  std::unordered_map<Compound, double> concentrations_;
  std::vector<Reaction> reactions_;
  double volume_; // liters
};

} // namespace chem
