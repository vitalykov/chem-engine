#pragma once

#include <initializer_list>
#include <vector>
#include <unordered_map>

#include "substance.hpp"
#include "reaction.hpp"

namespace chem {

class Flask {
public:
  Flask(std::initializer_list<Substance> substances,
        std::initializer_list<Reaction> reactions, double volume = 1.0);
  Flask() : volume_{1.0} {}

  inline void AddSubstrance(const Substance &substance) {
    concentrations_[substance.GetCompound()] += substance.Moles() / volume_;
  }

  inline void AddReaction(const Reaction &reaction) {
    reactions_.push_back(reaction);
  }

  void RunReactions(double dt = 0.1, double duration = 10.0);

private:
  std::unordered_map<Compound, double> concentrations_;
  std::vector<Reaction> reactions_;
  double volume_; // liters
};

} // namespace chem
