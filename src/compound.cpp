#include "compound.hpp"

#include <stack>
#include <utility>
#include <charconv>

namespace chem {

Compound::Compound(std::initializer_list<Atom> atoms, const std::string &name)
    : atoms_{atoms}, name_{name} {
  for (const auto &atom : atoms) {
    mol_mass_ += atom.Mass();
  }
}

Compound::Compound(std::string_view formula, std::string_view name) : name_{name} {
  std::stack<std::pair<std::vector<Atom>, int>> group_stack;
  group_stack.push({{}, 1}); // Main group with multiplier 1

  size_t i = 0;
  while (i < formula.length()) {
    if (formula[i] == '(') {
      // Start a new group
      group_stack.push({{}, 1});
      i++;
    } else if (formula[i] == ')') {
      // End current group and apply multiplier
      auto [group_atoms, multiplier] = group_stack.top();
      group_stack.pop();

      // Read multiplier after ')'
      i++;
      int group_multiplier = 1;
      if (i < formula.length() && std::isdigit(formula[i])) {
        size_t num_end = i;
        while (num_end < formula.length() && std::isdigit(formula[num_end])) {
          num_end++;
        }
        auto num_sv {formula.substr(i, num_end - i)};
        std::from_chars(num_sv.data(), num_sv.data() + num_sv.size(), group_multiplier);
        i = num_end;
      }

      // Add group atoms to parent group
      for (const auto &atom : group_atoms) {
        for (int j = 0; j < group_multiplier; j++) {
          group_stack.top().first.push_back(atom);
          // mol_mass_ += atom.Mass();
        }
      }
    } else if (std::isupper(formula[i])) {
      // Parse element symbol
      size_t symbol_end = i + 1;
      while (symbol_end < formula.length() &&
             std::islower(formula[symbol_end])) {
        symbol_end++;
      }

      auto symbol = formula.substr(i, symbol_end - i);

      // Parse count
      int count = 1;
      if (symbol_end < formula.length() && std::isdigit(formula[symbol_end])) {
        size_t num_end = symbol_end;
        while (num_end < formula.length() && std::isdigit(formula[num_end])) {
          num_end++;
        }
        auto num_sv {formula.substr(symbol_end, num_end - symbol_end)};
        std::from_chars(num_sv.data(), num_sv.data() + num_sv.size(), count);
        i = num_end;
      } else {
        i = symbol_end;
      }

      // Add atoms to current group
      for (int j = 0; j < count; j++) {
        Atom atom(symbol);
        group_stack.top().first.push_back(atom);
        // mol_mass_ += atom.Mass();
      }
    } else {
      i++; // Skip other characters (shouldn't happen in valid formulas)
    }
  }

  atoms_ = group_stack.top().first;
  mol_mass_ = 0.0;
  for (const auto& atom : atoms_) {
    mol_mass_ += atom.Mass();
  }
}

Compound::Compound(std::string_view formula) : Compound(formula, formula) {}

} // namespace chem
