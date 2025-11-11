#include "reaction.hpp"

#include <algorithm>
#include <iostream>
#include <cctype>
#include <charconv>

namespace chem {

Reaction::Reaction(std::initializer_list<Compound> reagents,
                   std::initializer_list<Compound> products, double k)
    : k_{k} {
  for (const auto &mol : reagents) {
    ++reagents_[mol];
  }
  for (const auto &mol : products) {
    ++products_[mol];
  }
}

Reaction::Reaction(std::string_view equation, double k) : k_{k} {
  auto eq_trimmed = std::string{};
  eq_trimmed.reserve(equation.size());
  std::copy_if(equation.begin(), equation.end(), std::back_inserter(eq_trimmed), [](auto ch) {
    return !std::isspace(ch);
  });

  constexpr auto kPlus{'+'};
  constexpr auto kEqual{'='};
  auto eq_pos {eq_trimmed.find(kEqual)};
  if (eq_pos == std::string::npos) {
    throw "Invalid equation: no '=' sign found.";
  }
  if (eq_pos + 1 == eq_trimmed.size()) {
    throw "Invalid equation: no right side of equation";
  }
  auto eq_sv {std::string_view(eq_trimmed)};
  auto left_side {eq_sv.substr(0, eq_pos)};
  auto right_side {eq_sv.substr(eq_pos + 1)};

  auto ParseEquationPart = [kPlus](std::string_view part, std::unordered_map<Compound, double>& mp) {
    size_t last {0};
    size_t pos;
    do {
      pos = part.find(kPlus, last);
      auto c_reag {part.substr(last, pos - last)};
      double coeff {1.0};
      auto formula {c_reag};
      auto [ptr, ec] = std::from_chars(c_reag.data(), c_reag.data() + c_reag.size(), coeff);
      if (ec == std::errc{} && ptr != c_reag.data()) {
        formula = c_reag.substr(ptr - c_reag.data());
      }
      auto compound {Compound(formula)};
      mp[compound] = coeff;

      last = pos + 1;
    } while (pos != std::string::npos);

  };
  ParseEquationPart(left_side, reagents_);
  ParseEquationPart(right_side, products_);
}

void Reaction::Print() const {
  std::cout << "Reagents" << '\n';
  for (const auto &[mol, count] : reagents_) {
    std::cout << count << ' ' << mol.Name() << ' ' << mol.Mass() << '\n';
  }
  std::cout << "Products" << '\n';
  for (const auto &[mol, count] : products_) {
    std::cout << count << ' ' << mol.Name() << ' ' << mol.Mass() << '\n';
  }
  std::cout.flush();
}

} // namespace chem