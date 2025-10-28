#include <iostream>
#include <vector>
#include <string>
#include <regex>

#include "atom.hpp"
#include "compound.hpp"
#include "reaction.hpp"
#include "flask.hpp"
#include "substance.hpp"

std::vector<std::string> ParseInput(std::string& input) {
  std::vector<std::string> result;
  std::string delimeters {" +="};
  size_t start {input.find_first_not_of(delimeters)};
  while (start != std::string::npos) {
    size_t sep {input.find_first_of(delimeters, start)};
    result.push_back(input.substr(start, sep - start));
    start = input.find_first_not_of(delimeters, sep);
  }

  return result;
}

struct ChemComponent {
  int num;
  int count;
};

// std::vector<ChemComponent> ParseComponent(const std::string& component, chem_engine::Elements& elements) {
//   static std::regex pattern {R"(([A-Z][a-z]?)(\d*))"};
//   constexpr size_t kSymbol {1};
//   constexpr size_t kIndex {2};
//   std::vector<ChemComponent> components;
//   for (std::sregex_iterator pit(component.begin(), component.end(), pattern); pit != std::sregex_iterator{}; ++pit) {
//     auto symbol {(*pit)[kSymbol]};
//     std::cout << symbol << ':' << elements.GetElementNumber(symbol) << ' ';
//     std::cout << (*pit)[kIndex] << ' ';
//   }
//   std::cout << '\n';

//   return components;
// }

int main() {
  auto C {chem::Atom("C")};
  auto H {chem::Atom("H")};
  auto O {chem::Atom("O")};
  auto Na {chem::Atom("Na")};
  auto NaOH {chem::Compound({Na, O, H}, "NaOH")};
  auto CO2 {chem::Compound({C, O, O}, "CO2")};
  auto Na2CO3 {chem::Compound({Na, Na, C, O, O, O}, "Na2CO3")};
  auto H2O {chem::Compound({H, H, O}, "H2O")};
  auto reaction {chem::Reaction({NaOH, NaOH, CO2}, {Na2CO3, H2O})};
  auto water {chem::Substance(H2O, 2.0)};
  auto carbon_dioxide {chem::Substance(CO2, 3.0)};
  auto hydroxide {chem::Substance(NaOH, 4.0)};
  auto soda {chem::Substance(Na2CO3, 0)};
  auto flask {chem::Flask({water, carbon_dioxide, hydroxide, soda}, {reaction})};
  flask.RunReactions();
  
  return 0;
}