#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include "atom.hpp"
#include "compound.hpp"
#include "flask.hpp"
#include "reaction.hpp"
#include "substance.hpp"

#include "elements.hpp"

std::vector<std::string> ParseInput(std::string &input) {
  std::vector<std::string> result;
  std::string delimeters{" +="};
  size_t start{input.find_first_not_of(delimeters)};
  while (start != std::string::npos) {
    size_t sep{input.find_first_of(delimeters, start)};
    result.push_back(input.substr(start, sep - start));
    start = input.find_first_not_of(delimeters, sep);
  }

  return result;
}

struct ChemComponent {
  int num;
  int count;
};

// std::vector<ChemComponent> ParseComponent(const std::string &component,
//                                           chem::Elements &elements) {
//   static std::regex pattern{R"((\d*[\.,]*\d*)([A-Z][a-z]?\d*[A-Za-z0-9]*))"};
//   constexpr size_t kSymbol{1};
//   constexpr size_t kIndex{2};
//   std::vector<ChemComponent> components;
//   for (std::sregex_iterator pit(component.begin(), component.end(), pattern);
//        pit != std::sregex_iterator{}; ++pit) {
//     auto symbol{(*pit)[0]};
//     std::cout << symbol << ':' << elements.Index(symbol) << ' ';
//     std::cout << (*pit)[kIndex] << ' ';
//   }
//   std::cout << '\n';

//   return components;
// }

int main() {
  // auto C {chem::Atom("C")};
  // auto H {chem::Atom("H")};
  // auto O {chem::Atom("O")};
  // auto Na {chem::Atom("Na")};
  // auto NaOH {chem::Compound("NaOH")};
  // auto CO2 {chem::Compound("CO2")};
  // auto Na2CO3 {chem::Compound("(NH4)2CO3")};
  // auto H2O {chem::Compound("H2O")};
  auto reaction {chem::Reaction("2 NaOH + CO2 = Na2CO3 + H2O")};
  reaction.Print();
  chem::Reaction("2 Al(OH)3 + 3 H2SO4 = Al2(SO4)3 + 6 H2O").Print();
  chem::Reaction("4 FeCl3 + 3 K4Fe(CN)6 = Fe4(Fe(CN)6)3 + 12 KCl").Print();
  // auto water {chem::Substance(H2O, 2.0)};
  // auto carbon_dioxide {chem::Substance(CO2, 3.0)};
  // auto hydroxide {chem::Substance(NaOH, 4.0)};
  // auto soda {chem::Substance(Na2CO3, 0)};
  // auto flask {chem::Flask({water, carbon_dioxide, hydroxide}, {reaction})};
  // flask.RunReactions();
  // auto elements = chem::Elements();

  // auto input = std::string{"2 NaOH + CO2 = Na2CO3 + H2O"};
  // auto components {ParseInput(input)};
  // for (const auto& com : components) {
  //   std::cout << com << '\n';
  //   ParseComponent(com, elements);
  // }


  return 0;
}