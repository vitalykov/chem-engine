// #include "atom.hpp"
// #include "compound.hpp"
// #include "substance.hpp"
#include "reaction.hpp"
#include "flask.hpp"

int main() {
  // auto reaction {chem::Reaction("2 NaOH + CO2 = Na2CO3 + H2O")};
  // reaction.Print();
  // chem::Reaction("2 Al(OH)3 + 3 H2SO4 = Al2(SO4)3 + 6 H2O").Print();
  // chem::Reaction("4 FeCl3 + 3 K4Fe(CN)6 = Fe4(Fe(CN)6)3 + 12 KCl").Print();
  // auto water {chem::Substance("H2O", 2.0)};
  // auto carbon_dioxide {chem::Substance("CO2", 3.0)};
  // auto hydroxide {chem::Substance("NaOH", 4.0)};
  auto flask {chem::Flask{}};
  flask.AddSubstance("H2", 6.0);
  flask.AddSubstance("N2", 2.0);
  // flask.AddSubstance("NH3", 10.0);
  flask.AddReaction("3H2 + N2 = 2NH3", 1.0e-3);
  flask.AddReaction("2NH3 = 3H2 + N2", 2.5e-7);
  flask.AddSubstance("O2", 3.0);
  flask.AddReaction("2 H2 + O2 = 2 H2O", 5.0e-3);
  flask.AddReaction("2 H2O = 2 H2 + O2", 2.0e-3);
  flask.RunReactions(0.1, 100.0);

  return 0;
}