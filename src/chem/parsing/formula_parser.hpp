#pragma once

#include <string_view>

#include "chem/core/molecular_graph.hpp"

namespace chem {

// Parses Hill-formula notation into an unconnected atom multiset.
// Throws ParseError on any invalid input.
[[nodiscard]] MolecularGraph parseFormula(std::string_view input);

} // namespace chem
