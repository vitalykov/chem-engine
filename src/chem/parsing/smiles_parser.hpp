#pragma once

#include <string_view>

#include "chem/core/molecular_graph.hpp"

namespace chem {

// Parses one connected SMILES molecule into a bonded graph.
// Throws ParseError on any invalid or unsupported input.
[[nodiscard]] MolecularGraph parseSmiles(std::string_view input);

} // namespace chem
