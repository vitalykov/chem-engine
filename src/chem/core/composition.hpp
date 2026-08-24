#pragma once

#include <map>

#include "chem/core/molecular_graph.hpp"

namespace chem {

using CompositionMap = std::map<int, int>;

[[nodiscard]] CompositionMap composition(const MolecularGraph& graph);
[[nodiscard]] double molar_mass(const MolecularGraph& graph);

} // namespace chem
