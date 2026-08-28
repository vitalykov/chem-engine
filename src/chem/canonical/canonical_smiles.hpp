#pragma once

#include <cstdint>
#include <string>

#include "chem/core/molecular_graph.hpp"

namespace chem {

// Version of the frozen canonicalization specification (docs/specs/m2b_canonicalization_spec.md).
// Any change to canonicalSmiles output requires a formal version bump.
inline constexpr std::uint32_t kCanonicalSpecVersion = 1;

// Pure function: graph -> canonical SMILES (spec version kCanonicalSpecVersion).
// Deterministic and side-effect free: equal graphs yield identical strings
// regardless of atom/bond insertion order.
// Throws ValidationError on an empty graph.
[[nodiscard]] std::string canonicalSmiles(const MolecularGraph& graph);

} // namespace chem