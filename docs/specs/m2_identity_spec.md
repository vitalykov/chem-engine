# Spec: M2 Identity

**Status**: Review

**Author**: Agent + Vitaly Kovalev

**Date**: 2026-08-25

**Related**: docs/design.md §3.2 (Molecule), §4 (identity contract), §7 (error handling), §8 (architecture), §9 (testing strategy), §10 (roadmap, milestone M2); docs/specs/m1_foundations_spec.md (delivered foundation)

## 1. Goal / Problem Statement

Milestone M2 (docs/design.md §10) delivers species identity: a SMILES parser
feeding the existing `MolecularGraph`, a canonicalizer whose output is frozen
by a versioned specification and a golden corpus, round-trip guarantees, and
the `chem::Molecule` type with equality/hash defined on canonical SMILES.

This milestone turns the engine from "parses inputs" into "identifies
species": structural isomers become distinct, equivalent spellings collapse to
one ID, and the species identifier comes into existence under the freezing
rules of docs/design.md §4.2.

Because identity spans three fairly independent subsystems, this umbrella spec
is accompanied by three sub-specs:

| Sub-spec | Subsystem |
|---|---|
| docs/specs/m2a_smiles_parser_spec.md | SMILES dialect subset → `MolecularGraph` |
| docs/specs/m2b_canonicalization_spec.md | Canonicalizer, golden corpus, round-trip |
| docs/specs/m2c_molecule_spec.md | `Molecule`, factories, formula→structure resolution |

## 2. Scope

### In Scope

- SMILES parser covering the v1 dialect subset (design.md §4.3)
- Canonical ranking (Weisfeiler-Lehman refinement) + canonical SMILES writer,
  pinned as canonicalization spec **version 1** (`kCanonicalSpecVersion = 1`)
- Golden corpus (~30 entries) locking canonical output on every build
- Round-trip invariant tests across the whole corpus
- `chem::Molecule` with explicit `fromFormula` / `fromSmiles` factories and
  equality/hash on canonical form
- Formula-to-connectivity resolution: unique valence-based inference, curated
  unambiguous-compound dictionary, loud ambiguity errors with alternatives
- `Atom` extension with an isotope field (required by the dialect)

### Out of Scope

- Stereochemistry (`@`, `/`, `\`) — parsed inputs are rejected loudly; no perception
- Aromaticity perception — kekulé alternating bonds only; lowercase aromatic
  atoms rejected loudly
- Dot-disconnected SMILES inputs (`[Na+].[Cl-]`) — rejected loudly in the
  parser; disconnected graphs arise only internally (formula fallback path,
  see m2c)
- Reactions, equation parsing, balance validation (M3)
- `KineticSystem`, `Vessel`, integrators (M3)
- Serialization of molecules (M4)
- Consteval compile-time molecule literals (deferred; Open Questions in m2c)
- Wildcards, reaction SMILES, atom mapping numbers, atom-class labels

## 3. User Stories or Use Cases

1. As a library user, I want `Molecule::fromSmiles("[OH2]") ==
   Molecule::fromSmiles("O")` so that equivalent spellings identify one species.
2. As a library user, I want ethanol (`CCO`) and dimethyl ether (`COC`) to be
   distinct molecules despite identical composition so that kinetics over
   mixtures stays chemically correct.
3. As a downstream contributor building the reaction layer (M3), I want
   `Molecule` values comparable and hashable so that reaction participants can
   live in maps keyed by identity.
4. As a library developer building the reaction layer (M3), I want equivalent
   inputs to yield equal `Molecule`s with equal hashes so that species
   collections contain no duplicates.

## 4. Functional Requirements

Requirements are distributed to the sub-specs; each is numbered there. This
umbrella carries only the cross-cutting rules:

### ⏳ FR-U1: Pipeline integrity

Every identifier MUST flow through the single pipeline of design.md §4.1:
input → parse → `MolecularGraph` → canonicalize → canonical SMILES string. No
format-specific shortcut to identity MAY exist (hard architectural rule 3).

### ⏳ FR-U2: Loud rejection

Unknown/unsupported syntax in any format MUST raise the typed exception named
in the relevant sub-spec, never silent acceptance (hard architectural rule 4).

### ⏳ FR-U3: Canonicalization spec version

The canonicalization output produced under these specs MUST be identified as
spec version 1 via `chem::kCanonicalSpecVersion` (defined in m2b FR-1). Any
future change to output requires a formal version bump per docs/design.md §4.2.

### ⏳ FR-U4: Identity guarantees

The published contract of design.md §4.4 MUST hold at the end of M2:

- Normalization: equivalent inputs yield identical canonical SMILES.
- Round-trip: canonical output re-parses to an identical graph.
- Comparison: `Molecule` equality/hashing are defined on canonical form.

For formula-derived molecules, normalization applies within a format and
across successfully resolved formulas (see m2c §6 for the documented caveat).

## 5. Non-Functional Requirements

- **Determinism**: canonicalization is a pure function of the graph; same
  input graph, same output string, always.
- **Correctness over speed**: no performance budget in M2 beyond "no obvious
  pathology"; refinement and traversal costs are acceptable at real-time
  scales for small molecules.
- **Zero dependencies**: only doctest (FetchContent) beyond the toolchain.
- **Testability**: golden corpus and round-trip suites run on every build via
  the existing `test-runner` target and CTest registration.
- **Safety**: RAII everywhere, ownership by value/containers, non-owning
  access via `std::span`/`std::string_view`, named casts only; debug/test runs
  verified under ASan+UBSan (AGENTS.md).

## 6. Technical Constraints & Architecture Notes

- Namespace `chem`; naming per AGENTS.md; `.clang-format` applied to touched files.
- Parsing and canonicalization stay separate stages with the graph between
  them (design.md §8); neither reaches into the other. The SMILES parser MUST
  NOT special-case canonical ordering; the canonicalizer MUST NOT contain
  dialect parsing logic beyond emitting strings that satisfy it.
- New code lands in `src/chem/parsing/`, `src/chem/canonical/`,
  `src/chem/core/molecule.*`; tests mirror this structure under `tests/`.
- Exceptions only in setup phase (all of M2 is setup phase); typed exceptions
  derived from `std::runtime_error`, never string literals (design.md §7).
- The proto branch is reference-only; none of its identity bugs (hashing by
  molar mass, composition-based equality) may carry over.

## 7. Data & Interfaces

See the Data & Interfaces sections of m2a (parser entry point, `Atom`
extension), m2b (canonical entry point, version constant, corpus file), and
m2c (`Molecule` type, compounds dictionary file).

New public headers added to the `chem` target's `FILE_SET HEADERS`:

```
src/chem/canonical/canonical_smiles.hpp
src/chem/core/molecule.hpp
```

## 8. Edge Cases & Error Handling

Distributed to the sub-specs: dialect edge cases in m2a §8, canonical edge
cases in m2b §8, formula-resolution ambiguity in m2c §8.

## 9. Acceptance Criteria / Definition of Done

- [ ] All functional requirements of m2a, m2b, m2c satisfied
- [ ] `cmake --build build --parallel` compiles warning-clean under
      `-Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wshadow -Werror`
- [ ] `ctest --test-dir build --output-on-failure` passes all suites, debug
      build running under ASan+UBSan
- [ ] Golden corpus suite passes unchanged twice in a row (determinism smoke)
      and round-trip suite passes over the whole corpus
- [ ] `clang-format -i` applied to touched sources; `clang-tidy` reports no
      new findings on changed files
- [ ] No new dependencies in CMakeLists.txt
- [ ] Public headers changed → README usage section updated
- [ ] Layout matches the sub-spec directory requirements; AGENTS.md conventions respected

## 10. Open Questions / Decisions Needed

None blocking. Deferred decisions recorded in m2c §10 (consteval literal,
compile-time checking strategy).

## 11. Implementation Notes

- Decisions locked during drafting (user-approved):
  - Explicit factories instead of an auto-detecting constructor: `Molecule`
    never sniffs input format; ambiguous strings like `"C"` or `"NO"` are
    routed explicitly.
  - DayLight-style valence model for implicit hydrogens and formula
    connectivity inference.
  - Canonicalization algorithm pinned precisely (invariants, refinement,
    tie-breaking) rather than specified contract-only.
  - Golden corpus starts at ~30 curated entries.
  - Compound dictionary restricted to unambiguous, time-proven formulas
    (e.g. `C2H5OH` for ethanol); ambiguous compositions throw with suggested
    alternative writings.
