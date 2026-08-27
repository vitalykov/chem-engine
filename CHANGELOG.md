# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to Semantic Versioning. During the 0.x series,
minor versions may still introduce breaking API changes; promotion to 1.0.0
requires the criteria in AGENTS.md ("Versioning").

## [Unreleased]

Canonicalization spec: v1 (`docs/design.md` §4.2-4.3), unchanged from the
previous release — the specification is defined but not yet implemented.

### Added

- SMILES parser (`chem::parseSmiles`) covering the v1 dialect subset: bare
  organic atoms (`B C N O P S F Cl Br I`), bracketed atoms (all 118 elements)
  with isotope labels and formal charges, single/double/triple bonds,
  branches, ring closures (single-digit and `%nn`), and Daylight-style
  implicit-hydrogen valence inference. All unsupported syntax (aromatic
  atoms, stereochemistry, disconnected components, wildcards, reaction
  SMILES) is rejected loudly with offset-quoting `ParseError`s.
- `Atom::isotope` field (mass number; 0 = unspecified), defaulting to 0 so
  existing formula-parser output and M1 behavior are unchanged.

## [0.1.0] - 2026-08-25

First release: M1 foundations — element data, molecular graphs,
composition and molar mass, Hill-formula parsing.

Canonicalization spec: v1 (`docs/design.md` §4.2-4.3), unchanged from
the previous release — the specification is defined but not yet
implemented in this release.

### Added

- Periodic table with flyweight `Element` handles: case-sensitive lookup
  by symbol or atomic number; unknown tokens throw `ParseError`; table
  data loads once per process from `data/elements.csv`.
- `MolecularGraph` type with atoms, typed bonds (`BondOrder`), and
  contiguous read-only views.
- Composition map (atomic number -> count, including implicit/explicit
  hydrogens) and molar mass computation over graphs.
- Hill-formula parser emitting `MolecularGraph`; all unsupported syntax
  is rejected loudly with position-quoting `ParseError`s.
- CMake build with doctest via FetchContent; strict warnings as errors;
  tests run under ASan+UBSan in Debug builds.
- Project documentation: design document, M1 foundations spec, agent
  guidelines (AGENTS.md), and README usage instructions.

### Changed

- Public API follows the documented LLVM-hybrid naming convention
  (camelCase functions, k-prefixed enum values, snake_case variables);
  this baseline applies to all future releases.
