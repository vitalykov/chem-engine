# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to Semantic Versioning. During the 0.x series,
minor versions may still introduce breaking API changes; promotion to 1.0.0
requires the criteria in AGENTS.md ("Versioning").

## [Unreleased]

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
