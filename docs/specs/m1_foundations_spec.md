# Spec: M1 Foundations

**Status**: Approved

**Author**: Agent + Vitaly Kovalev

**Date**: 2026-08-24

**Related**: docs/design.md §2 (key decisions), §3.1–3.2 (Element, molecule graph), §7 (error handling), §8 (architecture), §9 (testing strategy), §11 (roadmap, milestone M1)

## 1. Goal / Problem Statement

The repository on `main` is an empty slate — no build system, no sources. Milestone
M1 (docs/design.md §11) delivers the foundation every later milestone parses into:
a CMake skeleton with doctest wired up, the complete periodic table, the
`Element` access layer, the `MolecularGraph` type, and the formula parser.

The proto branch validated these concepts but violates several contracts that
are now binding: its parser silently skips invalid characters instead of raising
errors, hashes molecules by molar mass (double), carries integer pseudo-masses,
and covers only 26 elements. This milestone replaces those concepts with the
redesigned ones under `src/chem/`; nothing is merged from `proto`.

## 2. Scope

### In Scope

- CMake build skeleton: library target `chem`, test-runner target, C++20
- doctest via FetchContent, registered with CTest
- `data/elements.csv`: all 118 elements (symbol, atomic number, name, standard
  atomic weight)
- `chem::Element`: static lookup by symbol and by atomic number
- `chem::MolecularGraph`: full graph API (atoms with charge/H-count fields,
  bonds with order); formula parser populates atoms only, bonds remain empty
  until M2 introduces SMILES
- Composition map (element → count) derived from the graph
- Molar mass derived from composition × standard atomic weights
- Formula parser supporting element symbols, integer counts, and nested
  parenthesized groups with multipliers
- Loud rejection of all unsupported formula syntax as typed exceptions
- Unit tests mirroring the `src/chem/` structure

### Out of Scope

- Anything SMILES: parsing, canonicalization, golden corpus, dialect rules (M2)
- `chem::Molecule`, identity/equality/hash semantics (M2)
- Reactions, balance validation, reaction orders (M3)
- `KineticSystem`, `Vessel`, integrators (M3)
- Serialization / persistence (M4)
- Dot/hydrate notation in formulas (`CuSO4·5H2O`) — rejected loudly, not parsed
- Isotope-labeled formulas (e.g. `D2O`, `[13C]`) — rejected loudly
- Whitespace tolerance inside formulas — strict grammar, rejected loudly

## 3. User Stories or Use Cases

1. As a library developer, I want a working build-and-test loop so that all
   later milestones land incrementally with green CI-style checks.
2. As a downstream contributor implementing the SMILES parser (M2), I want a
   stable `MolecularGraph` type so that my parser emits into the same pipeline
   without touching core types.
3. As a library user, I want to parse `H2SO4` or `Fe(NO3)3` and obtain the
   element composition and molar mass so that I can validate inputs early.

## 4. Functional Requirements

### ⏳ FR-1: Build skeleton

The project MUST provide a CMake configuration (minimum 3.28) defining:

- FR-1a: A library target named `chem`, built as C++20 (`CXX_STANDARD 20`,
  `CXX_STANDARD_REQUIRED ON`), compiling everything under `src/chem/`.
- FR-1b: A test-runner executable linking `chem` and doctest, registering all
  test files with CTest via `add_test` (or `doctest_discover_tests` equivalent).
- FR-1c: Compiler warnings `-Wall -Wextra -Wpedantic -Wconversion
  -Wsign-conversion -Wshadow` enabled for both targets, with warnings treated
  as errors; builds MUST be warning-clean.
- FR-1d: Debug/test builds MUST compile with ASan+UBSan instrumentation;
  release builds MUST carry no sanitizer overhead.

### ⏳ FR-2: doctest wiring

The project MUST fetch doctest via FetchContent (pinned release tag) and
expose it to the test target. No other external dependency MAY be added
(docs/design.md §2).

### ⏳ FR-3: Periodic table data

The file `data/elements.csv` MUST contain exactly 118 data rows with columns
`symbol,atomic_number,name,standard_atomic_weight`, where:

- FR-3a: Weights are IUPAC/CIAAW conventional standard atomic weights in unified
  atomic mass units (u), e.g. H = 1.008, C = 12.011, Cl = 35.45.
- FR-3b: Elements without stable isotopes use the mass number of the
  longest-lived isotope (the bracketed tabulated value), e.g. Tc = 98,
  Pm = 145, Og = 294.
- FR-3c: The CSV MUST be embedded into the build (configured/copied next to the
  test binary or compiled in) so tests locate it deterministically at runtime.

### ⏳ FR-4: Element access

`src/chem/core/element.hpp/cpp` MUST define `chem::Element` as a lightweight
value type constructible from symbol or atomic number:

```cpp
namespace chem {
class Element {
public:
  explicit Element(std::string_view symbol);
  explicit Element(int atomic_number);
  std::string_view symbol() const noexcept;
  int atomic_number() const noexcept;
  double standard_weight() const noexcept;  // u
};
}
```

- FR-4a: Instances MUST be flyweight handles: each object holds only a pointer
  into the shared, once-loaded table; copying never duplicates data and is as
  cheap as copying a pointer. Construction performs one O(1) lookup; hot paths
  hoist constructed instances out of loops. The internal pointer is a
  non-owning view into function-local-static storage (justified bare `T*` per
  AGENTS.md); it MUST never be dereferenced-owned or freed.

- FR-4b: Lookup MUST be case-sensitive: `"Co"` resolves to cobalt, `"CO"` is an
  error.
- FR-4c: Unknown symbols and out-of-range atomic numbers MUST throw
  `chem::ParseError` with a message quoting the offending token.
- FR-4d: Data MUST load once (function-local static initialization); repeated
  lookups MUST not re-read the CSV.

### ⏳ FR-5: MolecularGraph type

`src/chem/core/molecular_graph.hpp/cpp` MUST define the graph that all parsers
emit (design.md §3.2, §8):

- FR-5a: An `Atom` record carrying a value-type `chem::Element element`
  (flyweight handle, see FR-4a), formal charge (`int charge`, default 0),
  implicit hydrogen count (`int implicit_h`, default 0), and explicit hydrogen
  count (`int explicit_h`, default 0).
- FR-5b: A bond representation with order enum
  `enum class BondOrder { Single, Double, Triple }`.
- FR-5c: Mutation API `std::ptrdiff_t add_atom(Atom)` returning the atom index
  and `void add_bond(std::ptrdiff_t a, std::ptrdiff_t b, BondOrder order)`
  (signed index arithmetic per AGENTS.md); read access via contiguous views
  (`std::span<const Atom> atoms()`, `std::span<const Bond> bonds()`).
- FR-5d: `add_bond` MUST debug-assert valid indices; the graph does not
  self-validate in release builds (setup-phase validation belongs to callers).

### ⏳ FR-6: Composition map

The composition (element → count) MUST be derivable from any graph:

- FR-6a: A free function `chem::CompositionMap composition(const MolecularGraph&)`
  where `CompositionMap` maps atomic numbers to total atom counts, including
  implicit and explicit hydrogens of each atom.
- FR-6b: Counts MUST aggregate multiplicities from repeated symbols and
  parenthesized groups correctly (no double counting; see proto bug fixed).

### ⏳ FR-7: Molar mass

A free function `double molar_mass(const MolecularGraph&)` MUST return the sum
of `element.standard_weight() * count` over the composition, in g/mol. For
water the result MUST equal approximately 18.015 g/mol (within 0.001).

### ⏳ FR-8: Formula parser

`src/chem/parsing/formula_parser.hpp/cpp` MUST expose:

```cpp
namespace chem {
// Throws ParseError on any invalid input.
MolecularGraph parse_formula(std::string_view input);
}
```

- FR-8a: Grammar: one or more element tokens; each token is an uppercase letter
  optionally followed by lowercase letters (full official symbols), optionally
  followed by a positive decimal integer count (default 1).
- FR-8b: Groups: parentheses with optional integer multiplier, nestable to
  arbitrary depth (`Ca(OH)2`, `Fe(NO3)3`, `K4[Fe...]` is NOT supported — square
  brackets are rejected; `Mg(NO2(OH))2`-style nesting IS valid).
- FR-8c: Output: a `MolecularGraph` containing one `Atom` per atom instance,
  zero bonds (unconnected multiset per design.md §3.2).
- FR-8d: Parsing MUST be pure: same input always yields the same graph; no
  global state.

### ⏳ FR-9: Loud rejection of invalid formulas

All of the following MUST throw `chem::ParseError` whose message quotes the
offending input substring (design.md §4.3 principle applied to formulas):

- FR-9a: Unknown element symbols (`Xx`, `Zz`, `Hx`).
- FR-9b: Characters outside the grammar at any position: whitespace, digits in
  invalid positions (`2H`), operators, dots (`.`), middle dots (`·`),
  square brackets, charges (`+`, `-`).
- FR-9c: Structural errors: empty input, unbalanced parentheses (either
  direction), a closing parenthesis without a preceding open, multiplier
  attached to nothing (`()2`), zero count (`H0`), count with leading zeros
  (`H007`).
- FR-9d: No silent skipping anywhere — the proto behavior of advancing past
  unrecognized characters MUST NOT reappear.

### ⏳ FR-10: Error types

`src/chem/core/errors.hpp` MUST define `chem::ParseError` and
`chem::ValidationError`, both derived from `std::runtime_error`. Nothing in the
codebase may throw string literals (design.md §7).

### ⏳ FR-11: Directory layout

Sources MUST follow design.md §8:

```
src/chem/
  core/          errors.hpp, element.hpp/cpp, molecular_graph.hpp/cpp
  parsing/       formula_parser.hpp/cpp
tests/
  core/          element_test.cpp, molecular_graph_test.cpp
  parsing/       formula_parser_test.cpp
data/
  elements.csv
```

### ⏳ FR-12: Test coverage

The doctest suite MUST cover:

- FR-12a: Element lookups for a representative sample (lightest, heaviest,
  two-letter symbols, bracket-weight elements like Tc), error paths FR-4c,
  case sensitivity FR-4b, flyweight copy semantics FR-4a.
- FR-12b: Parser acceptance cases: single element (`O`), repeated elements
  (`H2O`), groups and nesting (`Ca(OH)2`, `Fe(NO3)3`), deep nesting
  (≥3 levels), multi-digit counts.
- FR-12c: Every rejection class in FR-9 with exception-type assertions
  (doctest `REQUIRE_THROWS_AS`) and message-content checks quoting the input.
- FR-12d: Graph/composition/mass consistency: for each accepted test formula,
  `molar_mass` recomputed independently matches, and hydrogen counting through
  implicit/explicit fields sums correctly.
- FR-12e: Determinism smoke test: parse the same formula twice and compare
  graphs atom-by-atom.

## 5. Non-Functional Requirements

- **Determinism**: no hidden global state; CSV loads once, identically every run.
- **Correctness over speed**: lookup structures may be simple maps; performance
  budgets apply first in the numerics milestones, not here.
- **Zero dependencies**: only doctest (FetchContent) beyond the toolchain.
- **Testability**: every module reachable from tests via public headers only;
  no test-only backdoors in library code.
- **Portability**: builds with GCC and Clang in C++20 mode.
- **Safety**: RAII everywhere, no manual `new`/`delete`, ownership via value
  members or containers, non-owning access via references/`std::span`/
  `std::string_view`, named casts only (AGENTS.md "Safety and resource
  management"); debug/test runs verified under ASan+UBSan.

## 6. Technical Constraints & Architecture Notes

- Namespace `chem`; naming per AGENTS.md (`PascalCase` types, `snake_case`
  functions/variables, `k`-prefixed constants).
- Parsers produce `MolecularGraph` and nothing else (hard rule 3) — even though
  the formula parser could trivially emit counts directly, it MUST construct a
  graph and derive the composition from it, exercising the M2 path.
- Parsing/canonicalization separation does not bite yet (no canonicalizer);
  the parser MUST NOT reach into `core/` internals beyond public headers.
- Setup phase only at this stage: typed exceptions everywhere; the
  exception-free hot-loop contract has no code path yet.
- The proto branch is reference-only; no code is copied verbatim, and none of
  its bugs (silent skips, mass-based hashing, integer masses) may carry over.

## 7. Data & Interfaces

New public headers and their key signatures:

| Header | Interface |
|---|---|
| `core/errors.hpp` | `class ParseError : std::runtime_error`, `class ValidationError : std::runtime_error` |
| `core/element.hpp` | `Element(std::string_view)`, `Element(int)`, `.symbol()`, `.atomic_number()`, `.standard_weight()` (flyweight value type) |
| `core/molecular_graph.hpp` | `struct Atom`, `enum class BondOrder`, `class MolecularGraph` with `add_atom`, `add_bond`, `atoms()`, `bonds()` |
| `core/composition.hpp` | `CompositionMap composition(const MolecularGraph&)`, `double molar_mass(const MolecularGraph&)` |
| `parsing/formula_parser.hpp` | `MolecularGraph parse_formula(std::string_view)` |

Data file: `data/elements.csv`, header row
`symbol,atomic_number,name,standard_atomic_weight`, 118 rows, UTF-8, LF endings.

## 8. Edge Cases & Error Handling

| Input | Behavior |
|---|---|
| `""` (empty) | `ParseError`, message includes the empty-input context |
| `"H2O "` (trailing space) | `ParseError` at the space |
| `"h2O"` | `ParseError` — lowercase start is not an element token |
| `"Ca(OH)2"` | Valid: 1 Ca, 2 O, 2 H |
| `"Fe(NO3)3"` | Valid: 1 Fe, 3 N, 9 O |
| `"((H))"` | Valid: 1 H (nested group, multiplier defaults to 1) |
| `"H)"` | `ParseError` — close without open |
| `"(H"` | `ParseError` — unclosed group |
| `"H0"` | `ParseError` — zero count |
| `"C6H12O6"` | Valid: glucose composition, molar mass ≈ 180.156 g/mol |
| `"[OH-]"` | `ParseError` — brackets are not formula syntax (SMILES territory, M2) |

Unknown-element and malformed-syntax errors quote the offending substring and
its position where practical.

## 9. Acceptance Criteria / Definition of Done

- [ ] All functional requirements FR-1 … FR-12 satisfied
- [ ] `cmake -B build -S .` configures cleanly; FetchContent pins doctest
- [ ] `cmake --build build --parallel` compiles warning-clean under
      `-Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wshadow` with
      warnings as errors
- [ ] `ctest --test-dir build --output-on-failure` passes all suites, with the
      test build running under ASan+UBSan (FR-1d)
- [ ] Touched sources formatted via `clang-format -i`; `clang-tidy` reports no
      findings on changed files
- [ ] Tests exist for: element lookups + error paths, parser accept/reject
      classes, composition aggregation, molar mass values, determinism
- [ ] `data/elements.csv` contains exactly 118 rows; spot-checked weights match
      CIAAW tables (H, C, Cl, Tc, Og at minimum)
- [ ] No new dependencies beyond doctest in CMakeLists.txt
- [ ] No changes to canonicalization (nothing frozen exists yet — trivially met)
- [ ] Layout matches FR-11; AGENTS.md conventions respected throughout

## 10. Open Questions / Decisions Needed

None. Resolved during drafting:

- Graph scope: full API now, bonds unused until M2 (user decision).
- Mass convention: IUPAC standard atomic weights, bracketed values for
  monoisotopic-unstable elements (user decision).
- Composition map and molar mass included in M1, ahead of `Molecule` in M2
  (user decision).

## 11. Implementation Notes

- Embedding the CSV: prefer `configure_file(... COPYONLY)` into the build dir
  (proto precedent) plus a compile definition pointing tests at the copy, or
  generate a header at configure time — either is acceptable; loading must stay
  deterministic across runs (FR-4d).
- Parser shape: recursive descent mirroring the proto's group-stack approach is
  sufficient; replace silent-skip branches with error throws.
- Composition keyed by atomic number keeps it hashable and independent of
  string lifetime issues.
- Molar-mass test tolerances: 1e-3 g/mol absolute is loose enough for CSV
  rounding, tight enough to catch swapped weights.
