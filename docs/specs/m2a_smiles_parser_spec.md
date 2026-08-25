# Spec: M2a — SMILES Parser

**Status**: Review

**Author**: Agent + Vitaly Kovalev

**Date**: 2026-08-25

**Related**: docs/specs/m2_identity_spec.md (umbrella); docs/design.md §3.2, §4.1, §4.3, §7; docs/specs/m1_foundations_spec.md FR-5 (`MolecularGraph`)

## 1. Goal / Problem Statement

The engine's only input format so far is Hill-formula notation, which produces
unconnected atom multisets and therefore cannot distinguish structural
isomers. This sub-spec adds the SMILES parser: the second producer of
`MolecularGraph`, covering the v1 dialect subset of docs/design.md §4.3. It is
the connectivity-bearing input that makes graph-based identity meaningful.

The dialect is a strict subset of real-world SMILES. Anything outside the
subset — aromatic atoms, stereo markers, disconnected components — is rejected
loudly with `chem::ParseError` (hard architectural rule 4): silent acceptance
could collapse distinct molecules into one ID.

## 2. Scope

### In Scope

- `chem::parseSmiles(std::string_view)` in `src/chem/parsing/`
- Tokenizer + recursive-descent parser emitting into the existing
  `chem::MolecularGraph` (no new graph type)
- Bare organic-subset atoms and bracketed atoms (all 118 elements)
- Explicit (`-`, `=`, `#`) and implicit single bonds, branches, chains
- Ring closures with single digits and `%nn`
- Formal charges and isotope labels in brackets
- Implicit-hydrogen assignment via the DayLight-style valence model
- Loud rejection of all unsupported syntax as typed exceptions
- Extension of `Atom` with an isotope field
- Unit tests mirroring `src/chem/parsing/`

### Out of Scope

- Canonicalization (m2b) — the parser MUST NOT reorder or canonicalize anything
- Aromaticity perception; lowercase aromatic atoms (`c`, `o`, ...) are errors
- Stereochemistry: `@`, `@@`, `/`, `\` are errors
- Dot-disconnected SMILES (`.`) is an error in v1
- Reaction SMILES (`>`), wildcards (`*`, `%` as atom), atom maps (`:n`),
  atom-class labels (`:c`) are errors
- Bond orders beyond single/double/triple (quadruple, dative) are errors

## 3. User Stories or Use Cases

1. As a library user, I want to write `CC(=O)O` and get a bonded molecular
   graph so that acetic acid is distinguishable from its formula `C2H4O2`.
2. As a library user, I want `[13CH4]` or `[NH4+]` parsed faithfully so that
   isotopologues and charged species carry their labels through identity.
3. As a library developer implementing the canonicalizer (m2b), I want a
   parser that emits plain graphs so that canonicalization consumes the same
   `MolecularGraph` the formula parser already produces.

## 4. Functional Requirements

### ⏳ FR-1: Entry point

`src/chem/parsing/smiles_parser.hpp/cpp` MUST expose:

```cpp
namespace chem {
// Parses one connected SMILES molecule into a bonded graph.
// Throws ParseError on any invalid or unsupported input.
[[nodiscard]] MolecularGraph parseSmiles(std::string_view input);
}
```

- FR-1a: The result MUST contain exactly the atoms and bonds described by the
  input, with implicit hydrogens recorded per FR-7 and nothing else.
- FR-1b: Parsing MUST be pure: same input always yields an equal graph; no
  global state, no locale dependence.

### ⏳ FR-2: Bare organic atoms

The bare organic subset `B C N O P S F Cl Br I` MUST parse as unbracketed
atoms.

- FR-2a: Symbols MUST be matched case-sensitively, trying the two-character
  forms (`Cl`, `Br`) before their single-letter prefixes; `"Cl"` is chlorine.
  A lowercase letter that does not complete `Cl` or `Br` is an error (FR-10).
- FR-2b: Any other bare element symbol (e.g. `[Fe]` written as `Fe` bare) is
  an error; non-organic elements REQUIRE brackets.

### ⏳ FR-3: Bracket atoms

Bracket atoms `[...]` MUST support any element symbol from the periodic table,
with optional leading isotope mass number, optional explicit hydrogen count,
and optional formal charge:

```
'[' [isotope] symbol ['H' [count]] [charge] ']'
```

Examples: `[C]`, `[OH-]`, `[NH4+]`, `[13C]`, `[2H]`, `[Fe+3]`, `[O-2]`.

- FR-3a: Bracket atoms MUST set `explicit_h` to the given H count (default 0)
  and `implicit_h` to 0 — no valence inference inside brackets.
- FR-3b: Charges MUST accept `+`/`-` repeated up to a magnitude of 9
  (`[Fe+3]` and equivalent `[Fe++]`-style single-symbol runs limited to `++`,
  `--`; digit form required beyond that). Charge magnitude MUST fit `int`.
- FR-3c: Isotopes MUST be positive integers stored on the atom (FR-12); a
  value above 999 is an error.

### ⏳ FR-4: Bonds

Bond tokens `-`, `=`, `#` MUST map to `BondOrder::kSingle`, `kDouble`,
`kTriple`.

- FR-4a: A bond token between two atoms (or preceding a ring-closure digit)
  MUST add exactly one bond of the given order.
- FR-4b: Where no bond token appears between adjacent atoms, a single bond
  MUST be created (default bond).
- FR-4c: A bond token at any position where no following atom or ring closure
  exists is an error (e.g. trailing `C-`, `C=)`).

### ⏳ FR-5: Branches and chains

Parentheses MUST group branches: `( ... )` attaches to the current atom and
pops back to it on close.

- FR-5a: Branches nest to arbitrary depth.
- FR-5b: An empty branch `()` , an unclosed branch `(CC`, and a closing
  parenthesis without an open one `CC)` are errors.
- FR-5c: An atom MAY have multiple branches plus continuation bonds
  (`C(C)(C)C`).

### ⏳ FR-6: Ring closures

Ring-bond digits `0`–`9` and `%nn` (two decimal digits) MUST open and close
rings.

- FR-6a: A digit on an atom marks a ring-opening; the same digit on a later
  atom closes the ring by adding one bond between the two atoms.
- FR-6b: An optional bond token before the ring digit sets the ring bond order;
  if both ends specify orders they MUST match, otherwise error; if neither
  specifies, the bond is single.
- FR-6c: Reusing a currently-open digit, closing a digit that was never opened,
  re-closing an already-closed digit, and malformed `%` (non-two-digit
  follow-up, e.g. `%C` or `%7`) are errors.
- FR-6d: Ring closures MUST NOT create self-loops (same atom twice via one
  digit pair is impossible by construction, but a digit reused on the same
  atom immediately is an error).

### ⏳ FR-7: Implicit hydrogens (DayLight-style valence model)

For bare organic atoms, implicit hydrogens MUST be derived from the bond-order
sum using per-element allowed valences:

| Element | Allowed total valences |
|---|---|
| B | {3} |
| C | {4} |
| N | {3, 5} |
| O | {2} |
| P | {3, 5} |
| S | {2, 4, 6} |
| F, Cl, Br, I | {1} |

- FR-7a: Let `b` be the sum of bond orders over the atom's bonds (including
  ring-closure bonds). The element's smallest allowed valence `v ≥ b` is
  selected; `implicit_h = v − b`. If no allowed valence satisfies `v ≥ b`,
  the input is over-valent and MUST raise `ParseError`.
- FR-7b: Implicit hydrogens are computed after the whole molecule parses
  (ring closures may arrive late); they MUST NOT be assigned incrementally
  during traversal.
- FR-7c: Bracket atoms never receive implicit hydrogens (FR-3a).
- FR-7d: Examples: water `O` → O with implicit_h = 2; ammonia `N` → 3;
  methane `C` → 4; `[NH4+]` → explicit_h = 4, implicit_h = 0; `CS(=O)(=O)O`
  → S selects v = 6, implicit_h = 0.

### ⏳ FR-8: Charges

Formal charges exist only inside brackets (FR-3). Bare atoms always have
charge 0; `C+` outside brackets is an error.

### ⏳ FR-9: Isotope storage extension

`src/chem/core/molecular_graph.hpp` `struct Atom` MUST gain an isotope field:

```cpp
struct Atom {
  Element element;
  int charge = 0;
  int implicit_h = 0;
  int explicit_h = 0;
  int isotope = 0;   // mass number; 0 = unspecified
};
```

- FR-9a: Default 0 keeps every existing graph (formula parser output, M1
  tests) unchanged.
- FR-9b: Composition and molar mass (M1) remain keyed by atomic number and
  use standard weights regardless of isotope values; isotope-aware mass is out
  of scope (documented approximation).

### ⏳ FR-10: Loud rejection

All of the following MUST throw `ParseError` quoting the offending substring
and its offset (FR-11):

- FR-10a: Lowercase aromatic atoms anywhere (`c1ccccc1`, `co`).
- FR-10b: Stereo markers `@`, `@@`, directional bonds `/`, `\`.
- FR-10c: Disconnected components: `.` at top level or anywhere.
- FR-10d: Reaction SMILES: `>`.
- FR-10e: Wildcards `*`, atom maps (`:`), atom-class labels.
- FR-10f: Unknown symbols, invalid bracket contents (`[Hx]`, `[]`, `[123]`),
  unterminated bracket `[OH-`.
- FR-10g: Structural errors: empty input, whitespace anywhere, stray
  characters not covered by the grammar, unclosed ring digit at end of input
  (digit opened but never closed).
- FR-10h: No silent skipping anywhere; any unrecognized byte is an error.

### ⏳ FR-11: Error quality

Every `ParseError` message MUST include the zero-based byte offset of the
offending token and quote enough surrounding input to locate it
(e.g. `parse error at offset 4: '...' in "C1=CC=CC=C"`).

### ⏳ FR-12: Test coverage

`tests/parsing/smiles_parser_test.cpp` MUST cover:

- FR-12a: Acceptance: `C`, `O`, `[OH2]`, `CC(=O)O`, `C1=CC=CC=C1` (benzene
  kekulé), `C1CCC2CCCCC2C1`-style fused rings, `%10`-style multi-digit ring
  closures, `[13CH4]`, `[NH4+]`, `[Fe+3]`, `CS(=O)(=O)O`.
- FR-12b: Graph correctness for accepted inputs: atom counts, element
  identities, bond orders, implicit/explicit H totals, charges, isotopes.
- FR-12c: Valence model: each allowed-valence row exercised, including
  smallest-valence selection when several valences fit (e.g. S at bond sum 2
  picks v = 2, not 4 or 6), plus over-valence rejection (`CO(C)C`).
- FR-12d: Every rejection class FR-10a–g via `REQUIRE_THROWS_AS` plus
  message-content checks (offset + quoted input present).
- FR-12e: Purity smoke test: parse the same string twice, compare graphs
  field-by-field.
- FR-12f: Regression guard: existing formula-parser graphs still validate
  (isotope default 0 does not disturb M1 behavior).

## 5. Non-Functional Requirements

- **Determinism**: pure function of input bytes; no global state.
- **Fail-fast**: all validation happens during parsing (setup phase);
  no deferred validation state.
- **Zero dependencies**: doctest only.
- **Safety**: RAII/value semantics throughout the parser; indices as
  `std::ptrdiff_t`; no manual memory management; named casts only.
- **Portability**: GCC/Clang C++20.

## 6. Technical Constraints & Architecture Notes

- Namespace `chem`; naming per AGENTS.md; parser internals live in an
  anonymous namespace in the .cpp or `chem::detail`.
- The parser MUST emit into the public `MolecularGraph` API (`addAtom`,
  `addBond`) and MUST NOT reach into its private members or core internals.
- Implicit-hydrogen computation is part of parsing (the graph is complete when
  `parseSmiles` returns); it must not leak into canonicalization concerns.
- The tokenizer/parser shape (recursive descent mirroring branch depth) is
  implementer freedom; the observable contract is FR-1 … FR-11.

**References** — normative sources for the dialect:

- OpenSMILES specification (<https://opensmiles.org>) — grammar baseline. Our
  dialect is a strict subset; constructs outside FR-2 … FR-8 are rejected per
  FR-10 regardless of what either reference permits.
- Daylight SMILES theory manual
  (<https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html>) —
  background for the valence/implicit-hydrogen model (FR-7).

These references govern *parsing* conformance only. Canonical output ordering
is defined solely by this project (m2b); cross-toolkit canonical-form
compatibility is out of scope by design (docs/design.md §4).

## 7. Data & Interfaces

| Header | Interface |
|---|---|
| `parsing/smiles_parser.hpp` | `MolecularGraph parseSmiles(std::string_view)` |
| `core/molecular_graph.hpp` | `Atom` gains `int isotope = 0` |

No new data files. CMakeLists.txt gains `smiles_parser.cpp` in the `chem`
target sources, the header in `FILE_SET HEADERS`, and
`tests/parsing/smiles_parser_test.cpp` in `test-runner`.

Directory delta:

```
src/chem/parsing/smiles_parser.hpp/cpp   (new)
tests/parsing/smiles_parser_test.cpp     (new)
src/chem/core/molecular_graph.hpp        (extended)
```

## 8. Edge Cases & Error Handling

| Input | Behavior |
|---|---|
| `""` (empty) | `ParseError`, offset 0 |
| `"O"` | Valid: O with implicit_h = 2 (water) |
| `"[OH2]"` | Valid: O, explicit_h = 2, implicit_h = 0 — same composition as `O` |
| `"N"` | Valid: implicit_h = 3 |
| `"[NH4+]"` | Valid: charge +1, explicit_h = 4 |
| `"[Fe++]"` | Valid: charge +2 |
| `"[Fe+3]"` | Valid: charge +3 |
| `"[Fe+++]"` | `ParseError` — run form limited to `++`/`--`; use digit form |
| `"c1ccccc1"` | `ParseError` at offset 0 — aromatic lowercase atom |
| `"C/C=C/C"` | `ParseError` at first `/` — stereo |
| `"[Na+].[Cl-]"` | `ParseError` at `.` — disconnected components unsupported |
| `"C1CC"` | `ParseError` — ring digit 1 opened, never closed |
| `"C11"` | `ParseError` — digit reuse while open / immediate reclose |
| `"C%1CC%1"` | `ParseError` — `%` must be followed by two digits |
| `"C(C"` | `ParseError` — unclosed branch |
| `"C=O"` | Valid: double bond |
| `"C(F)(F)(F)F"` | Valid: tetrafluoromethane; F implicit_h = 0 each |
| `"CO(C)C"` | `ParseError` — O over-valent (bond sum 3 > max 2) |
| `"CS(=O)(=O)O"` | Valid: S bond sum 6, implicit_h = 0 |

## 9. Acceptance Criteria / Definition of Done

- [ ] All functional requirements FR-1 … FR-12 satisfied
- [ ] Warning-clean build; full `ctest` green under sanitizers
- [ ] Tests exist for: grammar acceptance/rejection classes, valence table,
      ring closures including `%nn`, isotopes/charges, purity, M1 regression
- [ ] Error messages verified to contain offset + offending substring
- [ ] `clang-format -i` applied; `clang-tidy` clean on changed files
- [ ] No new dependencies; layout matches §7

## 10. Open Questions / Decisions Needed

None. Resolved during drafting:

- Valence model: DayLight-style allowed-valence sets (user decision).
- Dots rejected in v1; disconnected graphs arise only via m2c's internal
  fallback path (user decision, umbrella scope).

## 11. Implementation Notes

- Two-pass shape works well: pass 1 tokenizes/traverses building atoms, raw
  bonds, and ring-bond bookkeeping; pass 2 resolves ring-closure pairs, then
  computes implicit hydrogens (FR-7b) in one sweep.
- Ring-digit state is a small fixed array indexed 0–99 (`%nn` range); "open"
  records the atom index and pending bond order.
- The valence lookup can be a `constexpr` table keyed by atomic number; bare
  atoms are restricted to the ten-element organic subset so the table is tiny.
- Keep offsets throughout: tracking `std::size_t pos` in the cursor makes
  FR-11 messages nearly free.
