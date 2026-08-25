# Spec: M2c — Molecule Type and Formula Resolution

**Status**: Review

**Author**: Agent + Vitaly Kovalev

**Date**: 2026-08-25

**Related**: docs/specs/m2_identity_spec.md (umbrella); docs/specs/m2a_smiles_parser_spec.md; docs/specs/m2b_canonicalization_spec.md; docs/design.md §3.2 (Molecule), §4.4

## 1. Goal / Problem Statement

M2a produces graphs and m2b freezes their canonical names. This sub-spec adds
the user-facing species type: `chem::Molecule`, which owns a graph, computes
its identity once at construction, and exposes composition and molar mass.
Equality and hashing are defined on canonical SMILES — never on composition,
mass, or raw input (hard architectural rule 2).

It also resolves the one structural gap left by M1: the formula parser emits
unconnected atom multisets, but identity needs connectivity. Formula inputs go
through a resolution pipeline — unique valence-based bond inference first, a
curated dictionary of unambiguous time-proven formulas second, and a loud
`ValidationError` with same-composition suggestions otherwise.

Construction is explicit per format (`fromFormula` / `fromSmiles`): no
auto-detection, so ambiguous strings like `"C"` or `"NO"` can never be silently
misrouted.

## 2. Scope

### In Scope

- `chem::Molecule` value type in `src/chem/core/molecule.hpp/cpp`
- Static factories `Molecule::fromFormula`, `Molecule::fromSmiles`,
  `Molecule::fromGraph`
- Cached canonical SMILES, composition map, molar mass
- `operator==` and `std::hash<Molecule>` on canonical form
- Formula-to-connectivity resolution: unique-inference → dictionary → error
  with alternatives
- `data/compounds.csv` (~150 curated entries: unambiguous formulas only)
- Unit tests mirroring `src/chem/core/`

### Out of Scope

- Parser or canonicalizer internals (m2a/m2b) — this type composes them
- Default-isomer guessing for ambiguous formulas (explicitly rejected:
  dictionary carries unambiguous writings only)
- Isotope-aware mass; charge-separated salts requiring disconnected SMILES
  input (dictionary entries must parse via `parseSmiles`, which rejects dots)
- Consteval compile-time literal validation (deferred, §10)
- Reactions referencing molecules (M3)

## 3. User Stories or Use Cases

1. As a library user, I want `Molecule::fromSmiles("CCO") ==
   Molecule::fromFormula("C2H5OH")` so that either spelling identifies ethanol.
2. As a library user, I want `Molecule::fromFormula("C2H6O")` to fail with a
   suggestion (`C2H5OH`) instead of silently picking an isomer.
3. As a library user building containers of species, I want `std::unordered_map<Molecule, ...>`
   to work correctly so that duplicates collapse by identity.

## 4. Functional Requirements

### ⏳ FR-1: Type shape

`src/chem/core/molecule.hpp/cpp` MUST define:

```cpp
namespace chem {
class Molecule {
public:
  static Molecule fromFormula(std::string_view input);
  static Molecule fromSmiles(std::string_view input);
  static Molecule fromGraph(const MolecularGraph& graph);

  [[nodiscard]] const MolecularGraph& graph() const noexcept;
  [[nodiscard]] const std::string& canonicalSmiles() const noexcept;
  [[nodiscard]] const CompositionMap& composition() const noexcept;
  [[nodiscard]] double molarMass() const noexcept;

private:
  explicit Molecule(MolecularGraph graph);

  MolecularGraph graph_;
  std::string canonical_smiles_;
  CompositionMap composition_;
  double molar_mass_ = 0.0;
};

bool operator==(const Molecule& a, const Molecule& b) noexcept;
}

template <>
struct std::hash<chem::Molecule> { ... };
```

- FR-1a: The private constructor MUST fully establish validity: it computes
  canonical SMILES, composition, and molar mass exactly once (design.md §3.2).
  No two-phase init; factories are the only construction path.
- FR-1b: All cached members MUST be computed in the constructor; accessors
  are trivial and `noexcept`.

### ⏳ FR-2: fromSmiles

`Molecule::fromFormula` aside, `fromSmiles(input)` MUST behave as
`parseSmiles(input)` followed by the shared constructor path.

- FR-2a: Parse failures propagate as the `ParseError`s of m2a unchanged.
- FR-2b: Canonicalization failure (`ValidationError` from m2b) propagates unchanged.

### ⏳ FR-3: fromGraph

`Molecule::fromGraph(graph)` MUST accept any programmatically built graph,
including disconnected ones.

- FR-3a: It MUST validate valences before accepting: every bare organic atom
  (m2a's ten-element subset) with bonds MUST satisfy the m2a valence model;
  violations throw `ValidationError`. Atoms outside the bare subset MUST have
  been bracket-style by construction (charge/isotope/H fields carried as set);
  no additional restriction applies beyond non-negativity of H counts.
- FR-3b: This is the sanctioned entry point for future input formats (e.g.
  drawer graphs): per design.md §8 they MUST produce a `MolecularGraph`, and
  this factory is where such graphs become molecules.

### ⏳ FR-4: Formula resolution — unique inference

`Molecule::fromFormula(input)` MUST first run `parseFormula`, then attempt
connectivity inference over the atom multiset:

- FR-4a: Candidate bond assignments are enumerated such that every atom's
  bond-order sum equals one of its allowed valences (the m2a table). Only
  assignments whose total bond count is minimal for satisfying all valences
  are candidates.
- FR-4b: If exactly one candidate exists up to graph isomorphism, its bonds
  are added to the graph and construction proceeds normally. Single-element
  multisets and noble-gas-like atoms with no allowed bonding resolve to
  disconnected singletons (canonicalized via m2b FR-8).
- FR-4c: If zero or multiple candidates exist, resolution falls through to
  FR-5. No guess is ever emitted.
- FR-4d: Examples: `H2O` → unique (H–O–H); `CH4` → unique; `NH3` → unique;
  `O2` → unique (double bond); `H2O2` → unique (H–O–O–H).

### ⏳ FR-5: Formula resolution — compound dictionary

When inference does not resolve (FR-4c), the pipeline consults
`data/compounds.csv`:

```
formula,name,smiles
```

- FR-5a: Entries are keyed by their exact conventional formula spelling. The
  lookup tries the user's input verbatim and accepts a hit only if the
  entry's composition equals the parsed input's composition (both sides
  computed via `parseFormula`, see FR-11 notes).
- FR-5b: A hit yields a `MolecularGraph` via `parseSmiles(smiles)`; the
  resolved graph's composition MUST equal the parsed formula's composition,
  else the data file is invalid — load-time validation asserts this for every
  row (fail fast at startup, not at lookup).
- FR-5c: Curation rule (user decision): entries carry **unambiguous,
  time-proven formulas only**. Compounds conventionally written ambiguously
  appear under their disambiguating spelling (`C2H5OH` for ethanol), never
  under ambiguous ones (`C2H6O`). Initial content ~150 common compounds
  (water, ammonia, major acids/bases/salts-as-covalent, simple organics,
  fuels, sugars); each entry hand-checked once at creation.
- FR-5d: The CSV loads once (function-local static, mirroring elements.csv),
  embedded into the build via `configure_file(... COPYONLY)` + compile
  definition like the other data files.

### ⏳ FR-6: Ambiguity rejection with alternatives

If inference fails and the dictionary has no entry for the input spelling:

- FR-6a: `fromFormula` MUST throw `chem::ValidationError`.
- FR-6b: The message MUST quote the input and list dictionary compounds whose
  composition equals the input's, formatted as suggested spellings with names
  (e.g. `"C2H6O" is ambiguous; known writings with this composition: C2H5OH
  (ethanol). Use SMILES to specify structure.`). If none exist, the message
  says so and advises SMILES.
- FR-6c: The error MUST NOT name a structure for the ambiguous input — no
  silent default isomer.

### ⏳ FR-7: Equality and hashing

- FR-7a: `operator==` MUST compare `canonicalSmiles()` strings only.
- FR-7b: `std::hash<Molecule>` MUST hash the same string; equal molecules
  hash equally; hash quality follows `std::hash<std::string>`.
- FR-7c: Neither equality nor hashing MAY consult composition, molar mass,
  source format, or raw input strings (hard architectural rule 2).
- FR-7d: Distinct structural isomers MUST compare unequal even when
  compositions and masses are identical (`CCO` vs `COC`).

### ⏳ FR-8: Directory layout

```
src/chem/core/molecule.hpp/cpp        (new)
tests/core/molecule_test.cpp          (new)
tests/core/formula_resolution_test.cpp (new)
data/compounds.csv                    (new)
```

CMakeLists.txt gains the sources, public header, test files, and the
compounds.csv copy step mirroring elements.csv.

### ⏳ FR-9: Test coverage

- FR-9a: Identity: cross-format equivalence pairs (`[OH2]`/`O`; `CCO` /
  `C2H5OH`; benzene kekulé / `C6H6`), isomer distinctness (`CCO` != `COC`),
  hash consistency with equality, cross-format **hash equality** (molecules
  equal via different input formats produce identical hashes), and use inside
  `std::unordered_set`.
- FR-9b: Inference: unique-resolution examples of FR-4d produce bonded graphs
  with correct bond counts and pass round-trip; fall-through cases reach FR-6.
- FR-9c: Dictionary: known entry resolves (`CH3COOH` if curated, else another
  listed writing); miss with same-composition suggestion raises
  `ValidationError` containing the alternative spelling; data-file integrity
  (composition agreement) verified by a dedicated suite iterating all rows.
- FR-9d: Error paths: empty input, parser errors propagated typed, ambiguity
  message content checks.
- FR-9e: Accessors: composition and mass match independently recomputed
  values from the stored graph.

## 5. Non-Functional Requirements

- **Determinism**: identical inputs yield identical molecules; dictionary
  loads identically every run.
- **Correctness over speed**: inference enumeration may be exponential in
  principle; acceptable for v1-scale formulas — a size guard (formulas above
  ~30 heavy atoms skip inference and go straight to dictionary/error, FR-4
  note) keeps worst cases bounded.
- **Zero dependencies**, **Safety**, **Portability**: per umbrella spec §5.
- **Fail-fast**: compound dictionary validated fully at first load.

## 6. Technical Constraints & Architecture Notes

- Namespace `chem`; naming per AGENTS.md (`snake_case` members with `_`
  suffix, factories camelCase).
- `Molecule` composes only public APIs of parsing and canonical modules —
  it owns no chemistry logic of its own beyond orchestration and caching.
- The inference engine lives beside the molecule code (`core/` or a
  `parsing/` helper); it MUST reuse the m2a valence table rather than define
  a parallel one — export that table from `smiles_parser.hpp` (or a small
  shared header) as part of this sub-spec if needed.
- Documented contract caveat (umbrella FR-U4): a formula input resolves to an
  identity only when inference or dictionary succeeds; until then there is no
  ID for that input. Normalization guarantee #1 therefore holds within each
  format and across successfully resolved formulas. This caveat goes into the
  README usage section with this milestone.

## 7. Data & Interfaces

| Header | Interface |
|---|---|
| `core/molecule.hpp` | `class Molecule` (factories, accessors), `operator==`, `std::hash` specialization |
| `parsing/smiles_parser.hpp` | (extended) exported valence table for inference reuse |

Data file: `data/compounds.csv`, columns `formula,name,smiles`, ~150 rows,
UTF-8, LF endings; every row's `smiles` parses clean under m2a rules and its
composition matches its `formula` column (load-time asserted, FR-5b).

## 8. Edge Cases & Error Handling

| Input | Behavior |
|---|---|
| `fromSmiles("O")` vs `fromFormula("H2O")` | Equal molecules (unique inference resolves water) |
| `fromFormula("H2O")` graph | Bonded H–O–H, not the disconnected fallback |
| `fromFormula("C2H6O")` | `ValidationError` suggesting `C2H5OH (ethanol)` |
| `fromFormula("C2H5OH")` | Resolves via dictionary to ethanol's graph |
| `fromFormula("Xx2O")` | `ParseError` (unknown element) — parser stage, before resolution |
| `fromFormula("NaCl")` | Inference finds no covalent assignment; dictionary miss (ionic, not representable) → `ValidationError`; documented limitation |
| `fromGraph(g)` with O having bond sum 3 | `ValidationError` (valence check, FR-3a) |
| `fromFormula("C30H60")`-scale input | Size guard: skips inference, dictionary/error path only (§5 note) |
| `Molecule` copied/moved | Trivially copyable semantics preserved; caches travel with it |

## 9. Acceptance Criteria / Definition of Done

- [ ] All functional requirements FR-1 … FR-9 satisfied
- [ ] Cross-format equivalence and isomer-distinctness tests green; hash/equality consistent
- [ ] Compound dictionary integrity suite passes over all ~150 rows
- [ ] Ambiguity errors verified to contain suggestions where available
- [ ] Warning-clean build; full `ctest` green under sanitizers
- [ ] `clang-format -i` applied; `clang-tidy` clean on changed files
- [ ] README usage section updated (public header added; formula caveat documented)
- [ ] No new dependencies; layout matches FR-8

## 10. Open Questions / Decisions Needed

None blocking. Deferred (recorded for future milestones):

- **Consteval molecule literal**: a `consteval operator""_mol` validating
  formulas/SMILES at compile time would move error discovery earlier but adds
  notable implementation cost (constexpr parsers). Deferred; revisit after
  the runtime pipeline stabilizes (user decision: defer).
- **Compile-time checking strategy** generally (static_assert-friendly
  factory wrappers) — same rationale.

## 11. Implementation Notes

- Unique-inference enumeration: atoms with allowed-valence sets bound the
  total bond count tightly (`Σ valence_i = 2·bonds`); enumerate degree
  sequences first, then pair up stubs, rejecting any partial assignment that
  exceeds an element's allowed valence. For v1-scale inputs (< 30 heavy atoms,
  guard above) exhaustive search with early pruning is fine.
- Graph isomorphism check for "exactly one candidate": compare canonical
  forms (m2b) of candidates — equal canonical strings mean the same molecule,
  collapsing isomorphic duplicates cheaply since canonicalization already exists.
- Dictionary normalization detail (FR-5a): store entries keyed by the exact
  conventional spelling; validate at load that `composition(parseFormula(row.formula))
  == composition(parseSmiles(row.smiles))`, and at lookup verify the same
  equality against the user input before accepting a hit.
- Keep `Molecule` copyable (value semantics, RAII); caches are plain values,
  so copies stay cheap enough for v1 and avoid lifetime questions entirely.
