# Spec: M2b — Canonicalization and Golden Corpus

**Status**: Review

**Author**: Agent + Vitaly Kovalev

**Date**: 2026-08-25

**Related**: docs/specs/m2_identity_spec.md (umbrella); docs/design.md §4 (identity contract), §8, §9; docs/specs/m1_foundations_spec.md FR-11 (layout)

## 1. Goal / Problem Statement

Identity is canonical SMILES. This sub-spec defines the canonicalizer: a pure
function `(MolecularGraph) -> std::string` whose output becomes the engine's
species identifier.

Because there is no universal canonical SMILES, ours is defined by our
implementation and frozen: the algorithm is pinned below as canonicalization
spec version 1, its output is locked by golden-file tests on every build, and
round-trip fidelity is asserted across the whole corpus (docs/design.md §4.2).

## 2. Scope

### In Scope

- `chem::canonicalSmiles(const MolecularGraph&)` in `src/chem/canonical/`
- Canonicalization spec version constant (`kCanonicalSpecVersion = 1`)
- Pinned algorithm: initial invariants, Weisfeiler-Lehman refinement,
  deterministic traversal, ring-closure numbering rules
- Handling of disconnected graphs (formula-derived fallback inputs)
- Round-trip invariant tests over the whole corpus
- Golden corpus file `data/golden/corpus.csv` (~30 entries) with build wiring
  mirroring `data/elements.csv`
- Unit tests mirroring `src/chem/canonical/`

### Out of Scope

- Parsing anything (m2a); the canonicalizer consumes graphs only
- Stereochemistry descriptors in output (inputs with stereo are rejected upstream)
- Aromaticity perception — output is kekulé only
- Isotope-aware molar mass or isotope-differentiated composition
- Canonical reaction SMILES (M3+)
- Any change to output after M2 lands without a formal spec version bump
  (docs/design.md §4.2; hard architectural rule 1)

## 3. User Stories or Use Cases

1. As a library user, I want `canonicalSmiles` of two equivalent graphs to be
   byte-identical so that I can compare molecules by string.
2. As a library developer building the reaction layer (M3), I want molecules
   parsed from different spellings of the same species to hash identically so
   that species collections and reaction participant maps stay duplicate-free.
3. As a library developer changing ranking internals, I want a failing golden
   test to tell me immediately that my change alters identity.

## 4. Functional Requirements

### ⏳ FR-1: Entry point and version constant

`src/chem/canonical/canonical_smiles.hpp/cpp` MUST expose:

```cpp
namespace chem {
inline constexpr std::uint32_t kCanonicalSpecVersion = 1;

// Pure function: graph -> canonical SMILES (spec version kCanonicalSpecVersion).
[[nodiscard]] std::string canonicalSmiles(const MolecularGraph& graph);
}
```

- FR-1a: The function MUST be deterministic and side-effect free: equal graphs
  yield identical strings regardless of atom/bond insertion order.
- FR-1b: An empty graph MUST throw `ValidationError` (nothing to identify).

### ⏳ FR-2: Initial invariants

Each atom MUST receive an initial invariant tuple:

```
(isotope, atomic_number, charge, total_h, degree)
```

where `total_h = implicit_h + explicit_h` and `degree` counts bonded
neighbors. Tuples compare lexicographically.

### ⏳ FR-3: Refinement

Invariants MUST be refined Weisfeiler-Lehman style:

- FR-3a: Each iteration replaces every atom's invariant with
  `(previous_invariant, sorted_multiset_of (neighbor_invariant, bond_order))`.
- FR-3b: Iteration repeats until the induced partition of atoms stops
  refining, with a hard cap of 64 iterations.
- FR-3c: Invariants are compared as tuples; no hashing collisions are
  possible because tuples are kept structural (no hash folding into integers).

### ⏳ FR-4: Canonical ranking and tie-breaking

- FR-4a: Atoms are ranked by final refined invariant, lexicographic order.
- FR-4b: Residual ties between non-equivalent atoms MUST NOT occur for the
  corpus (refinement separates them); if an implementation defect produces
  one, tie-break by atom index to stay deterministic, and this case MUST be
  surfaced in debug builds via assert.
- FR-4c: Residual ties between automorphic (symmetry-equivalent) atoms MAY be
  broken arbitrarily: any choice provably yields the same output string, which
  is why determinism survives.

### ⏳ FR-5: Traversal and writing

The writer MUST perform a depth-first emission:

- FR-5a: Start at the lowest-ranked atom not yet visited.
- FR-5b: Visit neighbors in increasing canonical rank order.
- FR-5c: Emit `(` before descending into all-but-the-first neighbor subtree
  at any atom, `)` after it returns; first neighbor continues inline.
- FR-5d: When traversal reaches a visited atom or re-crosses a ring-closure
  pair, emit a ring closure digit instead of recursing; digits are assigned
  sequentially in emission order: `1`–`9`, then `%10`–`%99`; exceeding 99
  simultaneous closures is a `ValidationError` (not reachable for v1-scale
  molecules).
- FR-5e: After finishing a component, continue with the next unvisited start
  atom until all atoms are emitted.

### ⏳ FR-6: Atom emission

- FR-6a: An atom MUST be emitted bare (e.g. `C`, `O`, `Cl`) if and only if it
  qualifies for the organic-subset bare form: element in {B, C, N, O, P, S, F,
  Cl, Br, I}, charge 0, isotope 0, and its hydrogens are exactly reproducible
  by the m2a valence model from the emitted bonds (implicit_h ≥ 0 and
  explicit_h == 0). Otherwise it MUST be bracketed per m2a FR-3 grammar:
  `[isotope?symbol[Hn][charge]]`.
- FR-6b: Bracketed form MUST include explicit H only when `explicit_h > 0`
  (as `[OH2]`, `[NH4+]`) — implicit H of bracket atoms is always 0 by FR-3a
  of m2a, so round-trip holds exactly.

### ⏳ FR-7: Bond emission

- FR-7a: Single bonds between consecutive emissions MUST be omitted (the
  parser's default bond reproduces them).
- FR-7b: Double and triple bonds MUST always be emitted explicitly (`=`, `#`)
  except where FR-7c suppresses them.
- FR-7c: Ring-closure digits carry their bond symbol when the bond is not
  single (`C=1...C=1` style); single ring closures omit it.
- FR-7d: Bond symbols preceding branch openings follow the same rules
  (`C(=O)O`).

### ⏳ FR-8: Disconnected graphs

Graphs with zero bonds or multiple components (formula-derived multisets,
internal fallback paths) MUST canonicalize:

- FR-8a: Each connected component is canonicalized independently under FR-2 …
  FR-7.
- FR-8b: Component strings are joined left-to-right in ascending lexicographic
  order of their canonical strings, separated by `.`.
- FR-8c: Note: m2a's parser rejects dot-inputs; dots in canonical output arise
  solely from disconnected input graphs constructed programmatically.

### ⏳ FR-9: Round-trip invariant

For every graph the engine can construct, the following MUST hold:

```
canonicalSmiles(parseSmiles(canonicalSmiles(g))) == canonicalSmiles(g)
```

and additionally the composition of the reparsed graph equals the original
composition. This fixed-point property is tested corpus-wide (FR-12b) and is
what makes IDs stable forever (design.md §4.2).

### ⏳ FR-10: Golden corpus data

`data/golden/corpus.csv` MUST exist with header row:

```
id,input_format,input,expected_canonical,spec_version
```

- FR-10a: Exactly 30 curated entries at M2 completion, UTF-8, LF endings.
- FR-10b: Coverage MUST include, at minimum:
  - symmetry/ring systems: benzene `C1=CC=CC=C1`, naphthalene kekulé,
    cyclohexane, fused bicyclics
  - symmetric acyclics: neopentane, dimethyl ether vs ethanol (isomer pair)
  - charged species: `[NH4+]`, `[OH-]`, acetate anion kekulé
  - isotope variants: `[13CH4]`, `[2H]O[2H]`-style labeled water
  - formula-vs-SMILES equivalence pairs (same expected_canonical reached via
    both `input_format` values)
  - single atoms (`[C]`, `O`), hypervalent sulfur (`CS(=O)(=O)O`),
    multi-digit ring closures (`%nn`), deep branching
- FR-10c: Every entry's `expected_canonical` was produced by this
  implementation and hand-reviewed for chemical sanity before freezing;
  entries record `spec_version=1`.

### ⏳ FR-11: Corpus build integration

- FR-11a: `configure_file(... COPYONLY)` copies `data/golden/corpus.csv` next
  to the build (mirroring the elements.csv pattern) with a compile definition
  exposing its path to tests.
- FR-11b: `tests/canonical/golden_corpus_test.cpp` loads the CSV, and for each
  entry parses via `parseFormula` or `parseSmiles` per `input_format`, asserts
  `canonicalSmiles(...) == expected_canonical`, and asserts the FR-9 round-trip.
- FR-11c: Corpus failures are treated as release blockers: a failing entry
  means the implementation changed identity, and per design.md §4.2 the change
  is wrong until the spec is formally re-versioned.

### ⏳ FR-12: Test coverage

- FR-12a: Unit tests: order-independence (same molecule built with different
  insertion orders yields identical strings), each FR-6 emission rule, FR-8
  disconnected handling, empty-graph rejection (FR-1b).
- FR-12b: Golden corpus suite per FR-11b, including round-trip assertions.
- FR-12c: Determinism smoke: run the full corpus twice per test session and
  require byte-equal results.
- FR-12d: Equivalence classes: `O` / `[OH2]`; `CCO` / ethanol formula path;
  benzene written with different starting atoms/rotations — all collapse to
  one string per class.

## 5. Non-Functional Requirements

- **Determinism**: pure function; no global state; no time/locale dependence.
- **Correctness over speed**: refinement capped at 64 iterations keeps worst
  cases bounded; no optimization pressure in M2.
- **Freeze discipline**: golden failures block merges; version bump requires
  migration story (design.md §4.2).
- **Zero dependencies**, **Safety**, **Portability**: per umbrella spec §5.

## 6. Technical Constraints & Architecture Notes

- Namespace `chem`; naming per AGENTS.md; internal helpers in `chem::detail`
  or anonymous namespaces.
- The canonicalizer MUST consume only public `MolecularGraph` accessors
  (`atoms()`, `bonds()`) — no friendship, no core internals.
- It MUST NOT include parsing headers; the round-trip tests live in `tests/`
  and are the only place parser and canonicalizer meet.
- Output strings MUST parse back through `parseSmiles` (or survive as
  disconnected forms for formula fallbacks); the writer never emits syntax
  outside the m2a dialect plus the FR-8 dot rule.

## 7. Data & Interfaces

| Header | Interface |
|---|---|
| `canonical/canonical_smiles.hpp` | `inline constexpr std::uint32_t kCanonicalSpecVersion`, `std::string canonicalSmiles(const MolecularGraph&)` |

Data file: `data/golden/corpus.csv`, columns
`id,input_format,input,expected_canonical,spec_version`, 30 rows + header.

Directory delta:

```
src/chem/canonical/canonical_smiles.hpp/cpp   (new)
tests/canonical/canonical_smiles_test.cpp     (new)
tests/canonical/golden_corpus_test.cpp        (new)
data/golden/corpus.csv                        (new)
```

## 8. Edge Cases & Error Handling

| Input graph | Behavior |
|---|---|
| Empty graph | `ValidationError` |
| Single bare carbon, no bonds (formula `C`) | `"C"` — writer emits bare form; parser re-derives implicit_h = 4; fixed point holds |
| Water built `O,H,H` unconnected (formula path) | Components canonicalized independently and lexicographically sorted: `"O.[H].[H]"` (hydrogen requires brackets, FR-6a) |
| `[13CH4]` | Bracketed emission (isotope forces brackets) |
| Benzene entered from any of its rotations | Byte-identical canonical string |
| Graph with >99 simultaneous ring closures | `ValidationError` (FR-5d guard) |

Note on bare-form emission: the writer relies on the m2a valence model to
reproduce hydrogens, so a bare atom's `implicit_h` is never written — it is
re-derived identically on re-parse (this is what makes FR-9 hold without
explicit H bookkeeping).

## 9. Acceptance Criteria / Definition of Done

- [ ] All functional requirements FR-1 … FR-12 satisfied
- [ ] Golden corpus passes unchanged on repeated runs; round-trip green corpus-wide
- [ ] Order-independence tests pass (insertion-order permutations)
- [ ] Warning-clean build; full `ctest` green under sanitizers
- [ ] `clang-format -i` applied; `clang-tidy` clean on changed files
- [ ] `data/golden/corpus.csv` contains exactly 30 entries, spot-checked for
      chemical sanity (benzene, naphthalene, charged, isotopes, isomer pairs)
- [ ] No new dependencies; layout matches §7

## 10. Open Questions / Decisions Needed

None. Resolved during drafting:

- Algorithm pinned precisely rather than contract-only (user decision): exact
  invariant tuple, tuple-based (hash-free) refinement, rank ordering, DFS
  emission rules are all normative above.
- Corpus size fixed at 30 entries (user decision).

## 11. Implementation Notes

- Ranking pass: represent invariants as `std::vector<std::tuple<...>>` and
  refine with sort-based partitioning; 64-iteration cap is generous (WL
  stabilizes in ≤ diameter iterations for molecular graphs).
- Automorphic ties (FR-4c): after stabilization, atoms sharing an invariant
  in a fully refined partition are automorphic for our purposes; picking the
  lowest index keeps code simple and output stable.
- Writer pass: recursive DFS with a `visited` array, `on_path` marker for
  ring detection, and a free-digit counter; branch parentheses track depth.
- Corpus curation tip: generate candidates with the implementation, then
  hand-verify each against known chemistry before committing; the review is
  the freeze, not the generator.
- External toolkits as curation aid: RDKit or Open Babel MAY be used while
  curating corpus entries to sanity-check that inputs parse as intended
  chemically. Their *canonical* forms necessarily differ from ours by design
  (design.md §4.2 non-goal: cross-toolkit compatibility) and MUST NOT be
  treated as expected output anywhere in this spec or its tests.
- Future-proofing: if a later milestone needs canonical ranks (e.g. reaction
  mapping), expose them as a separate function then — do not widen this API now.
