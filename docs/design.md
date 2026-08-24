# chem-engine design

Chemical reactions simulation engine. This document defines the architecture,
conventions, and contracts for the `main` line of development. The `proto`
branch is a throwaway prototype that validated the domain model; nothing on it
is binding.

## 1. Vision and goals

chem-engine is a **correctness-first chemical kinetics simulator**. It evolves
mixtures of chemical species over time according to reaction rate laws.

Primary goals, in priority order:

1. **Physical correctness** — species identity is exact, reactions conserve
   elements, rate laws are explicit and documented, results are reproducible.
2. **Determinism** — the same program state and inputs always produce the same
   trajectory. No hidden global state, no time-dependent behavior.
3. **Real-time capability** — a simulation step must be cheap and predictable
   enough for interactive use, but never at the cost of correctness.
4. **Embeddability** — a clean C++ library with no heavyweight dependencies,
   usable as a component inside larger applications (games, tools, education).
5. **Database readiness** — compounds and reactions are serializable entities
   identified by stable IDs, designed from day one to back a compound/reaction
   database.

Non-goals for v1:

- Quantum chemistry, molecular mechanics, any spatial/geometry modeling.
- Pressure/volume work, heat of reaction, adiabatic energy balances
  (thermodynamics beyond Arrhenius temperature dependence).
- Reaction network discovery or mechanism inference.
- Cross-toolkit canonical identifier compatibility (see §4).

## 2. Key decisions

| Decision | Choice | Rationale |
|---|---|---|
| Primary goal | Correctness-first simulator | Real-time is a constraint, not the driver |
| Language | C++20 | `std::span`, concepts, designated initializers; widely supported |
| Dependencies | None, except doctest via FetchContent | Lightweight builds; hand-rolled numerics are trivial at this scale |
| Solver strategy | `Integrator` interface; fixed-step RK4 default, adaptive RKF45 second | Deterministic core, extensibility without rewrites |
| Species identity | Own canonical SMILES over a molecular graph | Unique for structural isomers; string comparison stays cheap |
| Input formats | Plain formulas (`H2SO4`) *and* SMILES (`CS(=O)(=O)O`) | Formula parser proven in proto; SMILES needed for uniqueness |
| Periodic table | All 118 elements, precise masses | Cheap now, removes a class of "element not found" failures |
| Error handling | Exceptions at setup, exception-free hot loop | Parse/validation errors happen once; per-tick code stays lean |
| Tests | doctest | Single-header friction profile fits zero-dep philosophy |

## 3. Domain model

Namespace: `chem`. Clean slate relative to the prototype; concepts carry over,
names and semantics are redesigned.

### 3.1 Element

Static access to periodic table data: symbol, atomic number, precise atomic
mass. Loaded once from bundled data (see §8). `Element` is a lightweight
flyweight value type constructible from a symbol or atomic number; instances
hold only a view into the shared table, so copying never duplicates data.

### 3.2 Molecule

A distinct chemical species. Constructed from either format:

```cpp
auto water  = chem::Molecule("H2O");            // Hill-formula notation
auto acetic = chem::Molecule("CC(=O)O");        // SMILES notation
```

Internally a `Molecule` owns:

- **Molecular graph** — atoms (element, charge, implicit/explicit hydrogens)
  and bonds (order). Built directly by the SMILES parser; built as an
  unconnected multiset by the formula parser.
- **Composition map** — element → count, derived from the graph. Used for
  molar mass, conservation validation, and reaction balancing checks.
- **Canonical SMILES** — computed once at construction; the molecule's
  identity (§4).
- Molar mass, derived.

Equality and hashing operate on canonical SMILES. Structural isomers
(ethanol vs dimethyl ether, both `C2H6O`) remain distinct because identity is
graph-based, not composition-based.

### 3.3 Reaction

A stoichiometric transformation plus kinetic parameters:

- **Stoichiometry** — map of `Molecule` → coefficient on each side.
- **Reaction orders** per reactant — independent of stoichiometric
  coefficients by default (power-law kinetics `rate = k(T) · Π cᵢ^orderᵢ`).
  Default order for each reactant equals its coefficient only if the caller
  asks for that simplification; otherwise orders are explicit.
- **Rate parameters** — Arrhenius pair `(A, Ea)`; a bare constant `k` is
  sugar for `A = k, Ea = 0`.
- Reversibility is modeled as two `Reaction` objects sharing stoichiometry;
  no special-cased reverse machinery in v1.

Construction validates immediately: equation parses, every participant
resolves to a known `Molecule`, the equation is balanced (element counts
equal on both sides), coefficients positive.

### 3.4 Vessel

The simulated reactor — successor of the prototype's `Flask`:

- Species concentrations (mol/L), volume (L), temperature `T` (K).
- Isothermal in v1; `T` is a settable property consulted when computing
  `k(T)` at each step.
- API for adding substances/reactions (both object-based and string-based),
  stepping the system, and reading state.

### 3.5 KineticSystem

Pure ODE formulation decoupled from any container semantics:

```
dc/dt = Sᵀ · r(c, T)
```

where `S` is the stoichiometric matrix and `r` the vector of power-law rates.
Integrators consume a `KineticSystem`; `Vessel` delegates to one. This split
keeps numerics testable without chemistry plumbing and vice versa.

## 4. Identity: SMILES dialect and normalization contract

### 4.1 Pipeline

All identifiers flow through one path:

```
any input (formula | user SMILES | drawer graph)
    → parse → validate → MolecularGraph
    → canonicalize (Morgan/Weisfeiler-Lehman refinement, deterministic traversal)
    → canonical SMILES string
```

Users never need to write canonical SMILES. `Molecule("O") ==
Molecule("[OH2]")` holds because both normalize to the same canonical form.

### 4.2 Canonicalization is a frozen, versioned specification

There is no universal canonical SMILES: RDKit's differs from Open Babel's.
Ours is defined by our implementation and becomes part of our schema — the
future database keys on it. Therefore:

- The algorithm is a pure function `(MolecularGraph) -> std::string`.
- Its output is locked by **golden-file tests**: a committed corpus of
  molecules with their expected canonical strings, run on every build.
- The spec carries a version number. Any change to canonicalization output
  requires a migration story for stored IDs before merging.
- Invariant tested on the whole corpus: `parse(canonical(m)) == m`
  (round-trip fidelity). This is what makes IDs stable forever.

### 4.3 Dialect subset (v1)

Supported grammar: atoms (bracketed and bare organic subset), single/double/
triple bonds, branches, ring closures, formal charges, isotopes.
Explicitly **not** supported in v1: stereochemistry (`@`, `/`, `\`),
aromaticity perception (use kekulé alternating bonds), wildcards, reaction
SMILES. Encountering unsupported syntax is a loud `ParseError`, never silent
acceptance — silent drops could collapse distinct molecules into one ID.

### 4.4 Documented contract for users

Three short guarantees published wherever the engine is documented:

1. **Normalization**: equivalent inputs yield identical IDs.
2. **Round-trip**: canonical output re-parses to an identical graph.
3. **Comparison**: equality/hashing are defined on canonical form.

Hand-writers of SMILES read this; users of formulas or a future drawer need
not care.

## 5. Chemistry model

- Concentration-based power-law kinetics, isothermal.
- Temperature dependence via Arrhenius: `k(T) = A·exp(−Ea/RT)`.
- Equilibrium emerges naturally from forward/reverse reaction pairs;
  equilibrium constants are not input separately in v1.
- Conservation: element totals across all species are invariant during
  integration; violations indicate a bug and are asserted in debug builds.

## 6. Numerics

### 6.1 Integrator interface

```cpp
class Integrator {
public:
  virtual ~Integrator() = default;
  virtual void step(KineticSystem& system, std::span<double> concentrations,
                    double dt) = 0;
};
```

### 6.2 Implementations

- **`FixedStepRK4`** (default): deterministic cost per step, real-time
  friendly, accuracy appropriate for smooth kinetic curves.
- **`AdaptiveRKF45`**: error-controlled step size for offline/batch runs;
  non-deterministic cost per wall-clock tick, so not the default.
- Negative concentration guard: clamped at zero after each accepted step
  (documented approximation; adaptive solvers reduce how often it triggers).
- Future: BDF-type solver for stiff systems (fast vs slow reactions coexisting).
  The interface above must accommodate it; stiffness detection is out of scope
  until then.

Lessons from the proto folded in here: no double-buffering tricks that depend
on hash-map iteration order; state lives in flat contiguous storage
(`std::vector<double>` keyed by species index), which integrators view through
`std::span`.

## 7. Error handling

Two-phase policy:

- **Setup phase** (parsing, constructing molecules/reactions/vessels):
  throws typed exceptions — `ParseError`, `ValidationError` — derived from
  `std::runtime_error`. Fail fast with messages quoting the offending input.
- **Simulation phase** (stepping, integration): exception-free. No allocation
  after warm-up where feasible; precondition violations are debug asserts.

Never `throw` a string literal (a proto bug); always typed exception types.

## 8. Architecture and layout

```
src/chem/
  core/          element.hpp/cpp, molecule graph, composition, vessel
  parsing/       formula parser, SMILES parser, tokenizer, dialect rules
  canonical/     canonical ranking + SMILES writer (the frozen spec)
  numerics/      kinetic_system, integrator interface, rk4, rkf45
tests/           doctest suites mirroring src structure
data/            elements.csv (118 elements, precise masses), golden corpus
```

Rules:

- Library target `chem` plus a test runner target, mirroring the proto's
  CMake shape.
- Parsing and canonicalization are separate stages with the graph between
  them; neither reaches into the other.
- New input formats plug in beside the existing parsers; they must produce a
  `MolecularGraph` and nothing else.
- New integrators implement the interface; nothing outside `numerics/` knows
  integration details.

## 9. Testing strategy

Framework: doctest, fetched via FetchContent, registered with CTest.

Test categories:

1. **Unit tests** per module (parsers, canonicalizer, reaction validation).
2. **Golden corpus** — committed set of molecules with frozen canonical
   SMILES outputs, including tricky cases: ring systems, symmetric molecules
   (benzene, naphthalene), charged species, isotope variants, formula-vs-
   SMILES equivalence pairs.
3. **Analytic solutions** — first-order decay has an exact exponential
   solution; A ⇌ B approach to equilibrium has closed forms. Integrators are
   graded against these (error bounds per step size).
4. **Invariant tests** — mass/element conservation along trajectories;
   round-trip parse(canonical(m)) == m; determinism (two runs, same output).
5. **Error-path tests** — malformed formulas/SMILES raise the right
   exception type with useful messages.

## 10. Persistence (database forward-compatibility)

The engine is built so a compound/reaction database can be added later
without schema surgery:

- `Molecule` and `Reaction` serialize to a stable text representation
  (canonical-SMILES-keyed records; reactions reference participants by ID).
- Canonicalization spec version is stamped into serialized data.
- The planned molecule drawer targets the same `MolecularGraph` type the
  parsers emit, so drawn molecules enter the pipeline unchanged.

Schema sketch (for the future DB, not v1):

```
molecules(id: canonical_smiles, spec_version, composition_json, mol_mass, source_format)
reactions(id, stoichiometry_json /* id -> coeff */, orders_json, arrhenius_A, arrhenius_Ea)
```

## 11. Roadmap

- **M1 — Foundations**: CMake skeleton, doctest wiring, full periodic table,
  `Element`, molecular graph, formula parser.
- **M2 — Identity**: SMILES parser, canonicalizer + golden corpus, round-trip
  tests, `Molecule` with equality/hash on canonical form.
- **M3 — Reactions**: equation parsing, balance validation, orders,
  `KineticSystem`, `Vessel`, `FixedStepRK4`, analytic-solution tests.
- **M4 — Completeness**: `AdaptiveRKF45`, Arrhenius `k(T)` and vessel
  temperature, serialization, documentation of the SMILES dialect contract.
- **Beyond**: stiff solver (BDF), thermodynamics (heat of reaction), compound/
  reaction database, molecule drawer feeding `MolecularGraph`.
