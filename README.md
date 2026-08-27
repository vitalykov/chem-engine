# chem-engine

Chemical reactions simulation engine — a correctness-first kinetics simulation
library in C++20.

Status: early development (0.x). The current milestone provides element data,
molecular graphs, composition and molar-mass computation, a Hill-formula
parser, a SMILES parser covering the v1 dialect subset (atoms, bonds,
branches, ring closures, charges, isotopes), and a frozen canonical SMILES
canonicalizer with a golden corpus. Numerical integration arrives in a later
milestone; see `docs/design.md` for the full architecture.

## Build and test

Requires CMake >= 3.28 and a C++20 compiler (GCC or Clang):

```sh
cmake -B build -S .
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The only dependency is doctest (fetched automatically at configure time).
Debug builds run tests under ASan+UBSan by default.

## Usage

Parse a Hill formula into a molecular graph, then derive its elemental
composition and molar mass:

```cpp
#include <iostream>

#include "chem/core/composition.hpp"
#include "chem/core/errors.hpp"
#include "chem/parsing/formula_parser.hpp"

int main() {
  try {
    const chem::MolecularGraph glucose = chem::parseFormula("C6H12O6");
    // Ordered map of atomic number -> atom count: {{1,12},{6,6},{8,6}}
    const chem::CompositionMap counts = chem::composition(glucose);
    const double grams_per_mole = chem::molarMass(glucose);  // ~180.156

    std::cout << "atoms: " << glucose.atoms().size() << "\n";   // 24
    std::cout << "mass:  " << grams_per_mole << " g/mol\n";
  } catch (const chem::ParseError& e) {
    std::cerr << "invalid input: " << e.what() << "\n";
  }
}
```

Invalid input is rejected loudly — unsupported syntax throws `chem::ParseError`
with a message quoting the offending token and position; nothing is silently
ignored.

## SMILES

A connected SMILES molecule parses into a bonded `MolecularGraph` with the
v1 dialect subset — bare organic atoms (`B C N O P S F Cl Br I`), brackets
with all 118 elements, single/double/triple bonds, branches, ring closures
(including `%nn`), formal charges, and isotope labels. Unsupported constructs
(aromatic atoms, stereochemistry, disconnected components, wildcards,
reaction SMILES) throw `chem::ParseError`:

```cpp
#include <iostream>

#include "chem/core/errors.hpp"
#include "chem/parsing/smiles_parser.hpp"

int main() {
  try {
    const chem::MolecularGraph acetic = chem::parseSmiles("CC(=O)O");
    // 4 atoms, 3 bonds; carbon 0 has implicit_h 3, the carbonyl C has 0,
    // hydroxyl O has implicit_h 1.
    std::cout << "atoms: " << acetic.atoms().size() << "\n";  // 4
  } catch (const chem::ParseError& e) {
    std::cerr << "invalid SMILES: " << e.what() << "\n";
  }
}
```

Because every parser emits the same graph type, `parseSmiles` and
`parseFormula` compose with the same composition/molar-mass helpers.

## Canonical SMILES

`chem::canonicalSmiles` turns any `MolecularGraph` into a deterministic,
automorphism-invariant identifier string (spec version
`chem::kCanonicalSpecVersion`). Equivalent spellings collapse to one ID, so
`Molecule("O") == Molecule("[OH2]")` holds once the `Molecule` type lands; for
now the canonicalizer is a pure graph->string function:

```cpp
#include <iostream>

#include "chem/canonical/canonical_smiles.hpp"
#include "chem/parsing/smiles_parser.hpp"

int main() {
  // Two kekule spellings of benzene canonicalize identically.
  std::cout << chem::canonicalSmiles(chem::parseSmiles("C1=CC=CC=C1")) << "\n";  // C1=CC=CC=C1
  std::cout << chem::canonicalSmiles(chem::parseSmiles("C=1C=CC=CC=1")) << "\n"; // C1=CC=CC=C1
  // Equivalent bracket/bare spellings normalize.
  std::cout << chem::canonicalSmiles(chem::parseSmiles("[OH2]")) << "\n";        // O
}
```

The canonicalization specification is frozen: its output is locked by a golden
corpus (`data/golden/corpus.csv`) run on every build, and any change to the
output requires a formal spec-version bump (`docs/design.md` §4.2).

Elements are flyweight handles into a built-in periodic table, looked up once
by symbol or atomic number:

```cpp
#include "chem/core/element.hpp"

const chem::Element oxygen("O");
oxygen.symbol();         // "O"
oxygen.atomicNumber();   // 8
oxygen.standardWeight(); // 15.999 g/mol
```

Graphs can also be assembled directly; every parser emits the same type:

```cpp
auto graph = chem::MolecularGraph();
graph.addAtom(chem::Atom{.element = chem::Element("O")});
graph.addAtom(chem::Atom{.element = chem::Element("H")});
graph.addBond(0, 1, chem::BondOrder::kSingle);
```

## Linking against the library

Consume via CMake (`add_subdirectory`) — headers export automatically with the
`chem` target:

```cmake
add_subdirectory(chem-engine)
target_link_libraries(my_app PRIVATE chem)  # C++20 required
```

Element data is read from `data/elements.csv`; the path is baked in at
configure time, so binaries linked inside this build tree work out of the box.
Consumers embedding the sources elsewhere must provide the same compile
definition (`CHEM_ELEMENTS_CSV_PATH`).

## Layout

| Path | Contents |
|---|---|
| `src/chem/core/` | Elements, molecular graphs, composition, error types |
| `src/chem/parsing/` | Input-format parsers producing `MolecularGraph` (Hill formula, SMILES) |
| `src/chem/canonical/` | Canonical SMILES writer (frozen, versioned spec) |
| `tests/` | doctest suites mirroring `src/chem/` |
| `data/` | Periodic table (`elements.csv`), golden canonical corpus (`golden/corpus.csv`) |
| `docs/` | Architecture and per-milestone specifications |

Architecture and contracts live in `docs/design.md`; the current milestone
contract is `docs/specs/m2b_canonicalization_spec.md`.

## License

Licensed under the [Apache License 2.0](LICENSE).
