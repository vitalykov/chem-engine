# AGENTS.md

Guidance for AI coding agents working in this repository.

## Project

chem-engine — a correctness-first chemical kinetics simulation library in C++.
Read `docs/design.md` before making non-trivial changes; it is the source of truth
for architecture and contracts.

- Language: C++20
- Dependencies: none, except doctest (via FetchContent). Do not add libraries.
  Hand-roll numerics; the state vectors are small enough that this is cheap.
- Branches: work happens on feature branches off `main`. `proto` is a frozen
  prototype kept for reference only — never merge from it or build on it.

## Build and test

```sh
cmake -B build -S .
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Always build and run tests after changes. New code requires new tests.

## Code style

- Types, concepts: `PascalCase` (`MolecularGraph`, `KineticSystem`)
- Functions, variables: `snake_case` (`canonical_form`, `step_size`)
- Constants: `k` prefix + PascalCase (`kDefaultDt`)
- Namespace: `chem`
- Headers use `#pragma once`; one class per header where practical
- Include order: own header first, then project headers, then std
- Use `std::string_view` for non-owning string parameters
- No comments unless something is genuinely non-obvious. The design rationale
  lives in `docs/design.md`, not inline.
- No emojis anywhere (code, docs, commits)

## Error handling contract

- Setup phase (parsing, construction, validation): throw typed exceptions
  (`ParseError`, `ValidationError`) derived from `std::runtime_error`. Never
  throw string literals.
- Simulation hot loop: exception-free. Debug asserts for precondition
  violations only.

## Hard architectural rules

These protect contracts documented in `docs/design.md`; changing them requires a
deliberate design decision, not a drive-by edit:

1. **Canonicalization is frozen.** Never change the output of the canonical
   SMILES writer. Golden corpus tests lock it; if they fail, your change is
   wrong until the spec is formally re-versioned.
2. **Identity is canonical SMILES.** Equality/hash of molecules never compare
   composition, mass, or raw input strings.
3. **Parsers produce graphs.** Every input format (formula, SMILES, future
   formats) must go through `MolecularGraph`; no format-specific shortcuts to
   identity.
4. **Loud rejection of unsupported syntax.** Unknown/unsupported SMILES or
   formula constructs are errors, never silently ignored.
5. **Chemistry and numerics stay decoupled.** Integrators see
   `KineticSystem` + `std::span<double>`; nothing outside `numerics/`
   knows integration details.
6. **Conservation holds.** Element totals across species are invariant along
   any trajectory; debug builds assert this.

## Where to add things

| Adding... | Put it in | Notes |
|---|---|---|
| A new input format | `src/chem/parsing/` | Must emit `MolecularGraph` |
| A new integrator | `src/chem/numerics/` | Implement `Integrator`; add analytic-solution tests |
| Element data fixes | `data/elements.csv` | Precise masses; update tests if masses change |
| Golden corpus entries | `data/golden/` | Include the tricky case *and* its expected canonical form |
| Tests | `tests/` mirroring `src/chem/` structure | doctest |

## Commits

Follow Conventional Commits: `type(scope): imperative summary`. Types:
`feat`, `fix`, `refactor`, `test`, `docs`, `build`, `chore`. Scopes mirror
the module areas: `core`, `parsing`, `canonical`, `numerics`, `test`,
`build`, `docs`. Examples:

```
feat(parsing): Reject stereo markers with ParseError
feat(numerics): Add AdaptiveRKF45 step-size controller
fix(core): Fix composition map double-counting bracket hydrogens
```

Breaking changes require `!` after the type and a `BREAKING CHANGE:` footer.
Any change to canonical SMILES output is breaking by definition (it
invalidates stored identifiers) — see docs/design.md §4.2.

Commit only when explicitly asked by the user.

## Versioning

Semantic Versioning: `MAJOR.MINOR.PATCH`.

- Start at `0.x.y`; in 0.x, minor bumps may still break the API. Promote to
  `1.0.0` only once the golden corpus, round-trip invariant, and public API
  have survived real use.
- The library version and the canonicalization spec version are linked but
  distinct (docs/design.md §4.2). A canonical spec version bump always forces a
  SemVer major bump; never the reverse.
- Breaking changes must be visible in commit history via `!` or footer.

## Changelog

No `CHANGELOG.md` until 0.1.0 ships. From then on, follow Keep a Changelog:
sections `Added / Changed / Deprecated / Removed / Fixed / Security`,
newest first, one entry per release plus an `Unreleased` section.

Workflow: generate a draft from Conventional Commits (e.g. git-cliff), then
hand-curate before tagging. Every release entry must state the
canonicalization spec version and whether it changed since the previous
release.

## Verification checklist

Before finishing any task:

1. `cmake --build build --parallel` — clean compile, no new warnings
2. `ctest --test-dir build --output-on-failure` — all green
3. If parsing/canonicalization changed: golden corpus and round-trip tests
   specifically pass unchanged
4. No new dependencies in CMakeLists.txt
