#include "chem/canonical/canonical_smiles.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "chem/core/element.hpp"
#include "chem/core/errors.hpp"
#include "chem/core/molecular_graph.hpp"

namespace chem {
namespace detail {

// ---------------------------------------------------------------------------
// Valence model (local copy of the m2a table, needed for FR-6a bare-form tests)
// ---------------------------------------------------------------------------

// Returns the smallest allowed valence >= bond_sum, or -1 if none (over-valent).
int selectValence(int atomic_number, int bond_sum) {
  switch (atomic_number) {
  case 5: // B {3}
    return (bond_sum <= 3) ? 3 : -1;
  case 6: // C {4}
    return (bond_sum <= 4) ? 4 : -1;
  case 7: // N {3,5}
    if (bond_sum <= 3) {
      return 3;
    }
    return (bond_sum <= 5) ? 5 : -1;
  case 8: // O {2}
    return (bond_sum <= 2) ? 2 : -1;
  case 15: // P {3,5}
    if (bond_sum <= 3) {
      return 3;
    }
    return (bond_sum <= 5) ? 5 : -1;
  case 16: // S {2,4,6}
    if (bond_sum <= 2) {
      return 2;
    }
    if (bond_sum <= 4) {
      return 4;
    }
    return (bond_sum <= 6) ? 6 : -1;
  case 9:  // F {1}
  case 17: // Cl {1}
  case 35: // Br {1}
  case 53: // I {1}
    return (bond_sum <= 1) ? 1 : -1;
  default:
    return -1;
  }
}

bool isOrganicSubset(int atomic_number) {
  switch (atomic_number) {
  case 5:  // B
  case 6:  // C
  case 7:  // N
  case 8:  // O
  case 15: // P
  case 16: // S
  case 9:  // F
  case 17: // Cl
  case 35: // Br
  case 53: // I
    return true;
  default:
    return false;
  }
}

int bondOrderValue(BondOrder order) {
  switch (order) {
  case BondOrder::kSingle:
    return 1;
  case BondOrder::kDouble:
    return 2;
  case BondOrder::kTriple:
    return 3;
  }
  return 0;
}

// ---------------------------------------------------------------------------
// Adjacency list
// ---------------------------------------------------------------------------

struct Neighbor {
  std::ptrdiff_t atom;
  BondOrder order;
  std::ptrdiff_t bond; // index into graph.bonds()
};

struct Adjacency {
  explicit Adjacency(const MolecularGraph& graph)
      : adjacency(static_cast<std::size_t>(std::ssize(graph.atoms()))) {
    for (std::ptrdiff_t i = 0; i < std::ssize(graph.bonds()); ++i) {
      const Bond& bond = graph.bonds()[static_cast<std::size_t>(i)];
      const auto a = static_cast<std::size_t>(bond.a);
      const auto b = static_cast<std::size_t>(bond.b);
      adjacency[a].push_back({bond.b, bond.order, i});
      adjacency[b].push_back({bond.a, bond.order, i});
    }
  }

  std::vector<std::vector<Neighbor>> adjacency;
};

// ---------------------------------------------------------------------------
// FR-2: Initial invariants
// ---------------------------------------------------------------------------

struct InitialInvariant {
  int isotope;
  int atomic_number;
  int charge;
  int total_h;
  int degree;
};

// Lexicographic comparison of (isotope, atomic_number, charge, total_h, degree).
bool operator<(const InitialInvariant& lhs, const InitialInvariant& rhs) {
  return std::tie(lhs.isotope, lhs.atomic_number, lhs.charge, lhs.total_h, lhs.degree) <
         std::tie(rhs.isotope, rhs.atomic_number, rhs.charge, rhs.total_h, rhs.degree);
}

bool operator==(const InitialInvariant& lhs, const InitialInvariant& rhs) {
  return std::tie(lhs.isotope, lhs.atomic_number, lhs.charge, lhs.total_h, lhs.degree) ==
         std::tie(rhs.isotope, rhs.atomic_number, rhs.charge, rhs.total_h, rhs.degree);
}

std::vector<InitialInvariant> computeInitialInvariants(const MolecularGraph& graph,
                                                       const Adjacency& adj) {
  std::vector<InitialInvariant> invariants;
  invariants.reserve(static_cast<std::size_t>(std::ssize(graph.atoms())));
  for (std::size_t i = 0; i < graph.atoms().size(); ++i) {
    const Atom& atom = graph.atoms()[i];
    const int degree = static_cast<int>(std::ssize(adj.adjacency[i]));
    invariants.push_back({.isotope = atom.isotope,
                          .atomic_number = atom.element.atomicNumber(),
                          .charge = atom.charge,
                          .total_h = atom.implicit_h + atom.explicit_h,
                          .degree = degree});
  }
  return invariants;
}

// ---------------------------------------------------------------------------
// FR-3: Weisfeiler-Lehman refinement
// ---------------------------------------------------------------------------

struct RefinedInvariant {
  int prev_class;
  std::vector<std::pair<int, int>> neighbors; // (neighbor_class, bond_order_value), sorted
};

bool operator<(const RefinedInvariant& lhs, const RefinedInvariant& rhs) {
  if (lhs.prev_class != rhs.prev_class) {
    return lhs.prev_class < rhs.prev_class;
  }
  return lhs.neighbors < rhs.neighbors;
}

bool operator==(const RefinedInvariant& lhs, const RefinedInvariant& rhs) {
  return lhs.prev_class == rhs.prev_class && lhs.neighbors == rhs.neighbors;
}

// Assigns classes 0..k-1 by sorting invariants and relabeling with order preserved.
template <typename Inv>
int assignClasses(const std::vector<Inv>& invariants, std::vector<int>& classes) {
  const auto n = static_cast<std::ptrdiff_t>(std::ssize(invariants));
  std::vector<std::ptrdiff_t> order(static_cast<std::size_t>(n));
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](std::ptrdiff_t a, std::ptrdiff_t b) {
    const auto ua = static_cast<std::size_t>(a);
    const auto ub = static_cast<std::size_t>(b);
    return invariants[ua] < invariants[ub];
  });

  classes.assign(static_cast<std::size_t>(n), 0);
  int num_classes = 0;
  if (n > 0) {
    classes[static_cast<std::size_t>(order[0])] = 0;
    num_classes = 1;
    for (std::ptrdiff_t i = 1; i < n; ++i) {
      const auto iu = static_cast<std::size_t>(i);
      const auto prev_idx = static_cast<std::size_t>(order[iu - 1]);
      const auto cur_idx = static_cast<std::size_t>(order[iu]);
      if (!(invariants[cur_idx] == invariants[prev_idx])) {
        ++num_classes;
      }
      classes[cur_idx] = num_classes - 1;
    }
  }
  return num_classes;
}

// Refine classes until stable or 64 iterations reached (FR-3b).
// Returns final class indices. Asserts stabilization before the cap (FR-4b).
std::vector<int> refine(const Adjacency& adj, const std::vector<InitialInvariant>& initial) {
  const auto n = static_cast<std::ptrdiff_t>(adj.adjacency.size());

  std::vector<int> classes;
  static_cast<void>(assignClasses(initial, classes));

  constexpr int kMaxIterations = 64;
  for (int iter = 0; iter < kMaxIterations; ++iter) {
    std::vector<RefinedInvariant> refined(static_cast<std::size_t>(n));
    for (std::ptrdiff_t i = 0; i < n; ++i) {
      const auto ui = static_cast<std::size_t>(i);
      std::vector<std::pair<int, int>> neighbors;
      neighbors.reserve(adj.adjacency[ui].size());
      for (const Neighbor& nb : adj.adjacency[ui]) {
        const auto nb_idx = static_cast<std::size_t>(nb.atom);
        neighbors.emplace_back(classes[nb_idx], bondOrderValue(nb.order));
      }
      std::sort(neighbors.begin(), neighbors.end());
      refined[ui] = {.prev_class = classes[ui], .neighbors = std::move(neighbors)};
    }

    std::vector<int> new_classes;
    static_cast<void>(assignClasses(refined, new_classes));

    if (new_classes == classes) {
      return classes;
    }
    classes = std::move(new_classes);
  }
  // FR-4b: if we hit the cap without stabilizing, surface in debug.
  assert(false && "WL refinement did not stabilize within 64 iterations");
  return classes;
}

// ---------------------------------------------------------------------------
// FR-4: Canonical ranking (automorphism-invariant)
// ---------------------------------------------------------------------------

// WRITECOMPONENT requires a total rank order (no ties). Atoms of the same WL
// class are automorphic (spec §11); FR-4c allows any tie-break, but a plain
// index tie-break is NOT invariant under automorphism: two kekule spellings of
// benzene would canonicalize differently. To make the output invariant, the
// canonical string is the lexicographically smallest emitted string over all
// within-class rank permutations. The number of permutations is the product of
// the class factorials; pathological graphs (e.g. a 300-atom star used by the
// FR-5d stress test) are capped and fall back to index order, which is
// deterministic but not automorphism-minimal.
constexpr std::int64_t kMaxRankingCombinations = 100'000;

std::string writeComponent(const MolecularGraph& graph, const std::vector<int>& ranks);

std::string canonicalComponentString(const MolecularGraph& graph, const std::vector<int>& classes) {
  const auto n = static_cast<std::ptrdiff_t>(std::ssize(classes));
  const int num_classes =
      classes.empty() ? 0 : *std::max_element(classes.begin(), classes.end()) + 1;
  std::vector<std::vector<std::ptrdiff_t>> groups(static_cast<std::size_t>(num_classes));
  for (std::ptrdiff_t i = 0; i < n; ++i) {
    groups[static_cast<std::size_t>(classes[static_cast<std::size_t>(i)])].push_back(i);
  }
  for (auto& group : groups) {
    std::sort(group.begin(), group.end());
  }

  std::int64_t combos = 1;
  for (const auto& group : groups) {
    for (std::int64_t k = 2; k <= static_cast<std::int64_t>(group.size()); ++k) {
      combos *= k;
      if (combos > kMaxRankingCombinations) {
        break;
      }
    }
    if (combos > kMaxRankingCombinations) {
      break;
    }
  }

  std::vector<int> ranks(static_cast<std::size_t>(n), 0);
  if (combos > kMaxRankingCombinations) {
    // Fallback: index order within each class (deterministic).
    int next_rank = 0;
    for (const auto& group : groups) {
      for (const std::ptrdiff_t a : group) {
        ranks[static_cast<std::size_t>(a)] = next_rank++;
      }
    }
    return writeComponent(graph, ranks);
  }

  std::string best;
  bool have_best = false;

  std::function<void(std::size_t, int)> assign;
  assign = [&](std::size_t group_idx, int next_rank) {
    if (group_idx == groups.size()) {
      const std::string candidate = writeComponent(graph, ranks);
      if (!have_best || candidate < best) {
        best = candidate;
        have_best = true;
      }
      return;
    }
    const std::vector<std::ptrdiff_t> base = groups[group_idx];
    std::vector<std::ptrdiff_t> perm = base;
    do {
      for (std::size_t k = 0; k < perm.size(); ++k) {
        ranks[static_cast<std::size_t>(perm[k])] = next_rank + static_cast<int>(k);
      }
      assign(group_idx + 1, next_rank + static_cast<int>(perm.size()));
    } while (std::next_permutation(perm.begin(), perm.end()));
  };
  assign(0, 0);
  return best;
}

// ---------------------------------------------------------------------------
// FR-6: Atom emission
// ---------------------------------------------------------------------------

struct AtomEmission {
  bool bare; // organic-subset bare form (FR-6a)
};

// Bond sum of an atom within its component.
int bondSumOf(const MolecularGraph& graph, std::ptrdiff_t atom) {
  int sum = 0;
  for (const Bond& bond : graph.bonds()) {
    if (bond.a == atom || bond.b == atom) {
      sum += bondOrderValue(bond.order);
    }
  }
  return sum;
}

// Determines the emission form for an atom (FR-6a, FR-12d, design.md §4.1).
// Bare iff organic subset, charge 0, isotope 0, and the atom's stored
// hydrogens are compatible with the m2a valence model:
//   - total_h (implicit+explicit) equals the model re-derivation v - bond_sum
//     (equivalent spellings [OH2]/O, [CH4]/C normalize to the bare form), OR
//   - total_h == 0 (the formula path stores no hydrogens; bare form is used
//     as-if, re-deriving v - bond_sum on re-parse, per spec §8).
// Anything else (e.g. [OH3], charged, isotopically labelled) stays bracketed.
AtomEmission classifyAtom(const Atom& atom, int bond_sum) {
  const int atomic = atom.element.atomicNumber();
  if (!isOrganicSubset(atomic) || atom.charge != 0 || atom.isotope != 0) {
    return {.bare = false};
  }
  const int valence = selectValence(atomic, bond_sum);
  if (valence < 0) {
    return {.bare = false}; // over-valent: cannot be bare
  }
  const int total_h = atom.implicit_h + atom.explicit_h;
  if (total_h != valence - bond_sum && total_h != 0) {
    return {.bare = false}; // hydrogens not reproducible in bare form
  }
  return {.bare = true};
}

std::string formatCharge(int charge) {
  if (charge == 0) {
    return {};
  }
  const char sign = (charge > 0) ? '+' : '-';
  const int magnitude = (charge > 0) ? charge : -charge;
  if (magnitude == 1) {
    return {sign};
  }
  std::string result;
  result.push_back(sign);
  result.append(std::to_string(magnitude));
  return result;
}

char bondSymbol(BondOrder order) {
  switch (order) {
  case BondOrder::kSingle:
    return '\0'; // omitted (FR-7a)
  case BondOrder::kDouble:
    return '=';
  case BondOrder::kTriple:
    return '#';
  }
  return '\0';
}

std::string ringDigitString(int digit) {
  if (digit <= 9) {
    return {static_cast<char>('0' + digit)};
  }
  return "%" + std::to_string(digit);
}

// ---------------------------------------------------------------------------
// FR-5: Per-component traversal and writing
// ---------------------------------------------------------------------------

struct RingClosure {
  std::ptrdiff_t opener; // smaller dfs number (visited first)
  std::ptrdiff_t closer; // larger dfs number
  BondOrder order;
  int digit = -1;
};

struct DfsTree {
  std::vector<std::ptrdiff_t> parent;                // -1 for the root (FR-5a)
  std::vector<BondOrder> parent_order;               // bond order on the tree edge to parent
  std::vector<int> dfs_number;                       // visit index per atom (emission order)
  std::vector<std::vector<std::ptrdiff_t>> children; // tree children in rank order
  std::vector<RingClosure> closures;                 // back edges, each recorded once
};

// Single DFS pass: numbers atoms in visit order, identifies tree edges and
// back edges (ring closures). Each bond is classified exactly once: the bond
// that brings a vertex into the DFS tree is that vertex's tree bond (emitted
// inline or inside a branch); every other bond becomes a ring closure recorded
// at whichever endpoint encounters it first (guarded by a per-bond flag).
DfsTree buildTree(const MolecularGraph& graph, const Adjacency& adj, const std::vector<int>& ranks,
                  std::ptrdiff_t root) {
  const auto n = static_cast<std::ptrdiff_t>(std::ssize(graph.atoms()));
  DfsTree tree;
  tree.parent.assign(static_cast<std::size_t>(n), -1);
  tree.parent_order.assign(static_cast<std::size_t>(n), BondOrder::kSingle);
  tree.dfs_number.assign(static_cast<std::size_t>(n), -1);
  tree.children.resize(static_cast<std::size_t>(n));

  const auto num_bonds = static_cast<std::size_t>(std::ssize(graph.bonds()));
  std::vector<bool> bond_recorded(num_bonds, false);
  std::vector<std::ptrdiff_t> tree_bond(static_cast<std::size_t>(n), -1);

  struct Numberer {
    const Adjacency& adj;
    const std::vector<int>& ranks;
    DfsTree& tree;
    std::vector<bool>& bond_recorded;
    std::vector<std::ptrdiff_t>& tree_bond;
    int next_number = 0;

    Numberer(const Adjacency& a, const std::vector<int>& r, DfsTree& t, std::vector<bool>& recorded,
             std::vector<std::ptrdiff_t>& tbond)
        : adj(a), ranks(r), tree(t), bond_recorded(recorded), tree_bond(tbond) {}

    void visit(std::ptrdiff_t u) {
      const auto u_idx = static_cast<std::size_t>(u);
      tree.dfs_number[u_idx] = next_number++;
      std::vector<Neighbor> neighbors = adj.adjacency[u_idx];
      std::sort(neighbors.begin(), neighbors.end(), [&](const Neighbor& a, const Neighbor& b) {
        const auto ua = static_cast<std::size_t>(a.atom);
        const auto ub = static_cast<std::size_t>(b.atom);
        return ranks[ua] < ranks[ub];
      });
      for (const Neighbor& nb : neighbors) {
        const auto nb_idx = static_cast<std::size_t>(nb.atom);
        if (nb.bond == tree_bond[u_idx]) {
          continue; // the tree edge that brought u into the DFS tree
        }
        if (tree.dfs_number[nb_idx] < 0) {
          // Tree edge: u discovers v.
          tree.parent[nb_idx] = u;
          tree.parent_order[nb_idx] = nb.order;
          tree_bond[nb_idx] = nb.bond;
          tree.children[u_idx].push_back(nb.atom);
          visit(nb.atom);
        } else {
          // Back edge to an already-visited vertex: ring closure, recorded once.
          if (!bond_recorded[static_cast<std::size_t>(nb.bond)]) {
            bond_recorded[static_cast<std::size_t>(nb.bond)] = true;
            std::ptrdiff_t opener = u;
            std::ptrdiff_t closer = nb.atom;
            if (tree.dfs_number[u_idx] > tree.dfs_number[nb_idx]) {
              opener = nb.atom;
              closer = u;
            }
            tree.closures.push_back({.opener = opener, .closer = closer, .order = nb.order});
          }
        }
      }
    }
  };

  Numberer numberer(adj, ranks, tree, bond_recorded, tree_bond);
  numberer.visit(root);
  return tree;
}

// Assigns ring digits (FR-5d): closures are opened in emission order; each takes
// the lowest free digit 1-9 then %10-%99; the digit frees when the closure ends.
// More than 99 simultaneously open rings is a ValidationError.
void assignDigits(DfsTree& tree) {
  std::sort(tree.closures.begin(), tree.closures.end(),
            [&](const RingClosure& a, const RingClosure& b) {
              const auto a_opener = static_cast<std::size_t>(a.opener);
              const auto a_closer = static_cast<std::size_t>(a.closer);
              const auto b_opener = static_cast<std::size_t>(b.opener);
              const auto b_closer = static_cast<std::size_t>(b.closer);
              const std::pair<int, int> a_key{tree.dfs_number[a_opener], tree.dfs_number[a_closer]};
              const std::pair<int, int> b_key{tree.dfs_number[b_opener], tree.dfs_number[b_closer]};
              return a_key < b_key;
            });

  const auto n = static_cast<std::ptrdiff_t>(tree.dfs_number.size());
  std::array<bool, 100> in_use{};
  int open_count = 0;

  for (std::ptrdiff_t t = 0; t < n; ++t) {
    // Closures ending at atom with dfs number t free their digit first.
    for (RingClosure& c : tree.closures) {
      if (c.digit >= 0 && tree.dfs_number[static_cast<std::size_t>(c.closer)] == t) {
        in_use[static_cast<std::size_t>(c.digit)] = false;
        --open_count;
      }
    }
    // Closures opening at atom with dfs number t then take the lowest free digit.
    for (RingClosure& c : tree.closures) {
      if (tree.dfs_number[static_cast<std::size_t>(c.opener)] == t) {
        if (open_count >= 99) {
          throw ValidationError("canonicalSmiles: more than 99 simultaneous ring closures");
        }
        int digit = 1;
        while (digit <= 99 && in_use[static_cast<std::size_t>(digit)]) {
          ++digit;
        }
        c.digit = digit;
        in_use[static_cast<std::size_t>(digit)] = true;
        ++open_count;
      }
    }
  }
}

// Emits the canonical SMILES of one connected component (FR-5..FR-7).
std::string writeComponent(const MolecularGraph& graph, const std::vector<int>& ranks) {
  const auto n = static_cast<std::ptrdiff_t>(std::ssize(graph.atoms()));

  // Emission form of every atom (FR-6).
  std::vector<AtomEmission> emissions(static_cast<std::size_t>(n));
  for (std::ptrdiff_t i = 0; i < n; ++i) {
    const Atom& atom = graph.atoms()[static_cast<std::size_t>(i)];
    emissions[static_cast<std::size_t>(i)] = classifyAtom(atom, bondSumOf(graph, i));
  }

  // FR-5a: start at the lowest-ranked atom (rank 0; ranks are unique).
  auto root_it = std::min_element(ranks.begin(), ranks.end());
  const auto root = static_cast<std::ptrdiff_t>(std::distance(ranks.begin(), root_it));

  const Adjacency adj(graph);
  DfsTree tree = buildTree(graph, adj, ranks, root);
  assignDigits(tree);

  // Per-atom ring-closure digit tokens. Each closure contributes one token at
  // its opener (is_closing=false) and one at its closer (is_closing=true).
  // At a single atom, a digit may be freed by a close and immediately re-opened;
  // the close must be emitted before the open so the parser sees the slot free.
  struct DigitToken {
    int digit;
    BondOrder order;
    bool is_closing;
  };
  std::vector<std::vector<DigitToken>> atom_digits(static_cast<std::size_t>(n));
  for (const RingClosure& c : tree.closures) {
    atom_digits[static_cast<std::size_t>(c.opener)].push_back({c.digit, c.order, false});
    atom_digits[static_cast<std::size_t>(c.closer)].push_back({c.digit, c.order, true});
  }
  for (auto& list : atom_digits) {
    std::sort(list.begin(), list.end(), [](const DigitToken& a, const DigitToken& b) {
      if (a.digit != b.digit) {
        return a.digit < b.digit;
      }
      // Same digit: closes before opens.
      if (a.is_closing != b.is_closing) {
        return a.is_closing;
      }
      return bondOrderValue(a.order) < bondOrderValue(b.order);
    });
  }

  std::ostringstream out;

  struct Emitter {
    const MolecularGraph& graph;
    const DfsTree& tree;
    const std::vector<AtomEmission>& emissions;
    const std::vector<std::vector<DigitToken>>& atom_digits;
    std::ostringstream& out;

    void emitToken(std::ptrdiff_t u) {
      const Atom& atom = graph.atoms()[static_cast<std::size_t>(u)];
      const AtomEmission& em = emissions[static_cast<std::size_t>(u)];
      if (em.bare) {
        out << atom.element.symbol();
      } else {
        out << '[';
        if (atom.isotope != 0) {
          out << atom.isotope;
        }
        out << atom.element.symbol();
        const int h_count = atom.implicit_h + atom.explicit_h;
        if (h_count > 0) {
          out << 'H';
          if (h_count > 1) {
            out << h_count;
          }
        }
        if (atom.charge != 0) {
          out << formatCharge(atom.charge);
        }
        out << ']';
      }
    }

    void emitBond(BondOrder order) {
      const char symbol = bondSymbol(order);
      if (symbol != '\0') {
        out << symbol;
      }
    }

    // Emits the atom token, its ring-closure digits, and its descendants (FR-5c).
    // SMILES grammar attaches every "(" to the atom token immediately before it,
    // so all branch children (i >= 1) must be emitted right after the node's
    // token; only then does the first (inline) child continue the chain.
    void emitSubtree(std::ptrdiff_t u) {
      const auto u_idx = static_cast<std::size_t>(u);
      emitToken(u);
      for (const DigitToken& token : atom_digits[u_idx]) {
        emitBond(token.order); // FR-7c: symbol at both ends when non-single
        out << ringDigitString(token.digit);
      }
      const auto& children = tree.children[u_idx];
      for (std::size_t i = 1; i < children.size(); ++i) {
        const std::ptrdiff_t child = children[i];
        const auto child_idx = static_cast<std::size_t>(child);
        out << '(';
        emitBond(tree.parent_order[child_idx]); // FR-7d: branch bond inside parens
        emitSubtree(child);
        out << ')';
      }
      if (!children.empty()) {
        const std::ptrdiff_t child = children[0];
        const auto child_idx = static_cast<std::size_t>(child);
        emitBond(tree.parent_order[child_idx]); // first subtree inline
        emitSubtree(child);
      }
    }
  };

  Emitter emitter{graph, tree, emissions, atom_digits, out};
  emitter.emitSubtree(root);
  return out.str();
}

// ---------------------------------------------------------------------------
// FR-8: Connected components
// ---------------------------------------------------------------------------

std::vector<int> connectedComponents(const MolecularGraph& graph) {
  const auto n = static_cast<std::ptrdiff_t>(std::ssize(graph.atoms()));
  std::vector<int> component(static_cast<std::size_t>(n), -1);
  int current = 0;
  for (std::ptrdiff_t start = 0; start < n; ++start) {
    if (component[static_cast<std::size_t>(start)] != -1) {
      continue;
    }
    std::vector<std::ptrdiff_t> queue{start};
    component[static_cast<std::size_t>(start)] = current;
    std::size_t head = 0;
    while (head < queue.size()) {
      const std::ptrdiff_t u = queue[head];
      ++head;
      for (const Bond& bond : graph.bonds()) {
        std::ptrdiff_t other = -1;
        if (bond.a == u) {
          other = bond.b;
        } else if (bond.b == u) {
          other = bond.a;
        }
        if (other >= 0 && component[static_cast<std::size_t>(other)] == -1) {
          component[static_cast<std::size_t>(other)] = current;
          queue.push_back(other);
        }
      }
    }
    ++current;
  }
  return component;
}

} // namespace detail

// ---------------------------------------------------------------------------
// FR-1, FR-8: entry point
// ---------------------------------------------------------------------------

std::string canonicalSmiles(const MolecularGraph& graph) {
  if (graph.atoms().empty()) {
    throw ValidationError("canonicalSmiles: empty graph has no identity");
  }

  // FR-8a: canonicalize each connected component independently.
  const std::vector<int> components = detail::connectedComponents(graph);
  const int num_components =
      components.empty() ? 0 : *std::max_element(components.begin(), components.end()) + 1;

  std::vector<std::string> component_strings;
  component_strings.reserve(static_cast<std::size_t>(num_components));

  for (int comp = 0; comp < num_components; ++comp) {
    // Local index map for atoms of this component.
    std::vector<std::ptrdiff_t> local_index(graph.atoms().size(), -1);
    std::vector<std::ptrdiff_t> globals;
    for (std::size_t i = 0; i < components.size(); ++i) {
      if (components[i] == comp) {
        local_index[i] = static_cast<std::ptrdiff_t>(globals.size());
        globals.push_back(static_cast<std::ptrdiff_t>(i));
      }
    }

    MolecularGraph subgraph;
    for (const std::ptrdiff_t g : globals) {
      subgraph.addAtom(graph.atoms()[static_cast<std::size_t>(g)]);
    }
    for (const Bond& bond : graph.bonds()) {
      const std::ptrdiff_t local_a = local_index[static_cast<std::size_t>(bond.a)];
      const std::ptrdiff_t local_b = local_index[static_cast<std::size_t>(bond.b)];
      if (local_a >= 0 && local_b >= 0) {
        subgraph.addBond(local_a, local_b, bond.order);
      }
    }

    // Refine and rank within the component (FR-4: automorphism-minimal ordering).
    const detail::Adjacency sub_adj(subgraph);
    const std::vector<detail::InitialInvariant> initial =
        detail::computeInitialInvariants(subgraph, sub_adj);
    const std::vector<int> classes = detail::refine(sub_adj, initial);

    component_strings.push_back(detail::canonicalComponentString(subgraph, classes));
  }

  // FR-8b: components joined in ascending lexicographic order, '.'-separated.
  std::sort(component_strings.begin(), component_strings.end());
  std::string result;
  for (std::size_t i = 0; i < component_strings.size(); ++i) {
    if (i != 0) {
      result.push_back('.');
    }
    result += component_strings[i];
  }
  return result;
}

} // namespace chem
