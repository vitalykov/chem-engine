#pragma once

// Skeleton-key graph isomorphism for the golden corpus guard (FR-12e).
//
// Two MolecularGraphs are compared structurally: vertices carry a skeleton key
// (element, charge, isotope, per-bond-order degree) and edges are grouped by
// bond order. Isomorphism is decided by partitioning vertices into equal-key
// cells and backtracking over the within-cell permutations, checking edge
// preservation in BOTH directions at every step (a vertex bijection that
// preserves edges both ways is an isomorphism, not just a homomorphism).
//
// Simple graphs only: parallel bonds are not distinguished (a pair of atoms is
// either bonded or not per bond order). Every golden corpus entry is simple, so
// this is sufficient; a multigraph corpus entry would need per-edge counts.
//
// Stored hydrogens (implicit_h/explicit_h) are intentionally NOT part of the
// key: the canonicalizer's bare/bracket normalization ([OH2] -> O, [C] -> C)
// makes stored H counts a bookkeeping artifact, not a structural property.

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "chem/core/molecular_graph.hpp"

namespace chem::test::graph_iso {

struct SkeletonKey {
  int element; // atomic number
  int charge;
  int isotope;
  std::array<int, 3> degree; // [0]=single, [1]=double, [2]=triple
};

bool operator<(const SkeletonKey& lhs, const SkeletonKey& rhs) {
  return std::tie(lhs.element, lhs.charge, lhs.isotope, lhs.degree[0], lhs.degree[1],
                  lhs.degree[2]) < std::tie(rhs.element, rhs.charge, rhs.isotope, rhs.degree[0],
                                            rhs.degree[1], rhs.degree[2]);
}

struct SimpleGraph {
  int n = 0;
  std::vector<SkeletonKey> keys;
  std::array<std::vector<std::pair<int, int>>, 3> edges; // by bond order
};

SimpleGraph fromGraph(const chem::MolecularGraph& graph) {
  SimpleGraph result;
  result.n = static_cast<int>(graph.atoms().size());
  for (const auto& atom : graph.atoms()) {
    result.keys.push_back({atom.element.atomicNumber(), atom.charge, atom.isotope, {0, 0, 0}});
  }
  for (const auto& bond : graph.bonds()) {
    const int a = static_cast<int>(std::min(bond.a, bond.b));
    const int b = static_cast<int>(std::max(bond.a, bond.b));
    const auto order = static_cast<std::size_t>(bond.order);
    result.edges[order].emplace_back(a, b);
    ++result.keys[static_cast<std::size_t>(a)].degree[order];
    ++result.keys[static_cast<std::size_t>(b)].degree[order];
  }
  return result;
}

bool hasBond(const SimpleGraph& graph, int order, int a, int b) {
  for (const auto& edge : graph.edges[static_cast<std::size_t>(order)]) {
    if ((edge.first == a && edge.second == b) || (edge.first == b && edge.second == a)) {
      return true;
    }
  }
  return false;
}

bool isomorphic(const SimpleGraph& a, const SimpleGraph& b) {
  if (a.n != b.n) {
    return false;
  }

  // Partition vertices by skeleton key. Every key cell of A must exist in B
  // with the same size; cells then correspond 1:1 (std::map orders by key).
  std::map<SkeletonKey, std::vector<int>> cells_a;
  std::map<SkeletonKey, std::vector<int>> cells_b;
  for (int v = 0; v < a.n; ++v) {
    cells_a[a.keys[static_cast<std::size_t>(v)]].push_back(v);
  }
  for (int v = 0; v < b.n; ++v) {
    cells_b[b.keys[static_cast<std::size_t>(v)]].push_back(v);
  }
  std::vector<std::vector<int>> cell_vertices_a;
  std::vector<std::vector<int>> cell_vertices_b;
  for (const auto& [key, cell] : cells_a) {
    const auto it = cells_b.find(key);
    if (it == cells_b.end() || it->second.size() != cell.size()) {
      return false;
    }
    cell_vertices_a.push_back(cell);
    cell_vertices_b.push_back(it->second);
  }

  // Fixed assignment order on the A side (cell by cell); backtracking permutes
  // the B side within each cell only.
  const auto n = static_cast<std::size_t>(a.n);
  std::vector<int> order;
  std::vector<int> cell_of(n, 0);
  for (std::size_t c = 0; c < cell_vertices_a.size(); ++c) {
    for (const int v : cell_vertices_a[c]) {
      cell_of[static_cast<std::size_t>(v)] = static_cast<int>(c);
      order.push_back(v);
    }
  }

  std::vector<int> a_to_b(n, -1); // A vertex -> B vertex
  std::vector<int> b_to_a(n, -1); // inverse

  std::function<bool(std::size_t)> assign = [&](std::size_t idx) -> bool {
    if (idx == order.size()) {
      return true;
    }
    const int av = order[idx];
    const int cell = cell_of[static_cast<std::size_t>(av)];
    for (const int bv : cell_vertices_b[static_cast<std::size_t>(cell)]) {
      if (b_to_a[static_cast<std::size_t>(bv)] >= 0) {
        continue;
      }
      // Bidirectional edge check: every already-mapped A-neighbor of av must be
      // a B-neighbor of bv with the same bond order, and vice versa.
      bool ok = true;
      for (int k = 0; k < 3 && ok; ++k) {
        for (const auto& [u, w] : a.edges[static_cast<std::size_t>(k)]) {
          if (u != av && w != av) {
            continue;
          }
          if (u == w) { // self-loop (not expected in the simple corpus)
            if (!hasBond(b, k, bv, bv)) {
              ok = false;
            }
            continue;
          }
          const int other = (u == av) ? w : u;
          const int mapped = a_to_b[static_cast<std::size_t>(other)];
          if (mapped >= 0 && !hasBond(b, k, bv, mapped)) {
            ok = false;
            break;
          }
        }
        if (!ok) {
          break;
        }
        for (const auto& [u, w] : b.edges[static_cast<std::size_t>(k)]) {
          if (u != bv && w != bv) {
            continue;
          }
          if (u == w) { // self-loop
            if (!hasBond(a, k, av, av)) {
              ok = false;
            }
            continue;
          }
          const int other = (u == bv) ? w : u;
          const int mapped = b_to_a[static_cast<std::size_t>(other)];
          if (mapped >= 0 && !hasBond(a, k, av, mapped)) {
            ok = false;
            break;
          }
        }
      }
      if (!ok) {
        continue;
      }
      a_to_b[static_cast<std::size_t>(av)] = bv;
      b_to_a[static_cast<std::size_t>(bv)] = av;
      if (assign(idx + 1)) {
        return true;
      }
      a_to_b[static_cast<std::size_t>(av)] = -1;
      b_to_a[static_cast<std::size_t>(bv)] = -1;
    }
    return false;
  };

  return assign(0);
}

} // namespace chem::test::graph_iso
