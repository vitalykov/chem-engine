// SMILES parser — v1 dialect subset (docs/specs/m2a_smiles_parser_spec.md).
//
// Parser conformance is grounded in OpenSMILES (https://opensmiles.org) and the
// Daylight SMILES theory manual (https://www.daylight.com/dayhtml/doc/theory/
// theory.smiles.html). Three deliberate tightenings beyond those references:
//   1. charge magnitude limited to 9 (references permit ±15);
//   2. repeated-sign charge limited to ++/--, so [Fe+++] is an error;
//   3. over-valent bare atoms are errors (OpenSMILES falls back to zero implicit H).
// See docs/design.md §4.3 for the rationale.

#include "chem/parsing/smiles_parser.hpp"

#include <array>
#include <cstddef>
#include <iterator>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "chem/core/element.hpp"
#include "chem/core/errors.hpp"
#include "chem/core/molecular_graph.hpp"

namespace chem {

namespace {

// ---------------------------------------------------------------------------
// Part A — infrastructure (zero-based offsets, FR-11)
// ---------------------------------------------------------------------------

struct Cursor {
  std::string_view input;
  std::size_t pos;
};

bool isUpper(char c) { return c >= 'A' && c <= 'Z'; }

bool isLower(char c) { return c >= 'a' && c <= 'z'; }

bool isDigit(char c) { return c >= '0' && c <= '9'; }

std::string describeChar(char c) {
  const auto byte = static_cast<unsigned char>(c);
  if (byte >= 0x20 && byte < 0x7f) {
    return std::string("'") + c + "'";
  }
  std::ostringstream os;
  os << "(byte 0x" << std::hex << static_cast<int>(byte) << ")";
  return os.str();
}

std::string quotedAt(std::string_view input, std::size_t pos) {
  if (pos >= input.size()) {
    return "(end of input)";
  }
  return describeChar(input[pos]);
}

std::string errorAt(std::string_view message, std::string_view input, std::size_t pos) {
  std::ostringstream os;
  os << message << " at offset " << pos << ": " << quotedAt(input, pos) << " in \"" << input
     << "\"";
  return os.str();
}

// ---------------------------------------------------------------------------
// Part B — element resolution (FR-2)
// ---------------------------------------------------------------------------

struct OrganicEntry {
  std::string_view symbol;
  int atomic_number;
};

// The bare organic subset (FR-2). Two-char forms must be tried first.
constexpr std::array<OrganicEntry, 10> kOrganicSubset{{
    {"Cl", 17},
    {"Br", 35},
    {"B", 5},
    {"C", 6},
    {"N", 7},
    {"O", 8},
    {"P", 15},
    {"S", 16},
    {"F", 9},
    {"I", 53},
}};

// Returns the atomic number for a bare organic symbol, or 0 if not in the subset.
int organicAtomicNumber(std::string_view symbol) {
  for (const OrganicEntry& entry : kOrganicSubset) {
    if (entry.symbol == symbol) {
      return entry.atomic_number;
    }
  }
  return 0;
}

// Reads a bare organic atom at the cursor, advancing past it.
// Throws ParseError on any non-organic symbol or stray lowercase letter (FR-2b/FR-10a).
Atom parseBareAtom(Cursor& cursor) {
  const std::size_t start = cursor.pos;
  const char first = cursor.input[cursor.pos];
  ++cursor.pos;
  if (cursor.pos < cursor.input.size() && isLower(cursor.input[cursor.pos])) {
    const std::string_view two = cursor.input.substr(start, 2);
    if (organicAtomicNumber(two) != 0) {
      ++cursor.pos;
      return Atom{.element = Element(organicAtomicNumber(two))};
    }
    // Uppercase start (e.g. 'Fe', 'Co') means a non-organic element (FR-2b);
    // genuine aromatic atoms start lowercase and are rejected in rejectByte.
    throw ParseError(errorAt("non-organic element requires brackets", cursor.input, start));
  }
  const int atomic = organicAtomicNumber(std::string_view(&first, 1));
  if (atomic == 0) {
    throw ParseError(errorAt("non-organic element requires brackets", cursor.input, start));
  }
  return Atom{.element = Element(atomic)};
}

// ---------------------------------------------------------------------------
// Part C — bracket atoms (FR-3)
// ---------------------------------------------------------------------------

// Reads an unsigned integer from a digit run; empty run returns 0.
int readDigits(Cursor& cursor) {
  const std::size_t start = cursor.pos;
  while (cursor.pos < cursor.input.size() && isDigit(cursor.input[cursor.pos])) {
    ++cursor.pos;
  }
  int value = 0;
  for (std::size_t i = start; i < cursor.pos; ++i) {
    value = value * 10 + (cursor.input[i] - '0');
  }
  return value;
}

// Reads the element symbol inside brackets: one uppercase + optional lowercase.
std::string_view readBracketSymbol(Cursor& cursor) {
  const std::size_t start = cursor.pos;
  if (cursor.pos >= cursor.input.size() || !isUpper(cursor.input[cursor.pos])) {
    throw ParseError(errorAt("expected element symbol", cursor.input, cursor.pos));
  }
  ++cursor.pos;
  if (cursor.pos < cursor.input.size() && isLower(cursor.input[cursor.pos])) {
    ++cursor.pos;
  }
  return cursor.input.substr(start, cursor.pos - start);
}

Atom parseBracketAtom(Cursor& cursor) {
  const std::size_t bracket_pos = cursor.pos;
  ++cursor.pos; // consume '['

  Atom atom{.element = Element(1)}; // placeholder; overwritten below

  // Wildcards are rejected loudly (FR-10e).
  if (cursor.pos < cursor.input.size() && cursor.input[cursor.pos] == '*') {
    throw ParseError(errorAt("wildcards unsupported", cursor.input, cursor.pos));
  }

  // Isotope (optional leading digits, FR-3c).
  if (cursor.pos < cursor.input.size() && isDigit(cursor.input[cursor.pos])) {
    const std::size_t iso_start = cursor.pos;
    const int isotope = readDigits(cursor);
    if (isotope > 999) {
      throw ParseError(errorAt("isotope exceeds 999", cursor.input, iso_start));
    }
    atom.isotope = isotope;
  }

  // Element symbol (FR-3, FR-10a rejects lowercase aromatic).
  const std::string_view symbol = readBracketSymbol(cursor);
  try {
    atom.element = Element(symbol);
  } catch (const ParseError&) {
    throw ParseError(errorAt("unknown element symbol \"" + std::string(symbol) + "\"", cursor.input,
                             cursor.pos - symbol.size()));
  }

  // Explicit hydrogens (optional 'H' then optional digit, FR-3a).
  if (cursor.pos < cursor.input.size() && cursor.input[cursor.pos] == 'H') {
    ++cursor.pos;
    int hcount = 1;
    if (cursor.pos < cursor.input.size() && isDigit(cursor.input[cursor.pos])) {
      hcount = readDigits(cursor);
    }
    atom.explicit_h = hcount;
  }

  // Charge (optional, FR-3b).
  if (cursor.pos < cursor.input.size()) {
    const char c = cursor.input[cursor.pos];
    if (c == '+' || c == '-') {
      const std::size_t charge_pos = cursor.pos;
      ++cursor.pos;
      const int sign = (c == '+') ? 1 : -1;
      int magnitude = 1;
      if (cursor.pos < cursor.input.size() && cursor.input[cursor.pos] == c) {
        // Repeated-sign form: ++ or -- only (FR-3b, [Fe+++] rejected).
        ++cursor.pos;
        magnitude = 2;
        if (cursor.pos < cursor.input.size() && cursor.input[cursor.pos] == c) {
          throw ParseError(
              errorAt("repeated charge sign limited to two", cursor.input, cursor.pos));
        }
      } else if (cursor.pos < cursor.input.size() && isDigit(cursor.input[cursor.pos])) {
        magnitude = readDigits(cursor);
        if (magnitude < 1 || magnitude > 9) {
          throw ParseError(errorAt("charge magnitude must be 1..9", cursor.input, charge_pos));
        }
      }
      atom.charge = sign * magnitude;
    }
  }

  // Reject unsupported bracket contents (FR-10e/f/h).
  if (cursor.pos < cursor.input.size()) {
    const char c = cursor.input[cursor.pos];
    if (c == '@') {
      throw ParseError(errorAt("stereo markers unsupported", cursor.input, cursor.pos));
    }
    if (c == '*') {
      throw ParseError(errorAt("wildcards unsupported", cursor.input, cursor.pos));
    }
    if (c == ':') {
      throw ParseError(errorAt("atom-class labels unsupported", cursor.input, cursor.pos));
    }
  }

  // Closing ']'.
  if (cursor.pos >= cursor.input.size() || cursor.input[cursor.pos] != ']') {
    throw ParseError(errorAt("unterminated bracket atom", cursor.input, bracket_pos));
  }
  ++cursor.pos; // consume ']'

  atom.implicit_h = 0; // FR-3a: no valence inference inside brackets
  return atom;
}

// ---------------------------------------------------------------------------
// Part D, E — bonds, branches, chains, ring closures
// ---------------------------------------------------------------------------

std::optional<BondOrder> bondFromChar(char c) {
  switch (c) {
  case '-':
    return BondOrder::kSingle;
  case '=':
    return BondOrder::kDouble;
  case '#':
    return BondOrder::kTriple;
  default:
    return std::nullopt;
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

struct RingSlot {
  bool open = false;
  std::ptrdiff_t atom = -1;
  std::optional<BondOrder> pending_order;
};

constexpr std::size_t kRingSlotCount = 100;

struct Parser {
  explicit Parser(std::string_view input) : cursor{input, 0} {}

  Cursor cursor;
  MolecularGraph graph;
  std::array<RingSlot, kRingSlotCount> rings;
  // For each atom: whether it was a bare organic atom (eligible for implicit H).
  std::vector<bool> bare_organic;

  // Reads a ring-closure digit or %nn at the cursor (FR-6).
  void handleRingDigit(std::ptrdiff_t current_atom, std::optional<BondOrder> pending_bond) {
    const std::size_t digit_pos = cursor.pos;
    int rnum = 0;
    if (cursor.input[cursor.pos] == '%') {
      ++cursor.pos;
      if (cursor.pos + 1 >= cursor.input.size() || !isDigit(cursor.input[cursor.pos]) ||
          !isDigit(cursor.input[cursor.pos + 1])) {
        throw ParseError(errorAt("'%nn' requires two digits", cursor.input, digit_pos));
      }
      rnum = (cursor.input[cursor.pos] - '0') * 10 + (cursor.input[cursor.pos + 1] - '0');
      cursor.pos += 2;
    } else {
      rnum = cursor.input[cursor.pos] - '0';
      ++cursor.pos;
    }

    if (rnum < 0 || rnum >= static_cast<int>(kRingSlotCount)) {
      throw ParseError(errorAt("ring number out of range", cursor.input, digit_pos));
    }
    RingSlot& slot = rings[static_cast<std::size_t>(rnum)];
    if (!slot.open) {
      // Ring open.
      slot.open = true;
      slot.atom = current_atom;
      slot.pending_order = pending_bond;
    } else {
      // Ring close.
      if (slot.atom == current_atom) {
        throw ParseError(errorAt("ring closure on the same atom", cursor.input, digit_pos));
      }
      BondOrder order = BondOrder::kSingle;
      if (pending_bond && slot.pending_order) {
        if (*pending_bond != *slot.pending_order) {
          throw ParseError(errorAt("mismatched ring-closure bond orders", cursor.input, digit_pos));
        }
        order = *pending_bond;
      } else if (pending_bond) {
        order = *pending_bond;
      } else if (slot.pending_order) {
        order = *slot.pending_order;
      }
      graph.addBond(slot.atom, current_atom, order);
      slot.open = false;
      slot.atom = -1;
      slot.pending_order.reset();
    }
  }

  // Parses a chain; `prev_atom` is the atom to attach the first parsed atom to
  // (-1 for the top-level call); `pending_bond` is the bond order for attachment.
  void parseChain(std::ptrdiff_t prev_atom, std::optional<BondOrder> pending_bond) {
    std::ptrdiff_t current = prev_atom;
    std::optional<BondOrder> bond = pending_bond;

    while (cursor.pos < cursor.input.size()) {
      const char c = cursor.input[cursor.pos];

      if (c == '(') {
        if (current < 0) {
          throw ParseError(errorAt("branch with no preceding atom", cursor.input, cursor.pos));
        }
        const std::size_t open_pos = cursor.pos;
        ++cursor.pos;
        if (cursor.pos < cursor.input.size() && cursor.input[cursor.pos] == ')') {
          throw ParseError(errorAt("empty branch", cursor.input, cursor.pos));
        }
        parseChain(current, bond);
        if (cursor.pos >= cursor.input.size() || cursor.input[cursor.pos] != ')') {
          throw ParseError(errorAt("unclosed branch", cursor.input, open_pos));
        }
        ++cursor.pos;
        bond.reset(); // branch consumed the pending bond
        continue;
      }

      if (c == ')') {
        if (prev_atom < 0) {
          throw ParseError(errorAt("unbalanced closing parenthesis", cursor.input, cursor.pos));
        }
        if (bond) {
          throw ParseError(
              errorAt("trailing bond with no following atom", cursor.input, cursor.pos));
        }
        return; // end of this branch; caller consumes ')'
      }

      if (auto order = bondFromChar(c)) {
        if (current < 0) {
          throw ParseError(errorAt("bond with no preceding atom", cursor.input, cursor.pos));
        }
        if (bond) {
          throw ParseError(errorAt("consecutive bond tokens", cursor.input, cursor.pos));
        }
        bond = order;
        ++cursor.pos;
        continue;
      }

      if (isDigit(c) || c == '%') {
        if (current < 0) {
          throw ParseError(
              errorAt("ring closure with no preceding atom", cursor.input, cursor.pos));
        }
        handleRingDigit(current, bond);
        bond.reset();
        continue;
      }

      if (c == '[') {
        const Atom atom = parseBracketAtom(cursor);
        const std::ptrdiff_t index = graph.addAtom(atom);
        bare_organic.push_back(false);
        attachIfFirst(current, index, bond);
        current = index;
        bond.reset();
        continue;
      }

      if (isUpper(c)) {
        const Atom atom = parseBareAtom(cursor);
        const std::ptrdiff_t index = graph.addAtom(atom);
        bare_organic.push_back(true);
        attachIfFirst(current, index, bond);
        current = index;
        bond.reset();
        continue;
      }

      // Loud rejection of everything else (FR-10).
      rejectByte(c);
    }

    // End of input reached.
    if (bond) {
      throw ParseError(errorAt("trailing bond with no following atom", cursor.input, cursor.pos));
    }
  }

  void attachIfFirst(std::ptrdiff_t prev, std::ptrdiff_t current, std::optional<BondOrder> bond) {
    if (prev < 0) {
      return; // first atom of the whole molecule
    }
    graph.addBond(prev, current, bond.value_or(BondOrder::kSingle));
  }

  [[noreturn]] void rejectByte(char c) {
    if (isLower(c)) {
      throw ParseError(errorAt("unsupported aromatic atom", cursor.input, cursor.pos));
    }
    switch (c) {
    case '@':
      throw ParseError(errorAt("stereo markers unsupported", cursor.input, cursor.pos));
    case '/':
    case '\\':
      throw ParseError(errorAt("directional bonds unsupported", cursor.input, cursor.pos));
    case '.':
      throw ParseError(errorAt("disconnected components unsupported", cursor.input, cursor.pos));
    case '>':
      throw ParseError(errorAt("reaction SMILES unsupported", cursor.input, cursor.pos));
    case '*':
      throw ParseError(errorAt("wildcards unsupported", cursor.input, cursor.pos));
    case ':':
      throw ParseError(errorAt("atom-class labels unsupported", cursor.input, cursor.pos));
    default:
      if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
        throw ParseError(errorAt("whitespace not allowed", cursor.input, cursor.pos));
      }
      throw ParseError(
          errorAt("unexpected character " + describeChar(c), cursor.input, cursor.pos));
    }
  }
};

// ---------------------------------------------------------------------------
// Part F — valence table (FR-7)
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

} // namespace

// ---------------------------------------------------------------------------
// Part G — entry point (FR-1, FR-10)
// ---------------------------------------------------------------------------

MolecularGraph parseSmiles(std::string_view input) {
  if (input.empty()) {
    throw ParseError("empty SMILES input at offset 0");
  }

  Parser parser{input};
  parser.parseChain(-1, std::nullopt);

  // FR-10g: any still-open ring digit at end of input is an error.
  for (std::size_t i = 0; i < kRingSlotCount; ++i) {
    if (parser.rings[i].open) {
      throw ParseError(errorAt("unclosed ring digit " + std::to_string(i), input, input.size()));
    }
  }

  // Part F — assign implicit hydrogens after full traversal (FR-7b).
  const std::ptrdiff_t atom_count = std::ssize(parser.graph.atoms());
  std::vector<int> bond_sums(static_cast<std::size_t>(atom_count), 0);
  for (const Bond& bond : parser.graph.bonds()) {
    const int v = bondOrderValue(bond.order);
    bond_sums[static_cast<std::size_t>(bond.a)] += v;
    bond_sums[static_cast<std::size_t>(bond.b)] += v;
  }

  MolecularGraph result;
  for (std::ptrdiff_t i = 0; i < atom_count; ++i) {
    Atom atom = parser.graph.atoms()[static_cast<std::size_t>(i)];
    if (parser.bare_organic[static_cast<std::size_t>(i)]) {
      const int atomic = atom.element.atomicNumber();
      const int b = bond_sums[static_cast<std::size_t>(i)];
      const int v = selectValence(atomic, b);
      if (v < 0) {
        throw ParseError(errorAt("over-valent atom (element " + std::to_string(atomic) +
                                     ", bond sum " + std::to_string(b) + ")",
                                 input, input.size()));
      }
      atom.implicit_h = v - b;
    }
    result.addAtom(atom);
  }
  for (const Bond& bond : parser.graph.bonds()) {
    result.addBond(bond.a, bond.b, bond.order);
  }

  return result;
}

} // namespace chem
