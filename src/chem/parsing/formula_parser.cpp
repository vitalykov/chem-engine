#include "chem/parsing/formula_parser.hpp"

#include <cstdint>
#include <iterator>
#include <sstream>
#include <string>

#include "chem/core/element.hpp"
#include "chem/core/errors.hpp"

namespace chem {

namespace {

constexpr std::int64_t kMaxFormulaAtoms = 100'000;

struct Cursor {
  std::string_view input;
  std::size_t pos;
};

bool is_upper(char c) { return c >= 'A' && c <= 'Z'; }

bool is_lower(char c) { return c >= 'a' && c <= 'z'; }

bool is_digit(char c) { return c >= '0' && c <= '9'; }

std::string describe_char(char c) {
  const auto byte = static_cast<unsigned char>(c);
  if (byte >= 0x20 && byte < 0x7f) {
    return std::string("'") + c + "'";
  }
  std::ostringstream os;
  os << "(byte 0x" << std::hex << static_cast<int>(byte) << ")";
  return os.str();
}

std::string error_at(std::string_view message, std::string_view input, std::size_t pos) {
  std::ostringstream os;
  os << message << " at position " << (pos + 1) << " in \"" << input << "\"";
  return os.str();
}

Element resolve_element(std::string_view symbol, std::string_view input, std::size_t pos) {
  try {
    return Element{symbol};
  } catch (const ParseError&) {
    throw ParseError(
        error_at("unknown element symbol \"" + std::string(symbol) + "\"", input, pos));
  }
}

std::string_view read_symbol(Cursor& cursor) {
  const std::size_t start = cursor.pos;
  ++cursor.pos;
  while (cursor.pos < cursor.input.size() && is_lower(cursor.input[cursor.pos])) {
    ++cursor.pos;
  }
  return cursor.input.substr(start, cursor.pos - start);
}

int read_count(Cursor& cursor) {
  if (cursor.pos >= cursor.input.size() || !is_digit(cursor.input[cursor.pos])) {
    return 1;
  }
  const std::size_t start = cursor.pos;
  if (cursor.input[cursor.pos] == '0') {
    const bool leading_zeros =
        cursor.pos + 1 < cursor.input.size() && is_digit(cursor.input[cursor.pos + 1]);
    if (leading_zeros) {
      throw ParseError(error_at("count with leading zeros", cursor.input, start));
    }
    throw ParseError(error_at("zero count", cursor.input, start));
  }
  while (cursor.pos < cursor.input.size() && is_digit(cursor.input[cursor.pos])) {
    ++cursor.pos;
  }
  if (cursor.pos - start > 9) {
    throw ParseError(error_at("count too large", cursor.input, start));
  }
  int value = 0;
  for (std::size_t i = start; i < cursor.pos; ++i) {
    value = value * 10 + (cursor.input[i] - '0');
  }
  return value;
}

void append_multiplied(MolecularGraph& destination, const MolecularGraph& source, int multiplier) {
  for (int rep = 0; rep < multiplier; ++rep) {
    for (const Atom& atom : source.atoms()) {
      destination.add_atom(atom);
    }
  }
}

MolecularGraph parse_sequence(Cursor& cursor, bool inside_group) {
  MolecularGraph graph;
  while (cursor.pos < cursor.input.size()) {
    const char c = cursor.input[cursor.pos];
    if (c == '(') {
      ++cursor.pos;
      if (cursor.pos < cursor.input.size() && cursor.input[cursor.pos] == ')') {
        throw ParseError(error_at("empty group", cursor.input, cursor.pos));
      }
      MolecularGraph inner = parse_sequence(cursor, true);
      if (cursor.pos >= cursor.input.size() || cursor.input[cursor.pos] != ')') {
        throw ParseError(error_at("unclosed group, missing ')'", cursor.input, cursor.pos));
      }
      ++cursor.pos;
      const std::size_t group_end = cursor.pos;
      const int multiplier = read_count(cursor);
      const std::int64_t projected =
          std::ssize(graph.atoms()) + std::ssize(inner.atoms()) * multiplier;
      if (projected > kMaxFormulaAtoms) {
        throw ParseError(error_at("formula exceeds atom limit", cursor.input, group_end));
      }
      append_multiplied(graph, inner, multiplier);
    } else if (c == ')') {
      if (!inside_group) {
        throw ParseError(error_at("unbalanced closing parenthesis", cursor.input, cursor.pos));
      }
      return graph;
    } else if (is_upper(c)) {
      const std::size_t symbol_start = cursor.pos;
      const std::string_view symbol = read_symbol(cursor);
      Element element = resolve_element(symbol, cursor.input, symbol_start);
      const int count = read_count(cursor);
      if (std::ssize(graph.atoms()) + count > kMaxFormulaAtoms) {
        throw ParseError(error_at("formula exceeds atom limit", cursor.input, symbol_start));
      }
      for (int i = 0; i < count; ++i) {
        graph.add_atom(Atom{.element = element});
      }
    } else {
      throw ParseError(
          error_at("unexpected character " + describe_char(c), cursor.input, cursor.pos));
    }
  }
  return graph;
}

} // namespace

MolecularGraph parse_formula(std::string_view input) {
  if (input.empty()) {
    throw ParseError("empty formula input");
  }
  Cursor cursor{input, 0};
  return parse_sequence(cursor, false);
}

} // namespace chem
