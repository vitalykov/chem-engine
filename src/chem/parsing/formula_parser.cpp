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

std::string errorAt(std::string_view message, std::string_view input, std::size_t pos) {
  std::ostringstream os;
  os << message << " at position " << (pos + 1) << " in \"" << input << "\"";
  return os.str();
}

Element resolveElement(std::string_view symbol, std::string_view input, std::size_t pos) {
  try {
    return Element{symbol};
  } catch (const ParseError&) {
    throw ParseError(errorAt("unknown element symbol \"" + std::string(symbol) + "\"", input, pos));
  }
}

std::string_view readSymbol(Cursor& cursor) {
  const std::size_t start = cursor.pos;
  ++cursor.pos;
  while (cursor.pos < cursor.input.size() && isLower(cursor.input[cursor.pos])) {
    ++cursor.pos;
  }
  return cursor.input.substr(start, cursor.pos - start);
}

int readCount(Cursor& cursor) {
  if (cursor.pos >= cursor.input.size() || !isDigit(cursor.input[cursor.pos])) {
    return 1;
  }
  const std::size_t start = cursor.pos;
  if (cursor.input[cursor.pos] == '0') {
    const bool leading_zeros =
        cursor.pos + 1 < cursor.input.size() && isDigit(cursor.input[cursor.pos + 1]);
    if (leading_zeros) {
      throw ParseError(errorAt("count with leading zeros", cursor.input, start));
    }
    throw ParseError(errorAt("zero count", cursor.input, start));
  }
  while (cursor.pos < cursor.input.size() && isDigit(cursor.input[cursor.pos])) {
    ++cursor.pos;
  }
  if (cursor.pos - start > 9) {
    throw ParseError(errorAt("count too large", cursor.input, start));
  }
  int value = 0;
  for (const char c : cursor.input.substr(start, cursor.pos - start)) {
    value = value * 10 + (c - '0');
  }
  return value;
}

void appendMultiplied(MolecularGraph& destination, const MolecularGraph& source, int multiplier) {
  for (int rep = 0; rep < multiplier; ++rep) {
    for (const Atom& atom : source.atoms()) {
      destination.addAtom(atom);
    }
  }
}

MolecularGraph parseSequence(Cursor& cursor, bool inside_group) {
  MolecularGraph graph;
  while (cursor.pos < cursor.input.size()) {
    const char c = cursor.input[cursor.pos];
    if (c == '(') {
      ++cursor.pos;
      if (cursor.pos < cursor.input.size() && cursor.input[cursor.pos] == ')') {
        throw ParseError(errorAt("empty group", cursor.input, cursor.pos));
      }
      MolecularGraph inner = parseSequence(cursor, true);
      if (cursor.pos >= cursor.input.size() || cursor.input[cursor.pos] != ')') {
        throw ParseError(errorAt("unclosed group, missing ')'", cursor.input, cursor.pos));
      }
      ++cursor.pos;
      const std::size_t group_end = cursor.pos;
      const int multiplier = readCount(cursor);
      const std::int64_t projected =
          std::ssize(graph.atoms()) + std::ssize(inner.atoms()) * multiplier;
      if (projected > kMaxFormulaAtoms) {
        throw ParseError(errorAt("formula exceeds atom limit", cursor.input, group_end));
      }
      appendMultiplied(graph, inner, multiplier);
    } else if (c == ')') {
      if (!inside_group) {
        throw ParseError(errorAt("unbalanced closing parenthesis", cursor.input, cursor.pos));
      }
      return graph;
    } else if (isUpper(c)) {
      const std::size_t symbol_start = cursor.pos;
      const std::string_view symbol = readSymbol(cursor);
      Element element = resolveElement(symbol, cursor.input, symbol_start);
      const int count = readCount(cursor);
      if (std::ssize(graph.atoms()) + count > kMaxFormulaAtoms) {
        throw ParseError(errorAt("formula exceeds atom limit", cursor.input, symbol_start));
      }
      for (int i = 0; i < count; ++i) {
        graph.addAtom(Atom{.element = element});
      }
    } else {
      throw ParseError(
          errorAt("unexpected character " + describeChar(c), cursor.input, cursor.pos));
    }
  }
  return graph;
}

} // namespace

MolecularGraph parseFormula(std::string_view input) {
  if (input.empty()) {
    throw ParseError("empty formula input");
  }
  Cursor cursor{input, 0};
  return parseSequence(cursor, false);
}

} // namespace chem
