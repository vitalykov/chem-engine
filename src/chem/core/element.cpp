#include "chem/core/element.hpp"

#include <charconv>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>
#include <vector>

#include "chem/core/errors.hpp"

namespace chem {

struct ElementRecord {
  std::string_view symbol;
  int atomic_number;
  double standard_weight;
};

namespace {

constexpr std::size_t kExpectedElementCount = 118;

struct ElementTable {
  std::string storage;
  std::vector<ElementRecord> records;
  std::unordered_map<std::string_view, std::size_t> by_symbol;
};

std::vector<std::string_view> splitFields(std::string_view line, char delimiter) {
  std::vector<std::string_view> fields;
  std::size_t start = 0;
  while (true) {
    const std::size_t hit = line.find(delimiter, start);
    if (hit == std::string_view::npos) {
      fields.push_back(line.substr(start));
      return fields;
    }
    fields.push_back(line.substr(start, hit - start));
    start = hit + 1;
  }
}

template <typename T> T parseNumericField(std::string_view field) {
  if (field.empty()) {
    throw std::runtime_error("elements.csv: invalid numeric field ''");
  }
  T value{};
  const char* begin = field.data();
  const char* end = field.data() + field.size();
  const auto [ptr, ec] = std::from_chars(begin, end, value);
  if (ec != std::errc{} || ptr != end) {
    throw std::runtime_error("elements.csv: invalid numeric field '" + std::string(field) + "'");
  }
  return value;
}

ElementTable loadElementTable() {
  ElementTable table;

  constexpr const char* kCsvPath = CHEM_ELEMENTS_CSV_PATH;
  std::ifstream in(kCsvPath, std::ios::binary);
  if (!in) {
    throw std::runtime_error(std::string("cannot open element data file: ") + kCsvPath);
  }
  std::ostringstream buffer;
  buffer << in.rdbuf();

  table.storage = buffer.str();
  std::string_view remaining(table.storage);

  const std::string_view header = "symbol,atomic_number,name,standard_atomic_weight";
  if (remaining.substr(0, header.size()) != header) {
    throw std::runtime_error("elements.csv: unexpected header row");
  }
  remaining.remove_prefix(header.size());
  if (!remaining.empty() && remaining.front() == '\n') {
    remaining.remove_prefix(1);
  }

  std::size_t cursor = 0;
  while (cursor < remaining.size()) {
    const std::size_t newline = remaining.find('\n', cursor);
    const std::string_view line = remaining.substr(
        cursor, (newline == std::string_view::npos ? remaining.size() : newline) - cursor);
    if (newline == std::string_view::npos) {
      cursor = remaining.size();
    } else {
      cursor = newline + 1;
    }

    const std::vector<std::string_view> fields = splitFields(line, ',');
    if (fields.size() != 4 || fields[0].empty() || fields[2].empty()) {
      throw std::runtime_error("elements.csv: malformed row '" + std::string(line) + "'");
    }
    if (fields[1].empty() || fields[3].empty()) {
      throw std::runtime_error("elements.csv: malformed row '" + std::string(line) + "'");
    }

    const int atomic_number = parseNumericField<int>(fields[1]);
    const auto weight = parseNumericField<double>(fields[3]);
    if (atomic_number != static_cast<int>(table.records.size()) + 1) {
      throw std::runtime_error("elements.csv: non-sequential atomic number '" +
                               std::string(fields[1]) + "'");
    }
    if (!(weight > 0.0) || !std::isfinite(weight)) {
      throw std::runtime_error("elements.csv: invalid weight for '" + std::string(fields[0]) + "'");
    }

    table.by_symbol.emplace(fields[0], table.records.size());
    table.records.push_back(ElementRecord{fields[0], atomic_number, weight});
  }

  if (table.records.size() != kExpectedElementCount) {
    throw std::runtime_error("elements.csv: expected 118 data rows, found " +
                             std::to_string(table.records.size()));
  }
  if (table.by_symbol.size() != kExpectedElementCount) {
    throw std::runtime_error("elements.csv: duplicate symbols present");
  }
  return table;
}

const ElementTable& elementTable() {
  static const ElementTable table = loadElementTable();
  return table;
}

} // namespace

Element::Element(std::string_view symbol) {
  const ElementTable& table = elementTable();
  const auto it = table.by_symbol.find(symbol);
  if (it == table.by_symbol.end()) {
    throw ParseError("unknown element symbol \"" + std::string(symbol) + "\"");
  }
  record_ = &table.records[it->second];
}

Element::Element(int atomic_number) {
  const ElementTable& table = elementTable();
  if (atomic_number < 1 || atomic_number > static_cast<int>(table.records.size())) {
    throw ParseError("unknown element atomic number " + std::to_string(atomic_number));
  }
  record_ = &table.records[static_cast<std::size_t>(atomic_number) - 1];
}

std::string_view Element::symbol() const noexcept { return record_->symbol; }

int Element::atomicNumber() const noexcept { return record_->atomic_number; }

double Element::standardWeight() const noexcept { return record_->standard_weight; }

} // namespace chem
