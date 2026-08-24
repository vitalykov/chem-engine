#include <doctest/doctest.h>

#include <string_view>

#include "chem/core/element.hpp"
#include "chem/core/errors.hpp"

namespace {

using chem::Element;
using chem::ParseError;

TEST_CASE("Element lookup by symbol") {
  const Element hydrogen("H");
  CHECK(hydrogen.atomicNumber() == 1);
  CHECK(hydrogen.symbol() == "H");
  CHECK(hydrogen.standardWeight() == doctest::Approx(1.008));

  const Element oganesson("Og");
  CHECK(oganesson.atomicNumber() == 118);
  CHECK(oganesson.symbol() == "Og");

  const Element cobalt("Co");
  CHECK(cobalt.atomicNumber() == 27);
  CHECK(cobalt.symbol() == "Co");
}

TEST_CASE("Element lookup by atomic number") {
  const Element hydrogen(1);
  CHECK(hydrogen.symbol() == "H");
  CHECK(hydrogen.standardWeight() == doctest::Approx(1.008));

  const Element oganesson(118);
  CHECK(oganesson.symbol() == "Og");
  CHECK(oganesson.standardWeight() == doctest::Approx(294));
}

TEST_CASE("Both lookup paths agree") {
  const Element by_symbol("Fe");
  const Element by_number(26);
  CHECK(by_symbol.atomicNumber() == by_number.atomicNumber());
  CHECK(by_symbol.symbol() == by_number.symbol());
  CHECK(by_symbol.standardWeight() == doctest::Approx(by_number.standardWeight()));
}

TEST_CASE("Spot-checked standard atomic weights") {
  CHECK(Element("C").standardWeight() == doctest::Approx(12.011));
  CHECK(Element("Cl").standardWeight() == doctest::Approx(35.45));
  CHECK(Element("Tc").standardWeight() == doctest::Approx(98));
  CHECK(Element("Pm").standardWeight() == doctest::Approx(145));
}

TEST_CASE("Lookup is case-sensitive") {
  REQUIRE_THROWS_AS(Element("CO"), ParseError);
  REQUIRE_THROWS_AS(Element("co"), ParseError);
}

TEST_CASE("Unknown symbols throw ParseError quoting the token") {
  for (const std::string_view symbol : {"Xx", "Zz", "Hx"}) {
    try {
      const Element unknown(symbol);
      FAIL("expected ParseError for ", symbol);
    } catch (const ParseError& e) {
      const std::string message = e.what();
      CHECK_MESSAGE(message.find(symbol) != std::string::npos, message);
    }
  }
}

TEST_CASE("Out-of-range atomic numbers throw ParseError") {
  REQUIRE_THROWS_AS(Element(0), ParseError);
  REQUIRE_THROWS_AS(Element(-1), ParseError);
  REQUIRE_THROWS_AS(Element(119), ParseError);
}

TEST_CASE("Flyweight copies share data and compare equal") {
  const Element original("Au");
  const Element copy = original;
  CHECK(copy.atomicNumber() == original.atomicNumber());
  CHECK(copy.symbol() == original.symbol());
  CHECK(copy.standardWeight() == original.standardWeight());
  CHECK(sizeof(Element) == sizeof(void*));
}

TEST_CASE("Repeated lookups are stable") {
  for (int i = 0; i < 3; ++i) {
    CHECK(Element("Se").atomicNumber() == 34);
    CHECK(Element("Se").standardWeight() == doctest::Approx(78.971));
  }
}

} // namespace
