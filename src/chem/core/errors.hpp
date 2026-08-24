#pragma once

#include <stdexcept>

namespace chem {

class ParseError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class ValidationError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

} // namespace chem
