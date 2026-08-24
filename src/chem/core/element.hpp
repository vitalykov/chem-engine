#pragma once

#include <string_view>

namespace chem {

struct ElementRecord;

// Flyweight handle into the shared periodic table. Copying copies one pointer.
class Element {
public:
  explicit Element(std::string_view symbol);
  explicit Element(int atomic_number);

  std::string_view symbol() const noexcept;
  int atomic_number() const noexcept;
  double standard_weight() const noexcept;

private:
  // Non-owning view into function-local-static table storage; never freed
  // through this pointer (FR-4a).
  const ElementRecord* record_;
};

} // namespace chem
