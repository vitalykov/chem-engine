#pragma once

#include <unordered_map>
#include <string>

namespace chem {

class Elements {
public:
  struct ElementInfo {
    int index;
    double mass;
  };

  Elements();
  inline ElementInfo GetElementInfo(std::string_view symbol) const {
    return elements_.at(std::string(symbol));
  }
  inline int GetElementNumber(std::string_view symbol) const {
    return elements_.at(std::string(symbol)).index;
  }
  inline double GetElementMass(std::string_view symbol) const {
    return elements_.at(std::string(symbol)).mass;
  }

private:
  std::unordered_map<std::string, ElementInfo> elements_;
};

}  // namespace chem
