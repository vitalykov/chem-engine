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
  inline ElementInfo GetElementInfo(const std::string& symbol) const {
    return elements_.at(symbol);
  }
  inline int GetElementNumber(const std::string& symbol) const {
    return elements_.at(symbol).index;
  }
  inline double GetElementMass(const std::string& symbol) const {
    return elements_.at(symbol).mass;
  }

private:
  std::unordered_map<std::string, ElementInfo> elements_;
};

}  // namespace chem
