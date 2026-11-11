#include "SRothman/EECs/interface/OnePointEC.h"

using namespace eec;

struct WrongSizeFrame {
  static constexpr std::size_t D = 3;
  std::array<double, 2> operator()(const fastjet::PseudoJet&) const { return {}; }
};

int main() {
  OnePointEC<weight::Energy, norm::SumEnergy, WrongSizeFrame, output::Unbinned> c(1.0);
  return static_cast<int>(c({}, WrongSizeFrame{}).entries.size());
}
