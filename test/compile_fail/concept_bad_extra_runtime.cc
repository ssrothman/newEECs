#include "SRothman/EECs/interface/ProjectedENCRuntime.h"

using namespace eec;

// Valid only for pairs, so the runtime-N calculator (N = 2..6) must reject it.
struct PairOnlyExtra {
  static constexpr std::size_t K = 1;
  std::array<double, 1> operator()(const JetData&, const std::array<std::uint32_t, 2>&) const { return {{0.0}}; }
};

int main() {
  ProjectedENCRuntimeN<weight::Energy, norm::SumEnergy, angle::DeltaR, output::Unbinned, PairOnlyExtra> c(3, 1.0);
  return static_cast<int>(c({}).entries.size());
}
