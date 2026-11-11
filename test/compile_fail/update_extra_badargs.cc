#include "SRothman/EECs/interface/ProjectedENCRuntime.h"

using namespace eec;

struct OneArgExtra {
  static constexpr std::size_t K = 1;
  void update(double) {}
  template <std::size_t N>
  std::array<double, 1> operator()(const JetData&, const std::array<std::uint32_t, N>&) const {
    return {{0.0}};
  }
};

int main() {
  ProjectedENCRuntimeN<weight::Energy, norm::SumEnergy, angle::DeltaR, output::Unbinned, OneArgExtra> c(3, 1.0);
  c.update_extra(0.5, 2);
}
