#include "SRothman/EECs/interface/ProjectedENC.h"

using namespace eec;

int main() {
  ProjectedENC<2, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Unbinned, NoExtra, strategy::PairBins> c(1.0);
  return static_cast<int>(c({}).entries.size());
}
