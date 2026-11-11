#include "SRothman/EECs/interface/ProjectedENC.h"

using namespace eec;

struct NonStaticWeight {
  double value(const fastjet::PseudoJet& p) const { return p.E(); }
};

int main() {
  ProjectedENC<2, NonStaticWeight, norm::SumEnergy, angle::DeltaR, output::Unbinned> c(1.0);
  return static_cast<int>(c({}).entries.size());
}
