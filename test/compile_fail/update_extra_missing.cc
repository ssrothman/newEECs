#include "SRothman/EECs/interface/ProjectedENC.h"

using namespace eec;

int main() {
  ProjectedENC<2, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Unbinned> c(1.0);
  c.update_extra(0.5);
}
