#ifndef SRothman_EECs_test_TestJets_h
#define SRothman_EECs_test_TestJets_h

#include "SRothman/EECs/interface/ParticleInfo.h"

#include <fastjet/PseudoJet.hh>

#include <random>
#include <vector>

namespace eectest {

  // Massive-ish particles scattered around (eta0, phi0) with random charges in {-1, 0, 1}.
  inline std::vector<fastjet::PseudoJet> randomJet(std::size_t M, unsigned seed, double eta0 = 0.3, double phi0 = 1.0,
                                                   double spread = 0.4) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(-spread, spread);
    std::exponential_distribution<double> pt(0.2);
    std::uniform_int_distribution<int> q(-1, 1);
    std::vector<fastjet::PseudoJet> out;
    for (std::size_t i = 0; i < M; ++i) {
      auto p = fastjet::PtYPhiM(0.5 + pt(rng), eta0 + u(rng), phi0 + u(rng), 0.14);
      p.set_user_info(new eec::EECParticleInfo(q(rng)));
      out.push_back(p);
    }
    return out;
  }

}  // namespace eectest

#endif
