#ifndef SRothman_EECs_Policies_h
#define SRothman_EECs_Policies_h

#include "SRothman/EECs/interface/ParticleInfo.h"

#include <fastjet/PseudoJet.hh>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace eec {

  inline constexpr double kNoExternalNorm = std::numeric_limits<double>::quiet_NaN();

  // Per-particle quantity E_p: struct with `static double value(const PseudoJet&)`.
  namespace weight {
    struct Energy {
      static double value(const fastjet::PseudoJet& p) { return p.E(); }
    };
    struct Pt {
      static double value(const fastjet::PseudoJet& p) { return p.pt(); }
    };
    struct Charge {
      static double value(const fastjet::PseudoJet& p) { return particleCharge(p); }
    };
    struct ChargeEnergy {
      static double value(const fastjet::PseudoJet& p) { return particleCharge(p) * p.E(); }
    };
    struct ChargePt {
      static double value(const fastjet::PseudoJet& p) { return particleCharge(p) * p.pt(); }
    };
  }  // namespace weight

  // Jet normalisation E_Jet: struct with `static double value(particles, externalNorm)`.
  namespace norm {
    struct SumEnergy {
      static double value(const std::vector<fastjet::PseudoJet>& ps, double) {
        double s = 0;
        for (const auto& p : ps)
          s += p.E();
        return s;
      }
    };
    struct SumPt {
      static double value(const std::vector<fastjet::PseudoJet>& ps, double) {
        double s = 0;
        for (const auto& p : ps)
          s += p.pt();
        return s;
      }
    };
    struct External {
      static double value(const std::vector<fastjet::PseudoJet>&, double ext) {
        if (!std::isfinite(ext) || ext == 0)
          throw std::invalid_argument("eec: norm::External requires a finite, non-zero external normalisation");
        return ext;
      }
    };
    struct Unity {
      static double value(const std::vector<fastjet::PseudoJet>&, double) { return 1.0; }
    };
  }  // namespace norm

}  // namespace eec

#endif
