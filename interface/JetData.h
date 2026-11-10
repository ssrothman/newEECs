#ifndef SRothman_EECs_JetData_h
#define SRothman_EECs_JetData_h

#include <fastjet/PseudoJet.hh>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace eec {

  // Per-jet precomputed quantities shared by all correlators (structure of arrays).
  struct JetData {
    double E_jet = 0;
    double n = 1;
    const std::vector<fastjet::PseudoJet>* particles = nullptr;
    std::vector<double> w;  // (E_p / E_jet)^n
    std::vector<double> ux, uy, uz;
    std::vector<double> eta, phi;  // phi in [0, 2pi)

    std::size_t size() const { return w.size(); }
  };

  template <class Weight, class Norm>
  void fillJetData(JetData& jd, const std::vector<fastjet::PseudoJet>& ps, double n, double externalNorm) {
    const std::size_t M = ps.size();
    jd.particles = &ps;
    jd.n = n;
    jd.E_jet = Norm::value(ps, externalNorm);
    jd.w.resize(M);
    jd.ux.resize(M);
    jd.uy.resize(M);
    jd.uz.resize(M);
    jd.eta.resize(M);
    jd.phi.resize(M);
    const bool integerPower = (n == std::floor(n));
    for (std::size_t i = 0; i < M; ++i) {
      const auto& p = ps[i];
      double z = Weight::value(p) / jd.E_jet;
      if (n != 1.0) {
        // pow of a negative base (charge weights) is NaN for non-integer n.
        if (z < 0 && !integerPower)
          throw std::domain_error("eec: negative particle weight with non-integer power n");
        z = std::pow(z, n);
      }
      jd.w[i] = z;
      const double mod = p.modp();
      jd.ux[i] = p.px() / mod;
      jd.uy[i] = p.py() / mod;
      jd.uz[i] = p.pz() / mod;
      jd.eta[i] = p.pseudorapidity();
      jd.phi[i] = p.phi();
    }
  }

}  // namespace eec

#endif
