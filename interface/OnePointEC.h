#ifndef SRothman_EECs_OnePointEC_h
#define SRothman_EECs_OnePointEC_h

#include "SRothman/EECs/interface/Binning.h"
#include "SRothman/EECs/interface/JetData.h"
#include "SRothman/EECs/interface/Outputs.h"
#include "SRothman/EECs/interface/Policies.h"

#include <fastjet/PseudoJet.hh>

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace eec {

  // Frame policies: per-jet objects with `static constexpr std::size_t D` and
  // `std::array<double, D> operator()(const fastjet::PseudoJet&) const`.
  namespace frame {

    // Boosts particles into the rest frame of `rest`, then returns (theta, phi) with respect to
    // axes z and x given in that rest frame. x is orthogonalised against z; phi in [-pi, pi].
    class BoostRotate {
    public:
      static constexpr std::size_t D = 2;

      BoostRotate(const fastjet::PseudoJet& rest, std::array<double, 3> z, std::array<double, 3> x) : rest_(rest) {
        if (!(rest.m2() > 0))
          throw std::invalid_argument("eec::frame::BoostRotate: rest-frame vector must be timelike");
        z_ = unit(z);
        const double xz = dot(x, z_);
        for (int k = 0; k < 3; ++k)
          x[k] -= xz * z_[k];
        x_ = unit(x);
        y_ = {{z_[1] * x_[2] - z_[2] * x_[1], z_[2] * x_[0] - z_[0] * x_[2], z_[0] * x_[1] - z_[1] * x_[0]}};
      }

      std::array<double, 2> operator()(const fastjet::PseudoJet& p) const {
        fastjet::PseudoJet q(p);
        q.unboost(rest_);
        const std::array<double, 3> v{{q.px(), q.py(), q.pz()}};
        const double vx = dot(v, x_), vy = dot(v, y_), vz = dot(v, z_);
        return {{std::atan2(std::sqrt(vx * vx + vy * vy), vz), std::atan2(vy, vx)}};
      }

    private:
      static double dot(const std::array<double, 3>& a, const std::array<double, 3>& b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
      }
      static std::array<double, 3> unit(const std::array<double, 3>& a) {
        const double m = std::sqrt(dot(a, a));
        if (!(m > 1e-12))
          throw std::invalid_argument("eec::frame::BoostRotate: degenerate axis");
        return {{a[0] / m, a[1] / m, a[2] / m}};
      }

      fastjet::PseudoJet rest_;
      std::array<double, 3> z_, x_, y_;
    };

  }  // namespace frame

  // One-point correlator sum_i (E_i/E_jet)^n delta(coords_i - coords). Not thread-safe.
  template <class Weight, class Norm, class Frame, class Output>
  class OnePointEC {
  public:
    static constexpr std::size_t D = Frame::D;
    using Result = typename Output::template Result<D>;

    OnePointEC(double n, Binning binning)
      requires(Output::binned)
        : n_(n), binning_(std::move(binning)) {
      if (binning_.rank() != D)
        throw std::invalid_argument("eec::OnePointEC: binning rank must equal Frame::D");
    }

    explicit OnePointEC(double n)
      requires(!Output::binned)
        : n_(n) {}

    Result operator()(const std::vector<fastjet::PseudoJet>& particles,
                      const Frame& frame,
                      double externalNorm = kNoExternalNorm) const {
      fillJetData<Weight, Norm>(jd_, particles, n_, externalNorm);
      if constexpr (Output::binned)
        acc_.begin(binning_, particles.size());
      else
        acc_.begin(particles.size());
      for (std::size_t i = 0; i < particles.size(); ++i) {
        const std::array<double, D> c = frame(particles[i]);
        if constexpr (Output::binned)
          acc_.fill(binning_.linearIndex(c), jd_.w[i], 1);
        else
          acc_.fill(c, jd_.w[i], 1);
      }
      return acc_.finish(JetMeta{jd_.E_jet, n_, 1});
    }

  private:
    double n_;
    Binning binning_;
    mutable JetData jd_;
    mutable typename Output::template Accumulator<D> acc_;
  };

}  // namespace eec

#endif
