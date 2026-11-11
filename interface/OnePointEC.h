#ifndef SRothman_EECs_OnePointEC_h
#define SRothman_EECs_OnePointEC_h

#include "SRothman/EECs/interface/Binning.h"
#include "SRothman/EECs/interface/CalculatorBase.h"
#include "SRothman/EECs/interface/Concepts.h"
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

      BoostRotate(const fastjet::PseudoJet& rest, std::array<double, 3> z, std::array<double, 3> x) {
        if (!(rest.m2() > 0))
          throw std::invalid_argument("eec::frame::BoostRotate: rest-frame vector must be timelike");
        z_ = unit(z);
        const double xz = dot(x, z_);
        for (int k = 0; k < 3; ++k)
          x[k] -= xz * z_[k];
        x_ = unit(x);
        y_ = {{z_[1] * x_[2] - z_[2] * x_[1], z_[2] * x_[0] - z_[0] * x_[2], z_[0] * x_[1] - z_[1] * x_[0]}};
        restP_ = {{rest.px(), rest.py(), rest.pz()}};
        restE_ = rest.E();
        restM_ = rest.m();
        restProj_ = {{dot(restP_, x_), dot(restP_, y_), dot(restP_, z_)}};
      }

      std::array<double, 2> operator()(const fastjet::PseudoJet& p) const {
        // Same Lorentz transformation as fastjet::PseudoJet::unboost, without copying the PseudoJet.
        const std::array<double, 3> v{{p.px(), p.py(), p.pz()}};
        const double restFrameE = (p.E() * restE_ - dot(v, restP_)) / restM_;
        const double f = (restFrameE + p.E()) / (restE_ + restM_);
        const double vx = dot(v, x_) - f * restProj_[0];
        const double vy = dot(v, y_) - f * restProj_[1];
        const double vz = dot(v, z_) - f * restProj_[2];
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

      std::array<double, 3> z_, x_, y_;
      std::array<double, 3> restP_, restProj_;  // rest's three-momentum, and its components along x_, y_, z_
      double restE_, restM_;
    };

  }  // namespace frame

  // One-point correlator sum_i (E_i/E_jet)^n delta(coords_i - coords). Not thread-safe.
  template <WeightPolicy Weight, NormPolicy Norm, FramePolicy Frame, class Output>
    requires OutputPolicy<Output, Frame::D>
  class OnePointEC : private detail::CalculatorBase<Weight, Norm, Output, Frame::D> {
    using Base = detail::CalculatorBase<Weight, Norm, Output, Frame::D>;

  public:
    static constexpr std::size_t D = Frame::D;
    using typename Base::Result;

    OnePointEC(double n, Binning binning)
      requires(Output::binned)
        : Base("eec::OnePointEC", n, std::move(binning)) {}

    explicit OnePointEC(double n)
      requires(!Output::binned)
        : Base("eec::OnePointEC", n) {}

    Result operator()(const std::vector<fastjet::PseudoJet>& particles,
                      const Frame& frame,
                      double externalNorm = kNoExternalNorm) const {
      this->beginJet(particles, externalNorm, particles.size());
      for (std::size_t i = 0; i < particles.size(); ++i)
        this->fill(frame(particles[i]), this->jd_.w[i], 1);
      return this->finishJet(1);
    }
  };

}  // namespace eec

#endif
