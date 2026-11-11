#ifndef SRothman_EECs_CalculatorBase_h
#define SRothman_EECs_CalculatorBase_h

#include "SRothman/EECs/interface/Binning.h"
#include "SRothman/EECs/interface/JetData.h"
#include "SRothman/EECs/interface/Outputs.h"
#include "SRothman/EECs/interface/Policies.h"
#include "SRothman/EECs/interface/Warn.h"

#include <fastjet/PseudoJet.hh>

#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace eec::detail {

  // Warns once per flag when externalNorm is passed to a Norm that ignores it.
  template <class Norm>
  void warnIfExternalNormIgnored(bool& warned, double externalNorm) {
    if constexpr (ignoresExternalNorm<Norm>()) {
      if (!std::isnan(externalNorm) && !warned) {
        warned = true;
        warn("externalNorm was passed to a calculator whose Norm policy ignores it");
      }
    }
  }

  // Configuration and per-jet steps shared by the calculators: weights, binning and the output accumulator.
  template <class Weight, class Norm, class Output, std::size_t D>
  class CalculatorBase {
  public:
    using Result = typename Output::template Result<D>;

    const Binning& binning() const { return binning_; }

  protected:
    CalculatorBase(const char* name, double n, Binning binning)
      requires(Output::binned)
        : n_(checkedPower(name, n)), binning_(std::move(binning)) {
      if (binning_.rank() != D)
        throw std::invalid_argument(std::string(name) + ": binning rank must be " + std::to_string(D));
    }

    CalculatorBase(const char* name, double n)
      requires(!Output::binned)
        : n_(checkedPower(name, n)) {}

    // Fills jd_'s weights and starts the output for a jet.
    void beginJet(const std::vector<fastjet::PseudoJet>& ps, double externalNorm, std::size_t expectedFills) const {
      warnIfExternalNormIgnored<Norm>(warnedIgnoredExternalNorm_, externalNorm);
      fillWeights<Weight, Norm>(jd_, ps, n_, externalNorm);
      if constexpr (Output::binned)
        acc_.begin(binning_, expectedFills);
      else
        acc_.begin(expectedFills);
    }

    void fill(const std::array<double, D>& coord, double w, std::uint32_t sym) const {
      if constexpr (Output::binned)
        acc_.fill(binning_.linearIndex(coord), w, sym);
      else
        acc_.fill(coord, w, sym);
    }

    void fillBin(std::uint32_t bin, double w, std::uint32_t sym) const
      requires(Output::binned)
    {
      acc_.fill(bin, w, sym);
    }

    Result finishJet(std::size_t N) const { return acc_.finish(JetMeta{jd_.E_jet, n_, N}); }

    mutable JetData jd_;

  private:
    static double checkedPower(const char* name, double n) {
      if (!std::isfinite(n))
        throw std::invalid_argument(std::string(name) + ": n must be finite");
      return n;
    }

    double n_;
    Binning binning_;
    // Per calculator, so the ignored-externalNorm warning fires once per calculator.
    mutable bool warnedIgnoredExternalNorm_ = false;
    mutable typename Output::template Accumulator<D> acc_;
  };

}  // namespace eec::detail

#endif
