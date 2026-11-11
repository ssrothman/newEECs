#include "SRothman/EECs/interface/OnePointEC.h"
#include "SRothman/EECs/interface/ProjectedENC.h"
#include "TestHelpers.h"
#include "TestJets.h"

#include <catch.hpp>

#include <cmath>
#include <numbers>

using namespace eec;
using fastjet::PseudoJet;

// Sum rule: with n = 1 and E_jet = sum_i E_p_i, sum over ordered N-tuples of prod_i E_p_i / E_jet
// is (sum_i E_p_i / E_jet)^N = 1, for every N, output, strategy and weight definition.

namespace {

  constexpr double kTol = 1e-10;  // observed |sum-1| <= 2.2e-13; real defects shift the sum by O(1e-2)

  // E_jet = sum of the same per-particle quantity used as E_p.
  template <class W>
  struct SumOf {
    static double value(const std::vector<PseudoJet>& ps, double) {
      double s = 0;
      for (const auto& p : ps)
        s += W::value(p);
      return s;
    }
  };

  template <class R>
  double total(const R& r, const Binning& b) {
    return eectest::sum(eectest::toDense(r, b));
  }

  template <std::size_t N, class W, class Nrm>
  void checkAllConfigs(const std::vector<PseudoJet>& ps, const Binning& b) {
    eectest::forEachConfig<N, W, Nrm, angle::DeltaR>(ps, 1.0, b, [](const std::vector<double>& v) {
      const double s = eectest::sum(v);
      INFO("sum-1=" << s - 1);
      CHECK(std::abs(s - 1.0) < kTol);
    });
  }

  template <class W, class Nrm>
  void checkAllN(const std::vector<PseudoJet>& ps, const Binning& b) {
    checkAllConfigs<2, W, Nrm>(ps, b);
    checkAllConfigs<3, W, Nrm>(ps, b);
    checkAllConfigs<4, W, Nrm>(ps, b);
    checkAllConfigs<5, W, Nrm>(ps, b);
    checkAllConfigs<6, W, Nrm>(ps, b);
  }

}  // namespace

TEST_CASE("pENC weights sum to 1 for n=1 and E_jet = sum E_p", "[penc][sumrule]") {
  const auto ps = eectest::randomJet(10, 42);
  // Narrow range so many tuples land in flow bins, which must still be counted.
  Binning b({LogAxis(8, 0.05, 0.4)});

  checkAllN<weight::Energy, norm::SumEnergy>(ps, b);
  checkAllN<weight::Pt, norm::SumPt>(ps, b);
  checkAllN<weight::Energy, SumOf<weight::Energy>>(ps, b);

  // Signed weights: the rule holds whenever sum q_i E_i != 0.
  REQUIRE(std::abs(SumOf<weight::Charge>::value(ps, 0)) > 0.5);
  checkAllN<weight::Charge, SumOf<weight::Charge>>(ps, b);
  checkAllN<weight::ChargeEnergy, SumOf<weight::ChargeEnergy>>(ps, b);
  checkAllN<weight::ChargePt, SumOf<weight::ChargePt>>(ps, b);
}

TEST_CASE("pENC sum rule fails when n != 1 or E_jet is not the sum", "[penc][sumrule]") {
  const auto ps = eectest::randomJet(10, 42);
  Binning b({LogAxis(8, 0.05, 0.4)});
  using P2 = ProjectedENC<2, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Dense>;
  using P2pt = ProjectedENC<2, weight::Energy, norm::SumPt, angle::DeltaR, output::Dense>;
  CHECK(std::abs(total(P2(2.0, b)(ps), b) - 1.0) > 1e-3);
  CHECK(std::abs(total(P2pt(1.0, b)(ps), b) - 1.0) > 1e-3);
}

TEST_CASE("E1C weights sum to 1 for n=1 and E_jet = sum E_p", "[e1c][sumrule]") {
  const auto ps = eectest::randomJet(25, 9);
  const PseudoJet rest(0, 0, 50, std::sqrt(50 * 50 + 91.0 * 91.0));
  frame::BoostRotate f(rest, {{0, 0, 1}}, {{1, 0, 0}});
  Binning b({RegularAxis(6, 0.2, 1.0), RegularAxis(5, -1.0, 1.0)});
  CHECK(std::abs(total(OnePointEC<weight::Energy, norm::SumEnergy, frame::BoostRotate, output::Dense>(1.0, b)(ps, f), b) -
                 1.0) < kTol);
  CHECK(std::abs(total(OnePointEC<weight::Pt, norm::SumPt, frame::BoostRotate, output::SparseMerged>(1.0, b)(ps, f), b) -
                 1.0) < kTol);
  CHECK(std::abs(total(OnePointEC<weight::ChargePt, SumOf<weight::ChargePt>, frame::BoostRotate, output::Unbinned>(
                           1.0)(ps, f),
                       b) -
                 1.0) < kTol);
}
