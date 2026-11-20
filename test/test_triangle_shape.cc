#include "SRothman/EECs/interface/ProjectedENC.h"
#include "SRothman/EECs/interface/TriangleShape.h"
#include "TestHelpers.h"
#include "TestJets.h"

#include <catch.hpp>

#include <cmath>
#include <numbers>

using namespace eec;
using fastjet::PseudoJet;
using eectest::requireClose;
using eectest::toDense;

TEST_CASE("triangleShape analytic triangles", "[shape]") {
  const double halfPi = std::numbers::pi / 2;
  auto check = [](std::array<double, 2> got, double xi, double phi) {
    REQUIRE(got[0] == Approx(xi).margin(1e-12));
    REQUIRE(got[1] == Approx(phi).margin(1e-7));
  };
  check(triangleShape(1, 1, 1), 1.0, halfPi);                                 // equilateral
  check(triangleShape(3, 4, 5), 0.75, std::asin(std::sqrt(8.0 / 9.0)));       // generic
  check(triangleShape(5, 3, 4), 0.75, std::asin(std::sqrt(8.0 / 9.0)));       // any side order
  check(triangleShape(1, 2, 3), 0.5, 0.0);                                    // collinear
  check(triangleShape(1, 1, 0.01), 0.01, halfPi);                             // squeezed, isosceles
  check(triangleShape(2, 2.5, 1), 0.5, std::asin(std::sqrt(1 - 0.25)));      // R_L - R_M = R_S / 2
  check(triangleShape(0, 1, 1), 0.0, 0.0);                                    // coincident pair
  check(triangleShape(0, 0, 0), 0.0, 0.0);
  check(triangleShape(1, 2, 3 + 1e-15), 0.5, 0.0);                            // rounding past collinear
}

namespace {
  double refDeltaR(const PseudoJet& a, const PseudoJet& b) {
    return std::hypot(a.pseudorapidity() - b.pseudorapidity(), a.delta_phi_to(b));
  }

  // All ordered triples with fastjet angles: weight prod (E_i / sum E), filled at (R_L, xi, phi).
  std::vector<double> oracle(const std::vector<PseudoJet>& ps, const Binning& b) {
    std::vector<double> out(b.totalBins(), 0.0);
    double E = 0;
    for (const auto& p : ps)
      E += p.E();
    for (const auto& p : ps)
      for (const auto& q : ps)
        for (const auto& r : ps) {
          const double a = refDeltaR(p, q), c = refDeltaR(p, r), d = refDeltaR(q, r);
          const auto s = triangleShape(a, c, d);
          out[b.linearIndex(std::array<double, 3>{{std::max({a, c, d}), s[0], s[1]}})] +=
              p.E() / E * q.E() / E * r.E() / E;
        }
    return out;
  }

  template <class Out, class Strat>
  std::vector<double> resolved(const std::vector<PseudoJet>& ps, const Binning& b) {
    using Calc =
        ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, Out, TriangleShape<angle::DeltaR>, Strat>;
    if constexpr (Out::binned)
      return toDense(Calc(1.0, b)(ps), b);
    else
      return toDense(Calc(1.0)(ps), b);
  }
}  // namespace

TEST_CASE("resolved E3C matches the ordered-triple oracle", "[shape][oracle]") {
  const Binning b({LogAxis(12, 1e-3, 1.0), RegularAxis(5, 0.0, 1.0), RegularAxis(5, 0.0, std::numbers::pi / 2)});
  for (unsigned seed : {1u, 2u, 3u}) {
    const auto ps = eectest::randomJet(9, seed);
    const auto ref = oracle(ps, b);
    requireClose(resolved<output::Dense, strategy::PairBins>(ps, b), ref);
    requireClose(resolved<output::Dense, strategy::PairAngles>(ps, b), ref);
    requireClose(resolved<output::SparseMerged, strategy::PairBins>(ps, b), ref);
    requireClose(resolved<output::Unbinned, strategy::PairAngles>(ps, b), ref);
  }
}

TEST_CASE("resolved E3C summed over the shape is the projected E3C", "[shape]") {
  const Binning b3({LogAxis(12, 1e-3, 1.0), RegularAxis(5, 0.0, 1.0), RegularAxis(5, 0.0, std::numbers::pi / 2)});
  const Binning b1({LogAxis(12, 1e-3, 1.0)});
  const auto ps = eectest::randomJet(15, 4);
  const auto full = resolved<output::Dense, strategy::PairBins>(ps, b3);
  const auto proj =
      toDense(ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Dense>(1.0, b1)(ps), b1);
  std::vector<double> summed(b1.totalBins(), 0.0);
  for (std::size_t i = 0; i < full.size(); ++i)
    summed[i % b3.extent(0)] += full[i];
  requireClose(summed, proj);
}
