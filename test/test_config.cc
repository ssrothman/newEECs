#include "SRothman/EECs/interface/Config.h"

#include <catch.hpp>

#include <cmath>

using namespace eec;
using namespace eec::config;

namespace {

  constexpr double kPi = 3.14159265358979323846;

  // A particle with rest-frame direction (theta, phi) relative to rest-frame axes (x, y, z), boosted into the lab.
  fastjet::PseudoJet fromRestFrame(double theta,
                                   double phi,
                                   const std::array<double, 3>& x,
                                   const std::array<double, 3>& y,
                                   const std::array<double, 3>& z,
                                   const fastjet::PseudoJet& jet) {
    const double E = 7.0;
    double v[3];
    for (int k = 0; k < 3; ++k)
      v[k] = E * (std::sin(theta) * std::cos(phi) * x[k] + std::sin(theta) * std::sin(phi) * y[k] + std::cos(theta) * z[k]);
    fastjet::PseudoJet p(v[0], v[1], v[2], E);
    p.boost(jet);
    return p;
  }

}  // namespace

TEST_CASE("parseNorm and normValue", "[config]") {
  std::vector<fastjet::PseudoJet> ps{fastjet::PseudoJet(3, 4, 0, 6), fastjet::PseudoJet(0, 1, 2, 3)};
  const fastjet::PseudoJet jet(30, 40, 10, 100);
  CHECK(normValue(parseNorm("sumE"), ps, jet) == 9.0);
  CHECK(normValue(parseNorm("sumPt"), ps, jet) == Approx(6.0));
  CHECK(normValue(parseNorm("jetE"), ps, jet) == 100.0);
  CHECK(normValue(parseNorm("jetPt"), ps, jet) == Approx(50.0));
  CHECK(normValue(parseNorm("unity"), ps, jet) == 1.0);
  CHECK(normValue(NormSource::SumE, {}, jet) == 0.0);
  CHECK_THROWS_AS(parseNorm("SumE"), std::invalid_argument);
}

TEST_CASE("makeAxis, axisType and axisEdges", "[config]") {
  const Axis r = makeAxis({"regular", 4, 0.0, 2.0, {}});
  CHECK(axisType(r) == AxisType::Regular);
  CHECK(axisEdges(r) == std::vector<double>{0.0, 0.5, 1.0, 1.5, 2.0});

  const Axis l = makeAxis({"log", 2, 1e-2, 1.0, {}});
  CHECK(axisType(l) == AxisType::Log);
  const auto le = axisEdges(l);
  REQUIRE(le.size() == 3);
  CHECK(le[0] == Approx(1e-2));
  CHECK(le[1] == Approx(1e-1));
  CHECK(le[2] == Approx(1.0));

  const Axis v = makeAxis({"variable", 0, 0, 0, {0.0, 0.1, 0.5}});
  CHECK(axisType(v) == AxisType::Variable);
  CHECK(axisEdges(v) == std::vector<double>{0.0, 0.1, 0.5});

  CHECK_THROWS_AS(makeAxis({"regular", 0, 0.0, 1.0, {}}), std::invalid_argument);
  CHECK_THROWS_AS(makeAxis({"logarithmic", 3, 0.1, 1.0, {}}), std::invalid_argument);

  const Binning b = makeBinning({{"log", 10, 1e-3, 1.0, {}}, {"regular", 5, 0.0, 1.0, {}}});
  CHECK(b.rank() == 2);
  CHECK(b.totalBins() == 12 * 7);
}

TEST_CASE("Lab frame gives lab polar angle and azimuth", "[config][frame]") {
  const auto f = makeFrame(parseFrame("lab"), fastjet::PseudoJet(100, 0, 0, 200));
  for (const auto& p : {fastjet::PseudoJet(1, 2, 3, 5), fastjet::PseudoJet(-3, -0.5, -2, 4)}) {
    const auto c = f(p);
    CHECK(c[0] == Approx(p.theta()));
    CHECK(c[1] == Approx(p.phi_std()));
  }
}

TEST_CASE("Frames for a transverse jet along +x have the expected axes", "[config][frame]") {
  // Jet with pz = 0 moving along +x. Helicity: z = +x_lab, x = +z_lab, y = z cross x = -y_lab.
  // Collins-Soper: beams in the jet rest frame are (-gb, 0, +-1); z = +z_lab, x = +x_lab, y = +y_lab.
  const fastjet::PseudoJet jet(300, 0, 0, std::sqrt(300.0 * 300 + 90.0 * 90));
  const std::array<double, 3> ex{{1, 0, 0}}, ey{{0, 1, 0}}, ez{{0, 0, 1}}, mey{{0, -1, 0}};
  const auto hel = makeFrame(parseFrame("jetHelicity"), jet);
  const auto cs = makeFrame(parseFrame("jetCS"), jet);
  for (double theta : {0.3, 1.2, 2.5})
    for (double phi : {-2.8, -0.4, 1.1, 2.9}) {
      const auto ch = hel(fromRestFrame(theta, phi, ez, mey, ex, jet));
      CHECK(ch[0] == Approx(theta).margin(1e-9));
      CHECK(ch[1] == Approx(phi).margin(1e-9));
      const auto cc = cs(fromRestFrame(theta, phi, ex, ey, ez, jet));
      CHECK(cc[0] == Approx(theta).margin(1e-9));
      CHECK(cc[1] == Approx(phi).margin(1e-9));
    }
}

TEST_CASE("Collins-Soper z follows the sign of the jet pz", "[config][frame]") {
  // A forward and a mirrored backward jet: a particle mirrored in z must get the same CS theta.
  const fastjet::PseudoJet fwd(200, 50, 400, 500), bwd(200, 50, -400, 500);
  const auto csF = makeFrame(FrameKind::JetCS, fwd);
  const auto csB = makeFrame(FrameKind::JetCS, bwd);
  const fastjet::PseudoJet p(60, 10, 140, 160), pMirror(60, 10, -140, 160);
  CHECK(csF(p)[0] == Approx(csB(pMirror)[0]).margin(1e-9));
  // Without the sign flip the mirrored theta would be pi - theta; make sure that differs here.
  CHECK(std::abs(csF(p)[0] - (kPi - csF(p)[0])) > 0.1);
}

TEST_CASE("Helicity frame: a particle along the jet has theta 0", "[config][frame]") {
  const fastjet::PseudoJet jet(120, -80, 250, 330);
  const auto hel = makeFrame(FrameKind::JetHelicity, jet);
  fastjet::PseudoJet along = jet * 0.25;
  along.reset_momentum(along.px(), along.py(), along.pz(), along.modp());  // massless, same direction
  CHECK(hel(along)[0] == Approx(0.0).margin(1e-6));
  CHECK_THROWS_AS(parseFrame("helicity"), std::invalid_argument);
  CHECK_THROWS_AS(makeFrame(FrameKind::JetHelicity, fastjet::PseudoJet(1, 0, 0, 1)), std::invalid_argument);
}
