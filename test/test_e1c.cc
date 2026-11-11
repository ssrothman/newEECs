#include "SRothman/EECs/interface/OnePointEC.h"
#include "TestJets.h"

#include <catch.hpp>

#include <cmath>
#include <limits>
#include <numbers>
#include <random>

using namespace eec;
using fastjet::PseudoJet;

using std::numbers::pi;

TEST_CASE("BoostRotate without boost reproduces polar/azimuthal angles", "[e1c]") {
  const PseudoJet atRest(0, 0, 0, 91.0);
  frame::BoostRotate f(atRest, {{0, 0, 1}}, {{1, 0, 0}});
  CHECK(f(PseudoJet(0, 0, 5, 5))[0] == Approx(0).margin(1e-12));
  const auto cx = f(PseudoJet(5, 0, 0, 5));
  CHECK(cx[0] == Approx(pi / 2));
  CHECK(cx[1] == Approx(0).margin(1e-12));
  const auto cy = f(PseudoJet(0, 5, 0, 5));
  CHECK(cy[1] == Approx(pi / 2));
  const auto c = f(PseudoJet(1, -1, -std::sqrt(2.0), 2));
  CHECK(c[0] == Approx(3 * pi / 4));
  CHECK(c[1] == Approx(-pi / 4));
}

TEST_CASE("BoostRotate axes: rotated frame and orthogonalisation", "[e1c]") {
  const PseudoJet atRest(0, 0, 0, 91.0);
  // z along lab x; x given non-orthogonal (lab y plus some lab x), y = z cross x = lab z.
  frame::BoostRotate f(atRest, {{2, 0, 0}}, {{0.7, 3, 0}});
  CHECK(f(PseudoJet(4, 0, 0, 4))[0] == Approx(0).margin(1e-12));
  const auto a = f(PseudoJet(0, 4, 0, 4));
  CHECK(a[0] == Approx(pi / 2));
  CHECK(a[1] == Approx(0).margin(1e-12));
  CHECK(f(PseudoJet(0, 0, 4, 4))[1] == Approx(pi / 2));
  CHECK_THROWS_AS(frame::BoostRotate(atRest, {{1, 0, 0}}, {{2, 0, 0}}), std::invalid_argument);
  CHECK_THROWS_AS(frame::BoostRotate(PseudoJet(0, 0, 5, 5), {{1, 0, 0}}, {{0, 1, 0}}), std::invalid_argument);
}

TEST_CASE("BoostRotate boosts into the rest frame", "[e1c]") {
  // Two-body decay along rest-frame direction (theta0, phi0), then boosted along lab z.
  const double m = 91.0, theta0 = 0.7, phi0 = 0.4, E = m / 2;
  const PseudoJet d1(
      E * std::sin(theta0) * std::cos(phi0), E * std::sin(theta0) * std::sin(phi0), E * std::cos(theta0), E);
  const PseudoJet d2(-d1.px(), -d1.py(), -d1.pz(), E);
  const PseudoJet moving(0, 0, 200, std::sqrt(200 * 200 + m * m));
  PseudoJet l1(d1), l2(d2);
  l1.boost(moving);
  l2.boost(moving);
  REQUIRE(std::abs(l1.theta() - theta0) > 0.1);  // the boost really changes the lab angle

  frame::BoostRotate f(moving, {{0, 0, 1}}, {{1, 0, 0}});
  const auto c1 = f(l1), c2 = f(l2);
  CHECK(c1[0] == Approx(theta0));
  CHECK(c1[1] == Approx(phi0));
  CHECK(c1[0] + c2[0] == Approx(pi));
  CHECK(std::abs(c1[1] - c2[1]) == Approx(pi));
}

TEST_CASE("OnePointEC weights and outputs", "[e1c]") {
  const PseudoJet atRest(0, 0, 0, 91.0);
  frame::BoostRotate f(atRest, {{0, 0, 1}}, {{1, 0, 0}});
  std::vector<PseudoJet> ps{PseudoJet(0, 0, 1, 1), PseudoJet(3, 0, 0, 3), PseudoJet(0, 0, -4, 4)};
  ps[0].set_user_info(new EECParticleInfo(1));
  ps[1].set_user_info(new EECParticleInfo(-1));
  ps[2].set_user_info(new EECParticleInfo(0));
  Binning b({RegularAxis(4, 0.0, pi + 1e-9), RegularAxis(4, -pi, pi + 1e-9)});

  OnePointEC<weight::Energy, norm::SumEnergy, frame::BoostRotate, output::Dense> dense(2.0, b);
  const auto r = dense(ps, f);
  CHECK(r.meta.E_jet == 8);
  CHECK(r.meta.N == 1);
  CHECK(r.values[b.linearIndex(std::array<double, 2>{{0.0, 0.0}})] == Approx(1. / 64));
  CHECK(r.values[b.linearIndex(std::array<double, 2>{{pi / 2, 0.0}})] == Approx(9. / 64));
  double total = 0;
  for (double v : r.values)
    total += v;
  CHECK(total == Approx((1. + 9 + 16) / 64));

  OnePointEC<weight::ChargeEnergy, norm::External, frame::BoostRotate, output::Unbinned> unb(1.0);
  const auto u = unb(ps, f, 10.0);
  REQUIRE(u.entries.size() == 3);
  CHECK(u.entries[1].weight == Approx(-0.3));
  CHECK(u.entries[1].coord[0] == Approx(pi / 2));
  CHECK(u.entries[2].weight == 0);
  CHECK(u.entries[2].coord[0] == Approx(pi));
  CHECK(u.entries[0].symfactor == 1);

  OnePointEC<weight::Energy, norm::SumEnergy, frame::BoostRotate, output::SparseMerged> sp(1.0, b);
  const auto s = sp(ps, f);
  CHECK(s.entries.size() == 3);
}

TEST_CASE("OnePointEC accepts a zero-momentum particle", "[e1c]") {
  // At rest in the lab, it moves along -z in the rest frame of a vector moving along +z.
  const frame::BoostRotate f(PseudoJet(0, 0, 50, std::sqrt(50 * 50 + 91.0 * 91.0)), {{0, 0, 1}}, {{1, 0, 0}});
  OnePointEC<weight::Energy, norm::SumEnergy, frame::BoostRotate, output::Unbinned> e1c(1.0);
  const auto r = e1c({PseudoJet(1, 0, 0, 1), PseudoJet(0, 0, 0, 0.14)}, f);
  REQUIRE(r.entries.size() == 2);
  CHECK(r.entries[1].coord[0] == Approx(pi));
  CHECK(r.entries[1].weight == Approx(0.14 / 1.14));
}

TEST_CASE("BoostRotate matches fastjet unboost for general frames", "[e1c]") {
  std::mt19937 rng(5);
  std::uniform_real_distribution<double> u(-50, 50), um(1, 100);
  auto unit = [](std::array<double, 3> a) {
    const double m = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    return std::array<double, 3>{{a[0] / m, a[1] / m, a[2] / m}};
  };
  auto dot = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  };
  for (int trial = 0; trial < 20; ++trial) {
    const double rx = u(rng), ry = u(rng), rz = u(rng), rm = um(rng);
    const PseudoJet rest(rx, ry, rz, std::sqrt(rx * rx + ry * ry + rz * rz + rm * rm));
    const std::array<double, 3> z{{u(rng), u(rng), u(rng)}}, x{{u(rng), u(rng), u(rng)}};
    const frame::BoostRotate f(rest, z, x);

    const auto zn = unit(z);
    const double xz = dot(x, zn);
    const auto xn = unit({{x[0] - xz * zn[0], x[1] - xz * zn[1], x[2] - xz * zn[2]}});
    const std::array<double, 3> yn{
        {zn[1] * xn[2] - zn[2] * xn[1], zn[2] * xn[0] - zn[0] * xn[2], zn[0] * xn[1] - zn[1] * xn[0]}};
    for (int i = 0; i < 10; ++i) {
      const double px = u(rng), py = u(rng), pz = u(rng), pm = 0.1 * um(rng);
      const PseudoJet p(px, py, pz, std::sqrt(px * px + py * py + pz * pz + pm * pm));
      PseudoJet q(p);
      q.unboost(rest);
      const std::array<double, 3> v{{q.px(), q.py(), q.pz()}};
      const double vx = dot(v, xn), vy = dot(v, yn), vz = dot(v, zn);
      const double theta = std::atan2(std::sqrt(vx * vx + vy * vy), vz), phi = std::atan2(vy, vx);
      const auto c = f(p);
      INFO("trial " << trial << " particle " << i);
      CHECK(c[0] == Approx(theta).margin(1e-9));
      // cos/sin so that phi = -pi and +pi compare equal.
      CHECK(std::cos(c[1]) == Approx(std::cos(phi)).margin(1e-9));
      CHECK(std::sin(c[1]) == Approx(std::sin(phi)).margin(1e-9));
    }
  }
}

TEST_CASE("OnePointEC rejects a non-finite n", "[e1c]") {
  using Dense = OnePointEC<weight::Energy, norm::SumEnergy, frame::BoostRotate, output::Dense>;
  using Unb = OnePointEC<weight::Energy, norm::SumEnergy, frame::BoostRotate, output::Unbinned>;
  const Binning b({RegularAxis(4, 0.0, pi), RegularAxis(4, -pi, pi)});
  for (double bad : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                     -std::numeric_limits<double>::infinity()}) {
    CHECK_THROWS_AS(Dense(bad, b), std::invalid_argument);
    CHECK_THROWS_AS(Unb(bad), std::invalid_argument);
  }
  CHECK_NOTHROW(Dense(-2.0, b));
  CHECK_NOTHROW(Unb(0.0));
}
