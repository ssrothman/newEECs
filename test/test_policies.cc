#include "SRothman/EECs/interface/Angles.h"
#include "SRothman/EECs/interface/JetData.h"
#include "SRothman/EECs/interface/Policies.h"
#include "TestJets.h"

#include <catch.hpp>

#include <cmath>
#include <numbers>

using namespace eec;
using fastjet::PseudoJet;

namespace {
  struct OtherUserInfo : fastjet::PseudoJet::UserInfoBase {};
}  // namespace

TEST_CASE("weight policies", "[policies]") {
  PseudoJet p(3, 4, 12, 14);  // pt 5, E 14
  p.set_user_info(new EECParticleInfo(-1, 211));
  CHECK(weight::Energy::value(p) == 14);
  CHECK(weight::Pt::value(p) == Approx(5));
  CHECK(weight::Charge::value(p) == -1);
  CHECK(weight::ChargeEnergy::value(p) == -14);
  CHECK(weight::ChargePt::value(p) == Approx(-5));
  CHECK(p.user_info<EECParticleInfo>().pdgId() == 211);

  PseudoJet bare(3, 4, 12, 14);
  CHECK_THROWS_AS(weight::Charge::value(bare), std::invalid_argument);
  PseudoJet foreign(3, 4, 12, 14);
  foreign.set_user_info(new OtherUserInfo);
  CHECK_THROWS_AS(weight::Charge::value(foreign), std::invalid_argument);
  CHECK_THROWS_AS(weight::ChargePt::value(bare), std::invalid_argument);
  CHECK_NOTHROW(weight::Energy::value(bare));
}

TEST_CASE("norm policies", "[policies]") {
  std::vector<PseudoJet> ps{PseudoJet(3, 4, 0, 6), PseudoJet(0, 1, 0, 2)};
  CHECK(norm::SumEnergy::value(ps, kNoExternalNorm) == 8);
  CHECK(norm::SumPt::value(ps, kNoExternalNorm) == Approx(6));
  CHECK(norm::Unity::value(ps, 123) == 1);
  CHECK(norm::External::value(ps, 42.5) == 42.5);
  CHECK_THROWS_AS(norm::External::value(ps, kNoExternalNorm), std::invalid_argument);
  CHECK_THROWS_AS(norm::External::value(ps, 0), std::invalid_argument);
}

TEST_CASE("JetData weights and geometry", "[jetdata]") {
  std::vector<PseudoJet> ps{PseudoJet(3, 4, 0, 6), PseudoJet(0, 0, 2, 2)};
  JetData jd;
  fillWeights<weight::Energy, norm::SumEnergy>(jd, ps, 1.0, kNoExternalNorm);
  fillGeometry(jd, ps);
  CHECK(jd.E_jet == 8);
  CHECK(jd.w[0] == Approx(0.75));
  CHECK(jd.w[1] == Approx(0.25));
  CHECK(jd.ux[0] == Approx(0.6));
  CHECK(jd.uy[0] == Approx(0.8));
  CHECK(jd.uz[1] == Approx(1.0));

  fillWeights<weight::Energy, norm::SumEnergy>(jd, ps, 2.0, kNoExternalNorm);
  CHECK(jd.w[0] == Approx(0.5625));
  fillWeights<weight::Energy, norm::External>(jd, ps, 0.5, 2.0);
  CHECK(jd.E_jet == 2.0);
  CHECK(jd.w[0] == Approx(std::sqrt(3.0)));
}

TEST_CASE("JetData rejects negative weight with non-integer power", "[jetdata]") {
  std::vector<PseudoJet> ps{PseudoJet(3, 4, 0, 6)};
  ps[0].set_user_info(new EECParticleInfo(-1));
  JetData jd;
  CHECK_THROWS_AS((fillWeights<weight::Charge, norm::Unity>(jd, ps, 0.5, kNoExternalNorm)), std::domain_error);
  fillWeights<weight::Charge, norm::Unity>(jd, ps, 3.0, kNoExternalNorm);
  CHECK(jd.w[0] == -1);
}

TEST_CASE("Angle3D", "[angles]") {
  std::vector<PseudoJet> ps{PseudoJet(1, 0, 0, 1), PseudoJet(0, 1, 0, 1), PseudoJet(-1, 0, 0, 1),
                            PseudoJet(std::cos(1e-7), std::sin(1e-7), 0, 1)};
  JetData jd;
  fillGeometry(jd, ps);
  CHECK(angle::Angle3D::pair(jd, 0, 1) == Approx(std::numbers::pi / 2));
  CHECK(angle::Angle3D::pair(jd, 0, 2) == Approx(std::numbers::pi));
  CHECK(angle::Angle3D::pair(jd, 1, 1) == 0);
  // acos(dot) would lose this entirely in double precision.
  CHECK(angle::Angle3D::pair(jd, 0, 3) == Approx(1e-7).epsilon(1e-6));
}

TEST_CASE("DeltaR with phi wrap", "[angles]") {
  std::vector<PseudoJet> ps{fastjet::PtYPhiM(1, 0.5, 0.05), fastjet::PtYPhiM(1, 0.2, 2 * std::numbers::pi - 0.05),
                            fastjet::PtYPhiM(1, 0.5, 1.05)};
  JetData jd;
  fillGeometry(jd, ps);
  CHECK(angle::DeltaR::pair(jd, 0, 1) == Approx(std::hypot(0.1, 0.3)));
  CHECK(angle::DeltaR::pair(jd, 1, 0) == Approx(std::hypot(0.1, 0.3)));
  CHECK(angle::DeltaR::pair(jd, 0, 2) == Approx(1.0));
}

TEST_CASE("fillWeights leaves the geometry empty", "[jetdata]") {
  const auto ps = eectest::randomJet(5, 3);
  JetData jd;
  fillWeights<weight::Energy, norm::SumEnergy>(jd, ps, 1.0, kNoExternalNorm);
  CHECK(jd.size() == 5);
  CHECK(jd.ux.empty());
  CHECK(jd.eta.empty());
  CHECK(jd.phi.empty());
}

TEST_CASE("fillGeometry rejects a zero-momentum particle", "[jetdata]") {
  std::vector<PseudoJet> ps{PseudoJet(3, 4, 0, 6), PseudoJet(0, 0, 0, 0.14)};
  JetData jd;
  CHECK_NOTHROW(fillWeights<weight::Energy, norm::SumEnergy>(jd, ps, 1.0, kNoExternalNorm));
  CHECK_THROWS_AS(fillGeometry(jd, ps), std::invalid_argument);
}
