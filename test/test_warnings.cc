#include "SRothman/EECs/interface/EECs.h"
#include "TestJets.h"

#include <catch.hpp>

#include <iostream>
#include <sstream>
#include <string>

using namespace eec;
using fastjet::PseudoJet;

namespace {

  struct CerrCapture {
    std::ostringstream buf;
    std::streambuf* old = std::cerr.rdbuf(buf.rdbuf());
    ~CerrCapture() { std::cerr.rdbuf(old); }
    std::size_t count(const std::string& needle) const {
      const std::string s = buf.str();
      std::size_t n = 0;
      for (auto pos = s.find(needle); pos != std::string::npos; pos = s.find(needle, pos + 1))
        ++n;
      return n;
    }
  };

  const std::string kMsg = "externalNorm was passed to a calculator whose Norm policy ignores it";

  struct NormWithoutTrait {
    static double value(const std::vector<PseudoJet>&, double ext) { return ext; }
  };

}  // namespace

TEST_CASE("ignoresExternalNorm trait", "[warning]") {
  STATIC_REQUIRE(detail::ignoresExternalNorm<norm::SumEnergy>());
  STATIC_REQUIRE(detail::ignoresExternalNorm<norm::SumPt>());
  STATIC_REQUIRE(detail::ignoresExternalNorm<norm::Unity>());
  STATIC_REQUIRE(!detail::ignoresExternalNorm<norm::External>());
  STATIC_REQUIRE(!detail::ignoresExternalNorm<NormWithoutTrait>());
}

TEST_CASE("warnIfExternalNormIgnored sets its flag only for norms that ignore externalNorm", "[warning]") {
  bool warned = false;
  detail::warnIfExternalNormIgnored<norm::SumEnergy>(warned, kNoExternalNorm);
  CHECK(!warned);
  detail::warnIfExternalNormIgnored<norm::External>(warned, 5.0);
  CHECK(!warned);
  detail::warnIfExternalNormIgnored<NormWithoutTrait>(warned, 5.0);
  CHECK(!warned);
  detail::warnIfExternalNormIgnored<norm::SumEnergy>(warned, 5.0);
  CHECK(warned);
}

#ifndef CMSSW_GIT_HASH
TEST_CASE("calculators warn once each about an ignored externalNorm", "[warning]") {
  const auto ps = eectest::randomJet(6, 2);
  const Binning b({RegularAxis(10, 0.0, 1.2)});
  const frame::BoostRotate f(PseudoJet(0, 0, 0, 91.0), {{0, 0, 1}}, {{1, 0, 0}});
  CerrCapture cap;

  ProjectedENC<2, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Dense> a(1.0, b), a2(1.0, b);
  a(ps);
  CHECK(cap.count(kMsg) == 0);
  a(ps, 50.0);
  a(ps, 50.0);
  CHECK(cap.count(kMsg) == 1);
  a2(ps, 50.0);
  CHECK(cap.count(kMsg) == 2);

  ProjectedENC<2, weight::Energy, norm::External, angle::DeltaR, output::Dense> ext(1.0, b);
  ext(ps, 50.0);
  CHECK(cap.count(kMsg) == 2);

  OnePointEC<weight::Energy, norm::Unity, frame::BoostRotate, output::Unbinned> e1c(1.0);
  e1c(ps, f, 50.0);
  e1c(ps, f, 50.0);
  CHECK(cap.count(kMsg) == 3);

  ProjectedENCRuntimeN<weight::Energy, norm::SumPt, angle::DeltaR, output::Dense> rt(3, 1.0, b);
  rt(ps, 50.0);
  rt(ps, 50.0);
  CHECK(cap.count(kMsg) == 4);
}
#endif
