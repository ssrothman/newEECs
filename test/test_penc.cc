#include "SRothman/EECs/interface/ProjectedENC.h"
#include "SRothman/EECs/interface/ProjectedENCRuntime.h"
#include "TestJets.h"

#include <catch.hpp>

#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>

using namespace eec;
using fastjet::PseudoJet;

namespace {

  // ---- independent reference: ordered tuples, angles from fastjet directly ----
  double refDeltaR(const PseudoJet& a, const PseudoJet& b) {
    return std::hypot(a.pseudorapidity() - b.pseudorapidity(), a.delta_phi_to(b));
  }
  double ref3D(const PseudoJet& a, const PseudoJet& b) {
    const double c = (a.px() * b.px() + a.py() * b.py() + a.pz() * b.pz()) / (a.modp() * b.modp());
    return std::acos(std::clamp(c, -1.0, 1.0));
  }
  // Most tests use Energy weights, SumEnergy norm and DeltaR angles.
  template <int N, class Out, class Extra = NoExtra, class Strat = strategy::DefaultFor<Out>>
  using EDR = ProjectedENC<N, weight::Energy, norm::SumEnergy, angle::DeltaR, Out, Extra, Strat>;
  template <class Out, class Extra = NoExtra>
  using EDRRuntimeN = ProjectedENCRuntimeN<weight::Energy, norm::SumEnergy, angle::DeltaR, Out, Extra>;

  using RefAngle = std::function<double(const PseudoJet&, const PseudoJet&)>;
  using RefExtra = std::function<double(const std::vector<std::size_t>&)>;

  template <class W, class Nrm>
  std::vector<double> oracle(const std::vector<PseudoJet>& ps,
                             int N,
                             double n,
                             const RefAngle& ang,
                             const Binning& b,
                             double ext = kNoExternalNorm,
                             const RefExtra& extra = nullptr) {
    std::vector<double> out(b.totalBins(), 0.0);
    const std::size_t M = ps.size();
    if (M == 0)
      return out;
    const double E = Nrm::value(ps, ext);
    std::vector<std::size_t> t(N, 0);
    while (true) {
      double w = 1, th = 0;
      for (int a = 0; a < N; ++a) {
        w *= std::pow(W::value(ps[t[a]]) / E, n);
        for (int c = a + 1; c < N; ++c)
          th = std::max(th, ang(ps[t[a]], ps[t[c]]));
      }
      std::uint32_t bin = b.binIndex(0, th);
      if (extra)
        bin += b.stride(1) * b.binIndex(1, extra(t));
      out[bin] += w;
      int k = 0;
      while (k < N && ++t[k] == M)
        t[k++] = 0;
      if (k == N)
        break;
    }
    return out;
  }

  // ---- projections of every output type onto a dense histogram ----
  std::vector<double> toDense(const DenseResult& r, const Binning&) { return r.values; }
  std::vector<double> toDense(const SparseResult& r, const Binning& b) {
    std::vector<double> out(b.totalBins(), 0.0);
    for (std::size_t i = 0; i < r.entries.size(); ++i) {
      if (i > 0)
        REQUIRE(r.entries[i - 1].bin < r.entries[i].bin);
      out[r.entries[i].bin] += r.entries[i].weight;
    }
    return out;
  }
  std::vector<double> toDense(const SparseUnmergedResult& r, const Binning& b) {
    std::vector<double> out(b.totalBins(), 0.0);
    for (const auto& e : r.entries)
      out[e.bin] += e.weight * e.symfactor;
    return out;
  }
  template <std::size_t D>
  std::vector<double> toDense(const UnbinnedResult<D>& r, const Binning& b) {
    std::vector<double> out(b.totalBins(), 0.0);
    for (const auto& e : r.entries)
      out[b.linearIndex(e.coord)] += e.weight * e.symfactor;
    return out;
  }

  void requireClose(const std::vector<double>& a, const std::vector<double>& b) {
    REQUIRE(a.size() == b.size());
    double scale = 0;
    for (double x : b)
      scale = std::max(scale, std::abs(x));
    for (std::size_t i = 0; i < a.size(); ++i) {
      INFO("bin " << i);
      REQUIRE(a[i] == Approx(b[i]).margin(1e-12 * std::max(scale, 1.0)));
    }
  }

  template <int N, class W, class Nrm, class A, class Out, class Strat>
  void compareOne(const std::vector<PseudoJet>& ps, double n, const RefAngle& ang, const Binning& b) {
    INFO("N=" << N);
    const auto ref = oracle<W, Nrm>(ps, N, n, ang, b);
    if constexpr (Out::binned) {
      ProjectedENC<N, W, Nrm, A, Out, NoExtra, Strat> calc(n, b);
      requireClose(toDense(calc(ps), b), ref);
    } else {
      ProjectedENC<N, W, Nrm, A, Out, NoExtra, Strat> calc(n);
      requireClose(toDense(calc(ps), b), ref);
    }
  }

  template <int N, class W, class Nrm, class A, class Out, class Loop>
  void compareLoop(const std::vector<PseudoJet>& ps, double n, const RefAngle& ang, const Binning& b) {
    using strategy::Strategy;
    namespace pre = strategy::precompute;
    compareOne<N, W, Nrm, A, Out, Strategy<Loop, pre::None>>(ps, n, ang, b);
    compareOne<N, W, Nrm, A, Out, Strategy<Loop, pre::PairAngles>>(ps, n, ang, b);
    if constexpr (Out::binned)
      compareOne<N, W, Nrm, A, Out, Strategy<Loop, pre::PairBins>>(ps, n, ang, b);
  }

  template <int N, class W, class Nrm, class A, class Out>
  void compareStrategies(const std::vector<PseudoJet>& ps, double n, const RefAngle& ang, const Binning& b) {
    compareLoop<N, W, Nrm, A, Out, strategy::loop::Ordered>(ps, n, ang, b);
    compareLoop<N, W, Nrm, A, Out, strategy::loop::Multiset>(ps, n, ang, b);
  }

  template <int N, class W, class Nrm, class A>
  void compareAll(const std::vector<PseudoJet>& ps, double n, const RefAngle& ang, const Binning& b) {
    compareStrategies<N, W, Nrm, A, output::Dense>(ps, n, ang, b);
    compareStrategies<N, W, Nrm, A, output::SparseMerged>(ps, n, ang, b);
    compareStrategies<N, W, Nrm, A, output::SparseMergedSort>(ps, n, ang, b);
    compareStrategies<N, W, Nrm, A, output::SparseUnmerged>(ps, n, ang, b);
    compareStrategies<N, W, Nrm, A, output::Unbinned>(ps, n, ang, b);
  }

  double sum(const std::vector<double>& v) { return std::accumulate(v.begin(), v.end(), 0.0); }

  // Smallest pairwise DeltaR in the tuple; permutation-symmetric.
  struct MinPairDR {
    static constexpr std::size_t K = 1;
    template <std::size_t N>
    std::array<double, 1> operator()(const JetData& jd, const std::array<std::uint32_t, N>& idx) const {
      double m = std::numeric_limits<double>::infinity();
      for (std::size_t a = 0; a < N; ++a)
        for (std::size_t c = a + 1; c < N; ++c)
          if (idx[a] != idx[c])
            m = std::min(m, angle::DeltaR::pair(jd, idx[a], idx[c]));
      return {{std::isinf(m) ? 0.0 : m}};
    }
  };

  // Constant extra coordinate, reconfigurable through update().
  struct ConstExtra {
    static constexpr std::size_t K = 1;
    double value = 0;
    void update(double v) { value = v; }
    template <std::size_t N>
    std::array<double, 1> operator()(const JetData&, const std::array<std::uint32_t, N>&) const {
      return {{value}};
    }
  };

  // Throws on its throwAt-th call (0 = never).
  struct ThrowingExtra {
    static constexpr std::size_t K = 1;
    int throwAt = 0;
    mutable int calls = 0;
    template <std::size_t N>
    std::array<double, 1> operator()(const JetData&, const std::array<std::uint32_t, N>&) const {
      if (++calls == throwAt)
        throw std::runtime_error("ThrowingExtra");
      return {{0.5}};
    }
  };

}  // namespace

TEST_CASE("pE2C analytic two-particle jet", "[penc]") {
  // Same eta so E is proportional to pt; DeltaR = 0.205.
  std::vector<PseudoJet> ps{fastjet::PtYPhiM(1, 0.4, 1.0), fastjet::PtYPhiM(3, 0.4, 1.205)};
  Binning b({RegularAxis(30, 0.0, 0.3)});
  EDR<2, output::Dense> calc(1.0, b);
  const auto r = calc(ps);
  const double z1 = 0.25, z2 = 0.75;
  const auto b0 = b.binIndex(0, 0.0), b12 = b.binIndex(0, 0.205);
  CHECK(r.values[b0] == Approx(z1 * z1 + z2 * z2));
  CHECK(r.values[b12] == Approx(2 * z1 * z2));
  CHECK(sum(r.values) == Approx(1.0));
  CHECK(r.meta.N == 2);
  CHECK(r.meta.E_jet == Approx(ps[0].E() + ps[1].E()));
}

TEST_CASE("pE3C analytic three-particle jet", "[penc]") {
  // DeltaR12 = 0.105, DeltaR23 = 0.155, DeltaR13 = 0.26.
  std::vector<PseudoJet> ps{
      fastjet::PtYPhiM(1, 0.0, 0.0), fastjet::PtYPhiM(2, 0.0, 0.105), fastjet::PtYPhiM(3, 0.0, 0.26)};
  Binning b({RegularAxis(30, 0.0, 0.3)});
  const double z1 = 1. / 6, z2 = 2. / 6, z3 = 3. / 6;
  auto check = [&](const std::vector<double>& v) {
    CHECK(v[b.binIndex(0, 0.0)] == Approx(z1 * z1 * z1 + z2 * z2 * z2 + z3 * z3 * z3));
    CHECK(v[b.binIndex(0, 0.105)] == Approx(3 * (z1 * z1 * z2 + z1 * z2 * z2)));
    CHECK(v[b.binIndex(0, 0.155)] == Approx(3 * (z2 * z2 * z3 + z2 * z3 * z3)));
    CHECK(v[b.binIndex(0, 0.26)] == Approx(3 * (z1 * z1 * z3 + z1 * z3 * z3) + 6 * z1 * z2 * z3));
    CHECK(sum(v) == Approx(1.0));
  };
  EDR<3, output::Dense, NoExtra, strategy::PairBins> pb(1.0, b);
  check(pb(ps).values);
  EDR<3, output::Dense, NoExtra, strategy::Multiset> ms(1.0, b);
  check(ms(ps).values);

  EDR<3, output::SparseUnmerged> un(1.0, b);
  const auto r = un(ps);
  REQUIRE(r.entries.size() == 10);  // C(5, 3) multisets for M = 3
  for (const auto& e : r.entries) {
    const auto th = b.unravel(e.bin)[0];
    if (th == b.binIndex(0, 0.0))
      CHECK(e.symfactor == 1);
    else if (th == b.binIndex(0, 0.26) && std::abs(e.weight - z1 * z2 * z3) < 1e-12)
      CHECK(e.symfactor == 6);
    else
      CHECK(e.symfactor == 3);
  }
}

TEST_CASE("pENC matches ordered-tuple oracle", "[penc][oracle]") {
  const auto ps = eectest::randomJet(9, 7);
  Binning lin({RegularAxis(25, 0.0, 1.2)});
  Binning lg({LogAxis(20, 1e-2, 1.5)});
  Binning var({VariableAxis({0.0, 0.05, 0.1, 0.3, 0.6, 1.0})});

  compareAll<2, weight::Energy, norm::SumEnergy, angle::DeltaR>(ps, 1.0, refDeltaR, lin);
  compareAll<3, weight::Energy, norm::SumEnergy, angle::DeltaR>(ps, 1.0, refDeltaR, lg);
  compareAll<4, weight::Pt, norm::SumPt, angle::Angle3D>(ps, 1.5, ref3D, var);
  compareAll<3, weight::ChargePt, norm::SumPt, angle::DeltaR>(ps, 2.0, refDeltaR, lin);
  compareAll<2, weight::ChargeEnergy, norm::Unity, angle::Angle3D>(ps, 1.0, ref3D, lg);
  compareAll<5, weight::Energy, norm::SumEnergy, angle::Angle3D>(ps, 0.5, ref3D, lin);
  compareAll<2, weight::Charge, norm::Unity, angle::DeltaR>(ps, 1.0, refDeltaR, var);
}

TEST_CASE("pENC normalisation and multiset bookkeeping", "[penc]") {
  const auto ps = eectest::randomJet(14, 3);
  Binning b({LogAxis(10, 1e-2, 1.0)});
  // expected = C(M + N - 1, N) multisets for M = 14.
  auto checkN = [&](auto calc, int N, std::size_t expected) {
    INFO("N=" << N);
    const auto r = calc(ps);
    REQUIRE(r.entries.size() == expected);
    double ws = 0, syms = 0;
    for (const auto& e : r.entries) {
      ws += e.weight * e.symfactor;
      syms += e.symfactor;
    }
    CHECK(ws == Approx(1.0));
    CHECK(syms == Approx(std::pow(double(ps.size()), N)));
  };
  checkN(EDR<2, output::SparseUnmerged>(1.0, b), 2, 105);
  checkN(EDR<3, output::SparseUnmerged>(1.0, b), 3, 560);
  checkN(EDR<4, output::SparseUnmerged>(1.0, b), 4, 2380);
  checkN(EDR<5, output::SparseUnmerged>(1.0, b), 5, 8568);
  checkN(EDR<6, output::Unbinned>(1.0), 6, 27132);
}

TEST_CASE("pENC extra binning coordinate", "[penc][extra]") {
  const auto ps = eectest::randomJet(8, 11);
  Binning b({RegularAxis(12, 0.0, 1.2), RegularAxis(6, 0.0, 0.6)});
  auto refExtra = [&](const std::vector<std::size_t>& t) {
    double m = std::numeric_limits<double>::infinity();
    for (std::size_t a = 0; a < t.size(); ++a)
      for (std::size_t c = a + 1; c < t.size(); ++c)
        if (t[a] != t[c])
          m = std::min(m, refDeltaR(ps[t[a]], ps[t[c]]));
    return std::isinf(m) ? 0.0 : m;
  };
  const auto ref = oracle<weight::Energy, norm::SumEnergy>(ps, 3, 1.0, refDeltaR, b, kNoExternalNorm, refExtra);
  EDR<3, output::Dense, MinPairDR> dense(1.0, b);
  requireClose(dense(ps).values, ref);
  EDR<3, output::SparseMerged, MinPairDR, strategy::Multiset> sparse(1.0, b);
  requireClose(toDense(sparse(ps), b), ref);
  EDR<3, output::Unbinned, MinPairDR> unb(1.0);
  const auto u = unb(ps);
  requireClose(toDense(u, b), ref);
  for (const auto& e : u.entries)
    CHECK(e.coord[1] <= e.coord[0]);
  CHECK_THROWS_AS((EDR<3, output::Dense, MinPairDR>(1.0, Binning({RegularAxis(12, 0.0, 1.2)}))), std::invalid_argument);
}

TEST_CASE("pENC E_jet tracking and renormalisation", "[penc]") {
  const auto ps = eectest::randomJet(10, 5);
  Binning b({RegularAxis(20, 0.0, 1.0)});
  ProjectedENC<3, weight::Energy, norm::External, angle::DeltaR, output::Dense> ext(1.5, b);
  EDR<3, output::Dense> sum(1.5, b);
  const auto re = ext(ps, 100.0);
  const auto rs = sum(ps);
  CHECK(re.meta.E_jet == 100.0);
  CHECK(re.meta.n == 1.5);
  const double f = renormalization(re.meta, rs.meta.E_jet);
  std::vector<double> scaled(re.values);
  for (auto& x : scaled)
    x *= f;
  requireClose(scaled, rs.values);
  CHECK_THROWS_AS(ext(ps), std::invalid_argument);
}

TEST_CASE("pENC degenerate jets", "[penc]") {
  Binning b({RegularAxis(10, 0.0, 1.0)});
  EDR<3, output::SparseMerged> calc(1.0, b);
  CHECK(calc({}).entries.empty());
  const auto r = calc({fastjet::PtYPhiM(5, 0.1, 0.2)});
  REQUIRE(r.entries.size() == 1);
  CHECK(r.entries[0].bin == b.binIndex(0, 0.0));
  CHECK(r.entries[0].weight == Approx(1.0));
}

TEST_CASE("pENC rejects a zero-momentum particle", "[penc]") {
  Binning b({RegularAxis(10, 0.0, 1.2)});
  auto ps = eectest::randomJet(6, 8);
  auto bad = ps;
  bad.insert(bad.begin() + 2, PseudoJet(0, 0, 0, 0.14));
  EDR<3, output::Dense> dr(1.0, b);
  ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::Angle3D, output::Dense> a3(1.0, b);
  CHECK_THROWS_AS(dr(bad), std::invalid_argument);
  CHECK_THROWS_AS(a3(bad), std::invalid_argument);
  EDR<3, output::Dense> fresh(1.0, b);
  requireClose(dr(ps).values, fresh(ps).values);
}

TEST_CASE("pENC calculator reuse across jets", "[penc]") {
  Binning b({RegularAxis(15, 0.0, 1.2)});
  EDR<3, output::SparseMerged> reused(1.0, b);
  for (unsigned seed = 0; seed < 4; ++seed) {
    const auto ps = eectest::randomJet(5 + 3 * seed, 100 + seed);
    EDR<3, output::SparseMerged> fresh(1.0, b);
    requireClose(toDense(reused(ps), b), toDense(fresh(ps), b));
  }
}

TEST_CASE("runtime-N dispatcher matches compile-time N", "[penc][runtime]") {
  const auto ps = eectest::randomJet(10, 21);
  Binning b({LogAxis(15, 1e-2, 1.5)});
  EDRRuntimeN<output::Dense> r3(3, 1.0, b), r5(5, 1.0, b);
  EDR<3, output::Dense> c3(1.0, b);
  EDR<5, output::Dense> c5(1.0, b);
  requireClose(r3(ps).values, c3(ps).values);
  requireClose(r5(ps).values, c5(ps).values);
  CHECK(r5(ps).meta.N == 5);
  EDRRuntimeN<output::Unbinned> u4(4, 1.0);
  CHECK(u4(ps).entries.size() == 715);  // C(13, 4) multisets for M = 10
  CHECK_THROWS_AS((EDRRuntimeN<output::Dense>(7, 1.0, b)), std::invalid_argument);
}

TEST_CASE("PairBins rejects a decreasing axis 0", "[penc][strategy]") {
  const auto ps = eectest::randomJet(9, 7);
  for (const Axis& ax : {Axis(RegularAxis(10, 1.2, 0.0)), Axis(LogAxis(10, 1.5, 1e-2))}) {
    const Binning b({ax});
    CHECK_THROWS_AS((EDR<3, output::Dense>(1.0, b)), std::invalid_argument);
    CHECK_THROWS_AS((EDRRuntimeN<output::SparseMerged>(3, 1.0, b)), std::invalid_argument);
    CHECK_THROWS_AS(
        (EDR<3, output::Dense, NoExtra, strategy::Strategy<strategy::loop::Ordered, strategy::precompute::PairBins>>(
            1.0, b)),
        std::invalid_argument);
    // Strategies that bin each tuple's theta_max directly handle a decreasing axis correctly.
    compareOne<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Dense, strategy::PairAngles>(
        ps, 1.0, refDeltaR, b);
    compareOne<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Dense, strategy::Multiset>(
        ps, 1.0, refDeltaR, b);
  }
  // Only axis 0 is constrained.
  CHECK_NOTHROW(
      (EDR<3, output::Dense, MinPairDR>(1.0, Binning({RegularAxis(12, 0.0, 1.2), RegularAxis(6, 0.6, 0.0)}))));
}

TEST_CASE("pENC calculator is reusable after a policy throws mid-jet", "[penc][outputs]") {
  const auto ps = eectest::randomJet(10, 7);
  Binning b({RegularAxis(12, 0.0, 1.2), RegularAxis(4, 0.0, 1.0)});
  using Calc = EDR<3, output::SparseMerged, ThrowingExtra>;
  Calc reused(1.0, b, ThrowingExtra{5});
  CHECK_THROWS_AS(reused(ps), std::runtime_error);
  Calc fresh(1.0, b);
  requireClose(toDense(reused(ps), b), toDense(fresh(ps), b));
}

TEST_CASE("strategy aliases and output-dependent default", "[penc][strategy]") {
  using namespace strategy;
  STATIC_REQUIRE(std::is_same_v<Naive, Strategy<loop::Ordered, precompute::None>>);
  STATIC_REQUIRE(std::is_same_v<Multiset, Strategy<loop::Multiset, precompute::None>>);
  STATIC_REQUIRE(std::is_same_v<PairAngles, Strategy<loop::Multiset, precompute::PairAngles>>);
  STATIC_REQUIRE(std::is_same_v<PairBins, Strategy<loop::Multiset, precompute::PairBins>>);
  STATIC_REQUIRE(std::is_same_v<DefaultFor<output::Dense>, PairBins>);
  STATIC_REQUIRE(std::is_same_v<DefaultFor<output::SparseUnmerged>, PairBins>);
  STATIC_REQUIRE(std::is_same_v<DefaultFor<output::Unbinned>, PairAngles>);
}

TEST_CASE("updateExtra reconfigures Extra between jets", "[penc][extra]") {
  const auto ps = eectest::randomJet(8, 4);
  Binning b({RegularAxis(12, 0.0, 1.2), RegularAxis(4, 0.0, 1.0)});
  auto extraBinWeight = [&](const DenseResult& r, double x) {
    const auto want = b.binIndex(1, x);
    double s = 0;
    for (std::uint32_t lin = 0; lin < r.values.size(); ++lin)
      if (b.unravel(lin)[1] == want)
        s += r.values[lin];
    return s;
  };
  EDR<3, output::Dense, ConstExtra> calc(1.0, b, ConstExtra{0.1});
  CHECK(extraBinWeight(calc(ps), 0.1) == Approx(1.0));
  calc.updateExtra(0.8);
  const auto r = calc(ps);
  CHECK(extraBinWeight(r, 0.8) == Approx(1.0));
  CHECK(extraBinWeight(r, 0.1) == 0);

  EDRRuntimeN<output::Dense, ConstExtra> rt(4, 1.0, b);
  rt.updateExtra(0.6);
  CHECK(extraBinWeight(rt(ps), 0.6) == Approx(1.0));
}
