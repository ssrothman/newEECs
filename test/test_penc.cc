#include "SRothman/EECs/interface/ProjectedENC.h"
#include "SRothman/EECs/interface/ProjectedENCRuntime.h"
#include "TestJets.h"

#include <catch.hpp>

#include <cmath>
#include <functional>
#include <limits>
#include <numeric>

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
  using RefAngle = std::function<double(const PseudoJet&, const PseudoJet&)>;
  using RefExtra = std::function<double(const std::vector<std::size_t>&)>;

  template <class W, class Nrm>
  std::vector<double> oracle(const std::vector<PseudoJet>& ps, int N, double n, const RefAngle& ang, const Binning& b,
                             double ext = kNoExternalNorm, const RefExtra& extra = nullptr) {
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

  template <int N, class W, class Nrm, class A, class Out>
  void compareStrategies(const std::vector<PseudoJet>& ps, double n, const RefAngle& ang, const Binning& b) {
    compareOne<N, W, Nrm, A, Out, strategy::Naive>(ps, n, ang, b);
    compareOne<N, W, Nrm, A, Out, strategy::Multiset>(ps, n, ang, b);
    compareOne<N, W, Nrm, A, Out, strategy::PairAngles>(ps, n, ang, b);
    compareOne<N, W, Nrm, A, Out, strategy::PairBins>(ps, n, ang, b);
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

  using EnergyDR2 = ProjectedENC<2, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Dense>;

}  // namespace

TEST_CASE("pE2C analytic two-particle jet", "[penc]") {
  // Same eta so E is proportional to pt; DeltaR = 0.205.
  std::vector<PseudoJet> ps{fastjet::PtYPhiM(1, 0.4, 1.0), fastjet::PtYPhiM(3, 0.4, 1.205)};
  Binning b({RegularAxis(30, 0.0, 0.3)});
  EnergyDR2 calc(1.0, b);
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
  std::vector<PseudoJet> ps{fastjet::PtYPhiM(1, 0.0, 0.0), fastjet::PtYPhiM(2, 0.0, 0.105),
                            fastjet::PtYPhiM(3, 0.0, 0.26)};
  Binning b({RegularAxis(30, 0.0, 0.3)});
  const double z1 = 1. / 6, z2 = 2. / 6, z3 = 3. / 6;
  auto check = [&](const std::vector<double>& v) {
    CHECK(v[b.binIndex(0, 0.0)] == Approx(z1 * z1 * z1 + z2 * z2 * z2 + z3 * z3 * z3));
    CHECK(v[b.binIndex(0, 0.105)] == Approx(3 * (z1 * z1 * z2 + z1 * z2 * z2)));
    CHECK(v[b.binIndex(0, 0.155)] == Approx(3 * (z2 * z2 * z3 + z2 * z3 * z3)));
    CHECK(v[b.binIndex(0, 0.26)] == Approx(3 * (z1 * z1 * z3 + z1 * z3 * z3) + 6 * z1 * z2 * z3));
    CHECK(sum(v) == Approx(1.0));
  };
  ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Dense, NoExtra, strategy::PairBins> pb(1.0,
                                                                                                                  b);
  check(pb(ps).values);
  ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Dense, NoExtra, strategy::Multiset> ms(1.0,
                                                                                                                  b);
  check(ms(ps).values);

  ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::SparseUnmerged> un(1.0, b);
  const auto r = un(ps);
  REQUIRE(r.entries.size() == 10);  // multichoose(3, 3)
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
  auto checkN = [&](auto calc, int N) {
    INFO("N=" << N);
    const auto r = calc(ps);
    REQUIRE(r.entries.size() == detail::multichoose(ps.size(), N));
    double ws = 0, syms = 0;
    for (const auto& e : r.entries) {
      ws += e.weight * e.symfactor;
      syms += e.symfactor;
    }
    CHECK(ws == Approx(1.0));
    CHECK(syms == Approx(std::pow(double(ps.size()), N)));
  };
  checkN(ProjectedENC<2, weight::Energy, norm::SumEnergy, angle::DeltaR, output::SparseUnmerged>(1.0, b), 2);
  checkN(ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::SparseUnmerged>(1.0, b), 3);
  checkN(ProjectedENC<4, weight::Energy, norm::SumEnergy, angle::DeltaR, output::SparseUnmerged>(1.0, b), 4);
  checkN(ProjectedENC<5, weight::Energy, norm::SumEnergy, angle::DeltaR, output::SparseUnmerged>(1.0, b), 5);
  checkN(ProjectedENC<6, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Unbinned>(1.0), 6);
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
  ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Dense, MinPairDR> dense(1.0, b);
  requireClose(dense(ps).values, ref);
  ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::SparseMerged, MinPairDR, strategy::Multiset>
      sparse(1.0, b);
  requireClose(toDense(sparse(ps), b), ref);
  ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Unbinned, MinPairDR> unb(1.0);
  const auto u = unb(ps);
  requireClose(toDense(u, b), ref);
  for (const auto& e : u.entries)
    CHECK(e.coord[1] <= e.coord[0]);
  CHECK_THROWS_AS((ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Dense, MinPairDR>(
                      1.0, Binning({RegularAxis(12, 0.0, 1.2)}))),
                  std::invalid_argument);
}

TEST_CASE("pENC E_jet tracking and renormalisation", "[penc]") {
  const auto ps = eectest::randomJet(10, 5);
  Binning b({RegularAxis(20, 0.0, 1.0)});
  ProjectedENC<3, weight::Energy, norm::External, angle::DeltaR, output::Dense> ext(1.5, b);
  ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Dense> sum(1.5, b);
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
  ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::SparseMerged> calc(1.0, b);
  CHECK(calc({}).entries.empty());
  const auto r = calc({fastjet::PtYPhiM(5, 0.1, 0.2)});
  REQUIRE(r.entries.size() == 1);
  CHECK(r.entries[0].bin == b.binIndex(0, 0.0));
  CHECK(r.entries[0].weight == Approx(1.0));
}

TEST_CASE("pENC calculator reuse across jets", "[penc]") {
  Binning b({RegularAxis(15, 0.0, 1.2)});
  ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::SparseMerged> reused(1.0, b);
  for (unsigned seed = 0; seed < 4; ++seed) {
    const auto ps = eectest::randomJet(5 + 3 * seed, 100 + seed);
    ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::SparseMerged> fresh(1.0, b);
    requireClose(toDense(reused(ps), b), toDense(fresh(ps), b));
  }
}

TEST_CASE("runtime-N dispatcher matches compile-time N", "[penc][runtime]") {
  const auto ps = eectest::randomJet(10, 21);
  Binning b({LogAxis(15, 1e-2, 1.5)});
  ProjectedENCRuntimeN<weight::Energy, norm::SumEnergy, angle::DeltaR, output::Dense> r3(3, 1.0, b), r5(5, 1.0, b);
  ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Dense> c3(1.0, b);
  ProjectedENC<5, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Dense> c5(1.0, b);
  requireClose(r3(ps).values, c3(ps).values);
  requireClose(r5(ps).values, c5(ps).values);
  CHECK(r5(ps).meta.N == 5);
  ProjectedENCRuntimeN<weight::Energy, norm::SumEnergy, angle::DeltaR, output::Unbinned> u4(4, 1.0);
  CHECK(u4(ps).entries.size() == detail::multichoose(ps.size(), 4));
  CHECK_THROWS_AS((ProjectedENCRuntimeN<weight::Energy, norm::SumEnergy, angle::DeltaR, output::Dense>(7, 1.0, b)),
                  std::invalid_argument);
}
