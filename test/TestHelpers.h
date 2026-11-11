#ifndef SRothman_EECs_test_TestHelpers_h
#define SRothman_EECs_test_TestHelpers_h

#include "SRothman/EECs/interface/ProjectedENC.h"

#include <catch.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <typeinfo>
#include <vector>

namespace eectest {

  // ---- projections of every output type onto a dense, flow-inclusive histogram ----
  inline std::vector<double> toDense(const eec::DenseResult& r, const eec::Binning&) { return r.values; }
  inline std::vector<double> toDense(const eec::SparseResult& r, const eec::Binning& b) {
    std::vector<double> out(b.totalBins(), 0.0);
    for (std::size_t i = 0; i < r.entries.size(); ++i) {
      if (i > 0)
        REQUIRE(r.entries[i - 1].bin < r.entries[i].bin);
      out[r.entries[i].bin] += r.entries[i].weight;
    }
    return out;
  }
  inline std::vector<double> toDense(const eec::SparseUnmergedResult& r, const eec::Binning& b) {
    std::vector<double> out(b.totalBins(), 0.0);
    for (const auto& e : r.entries)
      out[e.bin] += e.weight * e.symfactor;
    return out;
  }
  template <std::size_t D>
  std::vector<double> toDense(const eec::UnbinnedResult<D>& r, const eec::Binning& b) {
    std::vector<double> out(b.totalBins(), 0.0);
    for (const auto& e : r.entries)
      out[b.linearIndex(e.coord)] += e.weight * e.symfactor;
    return out;
  }

  inline double sum(const std::vector<double>& v) { return std::accumulate(v.begin(), v.end(), 0.0); }

  inline void requireClose(const std::vector<double>& a, const std::vector<double>& b) {
    REQUIRE(a.size() == b.size());
    double scale = 0;
    for (double x : b)
      scale = std::max(scale, std::abs(x));
    for (std::size_t i = 0; i < a.size(); ++i) {
      INFO("bin " << i);
      REQUIRE(a[i] == Approx(b[i]).margin(1e-12 * std::max(scale, 1.0)));
    }
  }

  // Result of a fresh ProjectedENC on ps; b is used only by binned outputs.
  template <std::size_t N, class W, class Nrm, class A, class Out, class Strat>
  auto runPENC(const std::vector<fastjet::PseudoJet>& ps, double n, const eec::Binning& b) {
    using Calc = eec::ProjectedENC<N, W, Nrm, A, Out, eec::NoExtra, Strat>;
    if constexpr (Out::binned)
      return Calc(n, b)(ps);
    else
      return Calc(n)(ps);
  }

  // Calls f(toDense(result)) for every output x loop x precompute combination of ProjectedENC<N, W, Nrm, A>.
  template <std::size_t N, class W, class Nrm, class A, class F>
  void forEachConfig(const std::vector<fastjet::PseudoJet>& ps, double n, const eec::Binning& b, const F& f) {
    namespace st = eec::strategy;
    int calls = 0;
    auto counted = [&](const std::vector<double>& dense) {
      ++calls;
      f(dense);
    };
    auto one = [&]<class Out, class Strat>() {
      INFO("N=" << N << " output=" << typeid(Out).name() << " strategy=" << typeid(Strat).name());
      counted(toDense(runPENC<N, W, Nrm, A, Out, Strat>(ps, n, b), b));
    };
    auto loop = [&]<class Out, class Loop>() {
      one.template operator()<Out, st::Strategy<Loop, st::precompute::None>>();
      one.template operator()<Out, st::Strategy<Loop, st::precompute::PairAngles>>();
      if constexpr (Out::binned)
        one.template operator()<Out, st::Strategy<Loop, st::precompute::PairBins>>();
    };
    auto output = [&]<class Out>() {
      loop.template operator()<Out, st::loop::Ordered>();
      loop.template operator()<Out, st::loop::Multiset>();
    };
    output.template operator()<eec::output::Dense>();
    output.template operator()<eec::output::SparseMerged>();
    output.template operator()<eec::output::SparseMergedSort>();
    output.template operator()<eec::output::SparseUnmerged>();
    output.template operator()<eec::output::Unbinned>();
    // 4 binned outputs x 2 loops x 3 precomputes, plus Unbinned x 2 loops x 2 precomputes.
    REQUIRE(calls == 28);
  }

}  // namespace eectest

#endif
