#include "SRothman/EECs/interface/EECs.h"

#include <benchmark/benchmark.h>

#include <algorithm>
#include <array>
#include <algorithm>
#include <array>
#include <random>
#include <vector>
#include <vector>

using namespace eec;

namespace {

  std::vector<fastjet::PseudoJet> syntheticJet(std::size_t M, unsigned seed) {
    std::mt19937 rng(seed);
    std::normal_distribution<double> d(0.0, 0.15);
    std::exponential_distribution<double> pt(0.1);
    std::vector<fastjet::PseudoJet> out;
    for (std::size_t i = 0; i < M; ++i)
      out.push_back(fastjet::PtYPhiM(0.5 + pt(rng), 0.5 + d(rng), 1.0 + d(rng), 0.14));
    return out;
  }

  // Constant extra coordinate: isolates the cost of binning an Extra.
  struct ConstExtra {
    static constexpr std::size_t K = 1;
    template <std::size_t N>
    std::array<double, 1> operator()(const JetData&, const std::array<std::uint32_t, N>&) const {
      return {{0.5}};
    }
  };

  // Smallest DeltaR between distinct particles of the tuple: an Extra that depends on the indices.
  struct MinDeltaR {
    static constexpr std::size_t K = 1;
    template <std::size_t N>
    std::array<double, 1> operator()(const JetData& d, const std::array<std::uint32_t, N>& idx) const {
      double m = 0;
      bool any = false;
      for (std::size_t a = 0; a < N; ++a)
        for (std::size_t b = a + 1; b < N; ++b)
          if (idx[a] != idx[b]) {
            const double r = angle::DeltaR::pair(d, idx[a], idx[b]);
            m = any ? std::min(m, r) : r;
            any = true;
          }
      return {{m}};
    }
  };

  Binning thetaBinning(std::size_t K) {
    std::vector<Axis> axes{LogAxis(50, 1e-3, 1.0)};
    for (std::size_t k = 0; k < K; ++k)
      axes.push_back(RegularAxis(10, 0.0, 1.0));
    return Binning(std::move(axes));
  }

  constexpr std::size_t kJets = 32;

  template <std::size_t N, class Out, class Strat, class Extra = NoExtra>
  void BM_pENC(benchmark::State& state) {
    // Cycle through distinct jets: repeating one jet lets the branch predictor memorise data-dependent branches.
    std::vector<std::vector<fastjet::PseudoJet>> jets;
    for (unsigned k = 0; k < kJets; ++k)
      jets.push_back(syntheticJet(static_cast<std::size_t>(state.range(0)), k + 1));
    using Calc = ProjectedENC<N, weight::Energy, norm::SumEnergy, angle::DeltaR, Out, Extra, Strat>;
    auto make = [] {
      if constexpr (Out::binned)
        return Calc(1.0, thetaBinning(Extra::K));
      else
        return Calc(1.0);
    };
    const Calc calc = make();
    std::size_t i = 0;
    for (auto _ : state) {
      auto r = calc(jets[i++ % kJets]);
      benchmark::DoNotOptimize(r);
    }
    state.counters["tuples"] = static_cast<double>(detail::multichoose(jets[0].size(), N));
  }

  void Args2(benchmark::internal::Benchmark* b) {
    for (int M : {20, 50, 100, 150})
      b->Arg(M);
  }
  void Args3(benchmark::internal::Benchmark* b) {
    for (int M : {20, 50, 100})
      b->Arg(M);
  }
  void ArgsHigh(benchmark::internal::Benchmark* b) {
    for (int M : {20, 50})
      b->Arg(M);
  }

}  // namespace

#define EEC_BENCH_STRATS(N, OUT, ARGS)                                         \
  BENCHMARK(BM_pENC<N, OUT, strategy::Multiset>)->Apply(ARGS);                 \
  BENCHMARK(BM_pENC<N, OUT, strategy::PairAngles>)->Apply(ARGS);               \
  BENCHMARK(BM_pENC<N, OUT, strategy::PairBins>)->Apply(ARGS);

EEC_BENCH_STRATS(2, output::Dense, Args2)
EEC_BENCH_STRATS(3, output::Dense, Args3)
EEC_BENCH_STRATS(4, output::Dense, ArgsHigh)
EEC_BENCH_STRATS(5, output::Dense, ArgsHigh)
EEC_BENCH_STRATS(6, output::Dense, ArgsHigh)

using UnpeeledPairBins = strategy::Strategy<strategy::loop::Multiset, strategy::precompute::PairBins, strategy::innermost::Unpeeled>;
BENCHMARK(BM_pENC<3, output::Dense, UnpeeledPairBins>)->Apply(Args3);
BENCHMARK(BM_pENC<4, output::Dense, UnpeeledPairBins>)->Apply(ArgsHigh);
BENCHMARK(BM_pENC<5, output::Dense, UnpeeledPairBins>)->Apply(ArgsHigh);
BENCHMARK(BM_pENC<6, output::Dense, UnpeeledPairBins>)->Apply(ArgsHigh);

BENCHMARK(BM_pENC<2, output::Dense, strategy::Naive>)->Apply(Args2);
BENCHMARK(BM_pENC<3, output::Dense, strategy::Naive>)->Apply(Args3);

BENCHMARK(BM_pENC<3, output::SparseMerged, strategy::PairBins>)->Apply(Args3);
BENCHMARK(BM_pENC<3, output::SparseMergedSort, strategy::PairBins>)->Apply(Args3);
BENCHMARK(BM_pENC<3, output::SparseUnmerged, strategy::PairBins>)->Apply(Args3);
BENCHMARK(BM_pENC<3, output::Unbinned, strategy::PairAngles>)->Apply(Args3);
BENCHMARK(BM_pENC<4, output::SparseUnmerged, strategy::PairBins>)->Apply(ArgsHigh);
BENCHMARK(BM_pENC<4, output::Unbinned, strategy::PairAngles>)->Apply(ArgsHigh);

BENCHMARK(BM_pENC<3, output::Dense, strategy::PairBins, ConstExtra>)->Apply(Args3);
BENCHMARK(BM_pENC<3, output::Dense, strategy::PairBins, MinDeltaR>)->Apply(Args3);
BENCHMARK(BM_pENC<3, output::Dense, strategy::PairAngles, MinDeltaR>)->Apply(Args3);
BENCHMARK(BM_pENC<4, output::Dense, strategy::PairBins, MinDeltaR>)->Apply(ArgsHigh);
BENCHMARK(BM_pENC<3, output::Unbinned, strategy::PairAngles, MinDeltaR>)->Apply(Args3);

BENCHMARK_MAIN();
