#include "SRothman/EECs/interface/EECs.h"

#include <benchmark/benchmark.h>

#include <random>

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

  Binning thetaBinning() { return Binning({LogAxis(50, 1e-3, 1.0)}); }

  template <int N, class Out, class Strat>
  void BM_pENC(benchmark::State& state) {
    const auto ps = syntheticJet(static_cast<std::size_t>(state.range(0)), 1);
    using Calc = ProjectedENC<N, weight::Energy, norm::SumEnergy, angle::DeltaR, Out, NoExtra, Strat>;
    auto make = [] {
      if constexpr (Out::binned)
        return Calc(1.0, thetaBinning());
      else
        return Calc(1.0);
    };
    const Calc calc = make();
    for (auto _ : state) {
      auto r = calc(ps);
      benchmark::DoNotOptimize(r);
    }
    state.counters["tuples"] = static_cast<double>(detail::multichoose(ps.size(), N));
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

BENCHMARK(BM_pENC<2, output::Dense, strategy::Naive>)->Apply(Args2);
BENCHMARK(BM_pENC<3, output::Dense, strategy::Naive>)->Apply(Args3);

BENCHMARK(BM_pENC<3, output::SparseMerged, strategy::PairBins>)->Apply(Args3);
BENCHMARK(BM_pENC<3, output::SparseMergedSort, strategy::PairBins>)->Apply(Args3);
BENCHMARK(BM_pENC<3, output::SparseUnmerged, strategy::PairBins>)->Apply(Args3);
BENCHMARK(BM_pENC<3, output::Unbinned, strategy::PairAngles>)->Apply(Args3);

BENCHMARK_MAIN();
