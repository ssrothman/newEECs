#include "SRothman/EECs/interface/Covariance.h"

#include <benchmark/benchmark.h>

#include <random>
#include <vector>

namespace {

  // Args: total bins, distinct bins touched per event.
  void BM_CovarianceEvent(benchmark::State& state) {
    const std::size_t nBins = state.range(0);
    const std::size_t k = state.range(1);
    std::mt19937 rng(3);
    std::uniform_int_distribution<std::uint32_t> bin(0, nBins - 1);
    std::vector<std::vector<std::uint32_t>> events(64);
    for (auto& e : events)
      for (std::size_t a = 0; a < k; ++a)
        e.push_back(bin(rng));
    eec::CovarianceAccumulator acc(nBins, nBins);
    std::size_t i = 0;
    for (auto _ : state) {
      const auto& e = events[i++ % events.size()];
      acc.beginEvent();
      for (auto b : e)
        acc.add(b, 0.5);
      acc.endEvent(1.0);
    }
    const double* out = acc.values().data();
    benchmark::DoNotOptimize(out);
  }

}  // namespace

BENCHMARK(BM_CovarianceEvent)->Args({52, 52})->Args({2704, 50})->Args({2704, 300})->Args({2704, 2704});

BENCHMARK_MAIN();
