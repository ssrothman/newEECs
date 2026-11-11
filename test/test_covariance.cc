#include "SRothman/EECs/interface/Covariance.h"

#include <catch.hpp>

#include <random>
#include <vector>

using eec::CovarianceAccumulator;

namespace {

  struct Event {
    std::vector<std::pair<std::uint32_t, double>> adds;  // may repeat bins
    double w;
  };

  std::vector<Event> randomEvents(std::size_t nBins, std::size_t nEvents, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<std::uint32_t> bin(0, nBins - 1);
    std::uniform_int_distribution<int> nAdds(0, 12);
    std::uniform_real_distribution<double> val(-0.5, 2.0), wt(-1.5, 3.0);
    std::vector<Event> out(nEvents);
    for (auto& e : out) {
      const int k = nAdds(rng);
      for (int a = 0; a < k; ++a)
        e.adds.push_back({bin(rng), val(rng)});
      e.w = wt(rng);
    }
    return out;
  }

  // Brute force: dense per-event histogram, full outer product.
  void bruteForce(const std::vector<Event>& evs, std::size_t n, std::vector<double>& H, std::vector<double>& C) {
    H.assign(n, 0.0);
    C.assign(n * n, 0.0);
    for (const auto& e : evs) {
      std::vector<double> h(n, 0.0);
      for (const auto& [b, v] : e.adds)
        h[b] += v;
      for (std::size_t i = 0; i < n; ++i) {
        H[i] += e.w * h[i];
        for (std::size_t j = 0; j < n; ++j)
          C[i * n + j] += e.w * e.w * h[i] * h[j];
      }
    }
  }

  void feed(CovarianceAccumulator& acc, const std::vector<Event>& evs) {
    for (const auto& e : evs) {
      acc.beginEvent();
      for (const auto& [b, v] : e.adds)
        acc.add(b, v);
      acc.endEvent(e.w);
    }
  }

}  // namespace

TEST_CASE("CovarianceAccumulator matches a dense brute-force computation", "[covariance]") {
  for (std::size_t n : {1u, 2u, 7u, 40u}) {
    const auto evs = randomEvents(n, 300, 11 + n);
    std::vector<double> H, C;
    bruteForce(evs, n, H, C);
    CovarianceAccumulator acc(n, 1000);
    feed(acc, evs);
    REQUIRE(acc.values().size() == n);
    for (std::size_t i = 0; i < n; ++i)
      CHECK(acc.values()[i] == Approx(H[i]).margin(1e-12));
    const auto dense = acc.denseCovariance();
    REQUIRE(dense.size() == n * n);
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j) {
        CHECK(dense[i * n + j] == Approx(C[i * n + j]).margin(1e-12));
        CHECK(acc.covariance(i, j) == Approx(C[i * n + j]).margin(1e-12));
      }
  }
}

TEST_CASE("CovarianceAccumulator hand-computed two-event example", "[covariance]") {
  CovarianceAccumulator acc(3, 10);
  // Event 1, w = 2: h = (1, 0, 3) with bin 2 filled twice (1 + 2).
  acc.beginEvent();
  acc.add(0, 1.0);
  acc.add(2, 1.0);
  acc.add(2, 2.0);
  acc.endEvent(2.0);
  // Event 2, w = -1: empty.
  acc.beginEvent();
  acc.endEvent(-1.0);
  CHECK(acc.values() == std::vector<double>{2.0, 0.0, 6.0});
  // C = 4 * h h^T; bin 2 is (1+2)^2 * 4 = 36, not (1 + 4) * 4.
  CHECK(acc.covariance(0, 0) == 4.0);
  CHECK(acc.covariance(0, 2) == 12.0);
  CHECK(acc.covariance(2, 0) == 12.0);
  CHECK(acc.covariance(2, 2) == 36.0);
  CHECK(acc.covariance(1, 1) == 0.0);
  CHECK(acc.covariance(0, 1) == 0.0);
  CHECK(acc.sumW() == 1.0);
  CHECK(acc.sumW2() == 5.0);
  CHECK(acc.nEvents() == 2);
}

TEST_CASE("CovarianceAccumulator merge equals a single pass", "[covariance]") {
  const std::size_t n = 9;
  const auto evs = randomEvents(n, 200, 5);
  CovarianceAccumulator all(n, 100), a(n, 100), b(n, 100);
  feed(all, evs);
  feed(a, std::vector<Event>(evs.begin(), evs.begin() + 77));
  feed(b, std::vector<Event>(evs.begin() + 77, evs.end()));
  a.merge(b);
  for (std::size_t i = 0; i < n; ++i) {
    CHECK(a.values()[i] == Approx(all.values()[i]).margin(1e-12));
    for (std::size_t j = 0; j < n; ++j)
      CHECK(a.covariance(i, j) == Approx(all.covariance(i, j)).margin(1e-12));
  }
  CHECK(a.sumW() == Approx(all.sumW()));
  CHECK(a.sumW2() == Approx(all.sumW2()));
  CHECK(a.nEvents() == all.nEvents());
  CHECK_THROWS_AS(a.merge(CovarianceAccumulator(n + 1, 100)), std::invalid_argument);
}

TEST_CASE("CovarianceAccumulator discards an event abandoned before endEvent", "[covariance]") {
  CovarianceAccumulator acc(4, 10);
  acc.beginEvent();
  acc.add(1, 5.0);
  acc.add(3, 2.0);
  // No endEvent: the next beginEvent must clear these, and bin 1 must still be fillable.
  acc.beginEvent();
  acc.add(0, 1.0);
  acc.add(1, 2.0);
  acc.endEvent(1.0);
  CHECK(acc.values() == std::vector<double>{1.0, 2.0, 0.0, 0.0});
  CHECK(acc.covariance(1, 1) == 4.0);
  CHECK(acc.covariance(0, 1) == 2.0);
  CHECK(acc.covariance(3, 3) == 0.0);
  CHECK(acc.covariance(0, 0) == 1.0);
  CHECK(acc.nEvents() == 1);
}

TEST_CASE("CovarianceAccumulator argument checks", "[covariance]") {
  CHECK_THROWS_AS(CovarianceAccumulator(0, 10), std::invalid_argument);
  CHECK_THROWS_AS(CovarianceAccumulator(11, 10), std::invalid_argument);
  CHECK_NOTHROW(CovarianceAccumulator(10, 10));
  CovarianceAccumulator acc(3, 10);
  acc.beginEvent();
  CHECK_THROWS_AS(acc.add(3, 1.0), std::out_of_range);
}
