#include "SRothman/EECs/interface/Binning.h"
#include "SRothman/EECs/interface/Outputs.h"

#include <catch.hpp>

#include <random>

using namespace eec;

namespace {
  Binning twoD() { return Binning({RegularAxis(5, 0.0, 1.0), VariableAxis({-1.0, 0.0, 2.0, 5.0})}); }
}  // namespace

TEST_CASE("Binning flow-inclusive indices", "[binning]") {
  Binning b({RegularAxis(10, 0.0, 1.0)});
  CHECK(b.totalBins() == 12);
  CHECK(b.binIndex(0, -0.1) == 0);
  CHECK(b.binIndex(0, 0.0) == 1);
  CHECK(b.binIndex(0, 0.95) == 10);
  CHECK(b.binIndex(0, 1.0) == 11);

  Binning lg({LogAxis(4, 1e-3, 1.0)});
  CHECK(lg.binIndex(0, 0.0) == 0);  // log(0) underflows
  CHECK(lg.binIndex(0, 2e-3) == 1);
}

TEST_CASE("Binning linear index matches boost::histogram", "[binning]") {
  Binning b = twoD();
  CHECK(b.totalBins() == 7 * 5);
  CHECK(b.stride(0) == 1);
  CHECK(b.stride(1) == 7);
  auto ref = bh::make_histogram_with(bh::dense_storage<double>(), std::vector<Axis>{b.axis(0), b.axis(1)});
  std::vector<double> dense(b.totalBins(), 0.0);
  std::mt19937 rng(1);
  std::uniform_real_distribution<double> ux(-0.3, 1.3), uy(-2, 6), uw(0, 1);
  for (int i = 0; i < 2000; ++i) {
    const std::array<double, 2> x{{ux(rng), uy(rng)}};
    const double w = uw(rng);
    ref(x[0], x[1], bh::weight(w));
    dense[b.linearIndex(x)] += w;
    const auto idx = b.unravel(b.linearIndex(x));
    CHECK(idx[0] == b.binIndex(0, x[0]));
    CHECK(idx[1] == b.binIndex(1, x[1]));
    CHECK(b.linearIndex(b.binIndex(0, x[0]), std::array<double, 1>{{x[1]}}) == b.linearIndex(x));
  }
  const auto h = b.toHistogram(dense);
  for (auto&& x : bh::indexed(ref, bh::coverage::all))
    CHECK(h.at(x.index(0), x.index(1)) == Approx(*x));
}

TEST_CASE("Dense and sparse accumulators", "[outputs]") {
  Binning b({RegularAxis(4, 0.0, 1.0)});
  const JetMeta meta{10.0, 1.0, 2};
  auto run = [&](auto acc) {
    acc.begin(b, 4);
    acc.fill(3, 0.5, 2);
    acc.fill(1, 0.25, 1);
    acc.fill(3, 0.1, 1);
    acc.fill(2, 1.0, 1);
    acc.fill(2, -1.0, 1);  // cancels to an exact zero but the bin was filled
    return acc.finish(meta);
  };

  const auto dense = run(output::Dense::Accumulator<1>{});
  REQUIRE(dense.values.size() == 6);
  CHECK(dense.values[3] == Approx(1.1));
  CHECK(dense.values[1] == Approx(0.25));
  CHECK(dense.values[0] == 0);
  CHECK(dense.meta.E_jet == 10.0);

  for (const auto& sp : {run(output::SparseMerged::Accumulator<1>{}), run(output::SparseMergedSort::Accumulator<1>{})}) {
    REQUIRE(sp.entries.size() == 3);
    CHECK(sp.entries[0].bin == 1);
    CHECK(sp.entries[0].weight == Approx(0.25));
    CHECK(sp.entries[1].bin == 2);
    CHECK(sp.entries[1].weight == 0);
    CHECK(sp.entries[2].bin == 3);
    CHECK(sp.entries[2].weight == Approx(1.1));
    CHECK(sp.meta.N == 2);
  }

  const auto un = run(output::SparseUnmerged::Accumulator<1>{});
  REQUIRE(un.entries.size() == 5);
  CHECK(un.entries[0].bin == 3);
  CHECK(un.entries[0].weight == 0.5);
  CHECK(un.entries[0].symfactor == 2);
}

TEST_CASE("SparseMerged scratch is reset between jets", "[outputs]") {
  Binning b({RegularAxis(4, 0.0, 1.0)});
  output::SparseMerged::Accumulator<1> acc;
  acc.begin(b, 0);
  acc.fill(2, 1.0, 1);
  acc.finish({});
  acc.begin(b, 0);
  acc.fill(1, 0.5, 1);
  const auto r = acc.finish({});
  REQUIRE(r.entries.size() == 1);
  CHECK(r.entries[0].bin == 1);
  CHECK(r.entries[0].weight == 0.5);
}

TEST_CASE("SparseMerged recovers from a jet that never reached finish()", "[outputs]") {
  Binning b({RegularAxis(4, 0.0, 1.0)});
  output::SparseMerged::Accumulator<1> acc;
  acc.begin(b, 0);
  acc.fill(2, 1.0, 1);
  acc.fill(3, 7.0, 1);  // abandoned without finish(), as when a fill throws
  acc.begin(b, 0);
  acc.fill(2, 0.5, 1);
  const auto r = acc.finish({});
  REQUIRE(r.entries.size() == 1);
  CHECK(r.entries[0].bin == 2);
  CHECK(r.entries[0].weight == 0.5);
  acc.begin(b, 0);
  acc.fill(3, 0.25, 1);
  const auto r2 = acc.finish({});
  REQUIRE(r2.entries.size() == 1);
  CHECK(r2.entries[0].bin == 3);
  CHECK(r2.entries[0].weight == 0.25);
}

TEST_CASE("renormalization", "[outputs]") {
  const JetMeta m{50.0, 2.0, 3};
  CHECK(renormalization(m, 100.0) == Approx(std::pow(0.5, 6)));
}
