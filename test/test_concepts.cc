#include "SRothman/EECs/interface/EECs.h"

#include <catch.hpp>

#include <algorithm>
#include <cmath>

using namespace eec;
using fastjet::PseudoJet;

namespace {

  // ---- examples copied from docs/guide.md section 9 ----
  struct ChargedEnergy {
    static double value(const PseudoJet& p) { return particleCharge(p) != 0 ? p.E() : 0.0; }
  };
  struct VectorSumPt {
    static constexpr bool usesExternalNorm = false;
    static double value(const std::vector<PseudoJet>& ps, double) {
      PseudoJet s(0, 0, 0, 0);
      for (const auto& p : ps)
        s += p;
      return s.pt();
    }
  };
  struct DeltaEta {
    static double pair(const JetData& d, std::size_t i, std::size_t j) { return std::abs(d.eta[i] - d.eta[j]); }
  };
  struct EventScale {
    static constexpr std::size_t K = 1;
    double scale = 0;
    void update(double s) { scale = s; }
    template <std::size_t N>
    std::array<double, 1> operator()(const JetData&, const std::array<std::uint32_t, N>&) const {
      return {{scale}};
    }
  };
  struct LabPolar {
    static constexpr std::size_t D = 2;
    std::array<double, 2> operator()(const PseudoJet& p) const { return {{p.theta(), p.phi_std()}}; }
  };
  struct TotalWeight {
    static constexpr bool binned = false;
    template <std::size_t D>
    struct Result {
      JetMeta meta;
      double total = 0;
    };
    template <std::size_t D>
    class Accumulator {
    public:
      void begin(std::size_t) { r_ = {}; }
      void fill(const std::array<double, D>&, double w, std::uint32_t sym) { r_.total += w * sym; }
      Result<D> finish(const JetMeta& m) {
        r_.meta = m;
        return r_;
      }

    private:
      Result<D> r_;
    };
  };

  // ---- policies that each break one requirement ----
  struct WeightNonStatic {
    double value(const PseudoJet& p) const { return p.E(); }
  };
  struct WeightWrongArg {
    static double value(const std::vector<PseudoJet>&) { return 0; }
  };
  struct NormNoExternalArg {
    static double value(const std::vector<PseudoJet>&) { return 1; }
  };
  struct NormRuntimeFlag {
    static inline bool usesExternalNorm = false;
    static double value(const std::vector<PseudoJet>&, double) { return 1; }
  };
  struct AngleOnPseudoJets {
    static double pair(const PseudoJet&, const PseudoJet&) { return 0; }
  };
  struct ExtraWrongSize {
    static constexpr std::size_t K = 2;
    template <std::size_t N>
    std::array<double, 1> operator()(const JetData&, const std::array<std::uint32_t, N>&) const {
      return {{0.0}};
    }
  };
  struct ExtraNoK {
    template <std::size_t N>
    std::array<double, 1> operator()(const JetData&, const std::array<std::uint32_t, N>&) const {
      return {{0.0}};
    }
  };
  struct ExtraOnlyPairs {
    static constexpr std::size_t K = 1;
    std::array<double, 1> operator()(const JetData&, const std::array<std::uint32_t, 2>&) const { return {{0.0}}; }
  };
  struct FrameWrongSize {
    static constexpr std::size_t D = 3;
    std::array<double, 2> operator()(const PseudoJet&) const { return {}; }
  };
  struct FrameNonConstCall {
    static constexpr std::size_t D = 2;
    std::array<double, 2> operator()(const PseudoJet&) { return {}; }
  };
  // Binned flag set, but begin/fill have the unbinned signatures.
  struct OutputBinnedMismatch : TotalWeight {
    static constexpr bool binned = true;
  };
  struct OutputNoFinish {
    static constexpr bool binned = false;
    template <std::size_t D>
    using Result = double;
    template <std::size_t D>
    struct Accumulator {
      void begin(std::size_t) {}
      void fill(const std::array<double, D>&, double, std::uint32_t) {}
    };
  };

}  // namespace

TEST_CASE("library and guide policies satisfy their concepts", "[concepts]") {
  STATIC_REQUIRE(WeightPolicy<weight::Energy>);
  STATIC_REQUIRE(WeightPolicy<weight::Pt>);
  STATIC_REQUIRE(WeightPolicy<weight::Charge>);
  STATIC_REQUIRE(WeightPolicy<weight::ChargeEnergy>);
  STATIC_REQUIRE(WeightPolicy<weight::ChargePt>);
  STATIC_REQUIRE(WeightPolicy<ChargedEnergy>);

  STATIC_REQUIRE(NormPolicy<norm::SumEnergy>);
  STATIC_REQUIRE(NormPolicy<norm::SumPt>);
  STATIC_REQUIRE(NormPolicy<norm::External>);
  STATIC_REQUIRE(NormPolicy<norm::Unity>);
  STATIC_REQUIRE(NormPolicy<VectorSumPt>);

  STATIC_REQUIRE(AnglePolicy<angle::Angle3D>);
  STATIC_REQUIRE(AnglePolicy<angle::DeltaR>);
  STATIC_REQUIRE(AnglePolicy<DeltaEta>);

  STATIC_REQUIRE(ExtraPolicy<NoExtra, 1>);
  STATIC_REQUIRE(ExtraPolicy<NoExtra, 6>);
  STATIC_REQUIRE(ExtraPolicy<EventScale, 3>);

  STATIC_REQUIRE(FramePolicy<frame::BoostRotate>);
  STATIC_REQUIRE(FramePolicy<LabPolar>);

  STATIC_REQUIRE(OutputPolicy<output::Dense, 1>);
  STATIC_REQUIRE(OutputPolicy<output::SparseMerged, 2>);
  STATIC_REQUIRE(OutputPolicy<output::SparseMergedSort, 1>);
  STATIC_REQUIRE(OutputPolicy<output::SparseUnmerged, 3>);
  STATIC_REQUIRE(OutputPolicy<output::Unbinned, 2>);
  STATIC_REQUIRE(OutputPolicy<TotalWeight, 1>);
}

TEST_CASE("concepts reject policies that break one requirement", "[concepts]") {
  STATIC_REQUIRE(!WeightPolicy<WeightNonStatic>);
  STATIC_REQUIRE(!WeightPolicy<WeightWrongArg>);
  STATIC_REQUIRE(!NormPolicy<NormNoExternalArg>);
  STATIC_REQUIRE(!NormPolicy<NormRuntimeFlag>);
  STATIC_REQUIRE(!AnglePolicy<AngleOnPseudoJets>);
  STATIC_REQUIRE(!ExtraPolicy<ExtraWrongSize, 3>);
  STATIC_REQUIRE(!ExtraPolicy<ExtraNoK, 3>);
  STATIC_REQUIRE(ExtraPolicy<ExtraOnlyPairs, 2>);
  STATIC_REQUIRE(!ExtraPolicy<ExtraOnlyPairs, 3>);
  STATIC_REQUIRE(!FramePolicy<FrameWrongSize>);
  STATIC_REQUIRE(!FramePolicy<FrameNonConstCall>);
  STATIC_REQUIRE(!OutputPolicy<OutputBinnedMismatch, 1>);
  STATIC_REQUIRE(!OutputPolicy<OutputNoFinish, 1>);
}
