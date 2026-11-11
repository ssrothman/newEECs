#ifndef SRothman_EECs_ProjectedENC_h
#define SRothman_EECs_ProjectedENC_h

#include "SRothman/EECs/interface/Angles.h"
#include "SRothman/EECs/interface/Binning.h"
#include "SRothman/EECs/interface/Concepts.h"
#include "SRothman/EECs/interface/JetData.h"
#include "SRothman/EECs/interface/Outputs.h"
#include "SRothman/EECs/interface/Policies.h"

#include <fastjet/PseudoJet.hh>

#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace eec {

  // A strategy is a tuple loop combined with a per-jet precompute; any combination is valid,
  // except that precompute::PairBins needs a binned output.
  namespace strategy {
    namespace loop {
      // All M^N ordered tuples, one fill each with symfactor 1.
      struct Ordered {};
      // Multisets i1 <= ... <= iN, one fill each with its multinomial symfactor.
      struct Multiset {};
    }  // namespace loop

    namespace precompute {
      // Pair angles computed on the fly.
      struct None {};
      // M x M pair-angle matrix.
      struct PairAngles {};
      // M x M matrix of axis-0 bins, using bin(max theta) == max bin(theta) for increasing axes.
      struct PairBins {};
    }  // namespace precompute

    template <class Loop, class Precompute>
    struct Strategy {
      using loop_type = Loop;
      using precompute_type = Precompute;
    };

    using Naive = Strategy<loop::Ordered, precompute::None>;
    using Multiset = Strategy<loop::Multiset, precompute::None>;
    using PairAngles = Strategy<loop::Multiset, precompute::PairAngles>;
    using PairBins = Strategy<loop::Multiset, precompute::PairBins>;

    template <class Output>
    using DefaultFor = std::conditional_t<Output::binned, PairBins, PairAngles>;
  }  // namespace strategy

  // Extra binning coordinates: K values computed from the tuple's particle indices.
  // Must be symmetric under permutations of the tuple (entries are per multiset).
  struct NoExtra {
    static constexpr std::size_t K = 0;
    template <std::size_t N>
    std::array<double, 0> operator()(const JetData&, const std::array<std::uint32_t, N>&) const {
      return {};
    }
  };

  namespace detail {
    constexpr std::uint32_t factorial(int n) { return n <= 1 ? 1u : static_cast<std::uint32_t>(n) * factorial(n - 1); }

    template <class>
    inline constexpr bool isStrategy = false;
    template <class L, class P>
    inline constexpr bool isStrategy<strategy::Strategy<L, P>> = true;

    // Number of multisets of size N drawn from M elements.
    inline std::size_t multichoose(std::size_t M, int N) {
      if (M == 0)
        return 0;
      double r = 1;
      for (int k = 1; k <= N; ++k)
        r = r * static_cast<double>(M + k - 1) / k;
      return static_cast<std::size_t>(r + 0.5);
    }
  }  // namespace detail

  // Projected N-point correlator. Not thread-safe: use one instance per thread.
  template <int N,
            WeightPolicy Weight,
            NormPolicy Norm,
            AnglePolicy Angle,
            class Output,
            class Extra = NoExtra,
            class Strategy = strategy::DefaultFor<Output>>
    requires ExtraPolicy<Extra, N> && OutputPolicy<Output, 1 + Extra::K>
  class ProjectedENC {
    static_assert(N >= 1 && N <= 12, "symmetry factors are computed in 32-bit integers");
    static_assert(detail::isStrategy<Strategy>, "Strategy must be strategy::Strategy<Loop, Precompute>");

  public:
    static constexpr std::size_t K = Extra::K;
    static constexpr std::size_t D = 1 + K;
    using Result = typename Output::template Result<D>;

    ProjectedENC(double n, Binning binning, Extra extra = Extra{})
      requires(Output::binned)
        : n_(n), binning_(std::move(binning)), extra_(std::move(extra)) {
      if (binning_.rank() != D)
        throw std::invalid_argument("eec::ProjectedENC: binning rank must be 1 + Extra::K");
      if constexpr (kPairBins) {
        const auto& ax = binning_.axis(0);
        if (!(ax.value(0) < ax.value(ax.size())))
          throw std::invalid_argument("eec::ProjectedENC: precompute::PairBins requires an increasing axis 0");
      }
    }

    explicit ProjectedENC(double n, Extra extra = Extra{})
      requires(!Output::binned)
        : n_(n), extra_(std::move(extra)) {}

    const Binning& binning() const { return binning_; }

    // Forwards to Extra::update between jets.
    template <class... Args>
    void updateExtra(Args&&... args) {
      constexpr bool ok = requires(Extra& e, Args&&... a) { e.update(std::forward<Args>(a)...); };
      static_assert(ok, "eec::ProjectedENC::updateExtra: Extra has no update() accepting these arguments");
      if constexpr (ok)
        extra_.update(std::forward<Args>(args)...);
    }

    Result operator()(const std::vector<fastjet::PseudoJet>& particles, double externalNorm = kNoExternalNorm) const {
      fillJetData<Weight, Norm>(jd_, particles, n_, externalNorm);
      const std::size_t M = jd_.size();
      const std::size_t expected = kOrdered ? ipow(M) : detail::multichoose(M, N);
      if constexpr (Output::binned)
        acc_.begin(binning_, expected);
      else
        acc_.begin(expected);
      prepareMatrices();
      if (M > 0) {
        std::array<std::uint32_t, N> idx{};
        recurse<0>(0, 1.0, initialMax(), 1, 0, idx);
      }
      return acc_.finish(JetMeta{jd_.E_jet, n_, N});
    }

  private:
    using Loop = typename Strategy::loop_type;
    using Precompute = typename Strategy::precompute_type;
    static constexpr bool kOrdered = std::is_same_v<Loop, strategy::loop::Ordered>;
    static constexpr bool kPairBins = std::is_same_v<Precompute, strategy::precompute::PairBins>;
    static constexpr bool kAngleMatrix = std::is_same_v<Precompute, strategy::precompute::PairAngles>;
    static_assert(kOrdered || std::is_same_v<Loop, strategy::loop::Multiset>, "unknown strategy loop");
    static_assert(kPairBins || kAngleMatrix || std::is_same_v<Precompute, strategy::precompute::None>,
                  "unknown strategy precompute");
    static_assert(!kPairBins || Output::binned, "precompute::PairBins requires a binned output");
    static constexpr std::uint32_t kFactorial = detail::factorial(N);

    using MaxT = std::conditional_t<kPairBins, std::uint32_t, double>;

    static std::size_t ipow(std::size_t M) {
      std::size_t r = 1;
      for (int k = 0; k < N; ++k)
        r *= M;
      return r;
    }

    MaxT initialMax() const {
      if constexpr (kPairBins)
        return binning_.binIndex(0, 0.0);
      else
        return 0.0;
    }

    void prepareMatrices() const {
      const std::size_t M = jd_.size();
      if constexpr (kAngleMatrix) {
        angMat_.resize(M * M);
        for (std::size_t i = 0; i < M; ++i) {
          angMat_[i * M + i] = 0.0;
          for (std::size_t j = i + 1; j < M; ++j)
            angMat_[i * M + j] = angMat_[j * M + i] = Angle::pair(jd_, i, j);
        }
      } else if constexpr (kPairBins) {
        binMat_.resize(M * M);
        const std::uint32_t zero = binning_.binIndex(0, 0.0);
        for (std::size_t i = 0; i < M; ++i) {
          binMat_[i * M + i] = zero;
          for (std::size_t j = i + 1; j < M; ++j)
            binMat_[i * M + j] = binMat_[j * M + i] = binning_.binIndex(0, Angle::pair(jd_, i, j));
        }
      }
    }

    MaxT pairMax(std::uint32_t i, std::uint32_t j) const {
      if constexpr (kPairBins)
        return binMat_[i * jd_.size() + j];
      else if constexpr (kAngleMatrix)
        return angMat_[i * jd_.size() + j];
      else
        return Angle::pair(jd_, i, j);
    }

    // den is the running product of multiplicity factorials; run is the multiplicity of idx[d-1].
    template <int d>
    void recurse(std::uint32_t start,
                 double wprod,
                 MaxT mx,
                 std::uint32_t den,
                 std::uint32_t run,
                 std::array<std::uint32_t, N>& idx) const {
      const auto M = static_cast<std::uint32_t>(jd_.size());
      for (std::uint32_t j = start; j < M; ++j) {
        idx[d] = j;
        const double w = wprod * jd_.w[j];
        MaxT m = mx;
        for (int k = 0; k < d; ++k)
          m = std::max(m, pairMax(idx[k], j));
        std::uint32_t r = 1;
        if constexpr (d > 0 && !kOrdered) {
          if (j == idx[d - 1])
            r = run + 1;
        }
        const std::uint32_t dn = den * r;
        if constexpr (d + 1 == N)
          leaf(idx, w, m, dn);
        else
          recurse<d + 1>(kOrdered ? 0 : j, w, m, dn, r, idx);
      }
    }

    void leaf(const std::array<std::uint32_t, N>& idx, double w, MaxT m, std::uint32_t den) const {
      const std::uint32_t sym = kOrdered ? 1u : kFactorial / den;
      if constexpr (Output::binned) {
        std::uint32_t bin;
        if constexpr (kPairBins)
          bin = m;
        else
          bin = binning_.binIndex(0, m);
        if constexpr (K > 0) {
          const auto ex = extra_(jd_, idx);
          for (std::size_t k = 0; k < K; ++k)
            bin += static_cast<std::uint32_t>(binning_.stride(k + 1)) * binning_.binIndex(k + 1, ex[k]);
        }
        acc_.fill(bin, w, sym);
      } else {
        std::array<double, D> c;
        c[0] = m;
        if constexpr (K > 0) {
          const auto ex = extra_(jd_, idx);
          for (std::size_t k = 0; k < K; ++k)
            c[k + 1] = ex[k];
        }
        acc_.fill(c, w, sym);
      }
    }

    double n_;
    Binning binning_;
    Extra extra_;
    mutable JetData jd_;
    mutable std::vector<double> angMat_;
    mutable std::vector<std::uint32_t> binMat_;
    mutable typename Output::template Accumulator<D> acc_;
  };

}  // namespace eec

#endif
