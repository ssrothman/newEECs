#ifndef SRothman_EECs_ProjectedENC_h
#define SRothman_EECs_ProjectedENC_h

#include "SRothman/EECs/interface/Angles.h"
#include "SRothman/EECs/interface/Binning.h"
#include "SRothman/EECs/interface/JetData.h"
#include "SRothman/EECs/interface/Outputs.h"
#include "SRothman/EECs/interface/Policies.h"

#include <fastjet/PseudoJet.hh>

#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace eec {

  namespace strategy {
    // Ordered M^N loop, one fill per ordered tuple with symfactor 1. Reference only.
    struct Naive {};
    // Multiset loop i1 <= ... <= iN, pair angles computed on the fly.
    struct Multiset {};
    // Multiset loop with a precomputed M x M pair-angle matrix.
    struct PairAngles {};
    // Multiset loop with a precomputed M x M matrix of axis-0 bins, using
    // bin(max theta) == max bin(theta) for increasing axes. Falls back to PairAngles if unbinned.
    struct PairBins {};
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
            class Weight,
            class Norm,
            class Angle,
            class Output,
            class Extra = NoExtra,
            class Strategy = strategy::PairBins>
  class ProjectedENC {
    static_assert(N >= 1 && N <= 12, "symmetry factors are computed in 32-bit integers");

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
          throw std::invalid_argument("eec::ProjectedENC: strategy::PairBins requires an increasing axis 0");
      }
    }

    explicit ProjectedENC(double n, Extra extra = Extra{})
      requires(!Output::binned)
        : n_(n), extra_(std::move(extra)) {}

    const Binning& binning() const { return binning_; }

    Result operator()(const std::vector<fastjet::PseudoJet>& particles, double externalNorm = kNoExternalNorm) const {
      fillJetData<Weight, Norm>(jd_, particles, n_, externalNorm);
      const std::size_t M = jd_.size();
      const std::size_t expected = kNaive ? ipow(M) : detail::multichoose(M, N);
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
    static constexpr bool kNaive = std::is_same_v<Strategy, strategy::Naive>;
    static constexpr bool kPairBins = std::is_same_v<Strategy, strategy::PairBins> && Output::binned;
    static constexpr bool kAngleMatrix = std::is_same_v<Strategy, strategy::PairAngles> ||
                                         (std::is_same_v<Strategy, strategy::PairBins> && !Output::binned);
    static_assert(kNaive || kPairBins || kAngleMatrix || std::is_same_v<Strategy, strategy::Multiset>,
                  "unknown strategy");
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
        if constexpr (d > 0 && !kNaive) {
          if (j == idx[d - 1])
            r = run + 1;
        }
        const std::uint32_t dn = den * r;
        if constexpr (d + 1 == N)
          leaf(idx, w, m, dn);
        else
          recurse<d + 1>(kNaive ? 0 : j, w, m, dn, r, idx);
      }
    }

    void leaf(const std::array<std::uint32_t, N>& idx, double w, MaxT m, std::uint32_t den) const {
      const std::uint32_t sym = kNaive ? 1u : kFactorial / den;
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
