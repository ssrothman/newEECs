#ifndef SRothman_EECs_Outputs_h
#define SRothman_EECs_Outputs_h

#include "SRothman/EECs/interface/Binning.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

namespace eec {

  // Everything needed to change E_jet after the fact.
  struct JetMeta {
    double E_jet = 0;
    double n = 1;
    int N = 0;
  };

  // Factor that converts weights normalised by meta.E_jet into weights normalised by newEjet.
  inline double renormalization(const JetMeta& meta, double newEjet) {
    return std::pow(meta.E_jet / newEjet, meta.n * meta.N);
  }

  struct BinEntry {
    std::uint32_t bin;
    double weight;
  };
  struct UnmergedBinEntry {
    std::uint32_t bin;
    double weight;
    std::uint32_t symfactor;
  };
  template <std::size_t D>
  struct UnbinnedEntry {
    std::array<double, D> coord;
    double weight;
    std::uint32_t symfactor;
  };

  struct DenseResult {
    JetMeta meta;
    std::vector<double> values;
  };
  struct SparseResult {
    JetMeta meta;
    std::vector<BinEntry> entries;  // sorted by bin, unique bins
  };
  struct SparseUnmergedResult {
    JetMeta meta;
    std::vector<UnmergedBinEntry> entries;
  };
  template <std::size_t D>
  struct UnbinnedResult {
    JetMeta meta;
    std::vector<UnbinnedEntry<D>> entries;
  };

  // Output policies. Each provides `binned`, `Result<D>` and `Accumulator<D>` with
  // begin(...), fill(...), finish(meta). Accumulators keep scratch between jets.
  namespace output {

    struct Dense {
      static constexpr bool binned = true;
      template <std::size_t D>
      using Result = DenseResult;
      template <std::size_t D>
      class Accumulator {
      public:
        void begin(const Binning& b, std::size_t) { r_.values.assign(b.totalBins(), 0.0); }
        void fill(std::uint32_t bin, double w, std::uint32_t sym) { r_.values[bin] += w * sym; }
        DenseResult finish(const JetMeta& m) {
          r_.meta = m;
          return std::move(r_);
        }

      private:
        DenseResult r_;
      };
    };

    // Merges via a dense scratch array plus a touched-bin list.
    struct SparseMerged {
      static constexpr bool binned = true;
      template <std::size_t D>
      using Result = SparseResult;
      template <std::size_t D>
      class Accumulator {
      public:
        void begin(const Binning& b, std::size_t) {
          if (scratch_.size() != b.totalBins()) {
            scratch_.assign(b.totalBins(), 0.0);
            seen_.assign(b.totalBins(), 0);
          }
          touched_.clear();
        }
        void fill(std::uint32_t bin, double w, std::uint32_t sym) {
          if (!seen_[bin]) {
            seen_[bin] = 1;
            touched_.push_back(bin);
          }
          scratch_[bin] += w * sym;
        }
        SparseResult finish(const JetMeta& m) {
          std::sort(touched_.begin(), touched_.end());
          SparseResult r;
          r.meta = m;
          r.entries.reserve(touched_.size());
          for (auto b : touched_) {
            r.entries.push_back({b, scratch_[b]});
            scratch_[b] = 0;
            seen_[b] = 0;
          }
          touched_.clear();
          return r;
        }

      private:
        std::vector<double> scratch_;
        std::vector<char> seen_;
        std::vector<std::uint32_t> touched_;
      };
    };

    // Same result as SparseMerged, merging by sort + reduce instead of a dense scratch array.
    struct SparseMergedSort {
      static constexpr bool binned = true;
      template <std::size_t D>
      using Result = SparseResult;
      template <std::size_t D>
      class Accumulator {
      public:
        void begin(const Binning&, std::size_t expected) {
          buf_.clear();
          buf_.reserve(expected);
        }
        void fill(std::uint32_t bin, double w, std::uint32_t sym) { buf_.push_back({bin, w * sym}); }
        SparseResult finish(const JetMeta& m) {
          std::sort(buf_.begin(), buf_.end(), [](const BinEntry& a, const BinEntry& b) { return a.bin < b.bin; });
          SparseResult r;
          r.meta = m;
          for (const auto& e : buf_) {
            if (!r.entries.empty() && r.entries.back().bin == e.bin)
              r.entries.back().weight += e.weight;
            else
              r.entries.push_back(e);
          }
          buf_.clear();
          return r;
        }

      private:
        std::vector<BinEntry> buf_;
      };
    };

    struct SparseUnmerged {
      static constexpr bool binned = true;
      template <std::size_t D>
      using Result = SparseUnmergedResult;
      template <std::size_t D>
      class Accumulator {
      public:
        void begin(const Binning&, std::size_t expected) {
          r_.entries.clear();
          r_.entries.reserve(expected);
        }
        void fill(std::uint32_t bin, double w, std::uint32_t sym) { r_.entries.push_back({bin, w, sym}); }
        SparseUnmergedResult finish(const JetMeta& m) {
          r_.meta = m;
          return std::move(r_);
        }

      private:
        SparseUnmergedResult r_;
      };
    };

    struct Unbinned {
      static constexpr bool binned = false;
      template <std::size_t D>
      using Result = UnbinnedResult<D>;
      template <std::size_t D>
      class Accumulator {
      public:
        void begin(std::size_t expected) {
          r_.entries.clear();
          r_.entries.reserve(expected);
        }
        void fill(const std::array<double, D>& c, double w, std::uint32_t sym) { r_.entries.push_back({c, w, sym}); }
        UnbinnedResult<D> finish(const JetMeta& m) {
          r_.meta = m;
          return std::move(r_);
        }

      private:
        UnbinnedResult<D> r_;
      };
    };

  }  // namespace output

}  // namespace eec

#endif
