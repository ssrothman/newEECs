#ifndef SRothman_EECs_Covariance_h
#define SRothman_EECs_Covariance_h

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace eec {

  // Aggregates per-event histograms h_e (flow-inclusive linear bins) with event weights w_e:
  //   H = sum_e w_e h_e,   C = sum_e w_e^2 h_e h_e^T,   plus sum w, sum w^2 and the event count.
  // C is stored as its packed upper triangle and updated only over the bins touched in each event.
  // Not thread-safe.
  class CovarianceAccumulator {
  public:
    // Throws std::invalid_argument if nBins is 0 or exceeds maxBins (the packed C holds nBins(nBins+1)/2 doubles).
    CovarianceAccumulator(std::size_t nBins, std::size_t maxBins) : n_(nBins) {
      if (nBins == 0 || nBins > maxBins)
        throw std::invalid_argument("eec::CovarianceAccumulator: nBins = " + std::to_string(nBins) +
                                    " must be in [1, " + std::to_string(maxBins) + "]");
      H_.assign(n_, 0.0);
      C_.assign(n_ * (n_ + 1) / 2, 0.0);
      event_.assign(n_, 0.0);
      seen_.assign(n_, 0);
    }

    std::size_t nBins() const { return n_; }

    // Discards any partial event left by an exception between beginEvent and endEvent.
    void beginEvent() {
      for (auto b : touched_) {
        event_[b] = 0;
        seen_[b] = 0;
      }
      touched_.clear();
    }

    // Adds `value` to bin `bin` of the current event's histogram; repeated bins are summed.
    void add(std::uint32_t bin, double value) {
      if (bin >= n_)
        throw std::out_of_range("eec::CovarianceAccumulator: bin out of range");
      if (!seen_[bin]) {
        seen_[bin] = 1;
        touched_.push_back(bin);
      }
      event_[bin] += value;
    }

    // Folds the current event into H, C and the sums with weight w; an event with no add() still counts.
    void endEvent(double w) {
      std::sort(touched_.begin(), touched_.end());
      const double w2 = w * w;
      const std::size_t k = touched_.size();
      for (std::size_t a = 0; a < k; ++a) {
        const std::uint32_t i = touched_[a];
        const double hi = event_[i];
        H_[i] += w * hi;
        double* row = C_.data() + rowStart(i);
        const double whi = w2 * hi;
        for (std::size_t b = a; b < k; ++b) {
          const std::uint32_t j = touched_[b];
          row[j - i] += whi * event_[j];
        }
      }
      for (auto b : touched_) {
        event_[b] = 0;
        seen_[b] = 0;
      }
      touched_.clear();
      sumW_ += w;
      sumW2_ += w2;
      ++nEvents_;
    }

    // Adds another accumulator over the same bins (e.g. one per stream).
    void merge(const CovarianceAccumulator& o) {
      if (o.n_ != n_)
        throw std::invalid_argument("eec::CovarianceAccumulator::merge: different number of bins");
      for (std::size_t i = 0; i < H_.size(); ++i)
        H_[i] += o.H_[i];
      for (std::size_t i = 0; i < C_.size(); ++i)
        C_[i] += o.C_[i];
      sumW_ += o.sumW_;
      sumW2_ += o.sumW2_;
      nEvents_ += o.nEvents_;
    }

    const std::vector<double>& values() const { return H_; }
    double covariance(std::size_t i, std::size_t j) const {
      if (i > j)
        std::swap(i, j);
      return C_[rowStart(i) + (j - i)];
    }
    // Full symmetric C, row-major nBins x nBins.
    std::vector<double> denseCovariance() const {
      std::vector<double> out(n_ * n_);
      for (std::size_t i = 0; i < n_; ++i)
        for (std::size_t j = i; j < n_; ++j)
          out[i * n_ + j] = out[j * n_ + i] = C_[rowStart(i) + (j - i)];
      return out;
    }
    double sumW() const { return sumW_; }
    double sumW2() const { return sumW2_; }
    std::uint64_t nEvents() const { return nEvents_; }

  private:
    // Offset of element (i, i) in the packed row-major upper triangle.
    std::size_t rowStart(std::size_t i) const { return i * n_ - i * (i - 1) / 2; }

    std::size_t n_;
    std::vector<double> H_, C_;
    std::vector<double> event_;
    std::vector<char> seen_;
    std::vector<std::uint32_t> touched_;
    double sumW_ = 0, sumW2_ = 0;
    std::uint64_t nEvents_ = 0;
  };

}  // namespace eec

#endif
