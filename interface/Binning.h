#ifndef SRothman_EECs_Binning_h
#define SRothman_EECs_Binning_h

#include <boost/histogram.hpp>

#include <array>
#include <cstdint>
#include <vector>

namespace eec {

  namespace bh = boost::histogram;

  // All three axis types have underflow and overflow bins fixed by their type.
  using RegularAxis = bh::axis::regular<double>;
  using LogAxis = bh::axis::regular<double, bh::axis::transform::log>;
  using VariableAxis = bh::axis::variable<double>;
  using Axis = bh::axis::variant<RegularAxis, LogAxis, VariableAxis>;
  using Histogram = bh::histogram<std::vector<Axis>, bh::dense_storage<double>>;

  // Runtime binning. Bin indices include flow bins (0 = underflow, extent-1 = overflow);
  // the linear index has axis 0 fastest.
  class Binning {
  public:
    Binning() = default;
    explicit Binning(std::vector<Axis> axes);

    std::size_t rank() const { return axes_.size(); }
    std::size_t totalBins() const { return total_; }
    const Axis& axis(std::size_t k) const { return axes_[k]; }
    std::size_t extent(std::size_t k) const { return extents_[k]; }
    std::size_t stride(std::size_t k) const { return strides_[k]; }

    std::uint32_t binIndex(std::size_t k, double x) const {
      return static_cast<std::uint32_t>(axes_[k].index(x) + 1);
    }

    template <std::size_t D>
    std::uint32_t linearIndex(const std::array<double, D>& x) const {
      std::uint32_t lin = 0;
      for (std::size_t k = 0; k < D; ++k)
        lin += static_cast<std::uint32_t>(strides_[k]) * binIndex(k, x[k]);
      return lin;
    }

    // Linear index from an already computed axis-0 bin and the values on axes 1..K.
    template <std::size_t K>
    std::uint32_t linearIndex(std::uint32_t bin0, const std::array<double, K>& rest) const {
      std::uint32_t lin = bin0;
      for (std::size_t k = 0; k < K; ++k)
        lin += static_cast<std::uint32_t>(strides_[k + 1]) * binIndex(k + 1, rest[k]);
      return lin;
    }

    std::vector<std::uint32_t> unravel(std::uint32_t lin) const;

    // boost::histogram with the same axes, holding the given dense (flow-inclusive) values.
    Histogram toHistogram(const std::vector<double>& dense) const;

  private:
    std::vector<Axis> axes_;
    std::vector<std::size_t> extents_, strides_;
    std::size_t total_ = 0;
  };

}  // namespace eec

#endif
