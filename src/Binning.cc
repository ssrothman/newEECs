#include "SRothman/EECs/interface/Binning.h"

#include <limits>
#include <stdexcept>

namespace eec {

  Binning::Binning(std::vector<Axis> axes) : axes_(std::move(axes)) {
    if (axes_.empty())
      throw std::invalid_argument("eec::Binning: at least one axis required");
    std::size_t total = 1;
    for (const auto& ax : axes_) {
      const std::size_t ext = static_cast<std::size_t>(ax.size()) + 2;
      strides_.push_back(total);
      extents_.push_back(ext);
      total *= ext;
      if (total > std::numeric_limits<std::uint32_t>::max())
        throw std::invalid_argument("eec::Binning: total bin count exceeds 2^32-1");
    }
    total_ = total;
  }

  std::vector<std::uint32_t> Binning::unravel(std::uint32_t lin) const {
    std::vector<std::uint32_t> idx(rank());
    for (std::size_t k = 0; k < rank(); ++k) {
      idx[k] = lin % extents_[k];
      lin /= extents_[k];
    }
    return idx;
  }

  Histogram Binning::toHistogram(const std::vector<double>& dense) const {
    if (dense.size() != total_)
      throw std::invalid_argument("eec::Binning::toHistogram: size mismatch");
    auto h = bh::make_histogram_with(bh::dense_storage<double>(), axes_);
    for (auto&& x : bh::indexed(h, bh::coverage::all)) {
      std::size_t lin = 0;
      for (std::size_t k = 0; k < rank(); ++k)
        lin += strides_[k] * static_cast<std::size_t>(x.index(k) + 1);
      *x = dense[lin];
    }
    return h;
  }

}  // namespace eec
