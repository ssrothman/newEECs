#ifndef SRothman_EECs_ProjectedENCRuntime_h
#define SRothman_EECs_ProjectedENCRuntime_h

#include "SRothman/EECs/interface/ProjectedENC.h"

#include <stdexcept>
#include <utility>
#include <variant>

namespace eec {

  // ProjectedENC with N chosen at runtime (2 <= N <= 6); other parameters stay compile-time.
  template <class Weight,
            class Norm,
            class Angle,
            class Output,
            class Extra = NoExtra,
            class Strategy = strategy::DefaultFor<Output>>
  class ProjectedENCRuntimeN {
    template <int N>
    using Impl = ProjectedENC<N, Weight, Norm, Angle, Output, Extra, Strategy>;
    using Variant = std::variant<Impl<2>, Impl<3>, Impl<4>, Impl<5>, Impl<6>>;

  public:
    using Result = typename Impl<2>::Result;

    ProjectedENCRuntimeN(int N, double n, Binning binning, Extra extra = Extra{})
      requires(Output::binned)
        : impl_(make(N, [&]<int NN>() { return Impl<NN>(n, binning, extra); })) {}

    ProjectedENCRuntimeN(int N, double n, Extra extra = Extra{})
      requires(!Output::binned)
        : impl_(make(N, [&]<int NN>() { return Impl<NN>(n, extra); })) {}

    Result operator()(const std::vector<fastjet::PseudoJet>& particles, double externalNorm = kNoExternalNorm) const {
      return std::visit([&](const auto& c) { return c(particles, externalNorm); }, impl_);
    }

    template <class... Args>
    void update_extra(Args&&... args) {
      std::visit([&](auto& c) { c.update_extra(std::forward<Args>(args)...); }, impl_);
    }

  private:
    template <class F>
    static Variant make(int N, F&& f) {
      switch (N) {
        case 2:
          return Variant(f.template operator()<2>());
        case 3:
          return Variant(f.template operator()<3>());
        case 4:
          return Variant(f.template operator()<4>());
        case 5:
          return Variant(f.template operator()<5>());
        case 6:
          return Variant(f.template operator()<6>());
        default:
          throw std::invalid_argument("eec::ProjectedENCRuntimeN: N must be in [2, 6]");
      }
    }

    Variant impl_;
  };

}  // namespace eec

#endif
