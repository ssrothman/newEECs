#ifndef SRothman_EECs_Concepts_h
#define SRothman_EECs_Concepts_h

#include "SRothman/EECs/interface/Binning.h"
#include "SRothman/EECs/interface/JetData.h"
#include "SRothman/EECs/interface/Outputs.h"

#include <fastjet/PseudoJet.hh>

#include <array>
#include <concepts>
#include <cstdint>
#include <type_traits>
#include <vector>

namespace eec {

  // Per-particle quantity E_p.
  template <class W>
  concept WeightPolicy = requires(const fastjet::PseudoJet& p) {
    { W::value(p) } -> std::convertible_to<double>;
  };

  // Jet normalisation E_jet; an optional usesExternalNorm must be a constant bool.
  template <class Nm>
  concept NormPolicy = requires(const std::vector<fastjet::PseudoJet>& ps, double ext) {
    { Nm::value(ps, ext) } -> std::convertible_to<double>;
  } && (!requires { Nm::usesExternalNorm; } || requires { std::bool_constant<Nm::usesExternalNorm>{}; });

  // Pair angle theta_ij.
  template <class A>
  concept AnglePolicy = requires(const JetData& d, std::size_t i, std::size_t j) {
    { A::pair(d, i, j) } -> std::convertible_to<double>;
  };

  // K extra pENC coordinates for an N-tuple.
  template <class X, int N>
  concept ExtraPolicy = requires { std::integral_constant<std::size_t, X::K>{}; } &&
                        requires(const X& x, const JetData& d, const std::array<std::uint32_t, N>& idx) {
                          { x(d, idx) } -> std::convertible_to<std::array<double, X::K>>;
                        };

  // E1C frame: D coordinates per particle.
  template <class F>
  concept FramePolicy =
      requires { std::integral_constant<std::size_t, F::D>{}; } && requires(const F& f, const fastjet::PseudoJet& p) {
        { f(p) } -> std::convertible_to<std::array<double, F::D>>;
      };

  // Output with D coordinates per fill; begin/fill signatures depend on `binned`.
  template <class O, std::size_t D>
  concept OutputPolicy =
      requires {
        std::bool_constant<O::binned>{};
        typename O::template Result<D>;
        typename O::template Accumulator<D>;
      } && std::default_initializable<typename O::template Accumulator<D>> &&
      requires(typename O::template Accumulator<D>& acc, const JetMeta& m) {
        { acc.finish(m) } -> std::convertible_to<typename O::template Result<D>>;
      } &&
      (O::binned ? requires(typename O::template Accumulator<D>& acc, const Binning& b, std::size_t n, std::uint32_t bin, double w, std::uint32_t sym) {
        acc.begin(b, n);
        acc.fill(bin, w, sym);
      } : requires(typename O::template Accumulator<D>& acc, std::size_t n, const std::array<double, D>& c, double w, std::uint32_t sym) {
        acc.begin(n);
        acc.fill(c, w, sym);
      });

}  // namespace eec

#endif
