#ifndef SRothman_EECs_TriangleShape_h
#define SRothman_EECs_TriangleShape_h

#include "SRothman/EECs/interface/Concepts.h"
#include "SRothman/EECs/interface/JetData.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>

namespace eec {

  // Shape of a triangle with sides R_L >= R_M >= R_S, arXiv:2201.07800 Eq. (3): xi = R_S / R_M in [0, 1] and
  // phi = arcsin sqrt(1 - (R_L - R_M)^2 / R_S^2) in [0, pi/2]. R_S = 0 (coincident particles, e.g. a repeated index),
  // which the paper leaves undefined, gives (-1, -1): the underflow bin of axes starting at 0 holds exactly these.
  inline std::array<double, 2> triangleShape(double a, double b, double c) {
    std::array<double, 3> r{{a, b, c}};
    std::sort(r.begin(), r.end());
    const double rs = r[0], rm = r[1], rl = r[2];
    if (rs <= 0)
      return {{-1.0, -1.0}};
    const double d = (rl - rm) / rs;
    // The clamp absorbs rounding where the triangle inequality is saturated (collinear particles).
    return {{rs / rm, std::asin(std::sqrt(std::clamp(1.0 - d * d, 0.0, 1.0)))}};
  }

  // Extra coordinates (xi, phi) of a three-particle tuple; with ProjectedENC<3, ...> axis 0 is R_L, which makes the
  // resolved E3C. Angle must be the calculator's angle policy.
  template <AnglePolicy Angle>
  struct TriangleShape {
    static constexpr std::size_t K = 2;
    template <std::size_t N>
    std::array<double, 2> operator()(const JetData& d, const std::array<std::uint32_t, N>& idx) const {
      static_assert(N == 3, "eec::TriangleShape needs three-particle tuples");
      // By index: Angle3D's self-angle is rounding-sized, not 0.
      if (idx[0] == idx[1] || idx[0] == idx[2] || idx[1] == idx[2])
        return {{-1.0, -1.0}};
      return triangleShape(Angle::pair(d, idx[0], idx[1]), Angle::pair(d, idx[0], idx[2]), Angle::pair(d, idx[1], idx[2]));
    }
  };

}  // namespace eec

#endif
