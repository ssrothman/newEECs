#ifndef SRothman_EECs_Angles_h
#define SRothman_EECs_Angles_h

#include "SRothman/EECs/interface/JetData.h"

#include <cmath>
#include <numbers>

namespace eec {

  // Pair angle theta_ij: struct with `static double pair(const JetData&, i, j)`, >= 0 and 0 for i == j.
  namespace angle {
    struct Angle3D {
      static double pair(const JetData& d, std::size_t i, std::size_t j) {
        // atan2 form keeps precision at small angles where acos(dot) does not.
        const double cx = d.uy[i] * d.uz[j] - d.uz[i] * d.uy[j];
        const double cy = d.uz[i] * d.ux[j] - d.ux[i] * d.uz[j];
        const double cz = d.ux[i] * d.uy[j] - d.uy[i] * d.ux[j];
        const double dot = d.ux[i] * d.ux[j] + d.uy[i] * d.uy[j] + d.uz[i] * d.uz[j];
        return std::atan2(std::sqrt(cx * cx + cy * cy + cz * cz), dot);
      }
    };

    // Uses pseudorapidity, not rapidity.
    struct DeltaR {
      static double pair(const JetData& d, std::size_t i, std::size_t j) {
        double dphi = d.phi[i] - d.phi[j];
        if (dphi > std::numbers::pi)
          dphi -= 2 * std::numbers::pi;
        else if (dphi < -std::numbers::pi)
          dphi += 2 * std::numbers::pi;
        const double deta = d.eta[i] - d.eta[j];
        return std::sqrt(dphi * dphi + deta * deta);
      }
    };
  }  // namespace angle

}  // namespace eec

#endif
