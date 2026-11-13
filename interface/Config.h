#ifndef SRothman_EECs_Config_h
#define SRothman_EECs_Config_h

#include "SRothman/EECs/interface/Binning.h"
#include "SRothman/EECs/interface/OnePointEC.h"

#include <fastjet/PseudoJet.hh>

#include <cstdint>
#include <string>
#include <vector>

// Runtime choices made by framework configurations (CMSSW python) for calculators that use norm::External.
namespace eec::config {

  enum class NormSource { SumE, SumPt, JetE, JetPt, Unity };
  // Accepts "sumE", "sumPt", "jetE", "jetPt", "unity"; throws std::invalid_argument otherwise.
  NormSource parseNorm(const std::string& s);
  double normValue(NormSource src, const std::vector<fastjet::PseudoJet>& constituents, const fastjet::PseudoJet& jet);

  struct AxisSpec {
    std::string type;  // "regular", "log" or "variable"
    unsigned nbins = 0;
    double lo = 0, hi = 0;
    std::vector<double> edges;  // "variable" only
  };
  // Throws std::invalid_argument for an unknown type or nbins == 0 (regular/log).
  Axis makeAxis(const AxisSpec& spec);
  Binning makeBinning(const std::vector<AxisSpec>& specs);

  // Codes written next to the edges in output files.
  enum class AxisType : std::uint8_t { Regular = 0, Log = 1, Variable = 2 };
  AxisType axisType(const Axis& a);
  // Bin edges without flow bins: size() + 1 values.
  std::vector<double> axisEdges(const Axis& a);

  // E1C frames for a jet with four-momentum `jet`:
  //  Lab:         lab polar angle and azimuth.
  //  JetAxis:     lab frame (no boost); z along the jet direction, x in the plane of z and the beam (toward +z beam).
  //  JetHelicity: jet rest frame; z along the jet direction, x in the plane of z and the beam (toward +z beam).
  //  JetCS:       jet rest frame, Collins-Soper axes: z bisects p_A and -p_B, x opposite to p_A + p_B, where
  //               beam A is the beam moving along sign(jet pz) (beam +z if pz = 0).
  enum class FrameKind { Lab, JetAxis, JetHelicity, JetCS };
  // Accepts "lab", "jetAxis", "jetHelicity", "jetCS"; throws std::invalid_argument otherwise.
  FrameKind parseFrame(const std::string& s);
  // Throws std::invalid_argument (from frame::BoostRotate) for a non-timelike jet or degenerate axes.
  frame::BoostRotate makeFrame(FrameKind kind, const fastjet::PseudoJet& jet);

}  // namespace eec::config

#endif
