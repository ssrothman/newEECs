#include "SRothman/EECs/interface/Config.h"

#include <array>
#include <cmath>
#include <stdexcept>

namespace eec::config {

  NormSource parseNorm(const std::string& s) {
    if (s == "sumE")
      return NormSource::SumE;
    if (s == "sumPt")
      return NormSource::SumPt;
    if (s == "jetE")
      return NormSource::JetE;
    if (s == "jetPt")
      return NormSource::JetPt;
    if (s == "unity")
      return NormSource::Unity;
    throw std::invalid_argument("eec::config: unknown norm '" + s + "'");
  }

  double normValue(NormSource src, const std::vector<fastjet::PseudoJet>& ps, const fastjet::PseudoJet& jet) {
    double s = 0;
    switch (src) {
      case NormSource::SumE:
        for (const auto& p : ps)
          s += p.E();
        return s;
      case NormSource::SumPt:
        for (const auto& p : ps)
          s += p.pt();
        return s;
      case NormSource::JetE:
        return jet.E();
      case NormSource::JetPt:
        return jet.pt();
      case NormSource::Unity:
        return 1.0;
    }
    throw std::logic_error("eec::config::normValue: unhandled NormSource");
  }

  Axis makeAxis(const AxisSpec& a) {
    if (a.type == "variable")
      return VariableAxis(a.edges);
    if (a.nbins == 0)
      throw std::invalid_argument("eec::config: axis of type '" + a.type + "' needs nbins > 0");
    if (a.type == "regular")
      return RegularAxis(a.nbins, a.lo, a.hi);
    if (a.type == "log")
      return LogAxis(a.nbins, a.lo, a.hi);
    throw std::invalid_argument("eec::config: unknown axis type '" + a.type + "'");
  }

  Binning makeBinning(const std::vector<AxisSpec>& specs) {
    std::vector<Axis> axes;
    axes.reserve(specs.size());
    for (const auto& s : specs)
      axes.push_back(makeAxis(s));
    return Binning(std::move(axes));
  }

  AxisType axisType(const Axis& a) {
    if (bh::axis::get_if<RegularAxis>(&a))
      return AxisType::Regular;
    if (bh::axis::get_if<LogAxis>(&a))
      return AxisType::Log;
    return AxisType::Variable;
  }

  std::vector<double> axisEdges(const Axis& a) {
    std::vector<double> e(a.size() + 1);
    for (int i = 0; i <= a.size(); ++i)
      e[i] = a.value(i);
    return e;
  }

  FrameKind parseFrame(const std::string& s) {
    if (s == "lab")
      return FrameKind::Lab;
    if (s == "jetAxis")
      return FrameKind::JetAxis;
    if (s == "jetHelicity")
      return FrameKind::JetHelicity;
    if (s == "jetCS")
      return FrameKind::JetCS;
    throw std::invalid_argument("eec::config: unknown frame '" + s + "'");
  }

  namespace {
    std::array<double, 3> unit3(const fastjet::PseudoJet& p) {
      const double m = p.modp();
      return {{p.px() / m, p.py() / m, p.pz() / m}};
    }
  }  // namespace

  frame::BoostRotate makeFrame(FrameKind kind, const fastjet::PseudoJet& jet) {
    switch (kind) {
      case FrameKind::Lab:
        return frame::BoostRotate(fastjet::PseudoJet(0, 0, 0, 1), {{0, 0, 1}}, {{1, 0, 0}});
      case FrameKind::JetAxis:
        return frame::BoostRotate(fastjet::PseudoJet(0, 0, 0, 1), {{jet.px(), jet.py(), jet.pz()}}, {{0, 0, 1}});
      case FrameKind::JetHelicity: {
        // A boost along the jet direction leaves the plane of (jet, beam) unchanged, so lab vectors serve as rest-frame axes.
        return frame::BoostRotate(jet, {{jet.px(), jet.py(), jet.pz()}}, {{0, 0, 1}});
      }
      case FrameKind::JetCS: {
        const double sign = jet.pz() < 0 ? -1.0 : 1.0;
        fastjet::PseudoJet beamA(0, 0, sign, 1), beamB(0, 0, -sign, 1);
        beamA.unboost(jet);
        beamB.unboost(jet);
        const auto a = unit3(beamA), b = unit3(beamB);
        return frame::BoostRotate(jet,
                                  {{a[0] - b[0], a[1] - b[1], a[2] - b[2]}},
                                  {{-(a[0] + b[0]), -(a[1] + b[1]), -(a[2] + b[2])}});
      }
    }
    throw std::logic_error("eec::config::makeFrame: unhandled FrameKind");
  }

}  // namespace eec::config
