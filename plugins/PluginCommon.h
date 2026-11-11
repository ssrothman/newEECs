#ifndef SRothman_EECs_plugins_PluginCommon_h
#define SRothman_EECs_plugins_PluginCommon_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "SRothman/EECs/interface/Config.h"
#include "SRothman/EECs/interface/Outputs.h"

#include <exception>
#include <string>
#include <type_traits>
#include <vector>

namespace eecplugin {

  // Rethrows library exceptions raised while reading the configuration as cms::Exception("Configuration").
  template <class F>
  auto configGuard(const std::string& module, F&& f) {
    try {
      return f();
    } catch (const cms::Exception&) {
      throw;
    } catch (const std::exception& ex) {
      throw cms::Exception("Configuration") << module << ": " << ex.what();
    }
  }

  inline edm::ParameterSet axisPSet(const std::string& type, unsigned nbins, double lo, double hi) {
    edm::ParameterSet p;
    p.addParameter<std::string>("type", type);
    p.addParameter<unsigned>("nbins", nbins);
    p.addParameter<double>("lo", lo);
    p.addParameter<double>("hi", hi);
    p.addParameter<std::vector<double>>("edges", {});
    return p;
  }

  // Adds the "axes" VPSet: each axis has type ("regular" | "log" | "variable"), nbins, lo, hi, edges.
  inline void addAxes(edm::ParameterSetDescription& desc, const std::vector<edm::ParameterSet>& defaults) {
    edm::ParameterSetDescription axis;
    axis.add<std::string>("type", "log");
    axis.add<unsigned>("nbins", 50);
    axis.add<double>("lo", 1e-3);
    axis.add<double>("hi", 1.0);
    axis.add<std::vector<double>>("edges", {});
    desc.addVPSet("axes", axis, defaults);
  }

  inline eec::Binning binningFrom(const edm::ParameterSet& cfg) {
    std::vector<eec::config::AxisSpec> specs;
    for (const auto& a : cfg.getParameter<std::vector<edm::ParameterSet>>("axes"))
      specs.push_back({a.getParameter<std::string>("type"),
                       a.getParameter<unsigned>("nbins"),
                       a.getParameter<double>("lo"),
                       a.getParameter<double>("hi"),
                       a.getParameter<std::vector<double>>("edges")});
    return eec::config::makeBinning(specs);
  }

  // Result for a jet with no constituents: no entries (Dense: all bins zero).
  template <class R>
  R emptyResult(const eec::JetMeta& meta, const eec::Binning& b) {
    R r;
    r.meta = meta;
    if constexpr (std::is_same_v<R, eec::DenseResult>)
      r.values.assign(b.totalBins(), 0.0);
    return r;
  }

}  // namespace eecplugin

#endif
