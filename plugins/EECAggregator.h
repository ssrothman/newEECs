#ifndef SRothman_EECs_plugins_EECAggregator_h
#define SRothman_EECs_plugins_EECAggregator_h

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SRothman/EECs/interface/Covariance.h"
#include "SRothman/EECs/plugins/PluginCommon.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TNamed.h"
#include "TParameter.h"

#include <cmath>
#include <optional>
#include <string>
#include <vector>

// Sums the binned EEC results of all jets of an event into h_e and aggregates H = sum w h_e and
// C = sum w^2 h_e h_e^T with w the generator weight (1 if genEventInfo is empty). Writes to the module's
// TFileService directory: values (TH1D over linear bins, errors sqrt(C_ii)), covariance (TH2D), axis<k>
// (empty TH1D carrying the bin edges of axis k), axisTypes (TNamed, comma-separated), sumW, sumW2, nEvents (TParameter).
template <class Result>
class EECAggregator : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit EECAggregator(const edm::ParameterSet& cfg)
      : src_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        binning_(eecplugin::configGuard("EECAggregator", [&] { return eecplugin::binningFrom(cfg); })),
        acc_(eecplugin::configGuard("EECAggregator", [&] {
          return eec::CovarianceAccumulator(binning_.totalBins(), cfg.getParameter<unsigned>("maxBins"));
        })) {
    const auto gen = cfg.getParameter<edm::InputTag>("genEventInfo");
    if (!gen.label().empty())
      gen_ = consumes(gen);
    usesResource(TFileService::kSharedResource);
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("eec"));
    desc.add<edm::InputTag>("genEventInfo", edm::InputTag("generator"))->setComment("empty label: unit weights");
    desc.add<unsigned>("maxBins", 5000)->setComment("upper bound on the flow-inclusive bin count (memory guard)");
    eecplugin::addAxes(desc, {eecplugin::axisPSet("log", 50, 1e-3, 1.0)});
    descriptions.addWithDefaultLabel(desc);
  }

  void analyze(const edm::Event& event, const edm::EventSetup&) override {
    const double w = gen_ ? event.get(*gen_).weight() : 1.0;
    const auto& results = event.get(src_);
    acc_.beginEvent();
    for (const auto& r : results) {
      if constexpr (std::is_same_v<Result, eec::DenseResult>) {
        if (r.values.size() != acc_.nBins())
          throw cms::Exception("EECs") << "EECAggregator " << moduleDescription().moduleLabel() << ": dense result has " << r.values.size()
                                       << " bins, binning has " << acc_.nBins();
        for (std::size_t b = 0; b < r.values.size(); ++b)
          if (r.values[b] != 0)
            acc_.add(b, r.values[b]);
      } else {
        for (const auto& e : r.entries) {
          if (e.bin >= acc_.nBins())
            throw cms::Exception("EECs") << "EECAggregator " << moduleDescription().moduleLabel() << ": bin " << e.bin
                                         << " outside the binning (" << acc_.nBins() << " bins)";
          if constexpr (std::is_same_v<Result, eec::SparseUnmergedResult>)
            acc_.add(e.bin, e.weight * e.symfactor);
          else
            acc_.add(e.bin, e.weight);
        }
      }
    }
    acc_.endEvent(w);
  }

  void endJob() override {
    edm::Service<TFileService> fs;
    TFileService& dir = *fs;
    const int nb = acc_.nBins();
    auto* values = dir.make<TH1D>("values", "EEC values;linear bin", nb, -0.5, nb - 0.5);
    auto* cov = dir.make<TH2D>("covariance", "sum w^2 h h^T;linear bin;linear bin", nb, -0.5, nb - 0.5, nb, -0.5, nb - 0.5);
    for (int i = 0; i < nb; ++i) {
      values->SetBinContent(i + 1, acc_.values()[i]);
      values->SetBinError(i + 1, std::sqrt(std::max(0.0, acc_.covariance(i, i))));
      for (int j = 0; j < nb; ++j)
        cov->SetBinContent(i + 1, j + 1, acc_.covariance(i, j));
    }
    values->SetEntries(acc_.nEvents());
    cov->SetEntries(acc_.nEvents());
    std::string types;
    for (std::size_t k = 0; k < binning_.rank(); ++k) {
      const auto e = eec::config::axisEdges(binning_.axis(k));
      // TVectorD has no name of its own, so the edges are stored as the binning of an empty TH1D.
      const std::string hname = "axis" + std::to_string(k);
      dir.make<TH1D>(hname.c_str(), ("bin edges of axis " + std::to_string(k)).c_str(), int(e.size()) - 1, e.data());
      switch (eec::config::axisType(binning_.axis(k))) {
        case eec::config::AxisType::Regular:
          types += "regular";
          break;
        case eec::config::AxisType::Log:
          types += "log";
          break;
        case eec::config::AxisType::Variable:
          types += "variable";
          break;
      }
      if (k + 1 < binning_.rank())
        types += ",";
    }
    dir.make<TNamed>("axisTypes", types.c_str());
    dir.make<TParameter<double>>("sumW", acc_.sumW());
    dir.make<TParameter<double>>("sumW2", acc_.sumW2());
    dir.make<TParameter<Long64_t>>("nEvents", static_cast<Long64_t>(acc_.nEvents()));
  }

private:
  const edm::EDGetTokenT<std::vector<Result>> src_;
  std::optional<edm::EDGetTokenT<GenEventInfoProduct>> gen_;
  const eec::Binning binning_;
  eec::CovarianceAccumulator acc_;
};

#endif
