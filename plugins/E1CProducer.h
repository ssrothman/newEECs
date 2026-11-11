#ifndef SRothman_EECs_plugins_E1CProducer_h
#define SRothman_EECs_plugins_E1CProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "SRothman/DataFormats/interface/EECConversion.h"
#include "SRothman/DataFormats/interface/EECJet.h"
#include "SRothman/EECs/interface/OnePointEC.h"
#include "SRothman/EECs/plugins/PluginCommon.h"

#include <cmath>
#include <vector>

// One-point correlator (theta, phi) of every jet in `src`, in the per-jet frame chosen by `frame`.
template <class Weight, class Output>
class E1CProducer : public edm::stream::EDProducer<> {
  using Calc = eec::OnePointEC<Weight, eec::norm::External, eec::frame::BoostRotate, Output>;
  using Result = typename Calc::Result;

public:
  explicit E1CProducer(const edm::ParameterSet& cfg)
      : src_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        put_(produces()),
        n_(cfg.getParameter<double>("n")),
        norm_(eecplugin::configGuard("E1CProducer", [&] {
          return eec::config::parseNorm(cfg.getParameter<std::string>("norm"));
        })),
        frame_(eecplugin::configGuard("E1CProducer", [&] {
          return eec::config::parseFrame(cfg.getParameter<std::string>("frame"));
        })),
        binning_(eecplugin::configGuard("E1CProducer", [&] {
          if constexpr (Output::binned)
            return eecplugin::binningFrom(cfg);
          else
            return eec::Binning{};
        })),
        calc_(eecplugin::configGuard("E1CProducer", [&] {
          if constexpr (Output::binned)
            return Calc(n_, binning_);
          else
            return Calc(n_);
        })) {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("eecJets"));
    desc.add<double>("n", 1.0)->setComment("power of the energy weights");
    desc.add<std::string>("norm", "sumPt")->setComment("E_jet: sumE, sumPt, jetE, jetPt or unity");
    desc.add<std::string>("frame", "jetHelicity")->setComment("lab, jetHelicity or jetCS");
    if constexpr (Output::binned)
      eecplugin::addAxes(desc,
                         {eecplugin::axisPSet("regular", 20, 0.0, M_PI), eecplugin::axisPSet("regular", 20, -M_PI, M_PI)});
    descriptions.addWithDefaultLabel(desc);
  }

  void produce(edm::Event& event, const edm::EventSetup&) override {
    const auto& jets = event.get(src_);
    std::vector<Result> out;
    out.reserve(jets.size());
    for (std::size_t i = 0; i < jets.size(); ++i) {
      const auto ps = srothman::toPseudoJets(jets[i]);
      const auto jet = srothman::jetPseudoJet(jets[i]);
      const double E = eec::config::normValue(norm_, ps, jet);
      if (ps.empty()) {
        out.push_back(eecplugin::emptyResult<Result>({E, n_, 1}, binning_));
        continue;
      }
      try {
        out.push_back(calc_(ps, eec::config::makeFrame(frame_, jet), E));
      } catch (const std::exception& ex) {
        throw cms::Exception("EECs") << "E1CProducer: jet " << i << " (pt " << jets[i].pt << "): " << ex.what();
      }
    }
    event.emplace(put_, std::move(out));
  }

private:
  const edm::EDGetTokenT<srothman::EECJetCollection> src_;
  const edm::EDPutTokenT<std::vector<Result>> put_;
  const double n_;
  const eec::config::NormSource norm_;
  const eec::config::FrameKind frame_;
  const eec::Binning binning_;
  const Calc calc_;
};

#endif
