#ifndef SRothman_EECs_plugins_PENCProducer_h
#define SRothman_EECs_plugins_PENCProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "SRothman/DataFormats/interface/EECConversion.h"
#include "SRothman/DataFormats/interface/EECJet.h"
#include "SRothman/EECs/interface/ProjectedENCRuntime.h"
#include "SRothman/EECs/plugins/PluginCommon.h"

#include <vector>

// Projected N-point correlator of every jet in `src`. One output entry per jet, in input order.
// E_jet comes from the runtime `norm` choice (passed to norm::External).
template <class Weight, class Angle, class Output>
class PENCProducer : public edm::stream::EDProducer<> {
  using Calc = eec::ProjectedENCRuntimeN<Weight, eec::norm::External, Angle, Output>;
  using Result = typename Calc::Result;

public:
  explicit PENCProducer(const edm::ParameterSet& cfg)
      : src_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        put_(produces()),
        N_(cfg.getParameter<unsigned>("N")),
        n_(cfg.getParameter<double>("n")),
        norm_(eecplugin::configGuard("PENCProducer", [&] {
          return eec::config::parseNorm(cfg.getParameter<std::string>("norm"));
        })),
        binning_(eecplugin::configGuard("PENCProducer", [&] {
          if constexpr (Output::binned)
            return eecplugin::binningFrom(cfg);
          else
            return eec::Binning{};
        })),
        calc_(eecplugin::configGuard("PENCProducer", [&] {
          if constexpr (Output::binned)
            return Calc(N_, n_, binning_);
          else
            return Calc(N_, n_);
        })) {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("eecJets"));
    desc.add<unsigned>("N", 2)->setComment("number of particles per tuple, 2..6");
    desc.add<double>("n", 1.0)->setComment("power of the energy weights");
    desc.add<std::string>("norm", "sumPt")->setComment("E_jet: sumE, sumPt, jetE, jetPt or unity");
    if constexpr (Output::binned)
      eecplugin::addAxes(desc, {eecplugin::axisPSet("log", 50, 1e-3, 1.0)});
    descriptions.addWithDefaultLabel(desc);
  }

  void produce(edm::Event& event, const edm::EventSetup&) override {
    const auto& jets = event.get(src_);
    std::vector<Result> out;
    out.reserve(jets.size());
    for (std::size_t i = 0; i < jets.size(); ++i) {
      const auto ps = srothman::toPseudoJets(jets[i]);
      const double E = eec::config::normValue(norm_, ps, srothman::jetPseudoJet(jets[i]));
      if (ps.empty()) {
        out.push_back(eecplugin::emptyResult<Result>({E, n_, N_}, binning_));
        continue;
      }
      try {
        out.push_back(calc_(ps, E));
      } catch (const std::exception& ex) {
        throw cms::Exception("EECs") << "PENCProducer: jet " << i << " (pt " << jets[i].pt << "): " << ex.what();
      }
    }
    event.emplace(put_, std::move(out));
  }

private:
  const edm::EDGetTokenT<srothman::EECJetCollection> src_;
  const edm::EDPutTokenT<std::vector<Result>> put_;
  const std::size_t N_;
  const double n_;
  const eec::config::NormSource norm_;
  const eec::Binning binning_;
  const Calc calc_;
};

#endif
