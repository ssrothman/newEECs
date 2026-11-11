#ifndef SRothman_EECs_plugins_EECTableProducer_h
#define SRothman_EECs_plugins_EECTableProducer_h

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "SRothman/EECs/plugins/PluginCommon.h"

#include <memory>
#include <string>
#include <vector>

namespace eecplugin {

  template <class R>
  struct TableTraits;
  template <>
  struct TableTraits<eec::DenseResult> {
    static constexpr bool binned = true, symfactor = false;
    static constexpr std::size_t D = 0;
  };
  template <>
  struct TableTraits<eec::SparseResult> {
    static constexpr bool binned = true, symfactor = false;
    static constexpr std::size_t D = 0;
  };
  template <>
  struct TableTraits<eec::SparseUnmergedResult> {
    static constexpr bool binned = true, symfactor = true;
    static constexpr std::size_t D = 0;
  };
  template <std::size_t DD>
  struct TableTraits<eec::UnbinnedResult<DD>> {
    static constexpr bool binned = false, symfactor = true;
    static constexpr std::size_t D = DD;
  };

}  // namespace eecplugin

// NanoAOD tables for a vector of per-jet EEC results (index-aligned with the jets):
//  event  "<name>"        one row per jet: E_jet, n, N, nEntries
//  event  "<name>Entries" one row per entry: iJet, then bin (binned) or coord names (unbinned), weight[, symfactor]
//  run    "<name>Axis<k>" bin edges of axis k (binned only), and "<name>AxisType" with one type code per axis
// Dense results write every bin of every jet.
template <class Result>
class EECTableProducer : public edm::global::EDProducer<edm::BeginRunProducer> {
  using Traits = eecplugin::TableTraits<Result>;

public:
  explicit EECTableProducer(const edm::ParameterSet& cfg)
      : src_(consumes(cfg.getParameter<edm::InputTag>("src"))),
        name_(cfg.getParameter<std::string>("name")),
        doc_(cfg.getParameter<std::string>("doc")),
        weightPrecision_(cfg.getParameter<int>("weightPrecision")),
        jetPut_(produces()),
        entryPut_(produces("entries")) {
    if constexpr (Traits::binned) {
      const auto b = eecplugin::configGuard("EECTableProducer", [&] { return eecplugin::binningFrom(cfg); });
      for (std::size_t k = 0; k < b.rank(); ++k) {
        edges_.push_back(eec::config::axisEdges(b.axis(k)));
        types_.push_back(static_cast<std::uint8_t>(eec::config::axisType(b.axis(k))));
        axisPut_.push_back(produces<nanoaod::FlatTable, edm::Transition::BeginRun>("axis" + std::to_string(k)));
      }
      typePut_ = produces<nanoaod::FlatTable, edm::Transition::BeginRun>("axisTypes");
    } else {
      coordNames_ = cfg.getParameter<std::vector<std::string>>("coordNames");
      coordPrecision_ = cfg.getParameter<int>("coordPrecision");
      if (coordNames_.size() != Traits::D)
        throw cms::Exception("Configuration")
            << "EECTableProducer: coordNames needs " << Traits::D << " names, got " << coordNames_.size();
    }
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("eec"));
    desc.add<std::string>("name", "EEC")->setComment("table name prefix");
    desc.add<std::string>("doc", "");
    desc.add<int>("weightPrecision", -1)->setComment("mantissa bits of the float weights; -1 = full float");
    if constexpr (Traits::binned) {
      eecplugin::addAxes(desc, {eecplugin::axisPSet("log", 50, 1e-3, 1.0)});
    } else {
      desc.add<std::vector<std::string>>("coordNames",
                                         Traits::D == 1 ? std::vector<std::string>{"angle"}
                                                        : std::vector<std::string>{"theta", "phi"});
      desc.add<int>("coordPrecision", -1)->setComment("mantissa bits of the float coordinates; -1 = full float");
    }
    descriptions.addWithDefaultLabel(desc);
  }

  void produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const override {
    const auto& results = event.get(src_);
    std::size_t total = 0;
    for (const auto& r : results)
      total += size(r);

    std::vector<double> Ejet;
    std::vector<float> n;
    std::vector<std::uint8_t> N;
    std::vector<std::uint32_t> nEntries;
    std::vector<std::uint16_t> iJet;
    std::vector<std::uint32_t> bin, sym;
    std::vector<float> weight;
    std::vector<std::vector<float>> coords(Traits::D);
    iJet.reserve(total);
    weight.reserve(total);

    for (std::size_t j = 0; j < results.size(); ++j) {
      const auto& r = results[j];
      Ejet.push_back(r.meta.E_jet);
      n.push_back(r.meta.n);
      N.push_back(r.meta.N);
      nEntries.push_back(size(r));
      if constexpr (std::is_same_v<Result, eec::DenseResult>) {
        for (std::size_t b = 0; b < r.values.size(); ++b) {
          iJet.push_back(j);
          bin.push_back(b);
          weight.push_back(r.values[b]);
        }
      } else {
        for (const auto& e : r.entries) {
          iJet.push_back(j);
          weight.push_back(e.weight);
          if constexpr (Traits::binned)
            bin.push_back(e.bin);
          else
            for (std::size_t k = 0; k < Traits::D; ++k)
              coords[k].push_back(e.coord[k]);
          if constexpr (Traits::symfactor)
            sym.push_back(e.symfactor);
        }
      }
    }

    std::unique_ptr<nanoaod::FlatTable> jt = std::make_unique<nanoaod::FlatTable>(results.size(), name_, false);
    jt->setDoc(doc_);
    jt->addColumn<double>("E_jet", Ejet, "E_jet used to normalise the weights");
    jt->addColumn<float>("n", n, "power of the energy weights");
    jt->addColumn<std::uint8_t>("N", N, "number of particles per tuple (1 for E1C)");
    jt->addColumn<std::uint32_t>("nEntries", nEntries, "number of rows in " + name_ + "Entries for this jet");

    std::unique_ptr<nanoaod::FlatTable> et = std::make_unique<nanoaod::FlatTable>(iJet.size(), name_ + "Entries", false);
    et->setDoc(doc_);
    et->addColumn<std::uint16_t>("iJet", iJet, "index of the jet in " + name_);
    if constexpr (Traits::binned)
      et->addColumn<std::uint32_t>("bin", bin, "flow-inclusive linear bin index (axis 0 fastest)");
    else
      for (std::size_t k = 0; k < Traits::D; ++k)
        et->addColumn<float>(coordNames_[k], coords[k], "coordinate " + std::to_string(k), coordPrecision_);
    et->addColumn<float>("weight",
                         weight,
                         Traits::symfactor ? "product of (E_p/E_jet)^n, without the symmetry factor"
                                           : "summed weight including symmetry factors",
                         weightPrecision_);
    if constexpr (Traits::symfactor)
      et->addColumn<std::uint32_t>("symfactor", sym, "multinomial symmetry factor");

    event.put(jetPut_, std::move(jt));
    event.put(entryPut_, std::move(et));
  }

  void globalBeginRunProduce(edm::Run& run, const edm::EventSetup&) const override {
    if constexpr (Traits::binned) {
      for (std::size_t k = 0; k < edges_.size(); ++k) {
        std::unique_ptr<nanoaod::FlatTable> t = std::make_unique<nanoaod::FlatTable>(edges_[k].size(), name_ + "Axis" + std::to_string(k), false);
        t->addColumn<double>("edge", edges_[k], "bin edges of axis " + std::to_string(k) + " (without flow bins)");
        run.put(axisPut_[k], std::move(t));
      }
      std::unique_ptr<nanoaod::FlatTable> t = std::make_unique<nanoaod::FlatTable>(types_.size(), name_ + "AxisType", false);
      t->addColumn<std::uint8_t>("type", types_, "axis type: 0 regular, 1 log, 2 variable");
      run.put(typePut_, std::move(t));
    }
  }

private:
  static std::size_t size(const Result& r) {
    if constexpr (std::is_same_v<Result, eec::DenseResult>)
      return r.values.size();
    else
      return r.entries.size();
  }

  const edm::EDGetTokenT<std::vector<Result>> src_;
  const std::string name_, doc_;
  const int weightPrecision_;
  int coordPrecision_ = -1;
  std::vector<std::string> coordNames_;
  std::vector<std::vector<double>> edges_;
  std::vector<std::uint8_t> types_;
  const edm::EDPutTokenT<nanoaod::FlatTable> jetPut_, entryPut_;
  std::vector<edm::EDPutTokenT<nanoaod::FlatTable>> axisPut_;
  edm::EDPutTokenT<nanoaod::FlatTable> typePut_;
};

#endif
