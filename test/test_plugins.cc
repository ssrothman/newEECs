#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/TestProcessor/interface/TestProcessor.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SRothman/DataFormats/interface/EECConversion.h"
#include "SRothman/DataFormats/interface/EECJet.h"
#include "SRothman/EECs/interface/Config.h"
#include "SRothman/EECs/interface/EECs.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNamed.h"
#include "TParameter.h"

#include <cmath>
#include <cstdio>
#include <unistd.h>
#include <random>
#include <string>
#include <vector>

using srothman::EECJet;
using srothman::EECJetCollection;

namespace {

  const char* kHeader = R"_(from FWCore.TestProcessor.TestProcess import *
process = TestProcess()
)_";

  EECJet randomJet(unsigned seed, std::size_t M) {
    std::mt19937 rng(seed);
    std::normal_distribution<double> d(0.0, 0.2);
    std::exponential_distribution<double> pt(0.2);
    std::uniform_int_distribution<int> q(-1, 1);
    EECJet j;
    j.pt = 350;
    j.eta = 0.4;
    j.phi = 1.1;
    j.mass = 85;
    for (std::size_t i = 0; i < M; ++i) {
      const auto p = fastjet::PtYPhiM(0.5 + pt(rng), 0.4 + d(rng), 1.1 + d(rng), 0.0);
      j.constituents.push_back({p.px(), p.py(), p.pz(), p.E(), static_cast<std::int8_t>(q(rng)), 211, 1.0f, 0});
    }
    return j;
  }

  double sumPt(const EECJet& j) {
    double s = 0;
    for (const auto& c : j.constituents)
      s += std::hypot(c.px, c.py);
    return s;
  }

  // Jet energy from (pt, eta, m), independent of jetPseudoJet.
  double jetEnergy(const EECJet& j) {
    const double p = j.pt * std::cosh(j.eta);
    return std::sqrt(p * p + j.mass * j.mass);
  }

  const char* kLogAxis =
      "axes = cms.VPSet(cms.PSet(type = cms.string('log'), nbins = cms.uint32(20), lo = cms.double(1e-3), "
      "hi = cms.double(1.0), edges = cms.vdouble()))";

  double col(const nanoaod::FlatTable& t, const std::string& name, unsigned row) {
    const int c = t.columnIndex(name);
    REQUIRE(c >= 0);
    return t.getAnyValue(row, c);
  }

}  // namespace

TEST_CASE("PENCProducer Dense output equals a direct library call", "[plugins][penc]") {
  edm::test::TestProcessor::Config config{std::string(kHeader) +
                                          "process.toTest = cms.EDProducer('PENCProducer_Pt_DeltaR_Dense', src = "
                                          "cms.InputTag('jets'), N = cms.uint32(3), n = cms.double(1.0), norm = "
                                          "cms.string('sumPt'), " +
                                          kLogAxis + ")\nprocess.moduleToTest(process.toTest)\n"};
  auto put = config.produces<EECJetCollection>("jets");
  edm::test::TestProcessor tester(config);

  EECJetCollection jets{randomJet(1, 25), EECJet{}, randomJet(2, 12)};
  jets[1].pt = 200;
  auto ev = tester.test(std::make_pair(put, std::make_unique<EECJetCollection>(jets)));
  const auto& out = *ev.get<std::vector<eec::DenseResult>>();
  REQUIRE(out.size() == 3);

  const eec::Binning b({eec::LogAxis(20, 1e-3, 1.0)});
  eec::ProjectedENC<3, eec::weight::Pt, eec::norm::External, eec::angle::DeltaR, eec::output::Dense> ref(1.0, b);
  for (std::size_t i : {0u, 2u}) {
    const double E = sumPt(jets[i]);
    const auto r = ref(srothman::toPseudoJets(jets[i]), E);
    CHECK(out[i].meta.E_jet == Approx(E).epsilon(1e-12));
    CHECK(out[i].meta.N == 3);
    CHECK(out[i].meta.n == 1.0);
    REQUIRE(out[i].values.size() == r.values.size());
    double total = 0;
    for (std::size_t k = 0; k < r.values.size(); ++k) {
      CHECK(out[i].values[k] == Approx(r.values[k]).margin(1e-14));
      total += out[i].values[k];
    }
    CHECK(total == Approx(1.0));  // sum rule for n = 1 with E_jet = sum of the weights
  }
  // The empty jet: all bins zero, E_jet = sum pt = 0.
  CHECK(out[1].meta.E_jet == 0.0);
  CHECK(out[1].values == std::vector<double>(b.totalBins(), 0.0));
}

TEST_CASE("PENCProducer Unbinned output with norm jetE and n = 2", "[plugins][penc]") {
  edm::test::TestProcessor::Config config{std::string(kHeader) +
                                          "process.toTest = cms.EDProducer('PENCProducer_Energy_Angle3D_Unbinned', "
                                          "src = cms.InputTag('jets'), N = cms.uint32(2), n = cms.double(2.0), norm = "
                                          "cms.string('jetE'))\nprocess.moduleToTest(process.toTest)\n"};
  auto put = config.produces<EECJetCollection>("jets");
  edm::test::TestProcessor tester(config);
  EECJetCollection jets{randomJet(3, 9)};
  auto ev = tester.test(std::make_pair(put, std::make_unique<EECJetCollection>(jets)));
  const auto& out = *ev.get<std::vector<eec::UnbinnedResult<1>>>();
  REQUIRE(out.size() == 1);
  const double E = jetEnergy(jets[0]);
  CHECK(out[0].meta.E_jet == Approx(E).epsilon(1e-12));
  CHECK(out[0].meta.N == 2);
  CHECK(out[0].meta.n == 2.0);
  eec::ProjectedENC<2, eec::weight::Energy, eec::norm::External, eec::angle::Angle3D, eec::output::Unbinned> ref(2.0);
  const auto r = ref(srothman::toPseudoJets(jets[0]), E);
  REQUIRE(out[0].entries.size() == r.entries.size());
  REQUIRE(r.entries.size() == 9 * 10 / 2);
  for (std::size_t k = 0; k < r.entries.size(); ++k) {
    CHECK(out[0].entries[k].coord[0] == r.entries[k].coord[0]);
    CHECK(out[0].entries[k].weight == Approx(r.entries[k].weight).epsilon(1e-12));
    CHECK(out[0].entries[k].symfactor == r.entries[k].symfactor);
  }
}

TEST_CASE("E1CProducer lab frame gives lab angles and pt weights", "[plugins][e1c]") {
  edm::test::TestProcessor::Config config{std::string(kHeader) +
                                          "process.toTest = cms.EDProducer('E1CProducer_Pt_Unbinned', src = "
                                          "cms.InputTag('jets'), n = cms.double(1.0), norm = cms.string('sumPt'), "
                                          "frame = cms.string('lab'))\nprocess.moduleToTest(process.toTest)\n"};
  auto put = config.produces<EECJetCollection>("jets");
  edm::test::TestProcessor tester(config);
  EECJetCollection jets{randomJet(4, 7)};
  auto ev = tester.test(std::make_pair(put, std::make_unique<EECJetCollection>(jets)));
  const auto& out = *ev.get<std::vector<eec::UnbinnedResult<2>>>();
  REQUIRE(out.size() == 1);
  REQUIRE(out[0].entries.size() == 7);
  CHECK(out[0].meta.N == 1);
  const double S = sumPt(jets[0]);
  for (std::size_t i = 0; i < 7; ++i) {
    const auto& c = jets[0].constituents[i];
    const fastjet::PseudoJet p(c.px, c.py, c.pz, c.E);
    CHECK(out[0].entries[i].coord[0] == Approx(p.theta()).epsilon(1e-12));
    CHECK(out[0].entries[i].coord[1] == Approx(std::atan2(c.py, c.px)).epsilon(1e-12));
    CHECK(out[0].entries[i].weight == Approx(p.pt() / S).epsilon(1e-12));
  }
}

TEST_CASE("E1CProducer jetCS frame matches the library frame", "[plugins][e1c]") {
  edm::test::TestProcessor::Config config{
      std::string(kHeader) +
      "process.toTest = cms.EDProducer('E1CProducer_Energy_SparseMerged', src = cms.InputTag('jets'), n = "
      "cms.double(1.0), norm = cms.string('sumE'), frame = cms.string('jetCS'), axes = cms.VPSet("
      "cms.PSet(type = cms.string('regular'), nbins = cms.uint32(8), lo = cms.double(0), hi = cms.double(3.15), edges "
      "= cms.vdouble()), cms.PSet(type = cms.string('variable'), nbins = cms.uint32(0), lo = cms.double(0), hi = "
      "cms.double(0), edges = cms.vdouble(-3.15, -1, 0, 1, 3.15))))\nprocess.moduleToTest(process.toTest)\n"};
  auto put = config.produces<EECJetCollection>("jets");
  edm::test::TestProcessor tester(config);
  EECJetCollection jets{randomJet(5, 30)};
  auto ev = tester.test(std::make_pair(put, std::make_unique<EECJetCollection>(jets)));
  const auto& out = *ev.get<std::vector<eec::SparseResult>>();
  REQUIRE(out.size() == 1);
  const eec::Binning b({eec::RegularAxis(8, 0, 3.15), eec::VariableAxis(std::vector<double>{-3.15, -1, 0, 1, 3.15})});
  eec::OnePointEC<eec::weight::Energy, eec::norm::External, eec::frame::BoostRotate, eec::output::SparseMerged> ref(1.0, b);
  const auto ps = srothman::toPseudoJets(jets[0]);
  double S = 0;
  for (const auto& p : ps)
    S += p.E();
  const auto frame = eec::config::makeFrame(eec::config::FrameKind::JetCS, srothman::jetPseudoJet(jets[0]));
  const auto r = ref(ps, frame, S);
  REQUIRE(out[0].entries.size() == r.entries.size());
  for (std::size_t k = 0; k < r.entries.size(); ++k) {
    CHECK(out[0].entries[k].bin == r.entries[k].bin);
    CHECK(out[0].entries[k].weight == Approx(r.entries[k].weight).epsilon(1e-12));
  }
  // A different frame bins differently for this jet, so the frame choice is visible in the output.
  const auto rh = ref(ps, eec::config::makeFrame(eec::config::FrameKind::JetHelicity, srothman::jetPseudoJet(jets[0])), S);
  bool differs = rh.entries.size() != r.entries.size();
  for (std::size_t k = 0; !differs && k < r.entries.size(); ++k)
    differs = rh.entries[k].bin != r.entries[k].bin;
  CHECK(differs);
}

TEST_CASE("EEC producers reject bad configuration", "[plugins]") {
  const std::string badNorm = std::string(kHeader) +
                              "process.toTest = cms.EDProducer('PENCProducer_Pt_DeltaR_SparseMerged', src = "
                              "cms.InputTag('jets'), N = cms.uint32(2), n = cms.double(1.0), norm = cms.string('pt'), " +
                              kLogAxis + ")\nprocess.moduleToTest(process.toTest)\n";
  const std::string badN = std::string(kHeader) +
                           "process.toTest = cms.EDProducer('PENCProducer_Pt_DeltaR_SparseMerged', src = "
                           "cms.InputTag('jets'), N = cms.uint32(7), n = cms.double(1.0), norm = cms.string('sumPt'), " +
                           kLogAxis + ")\nprocess.moduleToTest(process.toTest)\n";
  const std::string badFrame = std::string(kHeader) +
                               "process.toTest = cms.EDProducer('E1CProducer_Pt_Unbinned', src = cms.InputTag('jets'), "
                               "n = cms.double(1.0), norm = cms.string('sumPt'), frame = cms.string('CS'))\n"
                               "process.moduleToTest(process.toTest)\n";
  for (const auto& cfg : {badNorm, badN, badFrame}) {
    edm::test::TestProcessor::Config config{cfg};
    config.produces<EECJetCollection>("jets");
    CHECK_THROWS_AS(edm::test::TestProcessor(config), cms::Exception);
  }
}

TEST_CASE("EECTableProducer writes per-jet, entry and run tables", "[plugins][table]") {
  edm::test::TestProcessor::Config config{
      std::string(kHeader) +
      "process.toTest = cms.EDProducer('EECTableProducer_SparseUnmerged', src = cms.InputTag('res'), name = "
      "cms.string('E3C'), doc = cms.string('test'), weightPrecision = cms.int32(-1), axes = cms.VPSet("
      "cms.PSet(type = cms.string('regular'), nbins = cms.uint32(4), lo = cms.double(0), hi = cms.double(2), edges = "
      "cms.vdouble()), cms.PSet(type = cms.string('log'), nbins = cms.uint32(2), lo = cms.double(0.01), hi = "
      "cms.double(1), edges = cms.vdouble())))\nprocess.moduleToTest(process.toTest)\n"};
  auto put = config.produces<std::vector<eec::SparseUnmergedResult>>("res");
  edm::test::TestProcessor tester(config);
  {
    // Run tables are produced at begin run.
    auto run = tester.testBeginRun(1);
    const auto& a0 = *run.get<nanoaod::FlatTable>("axis0");
    const auto& a1 = *run.get<nanoaod::FlatTable>("axis1");
    const auto& ty = *run.get<nanoaod::FlatTable>("axisTypes");
    CHECK(a0.name() == "E3CAxis0");
    REQUIRE(a0.nRows() == 5);
    CHECK(col(a0, "edge", 0) == 0.0);
    CHECK(col(a0, "edge", 4) == 2.0);
    REQUIRE(a1.nRows() == 3);
    CHECK(col(a1, "edge", 1) == Approx(0.1));
    REQUIRE(ty.nRows() == 2);
    CHECK(col(ty, "type", 0) == 0);
    CHECK(col(ty, "type", 1) == 1);
  }

  std::vector<eec::SparseUnmergedResult> res(2);
  res[0].meta = {120.5, 1.0, 3};
  res[0].entries = {{3, 0.25, 6}, {3, 0.125, 1}, {7, 0.5, 3}};
  res[1].meta = {80.0, 1.0, 3};
  auto ev = tester.test(std::make_pair(put, std::make_unique<std::vector<eec::SparseUnmergedResult>>(res)));
  const auto& jt = *ev.get<nanoaod::FlatTable>();
  const auto& et = *ev.get<nanoaod::FlatTable>("entries");
  CHECK(jt.name() == "E3C");
  REQUIRE(jt.nRows() == 2);
  CHECK(col(jt, "E_jet", 0) == 120.5);
  CHECK(col(jt, "E_jet", 1) == 80.0);
  CHECK(col(jt, "N", 0) == 3);
  CHECK(col(jt, "nEntries", 0) == 3);
  CHECK(col(jt, "nEntries", 1) == 0);
  CHECK(et.name() == "E3CEntries");
  REQUIRE(et.nRows() == 3);
  CHECK(col(et, "iJet", 2) == 0);
  CHECK(col(et, "bin", 2) == 7);
  CHECK(col(et, "weight", 1) == 0.125);
  CHECK(col(et, "symfactor", 0) == 6);
  CHECK(col(et, "symfactor", 1) == 1);

}

TEST_CASE("EECTableProducer for 2D unbinned results names the coordinates", "[plugins][table]") {
  edm::test::TestProcessor::Config config{std::string(kHeader) +
                                          "process.toTest = cms.EDProducer('EECTableProducer_Unbinned2', src = "
                                          "cms.InputTag('res'), name = cms.string('E1C'), coordNames = "
                                          "cms.vstring('theta', 'phi'))\nprocess.moduleToTest(process.toTest)\n"};
  auto put = config.produces<std::vector<eec::UnbinnedResult<2>>>("res");
  edm::test::TestProcessor tester(config);
  std::vector<eec::UnbinnedResult<2>> res(2);
  res[1].meta = {10.0, 1.0, 1};
  res[1].entries = {{{{0.5, -1.25}}, 0.75, 1}, {{{2.0, 3.0}}, 0.25, 1}};
  auto ev = tester.test(std::make_pair(put, std::make_unique<std::vector<eec::UnbinnedResult<2>>>(res)));
  const auto& et = *ev.get<nanoaod::FlatTable>("entries");
  REQUIRE(et.nRows() == 2);
  CHECK(col(et, "iJet", 0) == 1);
  CHECK(col(et, "theta", 0) == 0.5);
  CHECK(col(et, "phi", 0) == -1.25);
  CHECK(col(et, "theta", 1) == 2.0);
  CHECK(col(et, "weight", 1) == 0.25);
  CHECK(et.columnIndex("bin") == -1);
}

TEST_CASE("EECAggregator writes H, C and weight sums", "[plugins][aggregator]") {
  const std::string file = "testEECAggregator_" + std::to_string(::getpid()) + ".root";
  {
    edm::test::TestProcessor::Config config{
        std::string(kHeader) +
        "process.toTest = cms.EDAnalyzer('EECAggregator_SparseMerged', src = cms.InputTag('res'), "
        "genEventInfo = cms.InputTag('generator'), maxBins = cms.uint32(100), axes = cms.VPSet("
        "cms.PSet(type = cms.string('regular'), nbins = cms.uint32(2), lo = cms.double(0), hi = cms.double(1), edges "
        "= cms.vdouble())))\nprocess.moduleToTest(process.toTest)\nprocess.add_(cms.Service('JobReportService'))\nprocess.add_(cms.Service('TFileService', fileName "
        "= cms.string('" +
        file + "')))\n"};
    auto putRes = config.produces<std::vector<eec::SparseResult>>("res");
    auto putGen = config.produces<GenEventInfoProduct>("generator");
    edm::test::TestProcessor tester(config);

    auto gen = [](double w) {
      auto g = std::make_unique<GenEventInfoProduct>();
      g->setWeights({w});
      return g;
    };
    // 4 bins (2 + flow). Event 1 (w = 2): two jets, both hitting bin 1, so h = (0, 0.5 + 0.25, 0.5, 0).
    std::vector<eec::SparseResult> e1(2);
    e1[0].entries = {{1, 0.5}, {2, 0.5}};
    e1[1].entries = {{1, 0.25}};
    tester.test(std::make_pair(putRes, std::make_unique<std::vector<eec::SparseResult>>(e1)),
                std::make_pair(putGen, gen(2.0)));
    // Event 2 (w = -1): h = (0, 0, 0, 1).
    std::vector<eec::SparseResult> e2(1);
    e2[0].entries = {{3, 1.0}};
    tester.test(std::make_pair(putRes, std::make_unique<std::vector<eec::SparseResult>>(e2)),
                std::make_pair(putGen, gen(-1.0)));
    // Event 3 (w = 3): no jets.
    tester.test(std::make_pair(putRes, std::make_unique<std::vector<eec::SparseResult>>()),
                std::make_pair(putGen, gen(3.0)));
  }
  TFile f(file.c_str(), "READ");
  REQUIRE(f.IsOpen());
  auto* values = f.Get<TH1D>("toTest/values");
  auto* cov = f.Get<TH2D>("toTest/covariance");
  REQUIRE(values);
  REQUIRE(cov);
  REQUIRE(values->GetNbinsX() == 4);
  CHECK(values->GetBinContent(1) == 0.0);
  CHECK(values->GetBinContent(2) == Approx(1.5));  // 2 * 0.75
  CHECK(values->GetBinContent(3) == Approx(1.0));  // 2 * 0.5
  CHECK(values->GetBinContent(4) == Approx(-1.0));
  // C = 4 h1 h1^T + 1 h2 h2^T
  CHECK(cov->GetBinContent(2, 2) == Approx(4 * 0.75 * 0.75));
  CHECK(cov->GetBinContent(2, 3) == Approx(4 * 0.75 * 0.5));
  CHECK(cov->GetBinContent(3, 2) == Approx(4 * 0.75 * 0.5));
  CHECK(cov->GetBinContent(4, 4) == Approx(1.0));
  CHECK(cov->GetBinContent(2, 4) == 0.0);
  CHECK(values->GetBinError(2) == Approx(1.5));
  CHECK(f.Get<TParameter<double>>("toTest/sumW")->GetVal() == Approx(4.0));
  CHECK(f.Get<TParameter<double>>("toTest/sumW2")->GetVal() == Approx(14.0));
  CHECK(f.Get<TParameter<Long64_t>>("toTest/nEvents")->GetVal() == 3);
  auto* ax = f.Get<TH1D>("toTest/axis0");
  REQUIRE(ax);
  REQUIRE(ax->GetNbinsX() == 2);
  CHECK(ax->GetXaxis()->GetBinLowEdge(1) == 0.0);
  CHECK(ax->GetXaxis()->GetBinUpEdge(2) == 1.0);
  CHECK(std::string(f.Get<TNamed>("toTest/axisTypes")->GetTitle()) == "regular");
  f.Close();
  std::remove(file.c_str());
}

TEST_CASE("EECAggregator rejects a dense result of the wrong size", "[plugins][aggregator]") {
  const std::string file = "testEECAggregatorBad_" + std::to_string(::getpid()) + ".root";
  edm::test::TestProcessor::Config config{
      std::string(kHeader) +
      "process.toTest = cms.EDAnalyzer('EECAggregator_Dense', src = cms.InputTag('res'), "
      "genEventInfo = cms.InputTag(''), maxBins = cms.uint32(100), axes = cms.VPSet(cms.PSet(type = "
      "cms.string('regular'), nbins = cms.uint32(2), lo = cms.double(0), hi = cms.double(1), edges = "
      "cms.vdouble())))\nprocess.moduleToTest(process.toTest)\nprocess.add_(cms.Service('JobReportService'))\nprocess.add_(cms.Service('TFileService', fileName = "
      "cms.string('" +
      file + "')))\n"};
  auto put = config.produces<std::vector<eec::DenseResult>>("res");
  {
    edm::test::TestProcessor tester(config);
    std::vector<eec::DenseResult> r(1);
    r[0].values.assign(5, 0.1);
    CHECK_THROWS_AS(tester.test(std::make_pair(put, std::make_unique<std::vector<eec::DenseResult>>(r))),
                    cms::Exception);
  }
  std::remove(file.c_str());
}

TEST_CASE("EECAggregator multiplies unmerged weights by their symmetry factors", "[plugins][aggregator]") {
  const std::string file = "testEECAggregatorUnm_" + std::to_string(::getpid()) + ".root";
  {
    edm::test::TestProcessor::Config config{
        std::string(kHeader) +
        "process.toTest = cms.EDAnalyzer('EECAggregator_SparseUnmerged', src = cms.InputTag('res'), genEventInfo = "
        "cms.InputTag(''), maxBins = cms.uint32(100), axes = cms.VPSet(cms.PSet(type = cms.string('regular'), nbins "
        "= cms.uint32(2), lo = cms.double(0), hi = cms.double(1), edges = cms.vdouble())))\n"
        "process.moduleToTest(process.toTest)\nprocess.add_(cms.Service('JobReportService'))\n"
        "process.add_(cms.Service('TFileService', fileName = cms.string('" +
        file + "')))\n"};
    auto put = config.produces<std::vector<eec::SparseUnmergedResult>>("res");
    edm::test::TestProcessor tester(config);
    // h = (0, 0.1*6 + 0.2*3, 0, 0) = (0, 1.2, 0, 0); unit weight.
    std::vector<eec::SparseUnmergedResult> r(1);
    r[0].entries = {{1, 0.1, 6}, {1, 0.2, 3}};
    tester.test(std::make_pair(put, std::make_unique<std::vector<eec::SparseUnmergedResult>>(r)));
  }
  TFile f(file.c_str(), "READ");
  auto* values = f.Get<TH1D>("toTest/values");
  auto* cov = f.Get<TH2D>("toTest/covariance");
  REQUIRE(values);
  REQUIRE(cov);
  CHECK(values->GetBinContent(2) == Approx(1.2));
  CHECK(cov->GetBinContent(2, 2) == Approx(1.44));
  CHECK(f.Get<TParameter<double>>("toTest/sumW")->GetVal() == 1.0);
  f.Close();
  std::remove(file.c_str());
}

TEST_CASE("EECAggregator multiplies the extra weights into the generator weight", "[plugins][aggregator]") {
  const std::string file = "testEECAggregatorW_" + std::to_string(::getpid()) + ".root";
  {
    edm::test::TestProcessor::Config config{
        std::string(kHeader) +
        "process.toTest = cms.EDAnalyzer('EECAggregator_SparseMerged', src = cms.InputTag('res'), genEventInfo = "
        "cms.InputTag('generator'), weights = cms.VInputTag('sfA', 'sfB'), maxBins = cms.uint32(100), axes = "
        "cms.VPSet(cms.PSet(type = cms.string('regular'), nbins = cms.uint32(2), lo = cms.double(0), hi = "
        "cms.double(1), edges = cms.vdouble())))\nprocess.moduleToTest(process.toTest)\n"
        "process.add_(cms.Service('JobReportService'))\n"
        "process.add_(cms.Service('TFileService', fileName = cms.string('" +
        file + "')))\n"};
    auto putRes = config.produces<std::vector<eec::SparseResult>>("res");
    auto putGen = config.produces<GenEventInfoProduct>("generator");
    auto putA = config.produces<double>("sfA");
    auto putB = config.produces<double>("sfB");
    edm::test::TestProcessor tester(config);
    auto gen = std::make_unique<GenEventInfoProduct>();
    gen->setWeights({2.0});
    // w = 2 * 0.9 * 1.5 = 2.7; h = (0, 0.5, 0, 0).
    std::vector<eec::SparseResult> r(1);
    r[0].entries = {{1, 0.5}};
    tester.test(std::make_pair(putRes, std::make_unique<std::vector<eec::SparseResult>>(r)),
                std::make_pair(putGen, std::move(gen)), std::make_pair(putA, std::make_unique<double>(0.9)),
                std::make_pair(putB, std::make_unique<double>(1.5)));
  }
  TFile f(file.c_str(), "READ");
  auto* values = f.Get<TH1D>("toTest/values");
  auto* cov = f.Get<TH2D>("toTest/covariance");
  REQUIRE(values);
  REQUIRE(cov);
  CHECK(values->GetBinContent(2) == Approx(2.7 * 0.5));
  CHECK(cov->GetBinContent(2, 2) == Approx(2.7 * 2.7 * 0.25));
  CHECK(f.Get<TParameter<double>>("toTest/sumW")->GetVal() == Approx(2.7));
  f.Close();
  std::remove(file.c_str());
}
