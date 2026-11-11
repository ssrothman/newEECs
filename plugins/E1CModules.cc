#include "FWCore/Framework/interface/MakerMacros.h"
#include "SRothman/EECs/plugins/E1CProducer.h"

// Module type E1CProducer_<Weight>_<Output>; add a combination with one line.
#define EEC_DEFINE_E1C(W, O)                                                         \
  using E1CProducer_##W##_##O = E1CProducer<eec::weight::W, eec::output::O>; \
  DEFINE_FWK_MODULE(E1CProducer_##W##_##O)

EEC_DEFINE_E1C(Energy, Dense);
EEC_DEFINE_E1C(Energy, SparseMerged);
EEC_DEFINE_E1C(Energy, SparseUnmerged);
EEC_DEFINE_E1C(Energy, Unbinned);
EEC_DEFINE_E1C(Pt, Dense);
EEC_DEFINE_E1C(Pt, SparseMerged);
EEC_DEFINE_E1C(Pt, SparseUnmerged);
EEC_DEFINE_E1C(Pt, Unbinned);
