#include "FWCore/Framework/interface/MakerMacros.h"
#include "SRothman/EECs/plugins/PENCProducer.h"

// Module type PENCProducer_<Weight>_<Angle>_<Output>; add a combination with one line.
#define EEC_DEFINE_PENC(W, A, O)                                                                       \
  using PENCProducer_##W##_##A##_##O = PENCProducer<eec::weight::W, eec::angle::A, eec::output::O>; \
  DEFINE_FWK_MODULE(PENCProducer_##W##_##A##_##O)

EEC_DEFINE_PENC(Pt, DeltaR, Dense);
EEC_DEFINE_PENC(Pt, DeltaR, SparseMerged);
EEC_DEFINE_PENC(Pt, DeltaR, SparseUnmerged);
EEC_DEFINE_PENC(Pt, DeltaR, Unbinned);
EEC_DEFINE_PENC(Pt, Angle3D, Dense);
EEC_DEFINE_PENC(Pt, Angle3D, SparseMerged);
EEC_DEFINE_PENC(Pt, Angle3D, SparseUnmerged);
EEC_DEFINE_PENC(Pt, Angle3D, Unbinned);
