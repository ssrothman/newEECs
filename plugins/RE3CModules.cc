#include "FWCore/Framework/interface/MakerMacros.h"
#include "SRothman/EECs/plugins/RE3CProducer.h"

// Module type RE3CProducer_<Weight>_<Angle>_<Output>; add a combination with one line.
#define EEC_DEFINE_RE3C(W, A, O)                                                                       \
  using RE3CProducer_##W##_##A##_##O = RE3CProducer<eec::weight::W, eec::angle::A, eec::output::O>; \
  DEFINE_FWK_MODULE(RE3CProducer_##W##_##A##_##O)

EEC_DEFINE_RE3C(Pt, DeltaR, Dense);
EEC_DEFINE_RE3C(Pt, DeltaR, SparseMerged);
EEC_DEFINE_RE3C(Pt, DeltaR, SparseUnmerged);
EEC_DEFINE_RE3C(Pt, DeltaR, Unbinned);
EEC_DEFINE_RE3C(Pt, Angle3D, Dense);
EEC_DEFINE_RE3C(Pt, Angle3D, SparseMerged);
EEC_DEFINE_RE3C(Pt, Angle3D, SparseUnmerged);
EEC_DEFINE_RE3C(Pt, Angle3D, Unbinned);
EEC_DEFINE_RE3C(Energy, DeltaR, Dense);
EEC_DEFINE_RE3C(Energy, DeltaR, SparseMerged);
EEC_DEFINE_RE3C(Energy, DeltaR, SparseUnmerged);
EEC_DEFINE_RE3C(Energy, DeltaR, Unbinned);
EEC_DEFINE_RE3C(Energy, Angle3D, Dense);
EEC_DEFINE_RE3C(Energy, Angle3D, SparseMerged);
EEC_DEFINE_RE3C(Energy, Angle3D, SparseUnmerged);
EEC_DEFINE_RE3C(Energy, Angle3D, Unbinned);
