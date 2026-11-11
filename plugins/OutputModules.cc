#include "FWCore/Framework/interface/MakerMacros.h"
#include "SRothman/EECs/plugins/EECAggregator.h"
#include "SRothman/EECs/plugins/EECTableProducer.h"

// Names follow the output policy that produced the results (SparseMerged -> eec::SparseResult).
using EECTableProducer_Dense = EECTableProducer<eec::DenseResult>;
using EECTableProducer_SparseMerged = EECTableProducer<eec::SparseResult>;
using EECTableProducer_SparseUnmerged = EECTableProducer<eec::SparseUnmergedResult>;
using EECTableProducer_Unbinned1 = EECTableProducer<eec::UnbinnedResult<1>>;
using EECTableProducer_Unbinned2 = EECTableProducer<eec::UnbinnedResult<2>>;
DEFINE_FWK_MODULE(EECTableProducer_Dense);
DEFINE_FWK_MODULE(EECTableProducer_SparseMerged);
DEFINE_FWK_MODULE(EECTableProducer_SparseUnmerged);
DEFINE_FWK_MODULE(EECTableProducer_Unbinned1);
DEFINE_FWK_MODULE(EECTableProducer_Unbinned2);

using EECAggregator_Dense = EECAggregator<eec::DenseResult>;
using EECAggregator_SparseMerged = EECAggregator<eec::SparseResult>;
using EECAggregator_SparseUnmerged = EECAggregator<eec::SparseUnmergedResult>;
DEFINE_FWK_MODULE(EECAggregator_Dense);
DEFINE_FWK_MODULE(EECAggregator_SparseMerged);
DEFINE_FWK_MODULE(EECAggregator_SparseUnmerged);
