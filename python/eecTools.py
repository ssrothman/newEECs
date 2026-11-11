"""Helpers to add EEC producers, NanoAOD tables and aggregators to a cms.Process.

Module types are registered in plugins/*Modules.cc as
PENCProducer_<Weight>_<Angle>_<Output> and E1CProducer_<Weight>_<Output>.
"""
import FWCore.ParameterSet.Config as cms

OUTPUTS = ('Dense', 'SparseMerged', 'SparseUnmerged', 'Unbinned')


def regularAxis(nbins, lo, hi):
    return cms.PSet(type=cms.string('regular'), nbins=cms.uint32(nbins), lo=cms.double(lo),
                    hi=cms.double(hi), edges=cms.vdouble())


def logAxis(nbins, lo, hi):
    return cms.PSet(type=cms.string('log'), nbins=cms.uint32(nbins), lo=cms.double(lo),
                    hi=cms.double(hi), edges=cms.vdouble())


def variableAxis(edges):
    return cms.PSet(type=cms.string('variable'), nbins=cms.uint32(0), lo=cms.double(0),
                    hi=cms.double(0), edges=cms.vdouble(*edges))


def _finish(process, label, producer, output, axes, D, table, tableName, coordNames,
            aggregate, genEventInfo, task):
    """Attach producer, optional table (label 'Table' suffix, so the default NANO keep rule matches)
    and optional aggregator. Returns the aggregator, which must be put on a Path or EndPath."""
    binned = output != 'Unbinned'
    setattr(process, label, producer)
    if task is not None:
        task.add(producer)
    if table:
        if binned:
            tab = cms.EDProducer('EECTableProducer_' + output, src=cms.InputTag(label),
                                 name=cms.string(tableName or label), doc=cms.string(''),
                                 weightPrecision=cms.int32(-1), axes=cms.VPSet(*axes))
        else:
            tab = cms.EDProducer('EECTableProducer_Unbinned%d' % D, src=cms.InputTag(label),
                                 name=cms.string(tableName or label), doc=cms.string(''),
                                 weightPrecision=cms.int32(-1),
                                 coordNames=cms.vstring(*coordNames), coordPrecision=cms.int32(-1))
        setattr(process, label + 'Table', tab)
        if task is not None:
            task.add(tab)
    if not aggregate:
        return None
    agg = cms.EDAnalyzer('EECAggregator_' + output, src=cms.InputTag(label),
                         genEventInfo=cms.InputTag(genEventInfo), maxBins=cms.uint32(5000),
                         axes=cms.VPSet(*axes))
    setattr(process, label + 'Aggregator', agg)
    return agg


def _check(output, axes, rank, aggregate):
    if output not in OUTPUTS:
        raise ValueError('output must be one of %s, got %s' % (OUTPUTS, output))
    if output != 'Unbinned' and (axes is None or len(axes) != rank):
        raise ValueError('output %s needs %d axes' % (output, rank))
    if aggregate and output == 'Unbinned':
        raise ValueError('aggregation needs a binned output, got Unbinned')


def addPENC(process, label, src, N, n=1.0, weight='Pt', angle='DeltaR', output='SparseMerged',
            norm='sumPt', axes=None, table=True, tableName=None, aggregate=False,
            genEventInfo='generator', task=None):
    """Projected N-point correlator of the EECJetCollection `src`. axes: [one axis PSet] for binned outputs.
    genEventInfo='' aggregates with unit weights (data)."""
    _check(output, axes, 1, aggregate)
    p = cms.EDProducer('PENCProducer_%s_%s_%s' % (weight, angle, output), src=cms.InputTag(src),
                       N=cms.uint32(N), n=cms.double(n), norm=cms.string(norm))
    if output != 'Unbinned':
        p.axes = cms.VPSet(*axes)
    return _finish(process, label, p, output, axes, 1, table, tableName, ['angle'],
                   aggregate, genEventInfo, task)


def addE1C(process, label, src, n=1.0, weight='Pt', output='SparseMerged', norm='sumPt',
           frame='jetHelicity', axes=None, table=True, tableName=None, aggregate=False,
           genEventInfo='generator', task=None):
    """One-point correlator (theta, phi) in `frame` (lab, jetHelicity, jetCS). axes: [theta, phi] PSets."""
    _check(output, axes, 2, aggregate)
    p = cms.EDProducer('E1CProducer_%s_%s' % (weight, output), src=cms.InputTag(src),
                       n=cms.double(n), norm=cms.string(norm), frame=cms.string(frame))
    if output != 'Unbinned':
        p.axes = cms.VPSet(*axes)
    return _finish(process, label, p, output, axes, 2, table, tableName, ['theta', 'phi'],
                   aggregate, genEventInfo, task)
