"""Helpers to add EEC producers, NanoAOD tables and aggregators to a cms.Process from observable specs.

Every value comes from the caller's spec (see simulation/configs/default.yaml in boostedV); there are no physics
defaults here. Module types are registered in plugins/*Modules.cc as
PENCProducer_<Weight>_<Angle>_<Output> and E1CProducer_<Weight>_<Output>.

Observable spec keys:
  kind: penc | e1c;  n, weight, output, norm, axes (binned outputs only; a list of axis specs)
  penc: N, angle;  e1c: frame
  table: {weightPrecision, coordPrecision, coordNames (unbinned only)}
Axis spec: {type: regular|log, nbins, lo, hi} or {type: variable, edges}.
"""
import FWCore.ParameterSet.Config as cms

OUTPUTS = ('Dense', 'SparseMerged', 'SparseUnmerged', 'Unbinned')


def axisPSet(a):
    if a['type'] == 'variable':
        return cms.PSet(type=cms.string('variable'), nbins=cms.uint32(0), lo=cms.double(0), hi=cms.double(0),
                        edges=cms.vdouble(*a['edges']))
    return cms.PSet(type=cms.string(a['type']), nbins=cms.uint32(a['nbins']), lo=cms.double(a['lo']),
                    hi=cms.double(a['hi']), edges=cms.vdouble())


def _rank(spec):
    return 1 if spec['kind'] == 'penc' else 2


def _check(name, spec):
    if spec['output'] not in OUTPUTS:
        raise ValueError('%s: output must be one of %s, got %s' % (name, OUTPUTS, spec['output']))
    if spec['kind'] not in ('penc', 'e1c'):
        raise ValueError('%s: kind must be penc or e1c' % name)
    binned = spec['output'] != 'Unbinned'
    if binned and len(spec.get('axes', [])) != _rank(spec):
        raise ValueError('%s: output %s needs %d axes' % (name, spec['output'], _rank(spec)))
    if not binned and 'axes' in spec:
        raise ValueError('%s: unbinned output takes no axes' % name)


def addObservable(process, label, src, spec, tableName, task=None):
    """Producer `label` for the EECJetCollection `src` and its NanoAOD tables `<label>Table` (the suffix matches the
    default NANO keep rule), named tableName / tableName + 'Entries'."""
    _check(label, spec)
    binned = spec['output'] != 'Unbinned'
    if spec['kind'] == 'penc':
        p = cms.EDProducer('PENCProducer_%s_%s_%s' % (spec['weight'], spec['angle'], spec['output']),
                           src=cms.InputTag(src), N=cms.uint32(spec['N']), n=cms.double(spec['n']),
                           norm=cms.string(spec['norm']))
    else:
        p = cms.EDProducer('E1CProducer_%s_%s' % (spec['weight'], spec['output']), src=cms.InputTag(src),
                           n=cms.double(spec['n']), norm=cms.string(spec['norm']), frame=cms.string(spec['frame']))
    t = spec['table']
    if binned:
        axes = [axisPSet(a) for a in spec['axes']]
        p.axes = cms.VPSet(*axes)
        tab = cms.EDProducer('EECTableProducer_' + spec['output'], src=cms.InputTag(label), name=cms.string(tableName),
                             doc=cms.string(''), weightPrecision=cms.int32(t['weightPrecision']), axes=cms.VPSet(*axes))
    else:
        tab = cms.EDProducer('EECTableProducer_Unbinned%d' % _rank(spec), src=cms.InputTag(label),
                             name=cms.string(tableName), doc=cms.string(''),
                             weightPrecision=cms.int32(t['weightPrecision']), coordNames=cms.vstring(*t['coordNames']),
                             coordPrecision=cms.int32(t['coordPrecision']))
    for lab, mod in ((label, p), (label + 'Table', tab)):
        setattr(process, lab, mod)
        if task is not None:
            task.add(mod)
    return p


def addAggregator(process, label, src, spec, genEventInfo, maxBins, weights=()):
    """EECAggregator for the binned results `src` produced from `spec`. Put it on a Path (e.g. after a channel
    filter); its TFileService directory is `label`. genEventInfo '' means unit weights (data); `weights` are
    InputTags of extra per-event double weights."""
    _check(label, spec)
    if spec['output'] == 'Unbinned':
        raise ValueError('%s: aggregation needs a binned output' % label)
    agg = cms.EDAnalyzer('EECAggregator_' + spec['output'], src=cms.InputTag(src),
                         genEventInfo=cms.InputTag(genEventInfo), weights=cms.VInputTag(*weights),
                         maxBins=cms.uint32(maxBins),
                         axes=cms.VPSet(*[axisPSet(a) for a in spec['axes']]))
    setattr(process, label, agg)
    return agg
