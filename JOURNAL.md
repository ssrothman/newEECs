# EECs development journal

## Contents
- [2026-11-10 Initial implementation](#2026-11-10-initial-implementation)
- [2026-11-10 Benchmark: pENC strategies and outputs](#2026-11-10-benchmark-penc-strategies-and-outputs)
- [2026-11-10 Review fixes](#2026-11-10-review-fixes)
- [2026-11-10 Strategy factorisation and update_extra](#2026-11-10-strategy-factorisation-and-update_extra)
- [2026-11-11 Review round 2](#2026-11-11-review-round-2)

## 2026-11-10 Initial implementation
- Plan: `simulation/stepA.md` §A.1; decisions agreed with user are in README "Conventions".
- Tests: 24 Catch2 cases. Discrimination checked by 17 hand-applied mutations of the library
  (drop symmetry factor, sample without replacement, min instead of max, no phi wrap, acos angle,
  skip pow, skip boost, skip orthogonalisation, no scratch reset, merged/unmerged weight bugs,
  stride bug, histogram offset bug, wrong Unity norm, External check removed, wrong meta.N);
  every mutation made at least one assertion fail.

## 2026-11-10 Benchmark: pENC strategies and outputs
Setup: `bench/bench_penc.cc`, synthetic jets (Gaussian spread 0.15 in eta/phi, exponential pT),
Energy/SumEnergy/DeltaR, n=1, 50-bin log axis [1e-3, 1]. AMD EPYC, gcc 12.3.1 -O3, single thread,
`--benchmark_min_time=0.3s`. CPU time per jet:

| N | M | Naive | Multiset | PairAngles | PairBins |
|---|---|---|---|---|---|
| 2 | 50 | 48 us | 27 us | 20 us | 17 us |
| 2 | 150 | 423 us | 230 us | 166 us | 144 us |
| 3 | 50 | 3.2 ms | 0.74 ms | 0.43 ms | 0.060 ms |
| 3 | 100 | 25.6 ms | 5.8 ms | 3.3 ms | 0.41 ms |
| 4 | 50 | - | 11.6 ms | 6.8 ms | 0.74 ms |
| 5 | 50 | - | 126 ms | 72 ms | 8.8 ms |
| 6 | 50 | - | 1.21 s | 0.59 s | 0.13 s |

Outputs at N=3, PairBins (PairAngles for Unbinned), M=100: Dense 0.41 ms, SparseMerged 0.51 ms,
SparseMergedSort 4.4 ms, SparseUnmerged 0.40 ms, Unbinned 0.37 ms.

Findings:
- PairBins is ~8x faster than PairAngles for N>=3: per-leaf `axis::variant` index + log transform
  dominates PairAngles (Unbinned PairAngles, which skips that lookup, is 9x faster than Dense PairAngles).
- Scratch-array merging beats sort-reduce by ~8x at M=100.
- Default strategy is PairBins (binned) / PairAngles (unbinned).

## 2026-11-10 Review fixes
- boost accepts decreasing regular/log axes; with them PairBins was silently wrong (max |diff| 0.96 vs
  PairAngles, N=3, M=10). PairBins constructors now reject a decreasing axis 0; other strategies handle it.
- SparseMerged leaked state after a fill threw (next jet summed to 0.39, not 1); begin() now clears leftovers.
- Hoisting `kFactorial / den` out of the leaf loop: N=3 M=100 -16%, N=4 M=50 -10%, N=5 M=50 +8%. Not adopted.
- New tests checked by 7 mutations (each fix reverted or broken, plus over-broad axis check); each failed >= 1 test.

## 2026-11-10 Strategy factorisation and update_extra
- Decided with user: `Strategy<Loop, Precompute>` plus the four old names as aliases; PairBins with an unbinned
  output is a compile error (default strategy now `DefaultFor<Output>`); `update_extra` without a matching
  `Extra::update` is a compile error (static_assert), checked by `test/compile_fail` ctest cases.

## 2026-11-11 Review round 2
- Decided with user: `update_extra` renamed `updateExtra`; |p| = 0 particles throw in `fillJetData` (pENC and E1C);
  calculators stay `const` + `mutable` (callers must not share one across threads); stateless Weight/Norm/Angle kept.
- Row-max loop (incremental per-depth row max, one comparison per leaf), bitwise-identical output, median of 5,
  new/old time: PairBins N=3..5 1.55-1.83, N=6 0.90-1.02; Dense PairAngles N=3,4 1.00-1.20, N>=5 0.66-0.68;
  Unbinned PairAngles N=3 1.15, N>=4 0.93-0.98. PairAngles N=5 M=50 with row-max (40 ms) is still 5x PairBins (8 ms).
  Dropped: the leaf is not comparison-bound.
- Extra axes (K=1 constant Extra, Dense, PairBins): binning the Extra costs ~2.5x (N=3 M=100: ~1000 us vs 406 us).
  Tried `get_if` axis dispatch instead of `axis::variant::index`, for the Extra axes only (inlining it into the
  axis-0 lookup slowed N=3 PairAngles ~1.5x). Interleaved, core-pinned, min of 5: new/old 0.88-1.05 across
  N=2..5 and Regular/Log Extra axes, sign varying; rebuilding the same probe moved one ratio 0.85 -> 0.95, so code
  layout effects are as large as the change. Not adopted.
- Concepts (`interface/Concepts.h`) constrain all calculators; checked by `test/test_concepts.cc` and three
  compile_fail cases. Mutations: dropping each of 4 constraints/clauses fails >= 1 test. New zero-|p| and count tests:
  dropping the throw and a no-replacement loop each fail; a wrong `multichoose` (reserve only) passes, as intended.
