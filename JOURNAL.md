# EECs development journal

## Contents
- [2026-11-10 Initial implementation](#2026-11-10-initial-implementation)
- [2026-11-10 Benchmark: pENC strategies and outputs](#2026-11-10-benchmark-penc-strategies-and-outputs)
- [2026-11-10 Review fixes](#2026-11-10-review-fixes)
- [2026-11-10 Strategy factorisation and update_extra](#2026-11-10-strategy-factorisation-and-update_extra)
- [2026-11-11 Review round 2](#2026-11-11-review-round-2)
- [2026-11-11 Review round 3: de-duplication](#2026-11-11-review-round-3-de-duplication)
- [2026-11-11 Inner-loop optimisation study](#2026-11-11-inner-loop-optimisation-study)
- [2026-11-11 Peeled last level, rotating-jet benchmark](#2026-11-11-peeled-last-level-rotating-jet-benchmark)
- [2026-11-11 CovarianceAccumulator](#2026-11-11-covarianceaccumulator)
- [2026-11-11 CMSSW plugins](#2026-11-11-cmssw-plugins)

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

## 2026-11-11 Review round 3: de-duplication
- Decided with user: keep every strategy and output (an Extra may need `loop::Ordered`); no templated Norm/Weight policies;
  `Dense` keeps its once-per-jet allocation; E1C no longer throws on |p| = 0 (geometry is built only for pENC).
- `detail::CalculatorBase` now holds n (must be finite), binning, accumulator and the ignored-externalNorm flag for both
  calculators. `fillJetData` split into `fillWeights` + `fillGeometry`. N is `std::size_t` everywhere (std::array's N forces it).
  Template parameter `Strategy` renamed `Strat`.
- `BoostRotate` applies fastjet's unboost formula directly instead of copying the PseudoJet; with the geometry no longer
  built for E1C, E1C M=50 went 4.75 -> 2.64 us, M=150 14.2 -> 7.8 us (min of 3, core-pinned).
- pENC leaf, interleaved min of 3 new/old: 0.95-1.07 everywhere except Dense PairAngles N=3 in `eec_bench` (1.57). A
  binary with only that case and PairBins measured new/old ~1.0 (1.95-2.05 ms both), so this is the inlining sensitivity
  noted in round 2, not extra work. Routing binned PairAngles through `linearIndex(coord)` was 1.7x in `eec_bench`; not adopted.
- Compile-fail regexes accept GCC and clang 17 wording; verified both, and a compiling mutant matches neither.
- 13 mutations (n check, unboost formula, static_cast charge, bin0 stride, E1C geometry, fillWeights touching geometry,
  NaN externalNorm warning, shared warn flag, extras dropped in both leaf paths, test helper skipping an output or the
  callback, no symfactor): each fails >= 1 test.

## 2026-11-11 Inner-loop optimisation study
Prototypes in scratch copies only; library unchanged. Dense, Energy/SumEnergy/DeltaR, 50-bin log axis, core-pinned, min of 3-5.
- Cache locality is not limiting: ns/tuple is flat from L1-resident to beyond L3 (N=2 PairBins 12.1-13.5 ns for
  M=50..2500, the M=2500 matrix is 25 MB; N=3 PairBins 3.2 -> 2.1 ns for M=30..500; PairAngles likewise). Blocking or
  loop reordering not pursued.
- `eec_bench` reuses one jet, so the branch predictor can learn data-dependent branches: a run-length-merging leaf
  measured 0.66-0.76x there but 1.11-1.30x when the benchmark cycles through 32 distinct jets. Earlier loop decisions
  (row-max, sym hoisting) were measured on the single-jet benchmark.
- With 32 rotating jets, PairBins, new/old, reproduced on two cores:
  - `peel`: last tuple index in its own loop; first iteration (repeated index) peeled, so symfactor is one division
    per row, no per-leaf bookkeeping: N=3 0.77-0.80, N=4 0.84, N=5 0.89, N=2 0.99.
  - `il2`: Dense keeps 2 interleaved histogram copies, alternating per fill, summed in finish(): 0.94-0.96 alone.
  - `peel` + `il2`: N=3 0.69-0.73, N=4 0.67, N=5 0.70. Breaks the same-bin read-modify-write chain the peel exposes.
  - PairAngles and Multiset unchanged (0.99-1.04): their per-leaf axis lookup dominates.

## 2026-11-11 Peeled last level, rotating-jet benchmark
- Decided with user: adopt `peel` for every output/Extra; no 2-copy Dense (memory); `eec_bench` cycles 32 distinct jets
  and gains Extra (constant, MinDeltaR) and more output cases. The 2026-11-10 table used one repeated jet.
- `lastLevel` peels the first multiset iteration; `leaf` takes the symfactor. Ordered loop: no peel, symfactor 1.
  Mutations (peeled/constant symfactor wrong, peeled tuple dropped, peel applied to Ordered) each fail >= 1 test.
- Full `eec_bench`, new/HEAD, interleaved min of 3, three independent runs agree to ~0.03:
  - faster: Dense PairBins N=3 0.78-0.88, N=4 0.84-0.89; PairAngles N=4,5 0.80-0.93; Multiset 0.84-0.97;
    SparseUnmerged 0.77-0.90; Unbinned 0.81-0.97; SparseMerged 0.92-0.99; ConstExtra 0.74-0.88; MinDeltaR 0.94-0.99.
  - neutral: N=2 PairBins/PairAngles, Naive, SparseMergedSort, PairAngles N=3 (0.99-1.03).
  - slower: Dense PairBins N=5 1.08-1.15, N=6 1.20-1.36; PairAngles N=6 1.04-1.06.
- Tried against the N>=5 regression (all tests pass for each):
  - row pointers hoisted per row (with/without skipping the idx store when K=0): worse, up to 1.7x at N=6.
  - single loop with both symfactors hoisted (`j == start` select): fixes N=5,6 but loses most gains, N=4 PairBins 1.14.
  - `-fno-tree-vectorize`: N=6 1.05 but N=5 1.15-1.23. Not a clean vectorisation effect.
  - `[[gnu::noinline]]` lastLevel: 1.2-2.0x slower broadly (loop invariants no longer hoisted).
- Decided with user: keep `Peeled` as the default and expose the old loop as `strategy::innermost::Unpeeled`, a third
  defaulted `Strategy` parameter (per calculator, tested in the same binary). Oracle comparison runs every output x loop x
  precompute for both; mutations of each path's symfactor fail >= 1 test. Which path runs is not observable in tests;
  benchmark: Unpeeled/HEAD PairBins N=3..6 0.99-1.03 (it is HEAD's loop), Peeled/HEAD as above.

## 2026-11-11 CovarianceAccumulator
- `interface/Covariance.h`: packed upper triangle; sparse update over the touched bins (sorted, so row-wise writes).
- Bench (GCC 12, -O2): 2704 bins (52×52 E1C) costs 2.1 µs/event for 50 random adds, 68 µs for 300, and 1.07 ms with every bin touched; 52 bins, 52 adds: 0.6 µs.
- Tests compare against a dense brute force (n = 1, 2, 7, 40; signed weights; repeated bins). Mutations (w for w², no merge of repeated bins, no diagonal, no reset of an abandoned event, merge skipping C) all fail the tests now; the abandoned-event mutation needed a stronger test that re-fills an abandoned bin.

## 2026-11-11 CMSSW plugins
- 32 module types (16 pENC, 8 E1C, 5 tables, 3 aggregators); a full scram build of the plugins takes about 25 s of wall time with -j8.
- Traps: a `std::make_unique<FlatTable>(...)` with dependent arguments is itself dependent, so the type must be spelled out before calling `addColumn<T>`.
  TFileService already makes a directory per module label. `TVectorD` has no name, so edges are stored as an empty TH1D binning.
  `TestProcessor` needs `JobReportService` for TFileService. Calling `testBeginRun` after `test()` crashed in `clearPrincipal`, so run tables are tested first.
- Tests: TestProcessor results compared against direct library calls. The lab-frame E1C is checked against `PseudoJet::theta` and atan2. The aggregator is checked against a hand-computed H and C.
  Mutations caught: ignored norm, ignored frame, symfactor column, unit gen weight, dropped empty jets, wrong iJet, and unmerged symfactor in the aggregator (this last one needed an extra test).
- The library config helpers have their own mutation check: CS sign flip, CS x sign, helicity x axis, jetPt vs jetE, and a missing unboost all fail the tests.
