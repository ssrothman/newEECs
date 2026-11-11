# SRothman/EECs

Header-mostly C++20 library for projected N-point energy correlators (pENC) and one-point
correlators (E1C) on `fastjet::PseudoJet` inputs. Builds standalone with CMake or as the CMSSW
package `SRothman/EECs`.

## Build

Standalone (externals taken from the CMSSW toolbox of the `cmsenv`'d release):
```
cmsenv   # e.g. in CMSSW_15_0_22/src
cmake -S . -B build && cmake --build build -j8
ctest --test-dir build    # unit tests (./build/testEECs) and compile-fail checks
./build/eec_bench         # benchmarks
```
Without CMSSW, pass `-DFASTJET_ROOT=... -DBOOST_INCLUDE_DIR=... -DCATCH2_INCLUDE_DIR=... -DGBENCH_ROOT=...`.

CMSSW: place (clone or symlink) this repo at `$CMSSW_BASE/src/SRothman/EECs`, then
`scram b -j8 && scram b runtests`. Benchmarks are CMake-only.

## Usage
```cpp
#include "SRothman/EECs/interface/EECs.h"
using namespace eec;

Binning b({LogAxis(50, 1e-3, 1.0)});
ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::SparseMerged> e3c(/*n=*/1.0, b);
SparseResult r = e3c(constituents);           // r.meta holds E_jet, n, N

ProjectedENCRuntimeN<weight::Pt, norm::External, angle::Angle3D, output::Dense> enc(N, 1.0, b);
DenseResult d = enc(constituents, jetPt);

frame::BoostRotate f(bosonP4, zAxisInRestFrame, xAxisInRestFrame);
OnePointEC<weight::Energy, norm::SumEnergy, frame::BoostRotate, output::Unbinned> e1c(1.0);
UnbinnedResult<2> u = e1c(constituents, f);   // coords (theta, phi)
```
Charge-based weights need `p.set_user_info(new EECParticleInfo(charge, pdgId))`.
Calculators keep per-jet scratch space: use one instance per thread.

## Conventions
- pENC sums over ordered N-tuples with replacement; it is evaluated over multisets
  i1 <= ... <= iN with multinomial symmetry factor N!/prod(m_k!).
- Unmerged outputs (`SparseUnmerged`, `Unbinned`) hold one entry per multiset with `weight`
  (product of (E_p/E_jet)^n) and `symfactor` separately; `Dense` and `SparseMerged*` hold the product.
  With the `loop::Ordered` loop (e.g. `strategy::Naive`) they instead hold one entry per ordered tuple, each with `symfactor` 1.
- Bin indices include flow bins (0 = underflow); the linear index has axis 0 fastest.
  `Binning::toHistogram` converts dense values into a `boost::histogram`.
- `renormalization(meta, newEjet)` = (E_jet/newEjet)^(n N) rescales weights to another E_jet.
- `angle::DeltaR` uses pseudorapidity. `frame::BoostRotate` takes its z and x axes in the rest frame.
- Merged sparse outputs keep a bin that was filled even if its weights cancel to zero.

## Extending
See `docs/guide.md` for full documentation of policies, outputs, strategies and how to extend them.

All configuration is compile-time policy structs:
- particle weight E_p: `struct { static double value(const fastjet::PseudoJet&); }`
- jet norm E_jet: `struct { static double value(const std::vector<PseudoJet>&, double external); }`, optionally with
  `static constexpr bool usesExternalNorm`. If it is `false`, passing an `externalNorm` warns once per calculator
  (`edm::LogWarning` under scram, `std::cerr` otherwise); without it the norm is assumed to use the value.
- pair angle: `struct { static double pair(const JetData&, size_t i, size_t j); }` (>= 0, 0 on the diagonal)
- extra pENC binning coordinates: `struct { static constexpr size_t K; template <size_t N> std::array<double, K> operator()(const JetData&, const std::array<uint32_t, N>&) const; }`,
  permutation-symmetric; the binning then needs 1 + K axes (axis 0 is theta_max). An optional `update(...)` is
  reached through the calculator's `update_extra(...)`.
- E1C frame: class with `static constexpr size_t D` and `std::array<double, D> operator()(const PseudoJet&) const`.

Strategies are `strategy::Strategy<Loop, Precompute>` with `loop::{Ordered, Multiset}` and
`precompute::{None, PairAngles, PairBins}`; `Naive`, `Multiset`, `PairAngles`, `PairBins` are aliases, compared in `JOURNAL.md`.
`precompute::PairBins` needs a binned output and an increasing axis 0 (the constructor throws for a decreasing one).
