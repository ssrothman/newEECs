# EECs library guide

## Contents
- [1. Observables and conventions](#1-observables-and-conventions)
- [2. Calculators](#2-calculators)
- [3. Per-jet data flow](#3-per-jet-data-flow)
- [4. Policies](#4-policies)
- [5. Outputs](#5-outputs)
- [6. Strategies](#6-strategies)
- [7. Binning](#7-binning)
- [8. State, threading and errors](#8-state-threading-and-errors)
- [9. Extending the library](#9-extending-the-library)

Everything lives in namespace `eec`. Include `SRothman/EECs/interface/EECs.h` for all calculators.

## 1. Observables and conventions

**Projected N-point correlator (pENC)**

pENC(θ) = Σ over ordered N-tuples (p₁…p_N), drawn with replacement, of ∏ᵢ (E_pᵢ / E_jet)ⁿ · δ(θ_max − θ),

where θ_max is the largest pairwise angle θ_ij in the tuple (0 when all entries are the same particle).

A tuple's weight and θ_max do not depend on the order of its entries. So, except under `strategy::Naive`, the
library evaluates each multiset i₁ ≤ … ≤ i_N once and attaches the multinomial symmetry factor
N! / ∏ₖ mₖ!, where mₖ are the multiplicities of the distinct indices. Example for N = 3: (a,a,a) → 1,
(a,a,b) → 3, (a,b,c) → 6.

**One-point correlator (E1C)**

E1C(c) = Σᵢ (E_pᵢ / E_jet)ⁿ · δ(cᵢ − c),

where cᵢ = `frame(pᵢ)` is a D-dimensional coordinate of each particle in a caller-supplied frame.
E1C entries always have symmetry factor 1.

**E_jet tracking.** Every result carries `JetMeta{E_jet, n, N}`, where E1C has N = 1.
`renormalization(meta, newEjet)` = (E_jet / newEjet)^(n·N). Multiply by it to re-express the weights with another E_jet.

## 2. Calculators

| Class | Template parameters | Constructors | Call |
|---|---|---|---|
| `ProjectedENC` | `<int N, Weight, Norm, Angle, Output, Extra = NoExtra, Strategy = strategy::PairBins>`, with 1 ≤ N ≤ 12 | binned output: `(n, Binning, Extra = {})`; unbinned: `(n, Extra = {})` | `operator()(particles, externalNorm = kNoExternalNorm)` |
| `ProjectedENCRuntimeN` | same, without `N` | `(N, n, Binning, Extra = {})` / `(N, n, Extra = {})`, with 2 ≤ N ≤ 6 | same as above |
| `OnePointEC` | `<Weight, Norm, Frame, Output>` | `(n, Binning)` / `(n)` | `operator()(particles, frame, externalNorm = kNoExternalNorm)` |

- `particles` is a `std::vector<fastjet::PseudoJet>`.
- The binned constructor checks the binning rank: it must be 1 + `Extra::K` for pENC and `Frame::D` for E1C.
- The 12 upper bound on N exists because symmetry factors are computed in 32-bit integers.
- `ProjectedENCRuntimeN` holds a `std::variant` of `ProjectedENC<2..6>` and dispatches with `std::visit`. The result type does not depend on N.
- In `OnePointEC`, the frame object is passed on every call, so a per-jet frame (for example the boson rest frame) needs no rebuild of the calculator.

## 3. Per-jet data flow

1. `fillJetData<Weight, Norm>` fills a `JetData`, a structure of arrays:
   - `E_jet` = `Norm::value(particles, externalNorm)`
   - `w[i]` = (`Weight::value(pᵢ)` / E_jet)ⁿ. `pow` is skipped when n = 1.
   - `ux, uy, uz`: the unit three-momentum.
   - `eta`: pseudorapidity.
   - `phi`: azimuth in [0, 2π).
2. pENC only: the strategy may precompute an M×M matrix. It then walks the tuples, carrying the running weight
   product and θ_max, and calls `Output::Accumulator::fill` once per tuple or multiset.
   E1C: one fill per particle, at `frame(pᵢ)`.
3. `Accumulator::finish(meta)` returns the result.

`Angle` and `Extra` policies see the particles only through `JetData`. They do not get the `PseudoJet`s.

## 4. Policies

All policies are compile-time types, so there is no virtual dispatch on the hot path.

### Weight E_p — `eec::weight`
Interface: `static double value(const fastjet::PseudoJet&)`.

| Policy | E_p |
|---|---|
| `Energy` | `p.E()` |
| `Pt` | `p.pt()` |
| `Charge` | charge |
| `ChargeEnergy` | charge × E |
| `ChargePt` | charge × pt |

- Charge comes from `EECParticleInfo`, attached with `p.set_user_info(new EECParticleInfo(charge, pdgId))`. fastjet takes ownership.
- A particle without that info throws `std::invalid_argument`.
- Negative E_p is allowed. With a non-integer n, it throws `std::domain_error`, because `pow` would return NaN.

### Norm E_jet — `eec::norm`
Interface: `static double value(const std::vector<fastjet::PseudoJet>&, double externalNorm)`, plus an optional
`static constexpr bool usesExternalNorm`.

| Policy | E_jet | `usesExternalNorm` |
|---|---|---|
| `SumEnergy` | Σ E | false |
| `SumPt` | Σ pt | false |
| `External` | the `externalNorm` argument; throws `std::invalid_argument` unless it is finite and non-zero | true |
| `Unity` | 1 | false |

Passing an `externalNorm` (anything other than the NaN default `kNoExternalNorm`) to a calculator whose norm declares
`usesExternalNorm = false` warns once per calculator. The warning goes to `edm::LogWarning("EECs")` under scram and to
`std::cerr` otherwise. A norm that does not declare the flag is assumed to use the value and never warns.
The library does not check that E_jet ≠ 0, except for `External`.

### Angle θ_ij — `eec::angle`
Interface: `static double pair(const JetData&, std::size_t i, std::size_t j)`.

| Policy | θ_ij |
|---|---|
| `Angle3D` | the angle between three-momenta, computed as atan2(\|uᵢ×uⱼ\|, uᵢ·uⱼ). This keeps precision at small angles, unlike `acos`. |
| `DeltaR` | √(Δη² + Δφ²), with pseudorapidity η and Δφ wrapped into [−π, π] |

### Extra pENC coordinates
`NoExtra` (K = 0) is the default. An Extra adds K binning coordinates per tuple; see [9.4](#94-extra-coordinates).
Axis 0 is always θ_max, and axes 1…K are the Extra values.

### Frame (E1C)
`frame::BoostRotate(rest, z, x)`:
- It unboosts each particle into the rest frame of `rest`, which must be timelike.
- It returns (θ, φ) with respect to the axes `z` and `x`, both given in that rest frame. θ is in [0, π] and φ in [−π, π].
- `x` is orthogonalised against `z`, and y = z × x.
- It throws `std::invalid_argument` if `rest` is not timelike, or if `z`, or `x` after orthogonalisation, has norm ≤ 1e-12.
- To get plain lab angles, pass `rest = PseudoJet(0, 0, 0, m)` with m > 0.

## 5. Outputs
Policies are in `eec::output`.

| Policy | Result | Content |
|---|---|---|
| `Dense` | `DenseResult{meta, values}` | `values[linearIndex]` = Σ weight × symfactor over all flow-inclusive bins |
| `SparseMerged` | `SparseResult{meta, entries}` | `BinEntry{bin, weight}` sorted by bin, unique, weight = Σ weight × symfactor. Merged with a dense scratch array plus a touched-bin list. |
| `SparseMergedSort` | `SparseResult` | same result as `SparseMerged`, merged by sort + reduce; about 8× slower at N=3, M=100 (JOURNAL) |
| `SparseUnmerged` | `SparseUnmergedResult{meta, entries}` | `UnmergedBinEntry{bin, weight, symfactor}`, one per fill |
| `Unbinned` | `UnbinnedResult<D>{meta, entries}` | `UnbinnedEntry<D>{coord, weight, symfactor}`, one per fill |

- In unmerged outputs, `weight` is the bare product ∏(E_p/E_jet)ⁿ, and `symfactor` is stored separately.
- One fill means one multiset, giving multichoose(M, N) entries. Under `strategy::Naive` it means one ordered tuple, giving Mᴺ entries, each with `symfactor` = 1.
- Merged outputs keep a bin that was filled even if its weights cancel to zero. Signed (charge) weights can do this.
- `Unbinned` coordinates are (θ_max, extra₀, …) for pENC and `frame(p)` for E1C.

## 6. Strategies
These choose how pENC enumerates tuples; they live in `eec::strategy`. All strategies give the same histogram. The
tests compare each one against an ordered-tuple oracle.

| Strategy | Loop | Per-jet precompute | Notes |
|---|---|---|---|
| `Naive` | all Mᴺ ordered tuples, symfactor 1 | none | reference only |
| `Multiset` | i₁ ≤ … ≤ i_N | none | computes `Angle::pair` on the fly |
| `PairAngles` | multisets | M×M `double` angle matrix | |
| `PairBins` (default) | multisets | M×M `uint32` matrix of axis-0 bins | uses bin(max θ) = max bin(θ). Requires an increasing axis 0, and the constructor throws otherwise. With an unbinned output it behaves as `PairAngles`. |

- `PairBins` is 4.5–9× faster than `PairAngles` for N = 3–6 in the `JOURNAL.md` benchmarks, because it avoids an axis lookup per tuple.
- Cost grows as multichoose(M, N), and as Mᴺ for `Naive`.
- Unmerged and unbinned outputs store one entry per fill, so their memory grows the same way.
- The strategy set is closed: an unknown type fails a `static_assert`.

## 7. Binning
`Binning(std::vector<Axis>)` is a runtime binning built from boost::histogram axes:
- `RegularAxis`: `regular<double>`
- `LogAxis`: `regular<double, transform::log>`
- `VariableAxis`: `variable<double>`

- Every axis has underflow and overflow bins.
- `binIndex(k, x)` is flow-inclusive: 0 is underflow and `extent(k) − 1` is overflow. A `LogAxis` sends θ = 0 to underflow.
- `linearIndex(x)` = Σₖ `stride(k)` × `binIndex(k, xₖ)`, with axis 0 fastest. `unravel(lin)` inverts it.
- `toHistogram(dense)` returns a `boost::histogram` with the same axes, filled from a `Dense` values vector.
- The constructor throws `std::invalid_argument` for an empty axis list, or for more than 2³²−1 bins in total.

## 8. State, threading and errors
- Calculators keep per-jet scratch (`JetData`, matrices, accumulators) in `mutable` members and reuse it across jets.
  `operator()` is `const` but not thread-safe: use one instance per thread.
- If a policy throws mid-jet, the calculator can still be used for the next jet. The next call starts from clean state.

| Thrown by | Exception | Condition |
|---|---|---|
| `Binning` | `invalid_argument` | no axes; more than 2³²−1 bins |
| `Binning::toHistogram` | `invalid_argument` | size of the values vector ≠ `totalBins()` |
| `ProjectedENC` constructor | `invalid_argument` | rank ≠ 1 + K; `PairBins` with a decreasing axis 0 |
| `ProjectedENCRuntimeN` constructor | `invalid_argument` | N ∉ [2, 6] |
| `OnePointEC` constructor | `invalid_argument` | rank ≠ `Frame::D` |
| `frame::BoostRotate` | `invalid_argument` | `rest` not timelike; degenerate axis |
| `norm::External` | `invalid_argument` | `externalNorm` not finite, or zero |
| charge weights | `invalid_argument` | no `EECParticleInfo` |
| `fillJetData` | `domain_error` | negative E_p with non-integer n |

## 9. Extending the library
Each new policy is a plain type passed as a template argument. Nothing needs to be registered. The examples below
compile against the library. With n = 1 and a matching norm, they satisfy the sum rule Σ = 1.

### 9.1 Weight
```cpp
// Weight: energy of charged particles only.
struct ChargedEnergy {
  static double value(const fastjet::PseudoJet& p) { return eec::particleCharge(p) != 0 ? p.E() : 0.0; }
};
```

### 9.2 Norm
```cpp
// Norm: pt of the vector sum of the particles.
struct VectorSumPt {
  static constexpr bool usesExternalNorm = false;
  static double value(const std::vector<fastjet::PseudoJet>& ps, double) {
    fastjet::PseudoJet s(0, 0, 0, 0);
    for (const auto& p : ps)
      s += p;
    return s.pt();
  }
};
```
Declare `usesExternalNorm = false` if the norm ignores the second argument; this enables the warning described in
[section 4](#norm-e_jet--eecnorm).

### 9.3 Angle
```cpp
// Angle: pseudorapidity separation.
struct DeltaEta {
  static double pair(const eec::JetData& d, std::size_t i, std::size_t j) { return std::abs(d.eta[i] - d.eta[j]); }
};
```
Requirements:
- The angle must be ≥ 0, with `pair(d, i, i) == 0`. θ_max starts from 0, and the matrix diagonals are set to 0.
- It must be symmetric, `pair(i, j) == pair(j, i)`. The matrices fill only i < j, while `Naive` and `Multiset` call it in either order and with i == j.
- It can use only `JetData` fields. If you need another per-particle quantity (for example rapidity), add a field to `JetData` and fill it in `fillJetData`.

### 9.4 Extra coordinates
```cpp
// Extra coordinate: smallest DeltaR between distinct particles of the tuple (0 if all indices are equal).
struct MinDeltaR {
  static constexpr std::size_t K = 1;
  template <std::size_t N>
  std::array<double, 1> operator()(const eec::JetData& d, const std::array<std::uint32_t, N>& idx) const {
    double m = 0;
    bool any = false;
    for (std::size_t a = 0; a < N; ++a)
      for (std::size_t b = a + 1; b < N; ++b)
        if (idx[a] != idx[b]) {
          const double r = eec::angle::DeltaR::pair(d, idx[a], idx[b]);
          m = any ? std::min(m, r) : r;
          any = true;
        }
    return {{m}};
  }
};
// ProjectedENC<3, weight::Energy, norm::SumEnergy, angle::DeltaR, output::Dense, MinDeltaR>
//     calc(1.0, Binning({RegularAxis(20, 0.0, 1.0), RegularAxis(10, 0.0, 1.0)}));
```
Requirements:
- The result must be symmetric under permutations of `idx`. One fill stands for all orderings of a multiset, and `idx` arrives sorted non-decreasing except under `Naive`.
- Indices may repeat.
- The Extra is stored by value, constructed from the constructor's `extra` argument (default `Extra{}`), and called once per fill. It may hold configuration.
- With a binned output, the binning needs 1 + K axes.

### 9.5 Frame (E1C)
```cpp
// Frame: lab-frame polar and azimuthal angles.
struct LabPolar {
  static constexpr std::size_t D = 2;
  std::array<double, 2> operator()(const fastjet::PseudoJet& p) const { return {{p.theta(), p.phi_std()}}; }
};
// OnePointEC<weight::Pt, norm::SumPt, LabPolar, output::Unbinned> e1c(1.0);  e1c(particles, LabPolar{});
```
Any D works. A binned output needs D axes.

### 9.6 Output
An output policy provides:
- `static constexpr bool binned`
- `template <std::size_t D> Result`
- `template <std::size_t D> Accumulator`, a default-constructible class with these methods:
  - `begin(const Binning&, std::size_t expectedFills)` if `binned`, else `begin(std::size_t expectedFills)`
  - `fill(std::uint32_t bin, double weight, std::uint32_t symfactor)` if `binned`, else
    `fill(const std::array<double, D>& coord, double weight, std::uint32_t symfactor)`
  - `Result<D> finish(const JetMeta&)`

`weight` excludes the symmetry factor. `begin` must reset all per-jet state, including state left by a jet that threw
before `finish`.

```cpp
// Output: only the total weight of the jet.
struct TotalWeight {
  static constexpr bool binned = false;
  template <std::size_t D>
  struct Result {
    eec::JetMeta meta;
    double total = 0;
  };
  template <std::size_t D>
  class Accumulator {
  public:
    void begin(std::size_t) { r_ = {}; }
    void fill(const std::array<double, D>&, double w, std::uint32_t sym) { r_.total += w * sym; }
    Result<D> finish(const eec::JetMeta& m) {
      r_.meta = m;
      return r_;
    }

  private:
    Result<D> r_;
  };
};
```

### 9.7 Testing a new policy
- **pENC:** add it to `compareAll` in `test/test_penc.cc`. That compares every output and strategy against an ordered-tuple oracle. The oracle needs a reference angle, or a reference Extra, computed independently from the `PseudoJet`s.
- **Any weight or norm:** for n = 1 and E_jet = Σ E_p, the sum rule in `test/test_sum_rule.cc` checks the bookkeeping.
- **CMSSW:** new test files must also be listed in `test/BuildFile.xml`. CMake globs `test/*.cc`.
