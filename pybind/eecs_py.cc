// Python bindings (pybind11) of the projected N-point correlator, for offline computation from stored constituents.
// Build: pybind/build.sh. The calculator is the one of PENCProducer (eec::ProjectedENCRuntimeN, norm::External,
// output::Dense), so the results follow the conventions of the NANO tables (EECs/docs/guide.md).
#include "SRothman/EECs/interface/ProjectedENC.h"
#include "SRothman/EECs/interface/ProjectedENCRuntime.h"
#include "SRothman/EECs/interface/TriangleShape.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;

namespace {

  using Array = py::array_t<double, py::array::c_style | py::array::forcecast>;
  using Index = py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>;

  template <class Calc>
  void loop(const Calc& calc, const double* px, const double* py, const double* pz, const double* E,
            const std::int64_t* off, const double* norms, std::size_t nJets, std::size_t nb, double* out) {
    std::vector<fastjet::PseudoJet> ps;
    for (std::size_t j = 0; j < nJets; ++j) {
      ps.clear();
      for (std::int64_t k = off[j]; k < off[j + 1]; ++k)
        ps.emplace_back(px[k], py[k], pz[k], E[k]);
      if (ps.empty())
        continue;   // the row stays zero, as PENCProducer's empty result
      const auto r = calc(ps, norms[j]);
      for (std::size_t b = 0; b < nb; ++b)
        out[j * nb + b] = r.values[b];
    }
  }

  template <class Weight, class Angle>
  void run(const double* px, const double* py, const double* pz, const double* E, const std::int64_t* off,
           const double* norms, std::size_t nJets, std::size_t N, double n, const eec::Binning& binning, double* out) {
    const eec::ProjectedENCRuntimeN<Weight, eec::norm::External, Angle, eec::output::Dense> calc(N, n, binning);
    loop(calc, px, py, pz, E, off, norms, nJets, binning.totalBins(), out);
  }

  template <class Weight, class Angle>
  void runResolved(const double* px, const double* py, const double* pz, const double* E, const std::int64_t* off,
                   const double* norms, std::size_t nJets, double n, const eec::Binning& binning, double* out) {
    const eec::ProjectedENC<3, Weight, eec::norm::External, Angle, eec::output::Dense, eec::TriangleShape<Angle>> calc(
        n, binning);
    loop(calc, px, py, pz, E, off, norms, nJets, binning.totalBins(), out);
  }

  std::vector<double> checkedEdges(const Array& edges, const char* name) {
    if (edges.size() < 2)
      throw std::invalid_argument(std::string(name) + " needs at least two entries");
    std::vector<double> e(edges.data(), edges.data() + edges.size());
    for (std::size_t k = 1; k < e.size(); ++k)
      if (!(e[k] > e[k - 1]))
        throw std::invalid_argument(std::string(name) + " must increase");
    return e;
  }

  // Number of jets after checking the particle arrays, offsets and norms.
  std::size_t checkedJets(const Array& px, const Array& py, const Array& pz, const Array& E, const Index& offsets,
                          const Array& norms) {
    const std::size_t nPart = px.size();
    if (py.size() != static_cast<py::ssize_t>(nPart) || pz.size() != static_cast<py::ssize_t>(nPart) ||
        E.size() != static_cast<py::ssize_t>(nPart))
      throw std::invalid_argument("px, py, pz and E must have the same length");
    if (offsets.size() < 1)
      throw std::invalid_argument("offsets needs at least one entry");
    const std::size_t nJets = offsets.size() - 1;
    if (norms.size() != static_cast<py::ssize_t>(nJets))
      throw std::invalid_argument("norms needs one entry per jet (len(offsets) - 1)");
    const auto* off = offsets.data();
    if (off[0] != 0 || off[nJets] != static_cast<std::int64_t>(nPart))
      throw std::invalid_argument("offsets must start at 0 and end at the number of particles");
    for (std::size_t j = 0; j < nJets; ++j)
      if (off[j + 1] < off[j])
        throw std::invalid_argument("offsets must not decrease");
    for (std::size_t j = 0; j < nJets; ++j)
      if (!(std::isfinite(norms.data()[j]) && norms.data()[j] != 0))
        throw std::invalid_argument("norms must be finite and non-zero");
    return nJets;
  }

  py::array_t<double> projectedENC(Array px, Array py, Array pz, Array E, Index offsets, Array norms, std::size_t N,
                                   double n, Array edges, const std::string& weight, const std::string& angle) {
    const std::size_t nJets = checkedJets(px, py, pz, E, offsets, norms);
    const auto* off = offsets.data();
    const eec::Binning binning({eec::VariableAxis(checkedEdges(edges, "edges"))});
    const std::size_t nb = binning.totalBins();

    py::array_t<double> out({static_cast<py::ssize_t>(nJets), static_cast<py::ssize_t>(nb)});
    double* o = out.mutable_data();
    std::fill(o, o + nJets * nb, 0.0);
    const double *x = px.data(), *y = py.data(), *z = pz.data(), *en = E.data(), *nm = norms.data();
    {
      py::gil_scoped_release release;
      if (weight == "Pt" && angle == "DeltaR")
        run<eec::weight::Pt, eec::angle::DeltaR>(x, y, z, en, off, nm, nJets, N, n, binning, o);
      else if (weight == "Pt" && angle == "Angle3D")
        run<eec::weight::Pt, eec::angle::Angle3D>(x, y, z, en, off, nm, nJets, N, n, binning, o);
      else if (weight == "Energy" && angle == "DeltaR")
        run<eec::weight::Energy, eec::angle::DeltaR>(x, y, z, en, off, nm, nJets, N, n, binning, o);
      else if (weight == "Energy" && angle == "Angle3D")
        run<eec::weight::Energy, eec::angle::Angle3D>(x, y, z, en, off, nm, nJets, N, n, binning, o);
      else
        throw std::invalid_argument("weight must be Pt or Energy and angle DeltaR or Angle3D");
    }
    return out;
  }

  py::array_t<double> resolvedE3C(Array px, Array py, Array pz, Array E, Index offsets, Array norms, double n,
                                  Array edgesRL, Array edgesXi, Array edgesPhi, const std::string& weight,
                                  const std::string& angle) {
    const std::size_t nJets = checkedJets(px, py, pz, E, offsets, norms);
    const auto* off = offsets.data();
    const eec::Binning binning({eec::VariableAxis(checkedEdges(edgesRL, "edgesRL")),
                                eec::VariableAxis(checkedEdges(edgesXi, "edgesXi")),
                                eec::VariableAxis(checkedEdges(edgesPhi, "edgesPhi"))});
    const std::size_t nb = binning.totalBins();
    py::array_t<double> out({static_cast<py::ssize_t>(nJets), static_cast<py::ssize_t>(nb)});
    double* o = out.mutable_data();
    std::fill(o, o + nJets * nb, 0.0);
    const double *x = px.data(), *y = py.data(), *z = pz.data(), *en = E.data(), *nm = norms.data();
    {
      py::gil_scoped_release release;
      if (weight == "Pt" && angle == "DeltaR")
        runResolved<eec::weight::Pt, eec::angle::DeltaR>(x, y, z, en, off, nm, nJets, n, binning, o);
      else if (weight == "Pt" && angle == "Angle3D")
        runResolved<eec::weight::Pt, eec::angle::Angle3D>(x, y, z, en, off, nm, nJets, n, binning, o);
      else if (weight == "Energy" && angle == "DeltaR")
        runResolved<eec::weight::Energy, eec::angle::DeltaR>(x, y, z, en, off, nm, nJets, n, binning, o);
      else if (weight == "Energy" && angle == "Angle3D")
        runResolved<eec::weight::Energy, eec::angle::Angle3D>(x, y, z, en, off, nm, nJets, n, binning, o);
      else
        throw std::invalid_argument("weight must be Pt or Energy and angle DeltaR or Angle3D");
    }
    return out;
  }

}  // namespace

PYBIND11_MODULE(eecs_py, m) {
  m.doc() = "Projected N-point energy correlators of the SRothman/EECs library";
  m.def("projected_enc", &projectedENC, py::arg("px"), py::arg("py"), py::arg("pz"), py::arg("E"), py::arg("offsets"),
        py::arg("norms"), py::arg("N"), py::arg("n"), py::arg("edges"), py::arg("weight") = "Pt",
        py::arg("angle") = "DeltaR",
        R"doc(Binned projected EN C of each jet: an array (jets, len(edges) + 1) with the underflow first and the
overflow last. Jet j holds the particles offsets[j]:offsets[j+1] of px, py, pz, E; every tuple (self tuples
included, at angle 0) adds the product of (E_p / norms[j])^n with E_p the particle pt (weight Pt) or energy,
times its symmetry factor, at its largest pairwise angle (DeltaR in eta-phi, or Angle3D).)doc");
  m.def("resolved_e3c", &resolvedE3C, py::arg("px"), py::arg("py"), py::arg("pz"), py::arg("E"), py::arg("offsets"),
        py::arg("norms"), py::arg("n"), py::arg("edgesRL"), py::arg("edgesXi"), py::arg("edgesPhi"),
        py::arg("weight") = "Pt", py::arg("angle") = "DeltaR",
        R"doc(Binned resolved E3C of each jet: an array (jets, B) with B = (len(edgesRL) + 1) (len(edgesXi) + 1)
(len(edgesPhi) + 1), flow bins included on every axis, R_L varying fastest. Each triple (repeated particles
included) adds as projected_enc with N = 3 at (R_L, xi = R_S / R_M, phi = arcsin sqrt(1 - (R_L - R_M)^2 / R_S^2))
(arXiv:2201.07800); a triple with R_S = 0 is at xi = phi = 0.)doc");
}
