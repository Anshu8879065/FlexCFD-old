// main.cpp
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Numerics/casulliwrap/CasulliWrap.hpp"
#include "PDEs/SVE/SVECas.hpp"
#include "PDEs/SWE/SWE2dCas.hpp"

struct MyPDEOptions {};

namespace
{
using Real = double;

static std::string get_arg(int argc, char** argv, int k, std::string fallback = "")
{
  return (argc > k && argv[k]) ? std::string(argv[k]) : fallback;
}

static std::size_t clamp_stride(std::size_t n, std::size_t target_points = 12)
{
  if (n == 0) return 1;
  std::size_t stride = n / target_points;
  return std::max<std::size_t>(1, stride);
}

template <class T>
static void print_1d_sample(const std::string& name,
                            const std::vector<T>& v,
                            std::size_t stride,
                            std::size_t max_rows = 200)
{
  std::cout << "\n" << name << " (sampled), n=" << v.size() << ", stride=" << stride << "\n";
  std::size_t rows = 0;
  for (std::size_t i = 0; i < v.size(); i += stride) {
    std::cout << name << "[" << i << "]=" << std::setprecision(16) << v[i] << "\n";
    if (++rows >= max_rows) break;
  }
  if (!v.empty() && ((v.size() - 1) % stride != 0) && rows < max_rows) {
    std::cout << name << "[" << (v.size() - 1) << "]=" << std::setprecision(16) << v.back() << "\n";
  }
}

static void print_usage(char** argv)
{
  std::cout
    << "Usage:\n"
    << "  " << argv[0] << " N [ic]\n\n"
    << "Where:\n"
    << "  N   : grid size (1D: N cells, 2D: NxN)\n"
    << "  ic  : initial condition selector (optional)\n"
    << "       uniform  -> use PDE SetIC() values (uniform)\n"
    << "       step     -> 1D dam-break step / 2D gaussian bump (default)\n\n"
    << "Examples:\n"
    << "  " << argv[0] << " 41\n"
    << "  " << argv[0] << " 41 uniform\n";
}

} // namespace


// -------------------- 2D-only helpers (ONLY compiled in 2D) --------------------
#if defined(CASE_2D)
namespace
{
static inline std::size_t idx_eta(std::size_t i, std::size_t j, std::size_t nx) { return i + nx * j; }
static inline std::size_t idx_u  (std::size_t i, std::size_t j, std::size_t nx) { return i + (nx - 1) * j; } // i=0..nx-2
static inline std::size_t idx_v  (std::size_t i, std::size_t j, std::size_t nx) { return i + nx * j; }        // j=0..ny-2

static Real u_at_cell_center(std::size_t i, std::size_t j,
                             std::size_t nx, std::size_t /*ny*/,
                             std::vector<Real> const& uface)
{
  if (nx < 2) return Real(0);

  if (i == 0) {
    return uface[idx_u(0, j, nx)];
  }
  if (i >= nx - 1) {
    return uface[idx_u(nx - 2, j, nx)];
  }
  const Real ul = uface[idx_u(i - 1, j, nx)];
  const Real ur = uface[idx_u(i,     j, nx)];
  return Real(0.5) * (ul + ur);
}

static Real v_at_cell_center(std::size_t i, std::size_t j,
                             std::size_t nx, std::size_t ny,
                             std::vector<Real> const& vface)
{
  if (ny < 2) return Real(0);

  if (j == 0) {
    return vface[idx_v(i, 0, nx)];
  }
  if (j >= ny - 1) {
    return vface[idx_v(i, ny - 2, nx)];
  }
  const Real vb = vface[idx_v(i, j - 1, nx)];
  const Real vt = vface[idx_v(i, j,     nx)];
  return Real(0.5) * (vb + vt);
}

static void print_2d_sample_side_by_side(std::size_t nx, std::size_t ny,
                                         Real dx, Real dy,
                                         std::vector<Real> const& eta,
                                         std::vector<Real> const& uface,
                                         std::vector<Real> const& vface,
                                         std::size_t stride_i,
                                         std::size_t stride_j,
                                         std::size_t max_rows = 400)
{
  std::cout << "\n=== 2D sampled table (cell centers) ===\n";
  std::cout << "nx=" << nx << " ny=" << ny
            << "  stride_i=" << stride_i << " stride_j=" << stride_j << "\n\n";

  std::cout
    << std::left
    << std::setw(6)  << "i"
    << std::setw(6)  << "j"
    << std::setw(10) << "x"
    << std::setw(10) << "y"
    << std::setw(24) << "eta(i,j)"
    << std::setw(24) << "u_cc(i,j)"
    << std::setw(24) << "v_cc(i,j)"
    << "\n";

  std::cout << std::string(6+6+10+10+24+24+24, '-') << "\n";

  std::size_t rows = 0;
  for (std::size_t j = 0; j < ny; j += stride_j) {
    for (std::size_t i = 0; i < nx; i += stride_i) {
      const Real x = Real(i) * dx;
      const Real y = Real(j) * dy;

      const Real e = eta[idx_eta(i, j, nx)];
      const Real u = u_at_cell_center(i, j, nx, ny, uface);
      const Real v = v_at_cell_center(i, j, nx, ny, vface);

      std::cout
        << std::left
        << std::setw(6)  << i
        << std::setw(6)  << j
        << std::setw(10) << std::setprecision(6)  << std::fixed << x
        << std::setw(10) << std::setprecision(6)  << std::fixed << y
        << std::setw(24) << std::setprecision(16) << std::scientific << e
        << std::setw(24) << std::setprecision(16) << std::scientific << u
        << std::setw(24) << std::setprecision(16) << std::scientific << v
        << "\n";

      if (++rows >= max_rows) {
        std::cout << "... (truncated after " << max_rows << " rows)\n";
        return;
      }
    }
  }

  if (nx > 0 && ny > 0) {
    const std::size_t i = nx - 1, j = ny - 1;
    const Real x = Real(i) * dx, y = Real(j) * dy;
    const Real e = eta[idx_eta(i, j, nx)];
    const Real u = u_at_cell_center(i, j, nx, ny, uface);
    const Real v = v_at_cell_center(i, j, nx, ny, vface);

    std::cout << "\n(last cell)\n";
    std::cout << "i=" << i << " j=" << j
              << " x=" << std::fixed << std::setprecision(6) << x
              << " y=" << std::fixed << std::setprecision(6) << y
              << "  eta=" << std::scientific << std::setprecision(16) << e
              << "  u_cc=" << std::scientific << std::setprecision(16) << u
              << "  v_cc=" << std::scientific << std::setprecision(16) << v
              << "\n";
  }
}
} // namespace
#endif // CASE_2D

int main(int argc, char** argv)
{
  const std::string a1 = get_arg(argc, argv, 1, "");
  if (a1 == "--help" || a1 == "-h") {
    print_usage(argv);
    return 0;
  }

  const std::size_t N = (argc > 1) ? static_cast<std::size_t>(std::stoul(argv[1])) : 41;
  const std::string ic_mode = get_arg(argc, argv, 2, "step");
  const bool use_uniform_ic = (ic_mode == "uniform");

#if defined(CASE_1D)
  {
    std::cout << "=== 1D CASE (SVEcas) ===\n";
    std::cout << "N=" << N << "  ic=" << ic_mode << "\n";

    auto pde = std::make_unique<fcfd::pdemodel::SVEcas<Real, MyPDEOptions>>();

    // Model params (your SVEcas provides GetModelParams())
    pde->GetModelParams().g      = 9.81;
    pde->GetModelParams().nm     = 0.03;
    pde->GetModelParams().radi   = 1.0;
    pde->GetModelParams().gammat = 0.1;

    // Uniform IC values (used if ic=uniform)
    pde->SetIC(0.0, 0.0);      // (eta0, u0)
    pde->SetBed(0.0, 0.0);
    pde->SetWindSpeedU(0.0);

    // Example BCs (subcritical: specify 1 condition on each end)
    pde->SetLeftEta(0.0);
    pde->SetRightU(0.0);

    fcfd::pdenumerics::CasulliWrap<Real, MyPDEOptions> wrap;
    wrap.SetGrid1D(N);

    wrap.MethodProps().dt           = 0.1;
    wrap.MethodProps().dx           = 1.0;
    wrap.MethodProps().nSteps       = 5;
    wrap.MethodProps().dryTol       = 1e-12;
    wrap.MethodProps().newtonMaxIts = 20;
    wrap.MethodProps().newtonTol    = 1e-10;
    wrap.SetNewtonVerbose(true);

    // IC selection: if not uniform, provide nonuniform initial condition functions
    if (!use_uniform_ic) {
      wrap.SetInitialEta1D([&](std::size_t, double x){
        const double xmid = 0.5 * (N - 1) * wrap.MethodProps().dx;
        return (x < xmid) ? 0.2 : 0.0;
      });
      wrap.SetInitialUFace1D([](std::size_t, double){ return 0.0; });
    }

    wrap.SetPDE(std::move(pde));

    const int err = wrap.NumericalSolve();
    if (err) {
      std::cerr << "1D Solve failed with code " << err << "\n";
      return err;
    }

    const std::size_t s_eta = clamp_stride(wrap.Eta1D().size(), 14);
    const std::size_t s_u   = clamp_stride(wrap.UFace1D().size(), 14);

    print_1d_sample("eta", wrap.Eta1D(), s_eta);
    print_1d_sample("u",   wrap.UFace1D(), s_u);
  }
#endif

#if defined(CASE_2D)
  {
    std::cout << "=== 2D CASE (SWE2dcas) ===\n";
    const std::size_t nx = N, ny = N;
    std::cout << "nx=" << nx << " ny=" << ny << "  ic=" << ic_mode << "\n";

    auto pde = std::make_unique<fcfd::pdemodel::SWE2dcas<Real, MyPDEOptions>>();

    // Model params (your SWE2dcas provides ModelParams())
    pde->ModelParams().g      = 9.81;
    pde->ModelParams().nm     = 0.03;
    pde->ModelParams().radi   = 1.0;
    pde->ModelParams().gammat = 0.1;

    pde->SetIC(0.0, 0.0, 0.0);   // (eta0,u0,v0)
    pde->SetBed(1.0, 0.0, 0.0);
    pde->SetWind(0.0, 0.0);

    // BC
    fcfd::pdemodel::CasulliBC2D<Real> bc;
    bc.left   = fcfd::pdemodel::SubcriticalBC2D::StageEta;
    bc.right  = fcfd::pdemodel::SubcriticalBC2D::StageEta;
    bc.bottom = fcfd::pdemodel::SubcriticalBC2D::Velocity;
    bc.top    = fcfd::pdemodel::SubcriticalBC2D::Velocity;

    bc.left_value   = 0.0;
    bc.right_value  = 0.0;
    bc.bottom_value = 0.0;
    bc.top_value    = 0.0;
    pde->SetBC(bc);

    fcfd::pdenumerics::CasulliWrap<Real, MyPDEOptions> wrap;
    wrap.SetGrid2D(nx, ny);

    wrap.MethodProps().dt           = 0.01;
    wrap.MethodProps().dx           = 1.0;
    wrap.MethodProps().dy           = 1.0;
    wrap.MethodProps().nSteps       = 5;
    wrap.MethodProps().dryTol       = 1e-12;
    wrap.MethodProps().newtonMaxIts = 20;
    wrap.MethodProps().newtonTol    = 1e-10;
    wrap.MethodProps().cgMaxIts     = 200;
    wrap.MethodProps().cgTol        = 1e-10;
    wrap.SetNewtonVerbose(true);

    if (!use_uniform_ic) {
      wrap.SetInitialEta2D([&](std::size_t, std::size_t, double x, double y){
        const double xc = 0.5 * (nx - 1) * wrap.MethodProps().dx;
        const double yc = 0.5 * (ny - 1) * wrap.MethodProps().dy;
        const double r2 = (x - xc) * (x - xc) + (y - yc) * (y - yc);
        return 0.2 * std::exp(-r2 / 50.0);
      });
    }

    wrap.SetPDE(std::move(pde));

    const int err = wrap.NumericalSolve();
    if (err) {
      std::cerr << "2D Solve failed with code " << err << "\n";
      return err;
    }

    const std::size_t stride_i = clamp_stride(nx, 12);
    const std::size_t stride_j = clamp_stride(ny, 12);

    print_2d_sample_side_by_side(
      nx, ny,
      wrap.MethodProps().dx,
      wrap.MethodProps().dy,
      wrap.Eta2D(),
      wrap.UFace2D(),
      wrap.VFace2D(),
      stride_i,
      stride_j
    );

    std::cout << "\nTip: to print FULL fields, reduce N or change clamp_stride target_points.\n";
  }
#endif

  return 0;
}
