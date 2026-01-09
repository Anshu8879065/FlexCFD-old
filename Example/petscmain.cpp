#include <cstdlib>
#include <iostream>
#include <memory>
#include <optional>

#include "Numerics/petscwrap/PetscInitException.hpp"
#include "Numerics/petscwrap/PetscWrap.hpp"
#include "PDEs/ModelParams.hpp"
#include "PDEs/SVE/SVE.hpp"
#include "PDEs/SWE/SWE2d.hpp"

using fcfd::pdemodel::PDESystem;
using fcfd::pdemodel::SVE;
using fcfd::pdemodel::SWE2d;
using fcfd::pdenumerics::PetscWrap;

namespace
{
void PetscWrap1d(int const argc, char* argv[])
{
  using fcfd::pdemodel::Model1dParams;
  using fcfd::pdenumerics::PetscInitException;

  static constexpr char help[] = "Coupled Saint-Venant System.  Option prefix -ptn_.\n"
                                   "Incorporates form  F(t,Y,dot Y) = G(t,Y)  where F() is IFunction and G() is\n"
                                   "RHSFunction().  Implements IJacobian() and RHSJacobian().  Defaults to\n"
                                   "ARKIMEX (= adaptive Runge-Kutta implicit-explicit) TS type.\n\n";

  auto sve = std::make_unique<SVE<double, Model1dParams<double>>>();

  try {
    auto pwrap = PetscWrap<Model1dParams<double>>(argc, argv, help);
    pwrap.SetPDE(std::move(sve));

    PetscInt L = 10;
    PetscReal Nx = 200;
    PetscReal dx = L / (Nx - 1);

    // First-derivative stencil weights for [i-1, i, i+1]
    std::optional<std::vector<PetscReal>> stengrid =
      std::vector<PetscReal> {-PetscReal(0.5) / dx, PetscReal(0.0), PetscReal(0.5) / dx};
    std::vector<PetscReal> stengridd = {
      1.0 / (dx * dx),  // u_{i-1}
      -2.0 / (dx * dx),  // u_i
      1.0 / (dx * dx)  // u_{i+1}
    };

    // Not used in code path -> nullopt
    std::optional<std::vector<PetscReal>> stengridd3 = std::nullopt;

    pwrap.SetGrid(1,
      {DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED},
      {static_cast<PetscInt>(Nx), 2, 1},
      stengrid,
      stengridd,
      stengridd3);
    pwrap.SetFields({"h", "qx"});
    pwrap.InitializeTS();
    pwrap.MakeProbFuns();
    pwrap.SetSolOpts();
    pwrap.NumericalSolve();

    // ---- print solution ----
    auto interpol = [](std::vector<double> v) -> std::vector<double>
    {
      return v;  // identity: prints h and qx
    };
    pwrap.EvaluateSolution(interpol);
    // ------------------------
  }
  catch (PetscInitException const& e) {
    std::clog << "Exception: " << e.what() << '\n';
  }
}

void PetscWrap2d(int const argc, char* argv[])
{
  using fcfd::pdemodel::Model2dParams;
  using fcfd::pdenumerics::PetscInitException;

  static constexpr char help[] = "2D Shallow Water Equations.  Option prefix -ptn_.\n"
                                   "Incorporates form  F(t,Y,dot Y) = G(t,Y)  where F() is IFunction and G() is\n"
                                   "RHSFunction().  Implements IJacobian() and RHSJacobian().  Defaults to\n"
                                   "ARKIMEX (= adaptive Runge-Kutta implicit-explicit) TS type.\n\n";

  auto swe2d = std::make_unique<SWE2d<double, Model2dParams<double>>>();

  try {
    auto pwrap = PetscWrap<Model2dParams<double>>(argc, argv, help);
    pwrap.SetPDE(std::move(swe2d));

    PetscInt L = 10;
    PetscReal Nx = 200;  // Nx == Ny
    PetscReal dx = L / (Nx - 1);  // dx==dy

    std::optional<std::vector<PetscReal>> stengrid =
      std::vector<PetscReal> {-PetscReal(0.5) / dx, PetscReal(0.0), PetscReal(0.5) / dx};

    std::optional<std::vector<PetscReal>> stengridd =
      std::vector<PetscReal> {PetscReal(1.0) / (dx * dx), -PetscReal(2.0) / (dx * dx), PetscReal(1.0) / (dx * dx)};

    std::optional<std::vector<PetscReal>> stengridd3 = std::nullopt;

    pwrap.SetGrid(2,
      {DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_GHOSTED},
      {static_cast<PetscInt>(Nx), static_cast<PetscInt>(Nx), 3, 1},
      stengrid,
      stengridd,
      stengridd3);
    pwrap.SetFields({"h", "qx", "qy"});
    pwrap.InitializeTS();
    pwrap.MakeProbFuns();
    pwrap.SetSolOpts();
    pwrap.NumericalSolve();

    // ---- print solution ----
    auto interpol = [](std::vector<double> v) -> std::vector<double>
    {
      return v;  // identity: prints h, qx, qy
    };
    pwrap.EvaluateSolution(interpol);
    // ------------------------
  }
  catch (PetscInitException const& e) {
    std::clog << "Exception: " << e.what() << '\n';
  }
}

}  // namespace

auto main(int const argc, char* argv[]) -> int
{
  PetscWrap1d(argc, argv);

  return EXIT_SUCCESS;
}
