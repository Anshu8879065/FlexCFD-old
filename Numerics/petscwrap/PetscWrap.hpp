#pragma once

#include <array>
#include <functional>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include <petsc.h>

#include "Numerics/SolTools.hpp"
#include "PetscInitException.hpp"
#include "PetscOptions.hpp"

namespace fcfd::pdenumerics
{

///
/// \brief Check if a given PetscErrorCode is an error code
/// \param[in] code The code in question
/// \returns True if the code is an error code, and false otherwise
///
constexpr auto IsPetscError(PetscErrorCode const code) noexcept -> bool
{
  return code != 0;
}

struct CallBacks
{
  PetscBool iFuncCalled;
  PetscBool iJacobianCalled;
  PetscBool rhsFuncCalled;
  PetscBool rhsJacobianCalled;
};

///
/// \brief Fill in ghost cells in the x, y, and z directions
/// \param[in] info DMDALocalInfo struct containing grid indices and dimensions
/// \param[in] fieldArray Array of field variables
/// \param[in] numFields The number of field variable being processed
/// \param[in] netWidth The radius of the grid cell
///
void FillGhostCells(DMDALocalInfo const& info, std::span<PetscReal> fieldArray, std::size_t numFields, int netWidth);

///
/// \brief Fill in ghost cells in a 1d grid
/// \param[in] distributedMesh The DM object
/// \param[in] info The DMDALocalInfo struct containing grid indices and dimensions
/// \returns A PetscErrorCode indicating success or failure
///
auto FillGhostCells1d(DM const& distributedMesh, DMDALocalInfo const& info) -> PetscErrorCode;

///
/// \brief Fill in ghost cells in a 2d grid
/// \param[in] distributedMesh The DM object
/// \param[in] info The DMDALocalInfo struct containing grid indices and dimensions
/// \returns A PetscErrorCode indicating success or failure
///
auto FillGhostCells2d(DM const& distributedMesh, DMDALocalInfo const& info) -> PetscErrorCode;

///
/// \brief Fill in ghost cells in a 3d grid
/// \param[in] distributedMesh The DM object
/// \param[in] info The DMDALocalInfo struct containing grid indices and dimensions
/// \returns A PetscErrorCode indicating success or failure
///
auto FillGhostCells3d(DM const& distributedMesh, DMDALocalInfo const& info) -> PetscErrorCode;

///
/// \brief Create a 1d stencil grid
/// \param dmPtr The data manager object from which the grid is to be created
/// \param stencilData Stencil data (i.e., global dimensions in x direction, degrees of freedom per node, stencil width)
/// \param boundaryCond The type of boundary condition to apply on the grid
/// \return A PetscErrorCode indicating success or failure
///
auto Create1dStencilGrid(DM& dmPtr, std::span<PetscInt> stencilData, DMBoundaryType boundaryCond) -> PetscErrorCode;

///
/// \brief Create a 2d stencil grid
/// \param dmPtr The data manager object from which the grid is to be created
/// \param stencilData Stencil data (i.e., global dimensions in x and y directions, degrees of freedom per node, stencil
/// width)
/// \param boundaryCond The type of boundary condition to apply on the grid
/// \return A PetscErrorCode indicating success or failure
///
auto Create2dStencilGrid(DM& dmPtr, std::span<PetscInt> stencilData, DMBoundaryType boundaryCond) -> PetscErrorCode;

///
/// \brief Create a 3d stencil grid
/// \param dmPtr The data manager object from which the grid is to be created
/// \param stencilData Stencil data (i.e., global dimensions in x, y, and z directions, degrees of freedom per node,
/// stencil width)
/// \param boundaryCond The type of boundary condition to apply on the grid
/// \return A PetscErrorCode indicating success or failure
///
auto Create3dStencilGrid(DM& dmPtr, std::span<PetscInt> stencilData, DMBoundaryType boundaryCond) -> PetscErrorCode;

///
/// \brief Set the initial state of the PDE system
/// \param[in] distributedMesh The DM object
/// \param[in] vec The Vec object
/// \param[in] info The DMDALocalInfo object
/// \param[in] initFunc A function that accepts grid indices (i[,j[,k]]) as a vector of PetscReal and
/// returns a vector containing the initial DOF values at that location (vector length must equal the DMDA DOF)
/// \returns A PetscErrorCode indicating success or failure
///
auto SetInitialState(DM& distributedMesh,
  Vec& vec,
  DMDALocalInfo& info,
  std::function<std::vector<PetscReal>(std::vector<PetscReal> const&)> const& initFunc) -> PetscErrorCode;

///
/// \brief Set the initial state of the PDE system (1D)
/// \param[in] distributedMesh The DM object
/// \param[in] vec The Vec object
/// \param[in] info The DMDALocalInfo object
/// \param[in] initFunc A function that accepts grid index {i} (as a vector {i}, PetscReal)
/// and returns a vector containing the initial DOF values at that location (vector length must equal the DMDA DOF)
/// \returns A PetscErrorCode indicating success or failure
///
auto SetInitialState1d(DM const& distributedMesh,
  Vec const& vec,
  DMDALocalInfo const& info,
  std::function<std::vector<PetscReal>(std::vector<PetscReal> const&)> const& initFunc) -> PetscErrorCode;

///
/// \brief Set the initial state of the PDE system (2D)
/// \param[in] distributedMesh The DM object
/// \param[in] vec The Vec object
/// \param[in] info The DMDALocalInfo object
/// \param[in] initFunc A function that accepts grid indices {i,j} (PetscReal) and returns
/// a vector containing the initial DOF values at that location (vector length must equal the DMDA DOF)
/// \returns A PetscErrorCode indicating success or failure
///
auto SetInitialState2d(DM const& distributedMesh,
  Vec const& vec,
  DMDALocalInfo const& info,
  std::function<std::vector<PetscReal>(std::vector<PetscReal> const&)> const& initFunc) -> PetscErrorCode;

///
/// \brief Set the initial state of the PDE system (3D)
/// \param[in] distributedMesh The DM object
/// \param[in] vec The Vec object
/// \param[in] info The DMDALocalInfo object
/// \param[in] initFunc A function that accepts grid indices {i,j,k} (PetscReal) and returns
/// a vector containing the initial DOF values at that location (vector length must equal the DMDA DOF)
/// \returns A PetscErrorCode indicating success or failure
///
auto SetInitialState3d(DM const& distributedMesh,
  Vec const& vec,
  DMDALocalInfo const& info,
  std::function<std::vector<PetscReal>(std::vector<PetscReal> const&)> const& initFunc) -> PetscErrorCode;

///
/// \brief Evaluate the solution at a certain grid point using the given function
/// \param[in] distributedMesh The DM object
/// \param[in] globalVec The Vec object
/// \param[in] info The DMDALocalInfo object
/// \param[in] func The function which will be used to process the solution at a certain grid point
/// \returns A PetscErrorCode indicating success or failure
///
auto EvaluateSolution(DM const& distributedMesh,
  Vec const& globalVec,
  DMDALocalInfo const& info,
  std::function<std::vector<PetscReal>(std::vector<PetscReal>&)> const& func) -> PetscErrorCode;

///
/// \brief Evaluate the solution at a certain grid point using the given function
/// \param[in] distributedMesh The DM object
/// \param[in] globalVec The Vec object
/// \param[in] info The DMDALocalInfo object
/// \param[in] func The function which will be used to process the solution at a certain grid point
/// \returns A PetscErrorCode indicating success or failure
///
auto EvaluateSolution1d(DM const& distributedMesh,
  Vec const& globalVec,
  DMDALocalInfo const& info,
  std::function<std::vector<PetscReal>(std::vector<PetscReal>&)> const& func) -> PetscErrorCode;

///
/// \brief Evaluate the solution at a certain grid point using the given function
/// \param[in] distributedMesh The DM object
/// \param[in] globalVec The Vec object
/// \param[in] info The DMDALocalInfo object
/// \param[in] func The function which will be used to process the solution at a certain grid point
/// \returns A PetscErrorCode indicating success or failure
///
auto EvaluateSolution2d(DM const& distributedMesh,
  Vec const& globalVec,
  DMDALocalInfo const& info,
  std::function<std::vector<PetscReal>(std::vector<PetscReal>&)> const& func) -> PetscErrorCode;

///
/// \brief Evaluate the solution at a certain grid point using the given function
/// \param[in] distributedMesh The DM object
/// \param[in] globalVec The Vec object
/// \param[in] info The DMDALocalInfo object
/// \param[in] func The function which will be used to process the solution at a certain grid point
/// \returns A PetscErrorCode indicating success or failure
///
auto EvaluateSolution3d(DM const& distributedMesh,
  Vec const& globalVec,
  DMDALocalInfo const& info,
  std::function<std::vector<PetscReal>(std::vector<PetscReal>&)> const& func) -> PetscErrorCode;

template<std::semiregular PDEOptions>
class PetscWrap : public SolTool<PetscReal, PDEOptions, PetscErrorCode, PetscOptions, PetscSolveOpts>
{
  using BaseClass = SolTool<PetscReal, PDEOptions, PetscErrorCode, PetscOptions, PetscSolveOpts>;

public:
  ///
  /// \brief Constructor
  ///
  /// Calls `PetscInitialize()` and sets RHS and I functions
  /// \param[in] argc The number of command-line arguments passed in to the program
  /// \param[in] argv A string containing a list of command-line arguments
  /// \param[in] msg A help message string describing the possible options to set to modify the
  /// program's behaviour
  /// \throws PetscInitException on failure of PetscInitialize or the setting of RHS and I functions
  ///
  PetscWrap(int argc, char** argv, const char* msg)
  {
    if (IsPetscError(PetscInitialize(&argc, &argv, nullptr, msg))) {
      throw PetscInitException();
    }
  }

  ///
  /// \brief Constructor
  /// \param[in] argc The number of command-line arguments passed in to the program
  /// \param[in] argv A string containing a list of the command-line arguments passed to the program
  /// \param[in] msg A help message string describing the possible options to set to modify the
  /// program's behaviour
  /// \param[in] numFields The number of unknown fields in the PDE system
  /// \throws PetscInitException on failure
  ///
  PetscWrap(int argc, char** argv, const char* msg, int numFields)
    : PetscWrap(argc, argv, msg)
  {
    m_methodProps.numFields = numFields;
    BaseClass::SetInitFlag(0, FlagState::IsSet);
  }

  ///
  /// \brief Constructor
  /// \param[in] argc The number of command-line arguments passed to the program
  /// \param[in] argv A string containing a list of the command-line arguments passed to the program
  /// \param[in] msg A help message string describing the possible options to set to modify the
  /// program's behaviour
  /// \param[in] tsType A string with the name of a PETSc TS method (the time/ODE integrators that
  /// PETSc provides)
  /// \param[in] numFields The number of unknown fields in the PDE system
  /// \throws PetscInitException on failure
  ///
  PetscWrap(int argc, char** argv, const char* msg, TSType const tsType, int const numFields)
    : PetscWrap(argc, argv, msg)
  {
    m_methodProps.numFields = numFields;
    m_methodProps.tsType = tsType;
    BaseClass::SetInitFlag(0, FlagState::IsSet);
  }

  ///
  /// \brief Constructor
  /// \param[in] argc The number of command-line arguments passed to the program
  /// \param[in] argv A string containing a list of the command-line arguments passed to the program
  /// \param[in] msg A help message string describing the possible options to set to modify the
  /// program's behaviour
  /// \param[in] tsType A string with the name of a PETSc TS method (the time/ODE integrators that
  /// PETSc provides)
  /// \param[in] fieldNames The names of the unknown fields in the PDE system
  /// \throws PetscInitException on failure
  ///
  PetscWrap(int argc, char** argv, const char* msg, TSType const tsType, std::vector<std::string> fieldNames)
    : PetscWrap(argc, argv, msg)
  {
    m_methodProps.tsType = tsType;
    SetFields(std::move(fieldNames));
    BaseClass::SetInitFlag(0, FlagState::IsSet);
  }

  ///
  /// \brief Constructor
  /// \param[in] argc The number of command-line arguments passed to the program
  /// \param[in] argv A string containing a list of the command-line arguments passed to the program
  /// \param[in] msg A help message string describing the possible options to set to modify the
  /// program's behaviour
  /// \param[in] fieldNames The names of the unknown fields in the PDE system
  /// \throws PetscInitException on failure
  ///
  PetscWrap(int argc, char** argv, const char* msg, std::vector<std::string> fieldNames)
    : PetscWrap(argc, argv, msg)
  {
    SetFields(std::move(fieldNames));
    BaseClass::SetInitFlag(0, FlagState::IsSet);
  }

  ///
  /// \brief Copy constructor
  ///
  PetscWrap(PetscWrap const&) = default;

  ///
  /// \brief Copy-assignment operator
  ///
  auto operator=(PetscWrap const&) -> PetscWrap& = default;

  ///
  /// \brief Move constructor
  ///
  PetscWrap(PetscWrap&&) noexcept = default;

  ///
  /// \brief Move-assignment operator
  ///
  auto operator=(PetscWrap&&) noexcept -> PetscWrap& = default;

  ///
  /// \brief Destructor
  /// \details Calls PetscFinalize()
  ///
  ~PetscWrap()
  {
    if (m_vec != nullptr) {
      VecDestroy(&m_vec);
    }
    if (m_myTS != nullptr) {
      TSDestroy(&m_myTS);
    }
    if (m_da != nullptr) {
      DMDestroy(&m_da);
    }

    PetscFinalize();
  }

  ///
  /// \brief Write the names of the unknown fields in the PDE system to our DMDA instance
  /// \param[in] fieldNames The names of the unknown fields in the PDE system
  /// \returns A PetscErrorCode instance on success or failure
  ///
  auto SetFields(std::vector<std::string>&& fieldNames) -> PetscErrorCode
  {
    m_methodProps.fields = std::move(fieldNames);
    m_methodProps.numFields = std::size(m_methodProps.fields);

    return SetFields();
  }

  ///
  /// \brief Set the solution options (an instance of PetscSolveOpts) to be used
  /// This function is currently no-op because
  /// \returns PetscErrorCode indicating success or failure
  ///
  auto InitSolMethod() -> PetscErrorCode override
  {
    return 0;
  }

  ///
  /// \brief Set the PetscOptions for this current instance of PetscWrap
  /// \param options The options to be set
  ///
  void InitSolMethod(PetscOptions options)
  {
    m_methodProps = std::move(options);
    BaseClass::SetInitFlag(0, FlagState::IsSet);
  }

  ///
  /// \brief Creates a DMDA instance depending on the number of dimensions we're working with
  ///
  void SetGrid()
  {
  }

  ///
  /// \brief Creates a DMDA instance depending on the number of dimensions we're working with
  /// \param[in] ndim The number of dimensions
  /// \param[in] bdryConds The boundary conditions
  /// \param[in] stendata Stencil data (i.e., global dimensions in x,y,z directions; degrees of freedom per node,
  /// stencil width)
  /// \param[in] stengrid Stencil grid for first-order derivatives
  /// \param[in] stengridd Stencil grid for second-order derivatives
  /// \param[in] stengridd3 Stencil grid for third-order derivatives
  /// \returns A PetscErrorCode detailing success or failure
  ///
  auto SetGrid(int ndim,
    std::array<DMBoundaryType, 3> const& bdryConds,
    std::vector<PetscInt> stendata,
    std::optional<std::vector<PetscReal>> stengrid,
    std::optional<std::vector<PetscReal>> stengridd,
    std::optional<std::vector<PetscReal>> stengridd3) -> PetscErrorCode
  {
    BaseClass::m_nd = ndim;

    if (stengrid.has_value()) {
      m_myStencilGridd = std::move(stengrid);
    }

    if (stengridd.has_value()) {
      m_myStencilGridd = std::move(stengridd);
    }

    if (stengridd3.has_value()) {
      m_myStencilGridd = std::move(stengridd3);
    }

    m_inTmp = std::vector<PetscReal>(stendata[1], 0.0);
    m_outTmp = std::vector<PetscReal>(stendata[1], 0.0);

    if (ndim == 1) {
      Create1dStencilGrid(m_da, stendata, bdryConds[0]);
    }
    else if (ndim == 2) {
      Create2dStencilGrid(m_da, stendata, bdryConds[1]);
    }
    else if (ndim == 3) {
      Create3dStencilGrid(m_da, stendata, bdryConds[2]);
    }

    m_myStencilData = std::move(stendata);
    PetscCall(DMSetFromOptions(m_da));
    PetscCall(DMSetUp(m_da));

    BaseClass::SetInitFlag(2, FlagState::IsSet);

    return 0;
  }

  ///
  /// \brief Set the four local residual evaluation functions for use with the DMDA
  /// \returns A PetscErrorCode indicating success or failure
  ///
  auto MakeProbFuns() -> PetscErrorCode
  {
    assert(m_da != nullptr and "m_da is in an invalid state and must be initialized before use");

    PetscCall(DMDATSSetRHSFunctionLocal(m_da, INSERT_VALUES, (DMDATSRHSFunctionLocal)FormRHSFunctionLocal, this));
    PetscCall(DMDATSSetRHSJacobianLocal(m_da, (DMDATSRHSJacobianLocal)FormRHSJacobianLocal, this));
    PetscCall(DMDATSSetIFunctionLocal(m_da, INSERT_VALUES, (DMDATSIFunctionLocal)FormIFunctionLocal, this));
    PetscCall(DMDATSSetIJacobianLocal(m_da, (DMDATSIJacobianLocal)FormIJacobianLocal, this));

    return 0;
  }

  auto SetSolOpts() -> PetscErrorCode override
  {
    assert(m_myTS != nullptr and "m_myTS is in an invalid state and must be initialized before use");

    // if (BaseClass::GetInitFlag(3) == FlagState::IsSet) {
    PetscCall(TSSetType(m_myTS, m_algOpts.tsType));
    PetscCall(TSSetTime(m_myTS, m_algOpts.initTime));
    PetscCall(TSSetMaxTime(m_myTS, m_algOpts.maxTime));
    PetscCall(TSSetTimeStep(m_myTS, m_algOpts.timeStep));
    PetscCall(TSSetExactFinalTime(m_myTS, m_algOpts.finalTimeOption));
    PetscCall(TSSetFromOptions(m_myTS));
    // }

    return 0;
  }

  ///
  /// \brief Set solution options
  /// \param[in] opts The solution options to be used
  //
  void SetSolOpts(PetscSolveOpts opts) override
  {
    BaseClass::SetSolOpts(opts);
    SetSolOpts();
  }

  ///
  /// \brief Set the TS scheme to be used
  /// \param[in] tsScheme The TS scheme to set
  ///
  void SetTSScheme(TS tsScheme)
  {
    m_myTS = tsScheme;
  }

  ///
  /// \brief Set the initial state of the PDE system
  /// \param[in] distributedMesh The DM object
  /// \param[in] vec The Vec object
  /// \param[in] info The DMDALocalInfo object
  /// \returns A PetscErrorCode indicating success or failure
  ///
  auto InitialState(DM& distributedMesh, Vec& vec, DMDALocalInfo& info) -> PetscErrorCode
  {
    assert(BaseClass::m_myPDE != nullptr && "PDE has not been set");

    PetscFunctionBeginUser;

    auto initFunc = [this](std::vector<PetscReal> const& coords) -> std::vector<PetscReal>
    { return this->BaseClass::m_myPDE->InitCond(coords); };

    return SetInitialState(distributedMesh, vec, info, initFunc);
  }

  ///
  /// \brief Initialize the timestepper
  /// \returns A PetscErrorCode indicating success or failure
  ///
  auto InitializeTS() -> PetscErrorCode
  {
    assert(m_da and "The DM object is in an invalid state. SetGrid() must be called before InitializeTS()");

    PetscFunctionBeginUser;

    PetscCall(TSCreate(PETSC_COMM_WORLD, &m_myTS));
    PetscCall(TSSetProblemType(m_myTS, TS_NONLINEAR));
    PetscCall(TSSetDM(m_myTS, m_da));
    PetscCall(TSSetApplicationContext(m_myTS, this));

    PetscFunctionReturn(PETSC_SUCCESS);
  }

  ///
  /// \brief
  /// \param[in] dmPtr The DM object
  /// \param[in] vec A vector the same size as one obtained with DMCreateGlobalVector or DMCreateLocalVector
  /// \returns PetscErrorCode indicating success or failure
  auto NumericalSolve() -> PetscErrorCode override
  {
    assert(m_da and "The DM object is in an invalid state and must be initialized before use");
    assert(m_myTS and "The TS object is in an invalid state and must be initialized before use");

    PetscFunctionBeginUser;

    if (m_vec == nullptr) {
      PetscCall(DMCreateGlobalVector(m_da, &m_vec));
    }

    PetscCall(InitialState(m_da, m_vec, m_info));
    PetscCall(TSSolve(m_myTS, m_vec));

    PetscFunctionReturn(PETSC_SUCCESS);
  }

     ///
  /// \brief Calls a function on values of fields of a PDE and stores the evaluated result
  ///
  /// \details
  /// This function evaluates the solution stored in the global vector at each local grid
  /// point owned by the current MPI rank. The evaluation is performed by calling the
  /// dimension-specific free functions EvaluateSolution1d/2d/3d, which iterate over the
  /// DMDA local domain.
  ///
  /// The values returned by the user-provided function are stored in m_evalValues in the
  /// same order as the grid traversal. The corresponding grid indices (i,j,k) are stored
  /// in m_evalIJK. No assumptions are made about primitive variables; the stored values
  /// correspond directly to the conservative degrees of freedom defined by the PDE.
  ///
  /// \param[in] interpol The function used to evaluate the fields of a PDE at a grid point
  /// \returns PetscErrorCode indicating success or failure
  ///
  auto EvaluateSolution(std::function<std::vector<double>(std::vector<PetscReal>&)> const& interpol)
    -> PetscErrorCode override
  {
    assert(m_da != nullptr && "DM object is not initialized");
    assert(m_vec != nullptr && "Solution vector is null. Call NumericalSolve() before evaluation");

    PetscFunctionBeginUser;

    // Clear previously stored evaluation results
    m_evalValues.clear();
    m_evalIJK.clear();

    // Precompute grid indices in standard DMDA local ordering so that indices
    // correspond exactly to the order in which the solution is evaluated
    if (m_info.dim == 1) {
      m_evalIJK.reserve(static_cast<std::size_t>(m_info.xm));
      for (PetscInt i = m_info.xs; i < m_info.xs + m_info.xm; ++i) {
        m_evalIJK.push_back({i, 0, 0});
      }
    }

    if (m_info.dim == 2) {
      m_evalIJK.reserve(static_cast<std::size_t>(m_info.xm) * static_cast<std::size_t>(m_info.ym));
      for (PetscInt j = m_info.ys; j < m_info.ys + m_info.ym; ++j) {
        for (PetscInt i = m_info.xs; i < m_info.xs + m_info.xm; ++i) {
          m_evalIJK.push_back({i, j, 0});
        }
      }
    }

    if (m_info.dim == 3) {
      m_evalIJK.reserve(static_cast<std::size_t>(m_info.xm)
        * static_cast<std::size_t>(m_info.ym)
        * static_cast<std::size_t>(m_info.zm));
      for (PetscInt k = m_info.zs; k < m_info.zs + m_info.zm; ++k) {
        for (PetscInt j = m_info.ys; j < m_info.ys + m_info.ym; ++j) {
          for (PetscInt i = m_info.xs; i < m_info.xs + m_info.xm; ++i) {
            m_evalIJK.push_back({i, j, k});
          }
        }
      }
    }

    // Wrapper around user function to store evaluated values while preserving
    // the required function signature for EvaluateSolution1d/2d/3d
    auto storeFunc =
      [this, &interpol](std::vector<PetscReal>& dofs) -> std::vector<PetscReal>
    {
      std::vector<double> outD = interpol(dofs);

      std::vector<PetscReal> out;
      out.reserve(outD.size());
      for (double v : outD) {
        out.push_back(static_cast<PetscReal>(v));
      }

      m_evalValues.push_back(out);
      return out;
    };

    PetscErrorCode ierr = PETSC_SUCCESS;

    // Dispatch to the dimension-specific evaluation routine
    if (m_info.dim == 1) {
      ierr = fcfd::pdenumerics::EvaluateSolution1d(m_da, m_vec, m_info, storeFunc);
    }

    if (m_info.dim == 2) {
      ierr = fcfd::pdenumerics::EvaluateSolution2d(m_da, m_vec, m_info, storeFunc);
    }

    if (m_info.dim == 3) {
      ierr = fcfd::pdenumerics::EvaluateSolution3d(m_da, m_vec, m_info, storeFunc);
    }

    PetscFunctionReturn(ierr);
  }

private:
  PetscOptions m_methodProps {};
  PetscSolveOpts m_algOpts {};
  CallBacks m_user {.iFuncCalled = PETSC_FALSE,
    .iJacobianCalled = PETSC_FALSE,
    .rhsFuncCalled = PETSC_FALSE,
    .rhsJacobianCalled = PETSC_FALSE};

  //! The time stepper we're using
  TS m_myTS {};
  DM m_da {};
  DMDALocalInfo m_info {};
  Vec m_vec {};
  PetscBool m_noRhsJacobian {PETSC_FALSE};
  PetscBool m_noIJacobian {PETSC_FALSE};
  PetscBool m_callbackReport {PETSC_TRUE};

  //! This stores our stencil data (i.e., the global dimensions in x, y, and z directions, degrees of freedom per node,
  //! stencil width)
  std::vector<PetscInt> m_myStencilData {};

  std::optional<std::vector<PetscReal>> m_myStencilGridd {};
  std::optional<std::vector<PetscReal>> m_myStencilGriddd {};
  std::optional<std::vector<PetscReal>> m_myStencilGridd3 {};

    //! Stores the last evaluation output per visited grid point (local portion on this rank)
  std::vector<std::vector<PetscReal>> m_evalValues {};

  //! Stores the corresponding (i,j,k) indices (local portion on this rank)
  std::vector<std::array<PetscInt, 3>> m_evalIJK {};


  std::vector<PetscReal> m_inTmp {};

  //! This stores the first-order spatial derivatives
  std::vector<PetscReal> m_indTmp {};

  //! This stores the second-order spatial derivatives
  //! Optional because these derivatives may not exist depending on the PDE we're trying to solve
  std::optional<std::vector<PetscReal>> m_ind2Tmp {};

  //! This stores the third-order spatial derivatives
  //! Optional because these derivatives may not exist depending on the PDE we're trying to solve
  std::optional<std::vector<PetscReal>> m_ind3Tmp {};

  //! This is used to store the computed time derivatives
  std::vector<PetscReal> m_dTime {};

  std::vector<PetscReal> m_outTmp {};
  std::vector<PetscReal> m_outdTmp {};
  std::optional<std::vector<PetscReal>> m_outd2Tmp {};
  std::optional<std::vector<PetscReal>> m_outd3Tmp {};

  std::vector<PetscReal> m_vTmp {};
  std::vector<PetscReal> m_vTmpTmp {};
  std::vector<PetscReal> m_vdTmp {};
  std::optional<std::vector<PetscReal>> m_vd2Tmp {};
  std::optional<std::vector<PetscReal>> m_vd3Tmp {};

  std::vector<MatStencil> m_colTmp {};
  std::vector<MatStencil> m_coldTmp {};
  std::optional<std::vector<MatStencil>> m_cold2Tmp {};
  std::optional<std::vector<MatStencil>> m_cold3Tmp {};

  //!
  MatStencil m_rowTmp {};
  MatStencil m_rowdTmp {};
  std::optional<MatStencil> m_rowd2Tmp {};
  std::optional<MatStencil> m_rowd3Tmp {};

  ///
  /// \brief Set field names for the DMDA instance
  ///
  auto SetFields() -> PetscErrorCode
  {
    assert(m_da and "The DM object is in an invalid state. SetGrid() must be called before SetFields()");

    PetscFunctionBeginUser;

    for (PetscInt i = 0; auto const& field : m_methodProps.fields) {
      PetscCall(DMDASetFieldName(m_da, i, field.c_str()));
      ++i;
    }

    BaseClass::SetInitFlag(0, FlagState::IsSet);

    PetscFunctionReturn(PETSC_SUCCESS);
  }

  static auto FormRHSFunctionLocal(DMDALocalInfo* info, PetscReal time, PetscReal* aY, PetscReal* aG, void* pwp)
    -> PetscErrorCode;
  // ENDRHSFUNCTION

  static auto FormRHSJacobianLocal(
    DMDALocalInfo* info, PetscReal time, PetscReal* aY, Mat asmMat, Mat precondMat, void* pwp) -> PetscErrorCode;

  // in system form  F(t,Y,dot Y) = G(t,Y),  compute F():
  //     F^u(t,a,q,a_t,q_t) = a_t +  Dx[q]
  //     F^v(t,a,q,a_t,q_t) = q_t +  Dx[q^2/a+ga^2/2] = q_t + 2q Dx[q]/a-q^2/a^2
  //     Dx[a]+g a Dx[a]
  // STARTIFUNCTION
  static auto FormIFunctionLocal(
    DMDALocalInfo* info, PetscReal time, PetscReal* aY, PetscReal* aYdot, PetscReal* aF, void* pwp) -> PetscErrorCode;
  // ENDIFUNCTION

  // in system form  F(t,Y,dot Y) = G(t,Y),  compute combined/shifted
  // Jacobian of F():
  //     J = (shift) dF/d(dot Y) + dF/dY
  // STARTIJACOBIAN
  static auto FormIJacobianLocal(DMDALocalInfo* info,
    PetscReal time,
    PetscReal* aY,
    PetscReal* aYdot,
    PetscReal shift,
    Mat asmMat,
    Mat precondMat,
    void* pwp) -> PetscErrorCode;
};

///
/// \brief Computes the local right-hand side function for the shallow water equations
///
/// \param[in] info DMDALocalInfo structure containing grid dimensions and local indices for this
/// processor's domain
/// \param[in] time Current simulation time
/// \param[in] aY Array of flow variables at the current time step
/// \param[in] aG Array of computed time derivatives for flow variables
/// \param[in] pwp A void pointer to an instance of PetscWrap (the context)
/// \returns
///
template<std::semiregular PDEOptions>
auto PetscWrap<PDEOptions>::FormRHSFunctionLocal(
  DMDALocalInfo* info, [[maybe_unused]] PetscReal time, PetscReal* aY, [[maybe_unused]] PetscReal* aG, void* pwp)
  -> PetscErrorCode
{
  auto* petscWrap = static_cast<PetscWrap<PDEOptions>*>(pwp);
  petscWrap->m_user.rhsFuncCalled = PETSC_TRUE;

  assert(!petscWrap->m_methodProps.fields.empty() and "Fields have not been set");
  assert(petscWrap->BaseClass::m_myPDE != nullptr and "PDE has not been set");

  auto const numFields = petscWrap->m_methodProps.fields.size();

  // If the vectors aren't the right size, resize them once
  if (petscWrap->m_inTmp.size() != numFields) {
    petscWrap->m_inTmp.resize(numFields);
  }

  if (petscWrap->m_vTmp.size() != numFields) {
    petscWrap->m_vTmp.resize(numFields);
  }

  auto& inTmp = petscWrap->m_inTmp;
  auto& vTmp = petscWrap->m_vTmp;

  PetscInt const xGridPts = info->xm;
  PetscInt const yGridPts = info->ym;

  if (info->dim == 1) {
    for (PetscInt i = info->xs; i < info->xs + info->xm; ++i) {
      std::copy(&aY[i * numFields], &aY[(i + 1) * numFields], inTmp.begin());

      for (std::size_t field = 0; field < numFields; ++field) {
        petscWrap->BaseClass::m_myPDE->RHS(inTmp, vTmp, field);
      }

      for (std::size_t field = 0; field < numFields; ++field) {
        aG[i * numFields + field] = vTmp[field];
      }
    }
  }
  else if (info->dim == 2) {
    for (PetscInt j = info->ys; j < info->ys + info->ym; ++j) {
      for (PetscInt i = info->xs; i < info->xs + info->xm; ++i) {
        std::copy(&aY[(i + xGridPts * j) * numFields], &aY[(i + 1 + xGridPts * j) * numFields], inTmp.begin());

        for (std::size_t field = 0; field < numFields; ++field) {
          petscWrap->BaseClass::m_myPDE->RHS(inTmp, vTmp, field);
        }

        for (std::size_t field = 0; field < numFields; ++field) {
          aG[((i + xGridPts * j) * numFields) + field] = vTmp[field];
        }
      }
    }
  }
  else if (info->dim == 3) {
    for (PetscInt k = info->zs; k < info->zs + info->zm; ++k) {
      for (PetscInt j = info->ys; j < info->ys + info->ym; ++j) {
        for (PetscInt i = info->xs; i < info->xs + info->xm; ++i) {
          std::copy(&aY[(i + xGridPts * ((j + yGridPts * k))) * numFields],
            &aY[(i + 1 + xGridPts * (j + yGridPts * k)) * numFields],
            inTmp.begin());

          for (std::size_t field = 0; field < numFields; ++field) {
            petscWrap->BaseClass::m_myPDE->RHS(inTmp, vTmp, field);
          }

          for (std::size_t field = 0; field < numFields; ++field) {
            aG[(i + xGridPts * (j + yGridPts * k) * numFields) + field] = vTmp[field];
          }
        }
      }
    }
  }

  return 0;
}

///
/// \brief Computes the Jacobian of the right-hand side for a local DMDA domain partition
///
/// This function assembles the Jacobian matrix of the PDE system's right-hand side for a local
/// sub-domain managed by PETSc's Distributed Array (DMDA). It iterates over grid points in 1D, 2D,
/// or 3D depending on problem dimensionality, evaluates the Jacobian at each point using the
/// user-provided PetscWrap instance, and populates the PETSc matrix using stencil-based indexing.
///
/// \param[in] info A pointer to DMDA local sub-domain information (grid indices and sizes)
/// \param[in] time Current simulation time
/// \param[in] aY An array containing the initial values of the unknown fields in the system
/// \param[in] asmMat The Jacobian matrix to be assembled (for implicit solvers)
/// \param[in] precondMat The material from which PETSc can build a preconditioner; typically the same as P
/// but may differ
/// \param[in] pwp A void pointer to an instance of PetscWrap (the context)
/// \returns PetscErrorCode indicating success or failure
///
template<std::semiregular PDEOptions>
auto PetscWrap<PDEOptions>::FormRHSJacobianLocal(
  DMDALocalInfo* info, [[maybe_unused]] PetscReal time, PetscReal* aY, Mat asmMat, Mat precondMat, void* pwp)
  -> PetscErrorCode
{
  auto* petscWrap = static_cast<PetscWrap<PDEOptions>*>(pwp);
  petscWrap->m_user.rhsJacobianCalled = PETSC_TRUE;

  assert(!petscWrap->m_methodProps.fields.empty() and "Fields have not been set");
  assert(petscWrap->BaseClass::m_myPDE != nullptr and "PDE has not been set");

  auto const numFields = std::size(petscWrap->m_methodProps.fields);

  // If the vectors aren't the right size, resize them once
  if (petscWrap->m_inTmp.size() != numFields) {
    petscWrap->m_inTmp.resize(numFields);
  }

  if (petscWrap->m_colTmp.size() != numFields) {
    petscWrap->m_colTmp.resize(numFields);
  }
  if (petscWrap->m_vTmp.size() != numFields) {
    petscWrap->m_vTmp.resize(numFields);
  }

  auto& rowTmp = petscWrap->m_rowTmp;
  auto& colTmp = petscWrap->m_colTmp;
  auto& inTmp = petscWrap->m_inTmp;
  auto& vTmp = petscWrap->m_vTmp;

  PetscInt const xGridPts = info->xm;
  PetscInt const yGridPts = info->ym;

  if (info->dim == 1) {
    for (PetscInt i = info->xs; i < info->xs + info->xm; ++i) {
      rowTmp.i = i;
      std::copy(&aY[i * numFields], &aY[(i + 1) * numFields], inTmp.begin());

      for (std::size_t field = 0; field < numFields; ++field) {
        colTmp[field].i = i;
        colTmp[field].c = static_cast<PetscInt>(field);

        rowTmp.c = static_cast<PetscInt>(field);
        petscWrap->BaseClass::m_myPDE->JacRHS(inTmp, vTmp, static_cast<PetscInt>(field));
        PetscCall(MatSetValuesStencil(
          precondMat, 1, &rowTmp, static_cast<PetscInt>(numFields), colTmp.data(), vTmp.data(), INSERT_VALUES));
      }
    }
  }
  else if (info->dim == 2) {
    for (PetscInt j = info->ys; j < info->ys + info->ym; ++j) {
      rowTmp.j = j;

      for (PetscInt i = info->xs; i < info->xs + info->xm; ++i) {
        rowTmp.i = i;
        std::copy(&aY[(i + xGridPts * j) * numFields], &aY[(i + 1 + xGridPts * j) * numFields], inTmp.begin());

        for (std::size_t field = 0; field < numFields; ++field) {
          colTmp[field].i = i;
          colTmp[field].j = j;
          colTmp[field].c = static_cast<PetscInt>(field);

          rowTmp.c = static_cast<PetscInt>(field);

          petscWrap->BaseClass::m_myPDE->JacRHS(inTmp, vTmp, static_cast<PetscInt>(field));
          PetscCall(MatSetValuesStencil(
            precondMat, 1, &rowTmp, static_cast<PetscInt>(numFields), colTmp.data(), vTmp.data(), INSERT_VALUES));
        }
      }
    }
  }
  else if (info->dim == 3) {
    for (PetscInt k = info->zs; k < info->zs + info->zm; ++k) {
      rowTmp.k = k;

      for (PetscInt j = info->ys; j < info->ys + info->ym; ++j) {
        rowTmp.j = j;

        for (PetscInt i = info->xs; i < info->xs + info->xm; ++i) {
          rowTmp.i = i;
          std::copy(&aY[(i + xGridPts * (j + yGridPts * k)) * numFields],
            &aY[(i + 1 + xGridPts * (j + yGridPts * k)) * numFields],
            inTmp.begin());

          for (std::size_t field = 0; field < numFields; ++field) {
            colTmp[field].i = i;
            colTmp[field].j = j;
            colTmp[field].k = k;
            colTmp[field].c = static_cast<PetscInt>(field);

            rowTmp.c = static_cast<PetscInt>(field);

            petscWrap->BaseClass::m_myPDE->JacRHS(inTmp, vTmp, static_cast<PetscInt>(field));
            PetscCall(MatSetValuesStencil(
              precondMat, 1, &rowTmp, static_cast<PetscInt>(numFields), colTmp.data(), vTmp.data(), INSERT_VALUES));
          }
        }
      }
    }
  }

  PetscCall(MatAssemblyBegin(precondMat, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(precondMat, MAT_FINAL_ASSEMBLY));

  if (asmMat != precondMat) {
    PetscCall(MatAssemblyBegin(asmMat, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(asmMat, MAT_FINAL_ASSEMBLY));
  }

  return 0;
}

///
/// \brief FormIFunctionLocal - Compute the local implicit function residuals for the shallow water
/// equations
///
/// \param[in] info   DMDA local information structure containing grid indices and dimensions
/// \param[in] time   Current time value
/// \param[in] aY     Array of field variables
/// \param[out] aYdot Array of time derivatives of field variables
/// \param[in,out] aF Array to store computed residuals for the implicit function
/// \param[in] pwp    A pointer to an instance of PetscWrap
///
/// \returns PetscErrorCode indicating success (0) or error code
///
template<std::semiregular PDEOptions>
auto PetscWrap<PDEOptions>::FormIFunctionLocal(DMDALocalInfo* info,
  [[maybe_unused]] PetscReal time,
  PetscReal* aY,
  [[maybe_unused]] PetscReal* aYdot,
  PetscReal* aF,
  void* pwp) -> PetscErrorCode
{
  auto* petscWrap = static_cast<PetscWrap<PDEOptions>*>(pwp);
  petscWrap->m_user.iFuncCalled = PETSC_TRUE;

  assert(!petscWrap->m_methodProps.fields.empty() and "Fields have not been set");
  assert(petscWrap->BaseClass::m_myPDE != nullptr && "PDE has not been set");

  auto const& stencilData = petscWrap->m_myStencilData;
  auto const numFields = std::size(petscWrap->m_methodProps.fields);

  // If the vectors aren't the right size, resize them once
  if (petscWrap->m_inTmp.size() != numFields) {
    petscWrap->m_inTmp.resize(numFields);
  }

  if (petscWrap->m_outTmp.size() != numFields) {
    petscWrap->m_outTmp.resize(numFields);
  }

  if (petscWrap->m_outTmp.size() != numFields) {
    petscWrap->m_outTmp.resize(numFields);
  }

  if (petscWrap->m_vTmp.size() != numFields) {
    petscWrap->m_vTmp.resize(numFields);
  }

  if (petscWrap->m_dTime.size() != numFields) {
    petscWrap->m_dTime.resize(numFields);
  }

  if (petscWrap->m_indTmp.size() != numFields) {
    petscWrap->m_indTmp.resize(numFields);
  }

  if (petscWrap->m_myStencilGriddd.has_value()) {
    if (!petscWrap->m_ind2Tmp.has_value()) {
      petscWrap->m_ind2Tmp = std::vector<PetscReal>(numFields);
    }
    else if (petscWrap->m_ind2Tmp->size() != numFields) {
      petscWrap->m_ind2Tmp->resize(numFields);
    }
  }

  auto& inTmp = petscWrap->m_inTmp;
  auto& dTime = petscWrap->m_dTime;
  auto& vTmp = petscWrap->m_vTmp;

  auto& outTmp = petscWrap->m_outTmp;

  auto const& my1dStenGrid = petscWrap->m_myStencilGridd;
  auto const& my2dStenGrid = petscWrap->m_myStencilGriddd;

  auto& in1dTmp = petscWrap->m_indTmp;
  auto& in2dTmp = petscWrap->m_ind2Tmp;

  auto const stencilWidth = std::invoke(
    [&stencilData, &info]
    {
      if (info->dim == 1) {
        return stencilData[2];
      }

      if (info->dim == 2) {
        return stencilData[3];
      }

      return stencilData[4];
    });

  // Fill ghost cells at domain boundaries using nearest interior values
  auto const netWidth = static_cast<PetscInt>(std::floor(stencilWidth / 2));

  PetscInt const xGridPts = info->xm;
  PetscInt const yGridPts = info->ym;

  if (info->dim == 1) {
    if (info->bx == DM_BOUNDARY_GHOSTED) {
      FillGhostCells1d(petscWrap->m_da, petscWrap->m_info);
    }

    for (PetscInt i = info->xs; i < info->xs + info->xm; ++i) {
      std::copy(&aY[i * numFields], &aY[(i + 1) * numFields], inTmp.begin());

      if (my1dStenGrid.has_value()) {
        std::fill(in1dTmp.begin(), in1dTmp.end(), 0);

        for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
          for (std::size_t field = 0; field < numFields; ++field) {
            // First derivative in x direction
            in1dTmp[field] += my1dStenGrid->at(pos) * aY[(i - netWidth + pos) * numFields + field];
          }
        }
      }

      if (my2dStenGrid.has_value() and in2dTmp.has_value()) {
        std::fill(in2dTmp->begin(), in2dTmp->end(), 0);

        for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
          for (std::size_t field = 0; field < numFields; ++field) {
            // Second derivative in x direction
            in2dTmp->at(field) += my2dStenGrid->at(pos) * aY[(i - netWidth + pos) * numFields + field];
          }
        }
      }

      petscWrap->BaseClass::m_myPDE->LHSOpAsplit(inTmp, dTime, in1dTmp, in2dTmp, std::nullopt, vTmp);

      for (std::size_t field = 0; field < numFields; ++field) {
        aF[i * numFields + field] = vTmp[field];
      }
    }
  }
  else if (info->dim == 2) {
    if (info->bx == DM_BOUNDARY_GHOSTED and info->by == DM_BOUNDARY_GHOSTED) {
      FillGhostCells2d(petscWrap->m_da, petscWrap->m_info);
    }

    for (PetscInt j = info->ys; j < info->ys + info->ym; ++j) {
      for (PetscInt i = info->xs; i < info->xs + info->xm; ++i) {
        std::copy(&aY[(i + xGridPts * j) * numFields], &aY[(i + 1 + xGridPts * j) * numFields], inTmp.begin());
        std::copy(&aF[(i + xGridPts * j) * numFields], &aF[(i + 1 + xGridPts * j) * numFields], vTmp.begin());

        if (my1dStenGrid.has_value()) {
          std::fill(in1dTmp.begin(), in1dTmp.end(), 0);

          for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
            for (std::size_t field = 0; field < numFields; ++field) {
              // First derivative in x direction
              in1dTmp[field] += my1dStenGrid->at(pos) * aY[((i - netWidth + pos) + xGridPts * j) * numFields + field];

              // First derivative in y direction
              in1dTmp[field] += my1dStenGrid->at(pos) * aY[(i + xGridPts * (j - netWidth + pos)) * numFields + field];
            }
          }
        }

        if (my2dStenGrid.has_value() and in2dTmp.has_value()) {
          std::fill(in2dTmp->begin(), in2dTmp->end(), 0);

          for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
            for (std::size_t field = 0; field < numFields; ++field) {
              // Second derivative in x direction
              in2dTmp->at(field) +=
                my2dStenGrid->at(pos) * aY[((i - netWidth + pos) + xGridPts * j) * numFields + field];

              // Second derivative in y direction
              in2dTmp->at(field) +=
                my2dStenGrid->at(pos) * aY[(i + xGridPts * (j - netWidth + pos)) * numFields + field];
            }
          }
        }

        petscWrap->BaseClass::m_myPDE->LHSOpAsplit(inTmp, dTime, in1dTmp, in2dTmp, std::nullopt, vTmp);

        for (std::size_t field = 0; field < numFields; ++field) {
          aF[(i + xGridPts * j) * numFields + field] = aF[field];
        }
      }
    }
  }
  else if (info->dim == 3) {
    if (info->bx == DM_BOUNDARY_GHOSTED and info->by == DM_BOUNDARY_GHOSTED and info->bz == DM_BOUNDARY_GHOSTED) {
      FillGhostCells3d(petscWrap->m_da, petscWrap->m_info);
    }

    for (PetscInt k = info->zs; k < info->zs + info->zm; ++k) {
      for (PetscInt j = info->ys; j < info->ys + info->ym; ++j) {
        for (PetscInt i = info->xs; i < info->xs + info->xm; ++i) {
          std::copy(&aY[(i + xGridPts * (j + yGridPts * k)) * numFields],
            &aY[(i + 1 + xGridPts * (j + yGridPts * k)) * numFields],
            inTmp.begin());

          std::copy(&aF[(i + xGridPts * (j + yGridPts * k)) * numFields],
            &aF[(i + 1 + xGridPts * (j + yGridPts * k)) * numFields],
            vTmp.begin());

          if (my1dStenGrid.has_value()) {
            std::fill(in1dTmp.begin(), in1dTmp.end(), 0);

            for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
              for (std::size_t field = 0; field < numFields; ++field) {
                // First derivative in x direction
                in1dTmp[field] += my1dStenGrid->at(pos)
                  * aY[((i - netWidth + pos) + xGridPts * (j + yGridPts * k)) * numFields + field];

                // First derivative in y direction
                in1dTmp[field] += my1dStenGrid->at(pos)
                  * aY[(i + xGridPts * ((j - netWidth + pos) + yGridPts * k)) * numFields + field];

                // First derivative in z direction
                in1dTmp[field] += my1dStenGrid->at(pos)
                  * aY[(i + xGridPts * (j + yGridPts * (k - netWidth + pos))) * numFields + field];
              }
            }
          }

          if (my2dStenGrid.has_value() and in2dTmp.has_value()) {
            std::fill(in2dTmp->begin(), in2dTmp->end(), 0);

            for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
              for (std::size_t field = 0; field < numFields; ++field) {
                // Second derivative in x direction
                in2dTmp->at(field) += my2dStenGrid->at(pos)
                  * aY[((i - netWidth + pos) + xGridPts * (j + yGridPts * k)) * numFields + field];

                // Second derivative in y direction
                in2dTmp->at(field) += my2dStenGrid->at(pos)
                  * aY[(i + xGridPts * ((j - netWidth + pos) + yGridPts * k)) * numFields + field];

                // Second derivative in z direction
                in2dTmp->at(field) += my2dStenGrid->at(pos)
                  * aY[(i + xGridPts * (j + yGridPts * (k - netWidth + pos))) * numFields + field];
              }
            }
          }

          petscWrap->BaseClass::m_myPDE->LHSOpAsplit(inTmp, dTime, in1dTmp, in2dTmp, std::nullopt, vTmp);

          for (std::size_t field = 0; field < numFields; ++field) {
            aF[i + xGridPts * (j + yGridPts * k) * numFields + field] = vTmp[field];
          }
        }
      }
    }
  }

  return 0;
}

///
/// \brief Computes the Jacobian matrix for the implicit function F(t, Y, Ẏ) in the Shallow Water
/// equations
///
/// \param[in] info  A pointer to DMDA local sub-domain information (grid indices and sizes)
/// \param[in] time  Current simulation time
/// \param[in] aY    An array containing the initial values of the unknown fields in the system at
/// all grid points in the local sub-domain
/// \param[in] aYdot Array of time derivatives at all grid points
/// \param[in] shift Shift parameter (typically dt/theta) for implicit time integration; multiplies
/// the Ẏ contribution to the domain
/// \param[in,out] asmMat The Jacobian matrix to be assembled (for implicit solvers)
/// \param[in,out] precondMat The material from which PETSc can build a preconditioner matrix; typically the
/// same as P but may differ; filled with stencil values representing spatial discretization
/// \param[in] pwp   A void pointer to an instance of PetscWrap (the context)
/// \returns PetscErrorCode indicating success or failure
///
template<std::semiregular PDEOptions>
auto PetscWrap<PDEOptions>::FormIJacobianLocal(DMDALocalInfo* info,
  [[maybe_unused]] PetscReal time,
  PetscReal* aY,
  [[maybe_unused]] PetscReal* aYdot,
  PetscReal shift,
  Mat asmMat,
  Mat precondMat,
  void* pwp) -> PetscErrorCode
{
  PetscCall(MatZeroEntries(precondMat));

  auto* petscWrap = static_cast<PetscWrap<PDEOptions>*>(pwp);
  petscWrap->m_user.iJacobianCalled = PETSC_TRUE;

  assert(!petscWrap->m_methodProps.fields.empty() and "Fields have not been set");
  assert(petscWrap->BaseClass::m_myPDE != nullptr && "PDE has not been set");

  auto const& stenData = petscWrap->m_myStencilData;
  auto const numFields = std::size(petscWrap->m_methodProps.fields);

  // If the vectors aren't the right size, resize them once
  if (petscWrap->m_inTmp.size() != numFields) {
    petscWrap->m_inTmp.resize(numFields);
  }

  if (petscWrap->m_dTime.size() != numFields) {
    petscWrap->m_dTime.resize(numFields);
  }
  if (petscWrap->m_indTmp.size() != numFields) {
    petscWrap->m_indTmp.resize(numFields);
  }

  auto const stencilWidth = std::invoke(
    [&stenData, &info]
    {
      if (info->dim == 1) {
        return stenData[2];
      }

      if (info->dim == 2) {
        return stenData[3];
      }

      return stenData[4];
    });

  if (std::cmp_not_equal(petscWrap->m_vTmp.size(), stencilWidth)) {
    petscWrap->m_vTmp.resize((std::size_t)stencilWidth);
  }
  if (std::cmp_not_equal(petscWrap->m_vTmpTmp.size(), stencilWidth)) {
    petscWrap->m_vTmpTmp.resize((std::size_t)stencilWidth);
  }
  if (std::cmp_not_equal(petscWrap->m_colTmp.size(), stencilWidth)) {
    petscWrap->m_colTmp.resize((std::size_t)stencilWidth);
  }

  if (petscWrap->m_myStencilGriddd.has_value()) {
    if (!petscWrap->m_ind2Tmp.has_value()) {
      petscWrap->m_ind2Tmp = std::vector<PetscReal>(numFields);
    }
    else if (petscWrap->m_ind2Tmp->size() != numFields) {
      petscWrap->m_ind2Tmp->resize(numFields);
    }
  }

  auto& rowTmp = petscWrap->m_rowTmp;
  auto& inTmp = petscWrap->m_inTmp;
  auto& dTime = petscWrap->m_dTime;

  auto const& my1dStenGrid = petscWrap->m_myStencilGridd;
  auto const& my2dStenGrid = petscWrap->m_myStencilGriddd;

  auto& in1dTmp = petscWrap->m_indTmp;
  auto& in2dTmp = petscWrap->m_ind2Tmp;

  auto& vTmp = petscWrap->m_vTmp;
  auto& vTmpTmp = petscWrap->m_vTmpTmp;
  auto& colTmp = petscWrap->m_colTmp;

  PetscInt const xGridPts = info->xm;
  PetscInt const yGridPts = info->ym;

  auto const netWidth = static_cast<PetscInt>(std::floor(stencilWidth / 2));

  if (info->dim == 1) {
    if (info->bx == DM_BOUNDARY_GHOSTED) {
      FillGhostCells1d(petscWrap->m_da, petscWrap->m_info);
    }

    // Iterate over interior points and assemble Jacobian entries

    for (PetscInt i = info->xs; i < info->xs + info->xm; ++i) {
      std::copy(&aY[i * numFields], &aY[(i + 1) * numFields], inTmp.begin());

      // Compute first derivative stencil contributions if available
      if (my1dStenGrid.has_value()) {
        std::fill(in1dTmp.begin(), in1dTmp.end(), 0);

        for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
          for (std::size_t field = 0; field < numFields; ++field) {
            // First derivative in x-direction
            in1dTmp[field] += my1dStenGrid->at(pos) * aY[(i - netWidth + pos) * numFields];
          }
        }
      }

      // Set row index for current point
      rowTmp.i = i;

      // Loop over component row field components
      for (std::size_t compRow = 0; compRow < numFields; ++compRow) {
        rowTmp.c = compRow;

        // Loop over component column field components
        for (std::size_t compCol = 0; compCol < numFields; ++compCol) {
          std::fill(vTmp.begin(), vTmp.end(), 0);

          petscWrap->BaseClass::m_myPDE->JacLHSOpAsplit(
            inTmp, dTime, in1dTmp, std::nullopt, std::nullopt, vTmpTmp, compRow, compCol, 0);

          for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
            vTmp[pos] += vTmpTmp[pos];
          }

          // Add time derivative shift to diagonal
          if (compRow == compCol) {
            vTmp[netWidth] += shift;
          }

          // Compute diagonal Jacobian contribution
          if (my1dStenGrid.has_value()) {
            petscWrap->BaseClass::m_myPDE->JacLHSOpAsplit(
              inTmp, dTime, in1dTmp, std::nullopt, std::nullopt, vTmpTmp, compRow, compCol, 1);

            // Compute off-diagonal Jacobian contributions from spatial stencil
            for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
              vTmp[pos] += my1dStenGrid->at(pos) * vTmpTmp[pos];
            }
          }

          // Compute second derivative contributions if available
          if (my2dStenGrid.has_value() and in2dTmp.has_value()) {
            petscWrap->BaseClass::m_myPDE->JacLHSOpAsplit(
              inTmp, dTime, in1dTmp, in2dTmp, std::nullopt, vTmp, compRow, compCol, 0);

            for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
              vTmp[pos] += my2dStenGrid->at(pos) * vTmpTmp[pos];
            }
          }

          // Compute off-diagonal Jacobian contributions from spatial stencil
          for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
            colTmp[pos].i = i - netWidth + pos;
            colTmp[pos].c = compCol;
          }

          PetscCall(
            MatSetValuesStencil(precondMat, 1, &rowTmp, stencilWidth, colTmp.data(), vTmp.data(), INSERT_VALUES));
        }
      }
    }
  }
  else if (info->dim == 2) {
    if (info->bx == DM_BOUNDARY_GHOSTED and info->by == DM_BOUNDARY_GHOSTED) {
      FillGhostCells2d(petscWrap->m_da, petscWrap->m_info);
    }

    // Assemble Jacobian for 2D domain
    for (PetscInt j = info->ys; j < info->ys + info->ym; ++j) {
      rowTmp.j = j;

      for (PetscInt i = info->xs; i < info->xs + info->xm; ++i) {
        // Extract solution vector for current grid point
        std::copy(&aY[(i + xGridPts * j) * numFields], &aY[(i + 1 + xGridPts * j) * numFields], inTmp.begin());

        // Compute stencil contributions from spatial derivatives
        if (my1dStenGrid.has_value()) {
          std::fill(in1dTmp.begin(), in1dTmp.end(), 0);

          for (std::size_t pos = 0; pos < netWidth; ++pos) {
            for (std::size_t field = 0; field < numFields; ++field) {
              // First derivative in x-direction
              in1dTmp[field] += my1dStenGrid->at(pos) * aY[((i - netWidth + pos) + xGridPts * j) * numFields + field];

              // First derivative in y-direction
              in1dTmp[field] += my1dStenGrid->at(pos) * aY[(i + xGridPts * (j - netWidth + pos)) * numFields + field];
            }
          }
        }

        // Set row index for current 2D point
        rowTmp.i = i;

        // Assemble Jacobian entries for all field combinations
        for (std::size_t compRow = 0; compRow < numFields; ++compRow) {
          rowTmp.c = compRow;

          for (std::size_t compCol = 0; compCol < numFields; ++compCol) {
            std::fill(vTmp.begin(), vTmp.end(), 0);

            petscWrap->BaseClass::m_myPDE->JacLHSOpAsplit(
              inTmp, dTime, in1dTmp, std::nullopt, std::nullopt, vTmpTmp, compRow, compCol, 0);

            for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
              vTmp[pos] += vTmpTmp[pos];
            }

            if (compRow == compCol) {
              vTmp[netWidth] += shift;
            }

            // Compute diagonal term (time derivative + LHS operator)
            if (my1dStenGrid.has_value()) {
              petscWrap->BaseClass::m_myPDE->JacLHSOpAsplit(
                inTmp, dTime, in1dTmp, std::nullopt, std::nullopt, vTmpTmp, compRow, compCol, 1);

              // Compute stencil contributions from first derivative
              for (std::size_t pos = 0; pos < netWidth; ++pos) {
                vTmp[pos] += my1dStenGrid->at(pos) * vTmpTmp[pos];
              }
            }

            // Compute second derivative contributions (if available)
            if (my2dStenGrid.has_value()) {
              petscWrap->BaseClass::m_myPDE->JacLHSOpAsplit(
                inTmp, dTime, in1dTmp, in2dTmp, std::nullopt, vTmp, compRow, compCol, 0);

              // Compute off-diagonal Jacobian contributions from spatial stencil
              for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
                petscWrap->BaseClass::m_myPDE->JacLHSOpAsplit(
                  inTmp, dTime, in1dTmp, in2dTmp, std::nullopt, vTmpTmp, compRow, compCol, 1);

                vTmp[pos] += my2dStenGrid->at(pos) * vTmpTmp[pos];
              }
            }

            for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
              colTmp[pos].i = i - netWidth + pos;
              colTmp[pos].j = j;
              colTmp[pos].c = compCol;
            }

            // Insert row of Jacobian into matrix
            PetscCall(MatSetValuesStencil(precondMat, 1, &rowTmp, 3, colTmp.data(), vTmp.data(), INSERT_VALUES));
          }
        }
      }
    }
  }
  else if (info->dim == 3) {
    if (info->bx == DM_BOUNDARY_GHOSTED and info->by == DM_BOUNDARY_GHOSTED and info->bz == DM_BOUNDARY_GHOSTED) {
      FillGhostCells3d(petscWrap->m_da, petscWrap->m_info);
    }

    // Assemble Jacobian for 3D domain
    for (PetscInt k = info->zs; k < info->zs + info->zm; ++k) {
      rowTmp.k = k;

      for (PetscInt j = info->ys; j < info->ys + info->ym; ++j) {
        rowTmp.j = j;

        for (PetscInt i = info->xs; i < info->xs + info->xm; ++i) {
          // Extract solution vector for current grid point
          std::copy(&aY[(i + xGridPts * (j + yGridPts * k)) * numFields],
            &aY[(i + 1 + xGridPts * (j + yGridPts * k)) * numFields],
            inTmp.begin());

          // Compute stencil contributions from spatial derivatives
          if (my1dStenGrid.has_value()) {
            std::fill(in1dTmp.begin(), in1dTmp.end(), 0);

            for (std::size_t pos = 0; pos < netWidth; ++pos) {
              for (std::size_t field = 0; field < numFields; ++field) {
                // First derivative in x-direction
                in1dTmp[field] += my1dStenGrid->at(pos)
                  * aY[((i - netWidth + pos) + xGridPts * (j + yGridPts * k)) * numFields + field];

                // First derivative in y-direction
                in1dTmp[field] += my1dStenGrid->at(pos)
                  * aY[((i + xGridPts) * ((j - netWidth + pos) + yGridPts * k)) * numFields + field];

                // First derivative in z-direction
                in1dTmp[field] +=
                  my1dStenGrid->at(pos) * aY[(i + xGridPts * j + yGridPts * (k - netWidth + pos)) * numFields + field];
              }
            }
          }

          // Set row index for current 3D point
          rowTmp.i = i;

          // Assemble Jacobian entries for all field contributions
          for (std::size_t compRow = 0; compRow < numFields; ++compRow) {
            rowTmp.c = compRow;

            for (std::size_t compCol = 0; compCol < numFields; ++compCol) {
              std::fill(vTmp.begin(), vTmp.end(), 0.0);

              petscWrap->BaseClass::m_myPDE->JacLHSOpAsplit(
                inTmp, dTime, in1dTmp, std::nullopt, std::nullopt, vTmpTmp, compRow, compCol, 0);

              for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
                vTmp[pos] += vTmpTmp[pos];
              }

              if (compRow == compCol) {
                vTmp[netWidth] += shift;
              }

              // Compute diagonal term (time derivative + LHS operator)
              if (my1dStenGrid.has_value()) {
                petscWrap->BaseClass::m_myPDE->JacLHSOpAsplit(
                  inTmp, dTime, in1dTmp, std::nullopt, std::nullopt, vTmpTmp, compRow, compCol, 1);

                // Compute stencil contributions from first derivative
                for (std::size_t pos = 0; pos < netWidth; ++pos) {
                  vTmp[pos] += my1dStenGrid->at(pos) * vTmpTmp[pos];
                }
              }

              // Compute second derivative contributions (if available)
              if (my2dStenGrid.has_value()) {
                petscWrap->BaseClass::m_myPDE->JacLHSOpAsplit(
                  inTmp, dTime, in1dTmp, in2dTmp, std::nullopt, vTmpTmp, compRow, compCol, 0);

                // Compute off-diagonal Jacobian contributions from spatial stencil
                for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
                  petscWrap->BaseClass::m_myPDE->JacLHSOpAsplit(
                    inTmp, dTime, in1dTmp, in2dTmp, std::nullopt, vTmpTmp, compRow, compCol, 1);
                  vTmp[pos] += my2dStenGrid->at(pos) * vTmpTmp[pos];
                }
              }

              for (std::size_t pos = 0; std::cmp_less(pos, stencilWidth); ++pos) {
                colTmp[pos].i = i - netWidth + pos;
                colTmp[pos].j = j;
                colTmp[pos].k = k;
                colTmp[pos].c = compCol;
              }

              PetscCall(
                MatSetValuesStencil(precondMat, 1, &rowTmp, stencilWidth, colTmp.data(), vTmp.data(), INSERT_VALUES));
            }
          }
        }
      }
    }
  }

  PetscCall(MatAssemblyBegin(precondMat, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(precondMat, MAT_FINAL_ASSEMBLY));

  if (asmMat != precondMat) {
    PetscCall(MatAssemblyBegin(asmMat, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(asmMat, MAT_FINAL_ASSEMBLY));
  }

  return 0;
}

}  // namespace fcfd::pdenumerics
