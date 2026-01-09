cmake_minimum_required(VERSION 3.14)

include(cmake/prelude.cmake)
include(FindPkgConfig REQUIRED)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)

pkg_get_variable(CMAKE_C_COMPILER PETSc ccompiler)
pkg_get_variable(CMAKE_CXX_COMPILER PETSc cxxcompiler)

if (CXX_COMPILER)
  set(CMAKE_CXX_COMPILER ${CXX_COMPILER})
  enable_language(CXX)
endif (CXX_COMPILER)
pkg_get_variable(FORTRAN_COMPILER PETSc fcompiler)
if (FORTRAN_COMPILER)
  set(CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER})
  enable_language(Fortran)
endif (FORTRAN_COMPILER)

project(
    flexcfd
    VERSION 0.1.0
    DESCRIPTION "Short description"
    HOMEPAGE_URL "https://example.com/"
)

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Declare library --------

add_library(flexcfd_lib OBJECT)

set_target_properties(flexcfd_lib PROPERTIES LINKER_LANGUAGE CXX)

add_subdirectory(Utils/itertools)
add_subdirectory(PDEs)
add_subdirectory(Numerics)

target_link_libraries(flexcfd_lib
    PUBLIC
        flexcfd_itertools
        flexcfd_pdemodel
        flexcfd_pdenumerics
        $<TARGET_OBJECTS:flexcfd_pdenumerics>
)

target_compile_features(flexcfd_lib PUBLIC cxx_std_20)


# ---- Declare PETSc executable ----

add_executable(flexcfd_exe ${PROJECT_SOURCE_DIR}/Example/petscmain.cpp)
add_executable(flexcfd::exe ALIAS flexcfd_exe)

set_property(TARGET flexcfd_exe PROPERTY OUTPUT_NAME petsc)

target_compile_features(flexcfd_exe PRIVATE cxx_std_20)

target_link_libraries(flexcfd_exe PRIVATE flexcfd_lib)

# ---- Declare Casulli executable ----

add_executable(casulli_exe ${PROJECT_SOURCE_DIR}/Example/Casulliwrap.cpp)
add_executable(casulli::exe ALIAS casulli_exe)

set_property(TARGET casulli_exe PROPERTY OUTPUT_NAME casulli)

target_compile_features(casulli_exe PRIVATE cxx_std_20)

target_link_libraries(casulli_exe PRIVATE flexcfd_lib)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
  include(cmake/install-rules.cmake)
endif()

# ---- Developer mode ----

if(NOT flexcfd_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of flexcfd"
  )
endif()

include(cmake/dev-mode.cmake)
