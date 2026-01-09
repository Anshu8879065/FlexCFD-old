add_library(flexcfd_pdemodel INTERFACE)

target_sources(flexcfd_pdemodel
    INTERFACE
        BoundaryCondition.hpp
        GridParams.hpp
        ModelParams.hpp
        PDEParams.hpp
        PDESystem.hpp
        PDEType.hpp
        SVE/SVE.hpp
        SVE/SVECas.hpp
        SWE/SWE2d.hpp
        SWE/SWE2dCas.hpp
)

target_include_directories(flexcfd_pdemodel
    ${warning_guard}
    INTERFACE
        "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>"
)


target_compile_features(flexcfd_pdemodel
    INTERFACE
        cxx_std_20
)
