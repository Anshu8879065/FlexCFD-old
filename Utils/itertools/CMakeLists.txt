add_library(flexcfd_itertools INTERFACE)

target_sources(flexcfd_itertools
    INTERFACE
        IterTools.hpp
        TriangleRange.hpp
        SquareRange.hpp
)

target_include_directories(flexcfd_itertools
    ${warning_guard}
    INTERFACE
        "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>"
)

target_compile_features(flexcfd_itertools
    INTERFACE
        cxx_std_20
)
