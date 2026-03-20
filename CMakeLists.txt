cmake_minimum_required(VERSION 3.15)

if(POLICY CMP0091)
    cmake_policy(SET CMP0091 NEW)
endif()

project(DACEx)

include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    # Keep runtime ABI consistent with typical prebuilt DACE binaries on Windows.
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")
endif()

option(DACEX_FETCH_DACE "Automatically fetch/build DACE if not found locally" OFF)
option(DACEX_ALLOW_WIN32 "Allow Win32 builds (requires matching 32-bit DACE binaries)" OFF)

if(WIN32 AND CMAKE_SIZEOF_VOID_P EQUAL 4 AND NOT DACEX_ALLOW_WIN32)
    message(FATAL_ERROR
        "Win32 (32-bit) build detected. Typical DACE binaries are 64-bit.\n"
        "Please configure with x64 (e.g., -A x64), or set DACEX_ALLOW_WIN32=ON\n"
        "only if you have matching 32-bit DACE binaries."
    )
endif()

# Public-repo friendly DACE discovery: no machine-specific assumptions.
set(DACEX_DACE_ROOT "" CACHE PATH "Path to DACE install/build root (contains include/ and lib/)")
set(DACEX_DACE_INCLUDE_DIR "" CACHE PATH "Path to DACE include directory (contains dace/dace.h)")
set(DACEX_DACE_LIBRARY "" CACHE FILEPATH "Path to DACE library file (e.g., dace_s.lib, libdace_s.a, dace.lib)")

set(DACEX_DACE_HINTS)
if(DACEX_DACE_ROOT)
    list(APPEND DACEX_DACE_HINTS "${DACEX_DACE_ROOT}")
endif()
if(DEFINED ENV{DACE_ROOT})
    list(APPEND DACEX_DACE_HINTS "$ENV{DACE_ROOT}")
endif()

# Auto-discover common local DACE locations from the wiki workflow.
set(DACEX_USER_HOME "")
if(DEFINED ENV{USERPROFILE} AND NOT "$ENV{USERPROFILE}" STREQUAL "")
    set(DACEX_USER_HOME "$ENV{USERPROFILE}")
elseif(DEFINED ENV{HOMEDRIVE} AND DEFINED ENV{HOMEPATH} AND NOT "$ENV{HOMEDRIVE}$ENV{HOMEPATH}" STREQUAL "")
    set(DACEX_USER_HOME "$ENV{HOMEDRIVE}$ENV{HOMEPATH}")
elseif(DEFINED ENV{USERNAME} AND NOT "$ENV{USERNAME}" STREQUAL "")
    set(DACEX_USER_HOME "C:/Users/$ENV{USERNAME}")
endif()

if(DACEX_USER_HOME STREQUAL "")
    get_filename_component(DACEX_USER_HOME "~" ABSOLUTE)
endif()

set(DACEX_DEFAULT_DACE_PATHS)
if(NOT DACEX_USER_HOME STREQUAL "")
    list(APPEND DACEX_DEFAULT_DACE_PATHS
        "${DACEX_USER_HOME}/dace-install"
        "${DACEX_USER_HOME}/dace-build"
    )
endif()

foreach(_dace_path IN LISTS DACEX_DEFAULT_DACE_PATHS)
    if(EXISTS "${_dace_path}")
        list(APPEND DACEX_DACE_HINTS "${_dace_path}")
    endif()
endforeach()

# Preferred path: use DACE's exported CMake package.
# Keep explicit include/library overrides for fallback mode only.
set(DACEX_FORCE_FALLBACK OFF)
if(DACEX_DACE_INCLUDE_DIR OR DACEX_DACE_LIBRARY)
    set(DACEX_FORCE_FALLBACK ON)
endif()

if(NOT DACEX_FORCE_FALLBACK)
    set(DACEX_DACE_CONFIG_HINTS ${DACEX_DACE_HINTS})
    foreach(_dace_path IN LISTS DACEX_DACE_HINTS)
        list(APPEND DACEX_DACE_CONFIG_HINTS
            "${_dace_path}/lib/cmake/dace"
            "${_dace_path}/cmake"
        )
    endforeach()
    find_package(dace CONFIG QUIET HINTS ${DACEX_DACE_CONFIG_HINTS})
endif()

# Optional first-time bootstrap: fetch DACE from GitHub if not found.
if(NOT dace_FOUND AND DACEX_FETCH_DACE AND NOT DACEX_FORCE_FALLBACK)
    include(FetchContent)
    FetchContent_Declare(
        dace
        GIT_REPOSITORY https://github.com/dacelib/dace.git
        GIT_TAG master
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(dace)
endif()

# Fallback path: explicit include/library discovery.
if(NOT dace_FOUND)
    if(NOT DACEX_DACE_INCLUDE_DIR)
        foreach(_dace_root IN LISTS DACEX_DACE_HINTS)
            if(EXISTS "${_dace_root}/include/dace/dace.h")
                set(DACEX_DACE_INCLUDE_DIR "${_dace_root}/include")
                break()
            endif()
        endforeach()
    endif()

    if(NOT DACEX_DACE_LIBRARY)
        foreach(_dace_root IN LISTS DACEX_DACE_HINTS)
            foreach(_dace_lib_candidate
                "${_dace_root}/lib/dace_s.lib"
                "${_dace_root}/lib/dace.lib"
                "${_dace_root}/lib/libdace_s.a"
                "${_dace_root}/lib/libdace.a"
                "${_dace_root}/lib/libdace_s.so"
                "${_dace_root}/lib/libdace.so"
                "${_dace_root}/lib/libdace_s.dylib"
                "${_dace_root}/lib/libdace.dylib"
                "${_dace_root}/lib64/dace_s.lib"
                "${_dace_root}/lib64/dace.lib"
                "${_dace_root}/lib64/libdace_s.a"
                "${_dace_root}/lib64/libdace.a"
                "${_dace_root}/lib64/libdace_s.so"
                "${_dace_root}/lib64/libdace.so"
                "${_dace_root}/lib64/libdace_s.dylib"
                "${_dace_root}/lib64/libdace.dylib"
            )
                if(EXISTS "${_dace_lib_candidate}")
                    set(DACEX_DACE_LIBRARY "${_dace_lib_candidate}")
                    break()
                endif()
            endforeach()
            if(DACEX_DACE_LIBRARY)
                break()
            endif()
        endforeach()
    endif()

    if(NOT DACEX_DACE_INCLUDE_DIR)
        find_path(DACEX_DACE_INCLUDE_DIR
            NAMES dace/dace.h
            HINTS ${DACEX_DACE_HINTS}
            PATH_SUFFIXES include
        )
    endif()

    if(NOT DACEX_DACE_LIBRARY)
        find_library(DACEX_DACE_LIBRARY
            NAMES dace_s dace
            HINTS ${DACEX_DACE_HINTS}
            PATH_SUFFIXES lib lib64
        )
    endif()

    if(NOT DACEX_DACE_INCLUDE_DIR OR NOT DACEX_DACE_LIBRARY)
        message(FATAL_ERROR
            "Could not locate DACE.\n"
            "Use one of:\n"
            "  - Install DACE and set -Ddace_DIR=<path containing daceConfig.cmake>\n"
            "  - Set -DCMAKE_PREFIX_PATH=<DACE install prefix>\n"
            "  - Or provide fallback vars:\n"
            "      -DDACEX_DACE_ROOT=<dace root>\n"
            "      -DDACEX_DACE_INCLUDE_DIR=<dir containing dace/dace.h>\n"
            "      -DDACEX_DACE_LIBRARY=<path to dace library file>\n"
            "Typical workflow from the DACE wiki:\n"
            "  git clone https://github.com/dacelib/dace.git dace\n"
            "  cmake -S dace -B dace-build\n"
            "  cmake --build dace-build\n"
            "  cmake --install dace-build  (optional)\n"
        )
    endif()

    add_library(dace::dace UNKNOWN IMPORTED)
    set_target_properties(dace::dace PROPERTIES
        IMPORTED_LOCATION "${DACEX_DACE_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${DACEX_DACE_INCLUDE_DIR}"
    )
endif()

set(DACEX_DACE_TARGET "")
if(TARGET dace::dace_s)
    set(DACEX_DACE_TARGET dace::dace_s)
elseif(TARGET dace::dace)
    set(DACEX_DACE_TARGET dace::dace)
elseif(TARGET dace)
    set(DACEX_DACE_TARGET dace)
elseif(TARGET dace_s)
    set(DACEX_DACE_TARGET dace_s)
endif()

if(DACEX_DACE_TARGET STREQUAL "")
    message(FATAL_ERROR
        "DACE targets were not found after discovery/bootstrap.\n"
        "Expected one of: dace::dace, dace::dace_s, dace, dace_s"
    )
endif()

add_library(DACEx STATIC
    src/Chebyshev.cpp
    src/Taylor.cpp
)
target_include_directories(DACEx PUBLIC include)
target_link_libraries(DACEx PUBLIC ${DACEX_DACE_TARGET})

if(WIN32 AND (DACEX_DACE_TARGET STREQUAL "dace::dace_s" OR DACEX_DACE_TARGET STREQUAL "dace_s"))
    target_compile_definitions(DACEx PUBLIC DACE_API=)
    if(MSVC)
        target_compile_definitions(DACEx PUBLIC $<$<CONFIG:Debug>:_ITERATOR_DEBUG_LEVEL=0>)
        set_property(TARGET DACEx PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")
        message(STATUS "Using static DACE on MSVC. If your DACE libs were built in Release, build DACEx with --config Release.")
    endif()
endif()

add_executable(DACEx_Chebyshev_test tests/Chebyshev_test.cpp)
add_executable(DACEx_Gaussian_example examples/Gaussian_example.cpp)
add_executable(DACEx_Gaussian_2D_example examples/Gaussian_2D_example.cpp)

foreach(DACEX_EXE_TARGET
    DACEx_Chebyshev_test
    DACEx_Gaussian_example
    DACEx_Gaussian_2D_example)
    target_link_libraries(${DACEX_EXE_TARGET} PRIVATE DACEx)
    if(WIN32 AND MSVC AND (DACEX_DACE_TARGET STREQUAL "dace::dace_s" OR DACEX_DACE_TARGET STREQUAL "dace_s"))
        target_compile_definitions(${DACEX_EXE_TARGET} PRIVATE $<$<CONFIG:Debug>:_ITERATOR_DEBUG_LEVEL=0>)
        set_property(TARGET ${DACEX_EXE_TARGET} PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")
    endif()
    set_target_properties(${DACEX_EXE_TARGET} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endforeach()

if(BUILD_TESTING)
    add_test(NAME cheb_basic_arithmetic COMMAND DACEx_Chebyshev_test basic-arithmetic)
    add_test(NAME cheb_calculus COMMAND DACEx_Chebyshev_test calculus)
    add_test(NAME cheb_trigonometric COMMAND DACEx_Chebyshev_test trigonometric)
    add_test(NAME cheb_exponential_logarithmic COMMAND DACEx_Chebyshev_test exponential-logarithmic)
    add_test(NAME cheb_power_root COMMAND DACEx_Chebyshev_test power-root)
    add_test(NAME cheb_hyperbolic COMMAND DACEx_Chebyshev_test hyperbolic)
    add_test(NAME cheb_utility COMMAND DACEx_Chebyshev_test utility)
    add_test(NAME cheb_norm_analysis COMMAND DACEx_Chebyshev_test norm-analysis)
    add_test(NAME cheb_evaluation COMMAND DACEx_Chebyshev_test evaluation)
    add_test(NAME cheb_compound_assignment COMMAND DACEx_Chebyshev_test compound-assignment)
    add_test(NAME cheb_multivariate COMMAND DACEx_Chebyshev_test multivariate)
endif()

# Ensure dace.dll is available beside executables on Windows.
if(WIN32)
    if(DACEX_DACE_LIBRARY MATCHES "\\.dll$")
        set(DACEX_DACE_DLL "${DACEX_DACE_LIBRARY}")
    endif()
    if(DACEX_DACE_DLL AND EXISTS "${DACEX_DACE_DLL}")
        foreach(DACEX_EXE_TARGET
            DACEx_Chebyshev_test
            DACEx_Gaussian_example
            DACEx_Gaussian_2D_example)
            add_custom_command(TARGET ${DACEX_EXE_TARGET} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "${DACEX_DACE_DLL}"
                    "$<TARGET_FILE_DIR:${DACEX_EXE_TARGET}>"
            )
        endforeach()
    endif()
endif()
