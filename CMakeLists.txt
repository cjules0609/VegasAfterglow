cmake_minimum_required(VERSION 3.15)
project(vegasglow LANGUAGES CXX)

# Always verbose build output for easier debugging
set(CMAKE_VERBOSE_MAKEFILE ON)

# Set C++ standard cleanly
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default build type to Release if not set
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Let pybind11 use modern Python finding (important for scikit-build-core)
set(PYBIND11_FINDPYTHON ON)

# Find pybind11 and Python
find_package(pybind11 CONFIG REQUIRED)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# Gather all C++ source files
file(GLOB_RECURSE SRC_FILES src/*.cpp)

# Build static C++ library
add_library(vegasglow_static STATIC ${SRC_FILES})
target_include_directories(vegasglow_static PUBLIC include external)
set_target_properties(vegasglow_static PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Apply warnings and useful compile flags
target_compile_options(vegasglow_static PUBLIC
    -Wall -Wextra
    -Wno-unused-parameter
    -Wno-unused-but-set-variable
    -Wno-sign-compare
    -Wno-deprecated-declarations
)

# Build the Python module (only pybind files)
pybind11_add_module(vegasglow
    pybind/pybind.cpp
    pybind/mcmc.cpp
)

target_link_libraries(vegasglow PRIVATE vegasglow_static)
target_include_directories(vegasglow PRIVATE include external)

# Install for pip install
include(GNUInstallDirs)
install(TARGETS vegasglow
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)