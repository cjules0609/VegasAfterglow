cmake_minimum_required(VERSION 3.15)
project(VegasAfterglow LANGUAGES CXX)

# — 1) C++17 & PIC
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# — 2) Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-O3 -flto -DNDEBUG -ffast-math -pipe)
elseif(MSVC)
  add_compile_options(/O2 /DNDEBUG /fp:fast /MP /GL)
  add_link_options(/LTCG)
endif()

# — 3) Find pybind11 (must be installed by build-system)
find_package(pybind11 CONFIG REQUIRED)

# — 4) Find or fetch zlib
find_package(ZLIB QUIET)
if (NOT ZLIB_FOUND)
  message(STATUS "ZLIB not found on system → Fetching zlib via FetchContent")
  include(FetchContent)
  FetchContent_Declare(
    zlib
    URL https://zlib.net/zlib-1.2.13.tar.gz
  )
  FetchContent_MakeAvailable(zlib)
  # zlib’s CMake adds targets: zlibstatic (and zlibshared)
  # Create the alias CMake’s FindZLIB would normally create:
  add_library(ZLIB::ZLIB ALIAS zlibstatic)
endif()

# — 5) Collect all your C++ sources
file(GLOB_RECURSE ALL_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/pybind/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

# — 6) Build the Python extension as VegasAfterglowC
pybind11_add_module(VegasAfterglowC ${ALL_SOURCES})

# — 7) Include dirs & link zlib
target_include_directories(VegasAfterglowC PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/external"
)
target_link_libraries(VegasAfterglowC PRIVATE ZLIB::ZLIB)

# — 8) Default to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()