cmake_minimum_required(VERSION 3.15)
project(VegasAfterglow LANGUAGES CXX)

# — C++17 & PIC
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# — Compiler flags (mirror your setup.py)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-O3 -flto -DNDEBUG -ffast-math -pipe)
elseif(MSVC)
  add_compile_options(/O2 /DNDEBUG /fp:fast /MP /GL)
  add_link_options(/LTCG)
endif()

# — Dependencies
find_package(pybind11 CONFIG REQUIRED)
find_package(ZLIB    REQUIRED)

# — All your C++ sources
file(GLOB_RECURSE ALL_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/pybind/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

# — Build the Python extension as VegasAfterglowC
pybind11_add_module(VegasAfterglowC ${ALL_SOURCES})

# — Include dirs & zlib
target_include_directories(VegasAfterglowC PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/external"
)
target_link_libraries(VegasAfterglowC PRIVATE ZLIB::ZLIB)

# — Default to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()