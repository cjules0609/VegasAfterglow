cmake_minimum_required(VERSION 3.15)

project(vegasglow LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# 优化选项
if(MSVC)
    add_compile_options(/O2 /W4 /DNDEBUG)
else()
    add_compile_options(
        -O3 -g -Wall -Wextra -std=c++20 -lz
        -DNDEBUG
        -DXTENSOR_DISABLE_ASSERT
        -DXTENSOR_DISABLE_CHECK_DIMENSION
        -DXTENSOR_DISABLE_CHECK_SHAPE
    )
endif()

# 查找Python环境
find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)

# 自动安装 numpy
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import numpy"
    RESULT_VARIABLE numpy_check
)
if(NOT numpy_check EQUAL 0)
    message(STATUS "numpy not found, installing...")
    execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install numpy)
endif()

# 自动安装 pybind11
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import pybind11"
    RESULT_VARIABLE pybind11_check
)
if(NOT pybind11_check EQUAL 0)
    message(STATUS "pybind11 not found, installing...")
    execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install pybind11)
endif()

# 获取 NumPy include
execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NumPy_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# 获取 pybind11 include
execute_process(
    COMMAND "${Python3_EXECUTABLE}" -m pybind11 --includes
    OUTPUT_VARIABLE Pybind11_INCLUDE_FLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
separate_arguments(Pybind11_INCLUDE_FLAGS)

foreach(INC ${Pybind11_INCLUDE_FLAGS})
    string(REGEX REPLACE "^-I" "" CLEANED_INC "${INC}")
    list(APPEND EXTRA_INCLUDES ${CLEANED_INC})
endforeach()

# === 1️⃣ 编译静态库 ===

file(GLOB_RECURSE SOURCES
    src/*.cpp
)

add_library(vegasglow_static STATIC ${SOURCES})

target_include_directories(vegasglow_static PRIVATE
    include
    external
)

# === 2️⃣ 编译 Python 模块 (可pip install) ===

file(GLOB_RECURSE PYBIND_SOURCES
    pybind/*.cpp
)

add_library(vegasglow MODULE ${SOURCES} ${PYBIND_SOURCES})
set_target_properties(vegasglow PROPERTIES PREFIX "")

target_include_directories(vegasglow PRIVATE
    include
    external
    ${NumPy_INCLUDE_DIR}
    ${EXTRA_INCLUDES}
)

# Windows特有处理
if(Python3_LIBRARIES)
    target_link_libraries(vegasglow PRIVATE ${Python3_LIBRARIES})
endif()

# === 3️⃣ 编译测试程序 ===

file(GLOB TEST_FILES
    tests/*.cpp
)

foreach(test_file ${TEST_FILES})
    get_filename_component(test_name ${test_file} NAME_WE)
    add_executable(${test_name} ${test_file})

    target_include_directories(${test_name} PRIVATE
        include
        external
    )

    target_link_libraries(${test_name}
        vegasglow_static
    )
endforeach()