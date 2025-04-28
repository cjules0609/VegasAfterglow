# Compiler
CXX := g++
CXXFLAG := -std=c++20 -flto -Iinclude -Iexternal -g -lz -O3 -w -DNDEBUG \
           -DXTENSOR_DISABLE_ASSERT -DXTENSOR_DISABLE_CHECK_DIMENSION -DXTENSOR_DISABLE_CHECK_SHAPE
CXXFLAGS := $(CXXFLAG) -march=native

# Python and Pybind11
PYTHON_BIN := python3

PYTHON_INCLUDE := $(shell $(PYTHON_BIN)-config --includes)
PYTHON_LIBS := $(shell $(PYTHON_BIN)-config --ldflags)
NUMPY_INCLUDE := $(shell $(PYTHON_BIN) -c "import numpy; print(numpy.get_include())")
PYBIND11_INCLUDE := $(shell $(PYTHON_BIN) -m pybind11 --includes)

# Directories
SRC_DIR := src
TEST_DIR := tests
OBJ_DIR := build

# Source files
SRC_FILES := $(shell find $(SRC_DIR) -type f -name "*.cpp")
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))
TEST_SRCS := $(shell find $(TEST_DIR) -type f -name "*.cpp")
TEST_EXES := $(TEST_SRCS:%.cpp=%)

# Module
MODULE_NAME := vegasglow
MODULE_SRC := pybind/pybind.cpp pybind/mcmc.cpp
MODULE_OUT := $(MODULE_NAME).so

# Default build: static lib + test binaries
all: libvegasglow.a $(TEST_EXES)

# Build object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Build static library
libvegasglow.a: $(OBJ_FILES)
	@echo "Archiving static library: $@"
	ar rcs $@ $^

# Build test executables
$(TEST_DIR)/%: $(TEST_DIR)/%.cpp libvegasglow.a
	$(CXX) $(CXXFLAGS) $< -L. -lvegasglow -o $@

# Build python module
$(MODULE_OUT): $(MODULE_SRC) $(SRC_FILES)
	$(CXX) -shared -fPIC -arch x86_64 \
		$(PYBIND11_INCLUDE) $(PYTHON_INCLUDE) -I$(NUMPY_INCLUDE) \
		$(MODULE_SRC) $(SRC_FILES) -o $(MODULE_OUT) \
		$(PYTHON_LIBS) $(CXXFLAG) -undefined dynamic_lookup

# Shortcut target
pymodule: $(MODULE_OUT)
	@echo "Python module built at: $(MODULE_OUT)"

# Clean
clean:
	rm -rf $(OBJ_DIR) libvegasglow.a $(MODULE_OUT) $(TEST_EXES)

.PHONY: all clean pymodule