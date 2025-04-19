# Compiler and flags
CXX         := g++
CXXFLAGS    := -std=c++20 -Iinclude -Iexternal -O3 -w -march=native -DNDEBUG  #-DEXTREME_SPEED -flto

# Directories
SRC_DIR     := src
TEST_DIR    := tests
OBJ_DIR     := build

# Find all source files in src (recursively)
SRC_FILES   := $(shell find $(SRC_DIR) -type f -name "*.cpp")
# Create a list of object files corresponding to the source files (preserving directory structure)
OBJ_FILES   := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))

# Find all test source files in tests (recursively)
TEST_SRCS   := $(shell find $(TEST_DIR) -type f -name "*.cpp")
# Generate executable names: for each test source file in tests/XXX.cpp,
TEST_EXES   := $(TEST_SRCS:%.cpp=%)

# Default target builds the static library and all test executables
all: libtests.a $(TEST_EXES)

# Build object files from src files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Build static library from object files
libtests.a: $(OBJ_FILES)
	@echo "Archiving static library: $@"
	ar rcs $@ $^

# Compile each test source file into an executable, linking with the static library.
$(TEST_DIR)/%: $(TEST_DIR)/%.cpp libtests.a
	$(CXX) $(CXXFLAGS) $< -L. -ltests -o $@

MODULE_NAME := vegasglow
MODULE_SRC := pybind/pybind.cpp pybind/mcmc.cpp
MODULE_OUT := $(MODULE_NAME).so

OMP_PATH := $(shell brew --prefix libomp)
PYBIND11_FLAGS := $(shell python3 -m pybind11 --includes) -Xpreprocessor -fopenmp -I$(OMP_PATH)/include -L$(OMP_PATH)/lib -lomp
PYTHON_LIBS := $(shell python3-config --ldflags)

# Fixed rule (use TAB for indentation, not spaces)
$(MODULE_OUT): $(MODULE_SRC) $(SRC_FILES)
	$(CXX) -arch x86_64 -O3 -w -shared -std=c++20 -DNDEBUG -fPIC \
		$(PYBIND11_FLAGS) $(MODULE_SRC) $(SRC_FILES) -o $(MODULE_OUT) \
		$(PYTHON_LIBS) -Iinclude -Iexternal -undefined dynamic_lookup

pymodule: $(MODULE_OUT)
	@echo "Python module built at: $(MODULE_OUT)"

# Clean up the build and test binary directories and static library.
clean:
	rm -rf $(OBJ_DIR) libtests.a $(MODULE_OUT)

.PHONY: all clean pymodule