# Compiler and flags
CXX         := g++
CXXFLAGS    := -std=c++17 -Iinclude -Iexternal -O3  -w -DNDEBUG  -march=native -flto #-DEXTREME_SPEED

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

# Clean up the build and test binary directories and static library.
clean:
	rm -rf $(OBJ_DIR) libtests.a 

.PHONY: all clean