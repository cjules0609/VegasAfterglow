# === Platform detection ===
ifeq ($(OS),Windows_NT)
    PLATFORM := Windows
    EXE_EXT := .exe
    SHARED_EXT := .pyd
    SHARED_FLAGS := -shared
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Darwin)
        PLATFORM := macOS
        SHARED_EXT := .so
        SHARED_FLAGS := -shared -undefined dynamic_lookup
    else
        PLATFORM := Linux
        SHARED_EXT := .so
        SHARED_FLAGS := -shared -fPIC
    endif
    EXE_EXT :=
endif

# === Paths & sources ===
SRC_DIR := src
TEST_DIR := tests
BUILD_DIR := build
MODULE := vegasglow

# === Compiler & flags ===
CXX ?= g++
CXXFLAGS := -std=c++20 -O3 -march=native -flto -Iinclude -Iexternal -g3 -w  -DXTENSOR_USE_XSIMD=ON
LDFLAGS := -lz
AR := ar
ARFLAGS := rcs

# Core library sources and objects
SRCS := $(shell find $(SRC_DIR) -type f -name '*.cpp')
OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# Tests: each tests/foo/bar.cpp → tests/foo/bar(.exe)
TEST_SRCS := $(shell find $(TEST_DIR) -type f -name '*.cpp')
TEST_EXES := $(patsubst %.cpp,%$(EXE_EXT),$(TEST_SRCS))

# === Default target ===
.PHONY: all lib tests clean test

# Default to building library and tests
all: lib tests

# Individual targets
lib: libvegasglow.a
tests: $(TEST_EXES)

# Build rules
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Static library
libvegasglow.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

# Test executables
%$(EXE_EXT): %.cpp libvegasglow.a
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $< -L. -lvegasglow $(LDFLAGS) -o $@

# Run tests
test: $(TEST_EXES)
	@echo "Running tests..."
	@for test in $(TEST_EXES); do ./$$test || exit 1; done

# Clean up everything
clean:
	rm -rf $(BUILD_DIR) libvegasglow.a $(TEST_EXES)

# Show configuration
.PHONY: info
info:
	@echo "Platform: $(PLATFORM)"
	@echo "Compiler: $(CXX)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"