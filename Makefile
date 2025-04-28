# =====================================
# Compiler and Platform detection
# =====================================
ifeq ($(OS),Windows_NT)
    DETECTED_OS := windows
    WSL_CHECK := $(shell uname -a | grep -i microsoft 2>/dev/null || echo "")
    ifneq ($(WSL_CHECK),)
        DETECTED_OS := linux-wsl
    endif
    CXX := $(shell which clang++ 2>/dev/null || which g++ 2>/dev/null || echo g++)
    EXE_SUFFIX := .exe
    SHARED_EXT := dll
else
    DETECTED_OS := $(shell uname -s | tr '[:upper:]' '[:lower:]')
    CXX := $(shell which clang++ 2>/dev/null || which g++ 2>/dev/null || echo g++)
    EXE_SUFFIX :=
    ifeq ($(DETECTED_OS),darwin)
        PLATFORM := macos
        SHARED_EXT := dylib
    else ifeq ($(DETECTED_OS),linux)
        PLATFORM := linux
        SHARED_EXT := so
    else
        PLATFORM := unknown
        SHARED_EXT := so
    endif
endif

ARCH := $(shell uname -m 2>/dev/null || echo unknown)

# =====================================
# Directories
# =====================================
SRC_DIR := src
OBJ_DIR := build
TEST_DIR := tests
PYBIND_DIR := pybind

# =====================================
# Flags
# =====================================
CXXFLAGS_COMMON := -std=c++20 -O3 -g -Wall -Wextra -lz -Iinclude -Iexternal \
                   -DNDEBUG -DXTENSOR_DISABLE_ASSERT -DXTENSOR_DISABLE_CHECK_DIMENSION -DXTENSOR_DISABLE_CHECK_SHAPE

ifeq ($(ARCH),arm64)
    ifeq ($(findstring clang++,$(CXX)),clang++)
        CXXFLAGS_CPU := -mcpu=apple-m1
    endif
else
    MARCH_SUPPORTED := $(shell $(CXX) -march=native -x c++ -E -o /dev/null - </dev/null 2>/dev/null && echo yes || echo no)
    ifeq ($(MARCH_SUPPORTED),yes)
        CXXFLAGS_CPU := -march=native
    endif
endif

ifeq ($(PLATFORM),windows)
    CXXFLAGS_PLATFORM :=
    SHARED_FLAGS := -shared
    LIB_PREFIX :=
    AR := ar
else
    CXXFLAGS_PLATFORM := -fPIC
    ifeq ($(PLATFORM),macos)
        SHARED_FLAGS := -shared -undefined dynamic_lookup
    else
        SHARED_FLAGS := -shared
    endif
    LIB_PREFIX := lib
    AR := ar
endif

LTO_SUPPORTED := $(shell $(CXX) -flto -x c++ -E -o /dev/null - </dev/null 2>/dev/null && echo yes || echo no)
ifeq ($(LTO_SUPPORTED),yes)
    CXXFLAGS_COMMON += -flto
endif

CXXFLAGS := $(CXXFLAGS_COMMON) $(CXXFLAGS_CPU) $(CXXFLAGS_PLATFORM)

# =====================================
# Sources
# =====================================
SRC_FILES := $(shell find $(SRC_DIR) -type f -name "*.cpp" 2>/dev/null)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))
TEST_SRCS := $(shell find $(TEST_DIR) -type f -name "*.cpp" 2>/dev/null)
TEST_EXES := $(patsubst %.cpp,%$(EXE_SUFFIX),$(TEST_SRCS))

MODULE_NAME := vegasglow
MODULE_SRC := $(shell find $(PYBIND_DIR) -type f -name "*.cpp" 2>/dev/null)
MODULE_OUT := $(MODULE_NAME).$(SHARED_EXT)

# =====================================
# Targets
# =====================================
.PHONY: all clean pymodule detect-python install

all: $(LIB_PREFIX)vegasglow.a $(TEST_EXES)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(LIB_PREFIX)vegasglow.a: $(OBJ_FILES)
	@echo "Archiving static library..."
	$(AR) rcs $@ $^

$(TEST_DIR)/%$(EXE_SUFFIX): $(TEST_DIR)/%.cpp $(LIB_PREFIX)vegasglow.a
	@echo "Building test $@..."
	$(CXX) $(CXXFLAGS) $< -L. -lvegasglow -o $@

# =====================================
# Python Module
# =====================================
pymodule: detect-python $(MODULE_OUT)

detect-python:
	@echo "Detecting Python environment..."
	$(eval PYTHON_BIN := $(shell which python3 2>/dev/null || which python 2>/dev/null))
	@if [ -z "$(PYTHON_BIN)" ]; then echo "Error: Python not found!"; exit 1; fi

	$(eval PYTHON_VERSION := $(shell $(PYTHON_BIN) -c "import sys; print('{}.{}'.format(sys.version_info.major, sys.version_info.minor))"))
	$(eval PYTHON_ROOT := $(shell $(PYTHON_BIN) -c "import sys; print(sys.prefix)"))
	$(eval PYTHON_INCLUDE_DIR := $(shell $(PYTHON_BIN) -c "import sysconfig; print(sysconfig.get_path('include'))"))
	$(eval PYTHON_INCLUDE_PLATDIR := $(shell $(PYTHON_BIN) -c "import sysconfig; print(sysconfig.get_path('platinclude'))"))

	@echo "Python root: $(PYTHON_ROOT)"
	@echo "Python include: $(PYTHON_INCLUDE_DIR)"
	@echo "Python platinclude: $(PYTHON_INCLUDE_PLATDIR)"

	# macOS special fix
	@if [ "$(PLATFORM)" = "macos" ]; then \
		if [ ! -f "$(PYTHON_INCLUDE_DIR)/Python.h" ]; then \
			echo "Python.h not found, searching macOS frameworks..."; \
			if [ -f "$(PYTHON_ROOT)/Headers/Python.h" ]; then \
				echo "Found Python.h at $(PYTHON_ROOT)/Headers"; \
				$(eval PYTHON_INCLUDE_DIR := $(PYTHON_ROOT)/Headers) ; \
			else \
				echo "Error: Cannot find Python.h. Install Xcode CommandLineTools!"; exit 1; \
			fi ; \
		fi ; \
	fi

	# Final Python include flags
	$(eval PYTHON_INCLUDES := -I$(PYTHON_INCLUDE_DIR))
	@if [ "$(PYTHON_INCLUDE_DIR)" != "$(PYTHON_INCLUDE_PLATDIR)" ] && [ -n "$(PYTHON_INCLUDE_PLATDIR)" ]; then \
		$(eval PYTHON_INCLUDES := $(PYTHON_INCLUDES) -I$(PYTHON_INCLUDE_PLATDIR)) ; \
	fi

	# Python libs
	$(eval PYTHON_CONFIG := $(shell which python3-config 2>/dev/null || which python-config 2>/dev/null))
	ifneq ($(PYTHON_CONFIG),)
		$(eval PYTHON_LIBS := $(shell $(PYTHON_CONFIG) --ldflags 2>/dev/null))
	else
		$(eval PYTHON_LIBS := $(shell $(PYTHON_BIN) -c "import sysconfig; print('-L' + sysconfig.get_config_var('LIBDIR') + ' -lpython' + sysconfig.get_config_var('VERSION') + (sysconfig.get_config_var('ABIFLAGS') or ''))"))
	endif

	$(eval PYTHON_ARCH := $(shell $(PYTHON_BIN) -c "import platform; print(platform.machine())"))

	# pybind11 and numpy check
	$(eval PYBIND11_INCLUDE := $(shell $(PYTHON_BIN) -m pybind11 --includes 2>/dev/null))
	$(eval NUMPY_INCLUDE := $(shell $(PYTHON_BIN) -c "import numpy; print(numpy.get_include())"))

	@echo "Final Python settings:"
	@echo "  PYTHON_INCLUDES = $(PYTHON_INCLUDES)"
	@echo "  PYTHON_LIBS = $(PYTHON_LIBS)"
	@echo "  PYBIND11_INCLUDE = $(PYBIND11_INCLUDE)"
	@echo "  NUMPY_INCLUDE = $(NUMPY_INCLUDE)"

$(MODULE_OUT): $(MODULE_SRC) $(SRC_FILES)
	@echo "Building Python module..."
	$(CXX) $(SHARED_FLAGS) \
	$(CXXFLAGS_COMMON) $(PYTHON_INCLUDES) $(PYBIND11_INCLUDE) -I$(NUMPY_INCLUDE) \
	$(MODULE_SRC) $(SRC_FILES) \
	$(PYTHON_LIBS) -o $(MODULE_OUT) \
	$(if $(filter $(PLATFORM),macos),-arch $(ARCH))

# =====================================
# Clean
# =====================================
clean:
	@echo "Cleaning..."
	rm -rf $(OBJ_DIR) $(LIB_PREFIX)vegasglow.a $(MODULE_OUT) $(TEST_EXES)

# =====================================
# Install
# =====================================
install: $(LIB_PREFIX)vegasglow.a
	@echo "Installing..."
	install -d /usr/local/lib
	install -m 644 $(LIB_PREFIX)vegasglow.a /usr/local/lib
	install -d /usr/local/include/vegasglow
	find include -name "*.hpp" -exec install -D -m 644 {} /usr/local/include/{} \;