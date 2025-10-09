# Explicitly set default target - must be before all variable definitions and $(info)
.DEFAULT_GOAL := default

# Compiler
CC = gcc

# Compiler flags with parallel LTO
CFLAGS = -O3 -march=native -fopenmp -fPIC -flto=auto
# -fsanitize=address -static-libasan

# Linker flags with parallel LTO
LDFLAGS = -fopenmp -flto=auto

# Library paths
FLINT_LIB_PATH = /usr/local/lib
PML_LIB_PATH = /usr/local/lib

# Include paths - using paths from C_INCLUDE_PATH environment variable
# Extract paths from environment variable, use default if not set
FLINT_INCLUDE_PATH ?= /usr/local/include
PML_INCLUDE_PATH ?= /usr/local/include

# Check if C_INCLUDE_PATH environment variable is set
ifdef C_INCLUDE_PATH
$(info Using C_INCLUDE_PATH: $(C_INCLUDE_PATH))
endif

# Local directories
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build

# Get system library paths in a simple way
SYSTEM_LIB_PATHS := $(shell echo "$$LD_LIBRARY_PATH:/usr/lib:/usr/local/lib:/usr/lib64:/usr/local/lib64:/lib:/lib64:/usr/lib/x86_64-linux-gnu" | tr ':' ' ' | tr -s ' ')

# Smarter header file check - use compiler to find header files
# Include C_INCLUDE_PATH in the test if it's set
ifdef C_INCLUDE_PATH
HEADER_TEST_FLAGS := -I./$(INCLUDE_DIR) $(foreach path,$(subst :, ,$(C_INCLUDE_PATH)),-I$(path))
else
HEADER_TEST_FLAGS := -I./$(INCLUDE_DIR)
endif

# Use a simple test file approach for more reliable detection
TEMP_TEST_FILE := .header_test.c

# Helper function to test header availability
define test_header
$(shell echo '#include <$(1)>' > $(TEMP_TEST_FILE) && $(CC) -E $(HEADER_TEST_FLAGS) $(TEMP_TEST_FILE) >/dev/null 2>&1 && echo yes || echo no; rm -f $(TEMP_TEST_FILE))
endef

# Simple library search functions
define find_pml_so
$(shell for path in $(SYSTEM_LIB_PATHS); do \
	if [ -n "$$path" ] && [ -d "$$path" ]; then \
		if [ -f "$$path/libpml.so" ]; then echo "$$path/libpml.so"; exit 0; fi; \
		for f in "$$path"/libpml.so.*; do \
			if [ -f "$$f" ]; then echo "$$f"; exit 0; fi; \
		done; \
	fi; \
done)
endef

define find_pml_a
$(shell for path in $(SYSTEM_LIB_PATHS); do \
	if [ -n "$$path" ] && [ -f "$$path/libpml.a" ]; then \
		echo "$$path/libpml.a"; exit 0; \
	fi; \
done)
endef

define check_pml_so
$(shell for path in $(SYSTEM_LIB_PATHS); do \
	if [ -n "$$path" ] && [ -d "$$path" ]; then \
		if [ -f "$$path/libpml.so" ]; then echo yes; exit 0; fi; \
		for f in "$$path"/libpml.so.*; do \
			if [ -f "$$f" ]; then echo yes; exit 0; fi; \
		done; \
	fi; \
done; echo no)
endef

define check_pml_a
$(shell for path in $(SYSTEM_LIB_PATHS); do \
	if [ -n "$$path" ] && [ -f "$$path/libpml.a" ]; then \
		echo yes; exit 0; \
	fi; \
done; echo no)
endef

FLINT_HEADER_CHECK := $(call test_header,flint/flint.h)
PML_HEADER_CHECK := $(call test_header,pml.h)
NMOD_POLY_MAT_UTILS_CHECK := $(call test_header,nmod_poly_mat_utils.h)
NMOD_POLY_MAT_EXTRA_CHECK := $(call test_header,nmod_poly_mat_extra.h)

# Check for PML library files in all system library paths
PML_DYNAMIC_LIB_CHECK := $(call check_pml_so)
PML_STATIC_LIB_CHECK := $(call check_pml_a)

# Find actual PML library paths for linking
PML_SO_PATH := $(call find_pml_so)
PML_A_PATH := $(call find_pml_a)

# PML is only available if ALL required headers AND at least one library file are found
PML_AVAILABLE := $(shell if [ "$(PML_HEADER_CHECK)" = "yes" ] && [ "$(NMOD_POLY_MAT_UTILS_CHECK)" = "yes" ] && [ "$(NMOD_POLY_MAT_EXTRA_CHECK)" = "yes" ] && ([ "$(PML_DYNAMIC_LIB_CHECK)" = "yes" ] || [ "$(PML_STATIC_LIB_CHECK)" = "yes" ]); then echo yes; else echo no; fi)

# Old directory check (as fallback)
FLINT_DIR_EXISTS := $(shell if [ -d "$(FLINT_INCLUDE_PATH)" ]; then echo yes; else echo no; fi)
PML_DIR_EXISTS := $(shell if [ -d "$(PML_INCLUDE_PATH)" ]; then echo yes; else echo no; fi)

# Validate FLINT (required) - use header check instead of directory check
ifneq ($(FLINT_HEADER_CHECK),yes)
$(error FLINT headers not found! Compiler cannot find flint/flint.h. Check C_INCLUDE_PATH or install FLINT development packages)
endif

# System libraries
SYSTEM_LIBS = -lmpfr -lgmp -lm -lpthread -lstdc++

# Include flags (including local include directory)
# Add C_INCLUDE_PATH directories if set
INCLUDE_FLAGS = -I./$(INCLUDE_DIR)
ifdef C_INCLUDE_PATH
INCLUDE_FLAGS += $(foreach path,$(subst :, ,$(C_INCLUDE_PATH)),-I$(path))
endif

# Library flags
FLINT_FLAGS = -DHAVE_FLINT
PML_FLAGS = 
# Use our improved PML availability check
ifeq ($(PML_AVAILABLE),yes)
PML_FLAGS = -DHAVE_PML
endif

# Combined CFLAGS
ALL_CFLAGS = $(CFLAGS) $(INCLUDE_FLAGS) $(FLINT_FLAGS) $(PML_FLAGS)

# FLINT library linking
FLINT_LIBS = -L$(FLINT_LIB_PATH) -lflint
FLINT_STATIC_LIBS = $(FLINT_LIB_PATH)/libflint.a

# PML library linking (optional) - use found paths or fallback to -lpml
PML_LIBS = 
PML_STATIC_LIBS = 
# Use our improved PML availability check
ifeq ($(PML_AVAILABLE),yes)
ifneq ($(PML_SO_PATH),)
# Use specific path if found
PML_LIBS = $(PML_SO_PATH)
else
# Fallback to -lpml (let linker find it)
PML_LIBS = -lpml
endif
ifneq ($(PML_A_PATH),)
PML_STATIC_LIBS = $(PML_A_PATH)
else
PML_STATIC_LIBS = $(PML_LIB_PATH)/libpml.a
endif
endif

# Combined external libraries
EXTERNAL_LIBS = $(FLINT_LIBS) $(PML_LIBS) $(SYSTEM_LIBS)
EXTERNAL_STATIC_PML_LIBS = $(FLINT_LIBS) $(PML_STATIC_LIBS) $(SYSTEM_LIBS)
EXTERNAL_STATIC_ALL_LIBS = $(FLINT_STATIC_LIBS) $(PML_STATIC_LIBS) $(SYSTEM_LIBS)

# Runtime library path (rpath) flags
RPATH_FLAGS = -Wl,-rpath,$(FLINT_LIB_PATH) -Wl,-rpath,.
ifeq ($(PML_AVAILABLE),yes)
RPATH_FLAGS += -Wl,-rpath,$(PML_LIB_PATH)
endif

# Source files for the math library (in src directory)
MATH_SOURCES = $(SRC_DIR)/dixon_complexity.c \
               $(SRC_DIR)/dixon_flint.c \
               $(SRC_DIR)/dixon_interface_flint.c \
               $(SRC_DIR)/dixon_test.c \
               $(SRC_DIR)/dixon_with_ideal_reduction.c \
               $(SRC_DIR)/fq_mat_det.c \
               $(SRC_DIR)/fq_mpoly_mat_det.c \
               $(SRC_DIR)/fq_multivariate_interpolation.c \
               $(SRC_DIR)/fq_mvpoly.c \
               $(SRC_DIR)/fq_nmod_roots.c \
               $(SRC_DIR)/fq_poly_mat_det.c \
               $(SRC_DIR)/fq_sparse_interpolation.c \
               $(SRC_DIR)/fq_unified_interface.c \
               $(SRC_DIR)/gf2128_mpoly.c \
               $(SRC_DIR)/gf28_mpoly.c \
               $(SRC_DIR)/gf2n_field.c \
               $(SRC_DIR)/gf2n_poly.c \
               $(SRC_DIR)/polynomial_system_solver.c \
               $(SRC_DIR)/resultant_with_ideal_reduction.c \
               $(SRC_DIR)/unified_mpoly_det.c \
               $(SRC_DIR)/unified_mpoly_interface.c \
               $(SRC_DIR)/unified_mpoly_resultant.c

# Object files (in build directory)
MATH_OBJECTS = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(MATH_SOURCES))

# Main source file (in current directory)
DIXON_SRC = dixon.c

# All source files (for LTO compilation)
ALL_SOURCES = $(DIXON_SRC) $(MATH_SOURCES)

# Library names (in current directory)
DIXON_STATIC_LIB = libdixon.a
DIXON_SHARED_LIB = libdixon.so

# Output executable (in current directory)
DIXON_TARGET = dixon

# Attack programs directory and files
ATTACK_DIR = ../Attack
# Find all C files recursively in Attack directory and subdirectories
ATTACK_C_FILES := $(shell find $(ATTACK_DIR) -name "*.c" 2>/dev/null | grep -v ".ipynb_checkpoints" || echo "")
ATTACK_EXECUTABLES := $(patsubst %.c,%,$(ATTACK_C_FILES))

# Create build directory
$(BUILD_DIR):
	@echo "Creating build directory..."
	mkdir -p $(BUILD_DIR)

# Default target - first build libraries, then compile all sources with LTO for maximum inlining, then build attack programs
default: $(DIXON_STATIC_LIB) $(DIXON_SHARED_LIB)
	@echo "Building $(DIXON_TARGET) with LTO (Link Time Optimization)..."
	@echo "Libraries built, now compiling all sources together for maximum inlining..."
	$(CC) $(ALL_CFLAGS) -o $(DIXON_TARGET) $(ALL_SOURCES) $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS)
	@echo "Build complete: $(DIXON_TARGET) (LTO optimized with libraries available)"
	@echo ""
	@echo "Now building Attack programs..."
	@$(MAKE) attack-programs-verbose

# Also build libraries with LTO for better performance
all: default
	@echo "Built dixon executable, libraries, and attack programs with LTO optimization"

# LTO target - compile all sources together for maximum optimization (same as default now)
$(DIXON_TARGET)-lto: $(DIXON_STATIC_LIB) $(DIXON_SHARED_LIB)
	@echo "Building $(DIXON_TARGET) with LTO (Link Time Optimization)..."
	@echo "Libraries built, now compiling all sources together for maximum inlining..."
	$(CC) $(ALL_CFLAGS) -o $(DIXON_TARGET) $(ALL_SOURCES) $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS)
	@echo "Build complete: $(DIXON_TARGET) (LTO optimized)"

# Traditional dynamic library target
$(DIXON_TARGET)-dynamic: $(DIXON_SRC) $(DIXON_SHARED_LIB)
	@echo "Building $(DIXON_TARGET) with dynamic dixon library..."
	$(CC) $(ALL_CFLAGS) -o $(DIXON_TARGET) $< -L. -ldixon $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS)
	@echo "Build complete: $(DIXON_TARGET) (dynamic dixon, dynamic FLINT/PML)"

# Build dynamic dixon library
dynamic-lib: $(DIXON_SHARED_LIB)

$(DIXON_SHARED_LIB): $(MATH_OBJECTS)
	@echo "Building dynamic dixon library..."
	$(CC) -shared -o $@ $^ $(EXTERNAL_LIBS) $(LDFLAGS)
	@echo "Dynamic library built: $(DIXON_SHARED_LIB)"

# Build static dixon library
static-lib: $(DIXON_STATIC_LIB)

$(DIXON_STATIC_LIB): $(MATH_OBJECTS)
	@echo "Building static dixon library..."
	ar rcs $@ $^
	@echo "Static library built: $(DIXON_STATIC_LIB)"

# Build with static dixon library (but dynamic FLINT/PML)
static: $(DIXON_TARGET)-static
	@echo "Built dixon with static library"

$(DIXON_TARGET)-static: $(DIXON_SRC) $(DIXON_STATIC_LIB)
	@echo "Building $(DIXON_TARGET) with static dixon library (dynamic FLINT/PML)..."
	$(CC) $(ALL_CFLAGS) -o $(DIXON_TARGET) $< $(DIXON_STATIC_LIB) $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS)
	@echo "Build complete: $(DIXON_TARGET) (static dixon, dynamic FLINT/PML)"

# Build with static dixon + static PML (but dynamic FLINT)
static-pml: static-lib $(DIXON_TARGET)-static-pml

$(DIXON_TARGET)-static-pml: $(DIXON_SRC) $(DIXON_STATIC_LIB)
	@echo "Building $(DIXON_TARGET) with static dixon + PML libraries (dynamic FLINT)..."
	$(CC) $(ALL_CFLAGS) -o $(DIXON_TARGET) $< $(DIXON_STATIC_LIB) $(EXTERNAL_STATIC_PML_LIBS) $(RPATH_FLAGS) $(LDFLAGS)
	@echo "Build complete: $(DIXON_TARGET) (static dixon+PML, dynamic FLINT)"

# Build with all static libraries (dixon + PML + FLINT)
static-all: static-lib $(DIXON_TARGET)-static-all

$(DIXON_TARGET)-static-all: $(DIXON_SRC) $(DIXON_STATIC_LIB)
	@echo "Building $(DIXON_TARGET) with all static libraries..."
	$(CC) $(ALL_CFLAGS) -o $(DIXON_TARGET) $< $(DIXON_STATIC_LIB) $(EXTERNAL_STATIC_ALL_LIBS) $(LDFLAGS)
	@echo "Build complete: $(DIXON_TARGET) (fully static)"

# Attack programs compilation with verbose output
attack-programs-verbose: $(DIXON_STATIC_LIB) $(DIXON_SHARED_LIB)
	@echo "Building Attack programs with verbose output..."
	@if [ -z "$(ATTACK_C_FILES)" ]; then \
		echo "No C files found in $(ATTACK_DIR)"; \
		echo "Checked directory: $(ATTACK_DIR)"; \
		if [ ! -d "$(ATTACK_DIR)" ]; then \
			echo "Directory $(ATTACK_DIR) does not exist!"; \
		fi; \
	else \
		echo "Found $(words $(ATTACK_C_FILES)) C files:"; \
		for cfile in $(ATTACK_C_FILES); do \
			echo "  $$cfile"; \
		done; \
		echo ""; \
		success_count=0; \
		fail_count=0; \
		for cfile in $(ATTACK_C_FILES); do \
			if [ -f "$$cfile" ]; then \
				executable="$${cfile%.c}"; \
				echo "Compiling $$cfile -> $$executable"; \
				echo "Command: $(CC) $(ALL_CFLAGS) -o \"$$executable\" \"$$cfile\" -L. -ldixon $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS)"; \
				if $(CC) $(ALL_CFLAGS) -o "$$executable" "$$cfile" -L. -ldixon $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS) 2>&1; then \
					echo "Successfully compiled: $$executable"; \
					ls -la "$$executable" | sed 's/^/  /'; \
					success_count=$$((success_count + 1)); \
				else \
					echo "Failed to compile $$cfile"; \
					fail_count=$$((fail_count + 1)); \
				fi; \
				echo ""; \
			else \
				echo "File not found: $$cfile"; \
				fail_count=$$((fail_count + 1)); \
			fi; \
		done; \
		echo "Attack programs compilation summary:"; \
		echo "  Success: $$success_count"; \
		echo "  Failed: $$fail_count"; \
		echo "  Total: $(words $(ATTACK_C_FILES))"; \
	fi

# Attack programs compilation (silent version for default target)
attack-programs: $(DIXON_STATIC_LIB) $(DIXON_SHARED_LIB)
	@echo "Building Attack programs..."
	@if [ -z "$(ATTACK_C_FILES)" ]; then \
		echo "No C files found in $(ATTACK_DIR)"; \
	else \
		echo "Found C files: $(ATTACK_C_FILES)"; \
		for cfile in $(ATTACK_C_FILES); do \
			if [ -f "$$cfile" ]; then \
				executable="$${cfile%.c}"; \
				echo ""; \
				echo "Compiling $$cfile -> $$executable"; \
				echo "Command: $(CC) $(ALL_CFLAGS) -o \"$$executable\" \"$$cfile\" -L. -ldixon $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS)"; \
				if $(CC) $(ALL_CFLAGS) -o "$$executable" "$$cfile" -L. -ldixon $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS); then \
					echo "Successfully compiled: $$executable"; \
					ls -la "$$executable" | sed 's/^/  /'; \
				else \
					echo "Failed to compile $$cfile"; \
				fi; \
			fi; \
		done; \
		echo ""; \
		echo "Attack programs compilation complete"; \
	fi

# Attack programs with static dixon library
attack-static: $(DIXON_STATIC_LIB)
	@echo "Building Attack programs with static dixon library..."
	@if [ -z "$(ATTACK_C_FILES)" ]; then \
		echo "No C files found in $(ATTACK_DIR)"; \
	else \
		echo "Found C files: $(ATTACK_C_FILES)"; \
		for cfile in $(ATTACK_C_FILES); do \
			if [ -f "$$cfile" ]; then \
				executable="$${cfile%.c}"; \
				echo "Compiling $$cfile -> $$executable (static dixon)"; \
				$(CC) $(ALL_CFLAGS) -o "$$executable" "$$cfile" $(DIXON_STATIC_LIB) $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS) || echo "Failed to compile $$cfile"; \
			fi; \
		done; \
		echo "Attack programs compilation complete (static dixon)"; \
	fi

# Clean attack programs
clean-attack:
	@echo "Cleaning Attack programs..."
	@if [ -n "$(ATTACK_EXECUTABLES)" ]; then \
		for exe in $(ATTACK_EXECUTABLES); do \
			if [ -f "$$exe" ]; then \
				echo "Removing $$exe"; \
				rm -f "$$exe"; \
			fi; \
		done; \
	fi
	@echo "Attack programs cleaned"

# Object file compilation (src/*.c -> build/*.o)
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	@echo "Compiling $<..."
	$(CC) $(ALL_CFLAGS) -c -o $@ $<

# Clean
clean: clean-attack
	rm -f $(DIXON_TARGET) $(DIXON_STATIC_LIB) $(DIXON_SHARED_LIB)
	rm -rf $(BUILD_DIR)
	@echo "Cleaned all build artifacts"

# Clean only build directory (keep executables and libraries)
clean-build:
	rm -rf $(BUILD_DIR)
	@echo "Cleaned build directory"

# Test library detection
test-paths:
	@echo "Testing library path detection..."
	@echo "LD_LIBRARY_PATH: $$LD_LIBRARY_PATH"
	@echo "SYSTEM_LIB_PATHS: $(SYSTEM_LIB_PATHS)"
	@echo "PML_DYNAMIC_LIB_CHECK: $(PML_DYNAMIC_LIB_CHECK)"
	@echo "PML_STATIC_LIB_CHECK: $(PML_STATIC_LIB_CHECK)"
	@echo "PML_SO_PATH: $(PML_SO_PATH)"
	@echo "PML_A_PATH: $(PML_A_PATH)"
	@echo "PML_AVAILABLE: $(PML_AVAILABLE)"

# Test Attack directory detection
test-attack:
	@echo "Testing attack detection..."
	@echo "ATTACK_DIR: $(ATTACK_DIR)"
	@echo "ATTACK_C_FILES: $(ATTACK_C_FILES)"
	@echo "ATTACK_EXECUTABLES: $(ATTACK_EXECUTABLES)"
	@echo ""
	@echo "Manual find test:"
	@find $(ATTACK_DIR) -name "*.c" 2>/dev/null | grep -v ".ipynb_checkpoints" || echo "No files found or directory doesn't exist"

# Show configuration
info:
	@echo "=== Build Configuration ==="
	@echo "CC: $(CC)"
	@echo "CFLAGS: $(ALL_CFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
	@echo "RPATH_FLAGS: $(RPATH_FLAGS)"
	@echo ""
	@echo "=== Directory Structure ==="
	@echo "Source directory: $(SRC_DIR)/"
	@echo "Include directory: $(INCLUDE_DIR)/"
	@echo "Build directory: $(BUILD_DIR)/"
	@echo "Output directory: ./"
	@echo ""
	@echo "=== System Library Paths ==="
	@echo "LD_LIBRARY_PATH: $$LD_LIBRARY_PATH"
	@echo "Library search paths: $(SYSTEM_LIB_PATHS)"
	@echo ""
	@echo "=== Library Status ==="
	@echo "FLINT headers found: $(FLINT_HEADER_CHECK)"
	@echo "FLINT directory exists: $(FLINT_DIR_EXISTS) at $(FLINT_INCLUDE_PATH)"
	@echo "PML headers found: $(PML_HEADER_CHECK)"
	@echo "nmod_poly_mat_utils.h found: $(NMOD_POLY_MAT_UTILS_CHECK)"
	@echo "nmod_poly_mat_extra.h found: $(NMOD_POLY_MAT_EXTRA_CHECK)"
	@echo "PML dynamic library found: $(PML_DYNAMIC_LIB_CHECK)"
	@echo "PML static library found: $(PML_STATIC_LIB_CHECK)"
	@echo "PML available (all headers + libraries): $(PML_AVAILABLE)"
	@echo "PML directory exists: $(PML_DIR_EXISTS) at $(PML_INCLUDE_PATH)"
	@echo ""
	@echo "=== Found Library Paths ==="
	@echo "PML dynamic library path: $(PML_SO_PATH)"
	@echo "PML static library path: $(PML_A_PATH)"
	@echo ""
	@echo "=== Library Paths ==="
	@echo "FLINT lib: $(FLINT_LIB_PATH)"
	@echo "PML lib: $(PML_LIB_PATH)"
	@echo ""
	@echo "=== Library Files ==="
	@echo "FLINT dynamic: $(FLINT_LIBS)"
	@echo "FLINT static: $(FLINT_STATIC_LIBS)"
	@echo "PML dynamic: $(PML_LIBS)"
	@echo "PML static: $(PML_STATIC_LIBS)"
	@echo ""
	@echo "=== Attack Directory Structure ==="
	@echo "Attack directory: $(ATTACK_DIR)"
	@echo -n "Directory exists: "
	@if [ -d "$(ATTACK_DIR)" ]; then \
		echo "YES"; \
		echo "C files found:"; \
		if [ -n "$(ATTACK_C_FILES)" ]; then \
			for cfile in $(ATTACK_C_FILES); do \
				echo "  $$cfile"; \
			done; \
		else \
			echo "  No C files found"; \
		fi; \
	else \
		echo "NO"; \
	fi
	@echo "EXTERNAL_LIBS (dynamic): $(EXTERNAL_LIBS)"
	@echo "EXTERNAL_STATIC_PML_LIBS: $(EXTERNAL_STATIC_PML_LIBS)"
	@echo "EXTERNAL_STATIC_ALL_LIBS: $(EXTERNAL_STATIC_ALL_LIBS)"

# Debug header file detection
debug-headers:
	@echo "=== Header File Detection Debug ==="
	@echo ""
	@echo "=== Compiler Search Paths ==="
	@echo "Getting GCC include search paths..."
	@$(CC) -E -v -x c /dev/null 2>&1 | sed -n '/#include <...> search starts here:/,/End of search list./p' | sed 's/^/ /'
	@echo ""
	@echo "=== Environment Variables ==="
	@echo "C_INCLUDE_PATH: $(C_INCLUDE_PATH)"
	@echo "CPLUS_INCLUDE_PATH: $(CPLUS_INCLUDE_PATH)"
	@echo ""
	@echo "=== Header File Tests ==="
	@echo -n "FLINT headers (flint/flint.h): "
	@echo '#include <flint/flint.h>' > .header_test.c && \
		if $(CC) -E $(HEADER_TEST_FLAGS) .header_test.c >/dev/null 2>&1; then \
			echo "FOUND"; \
		else \
			echo "NOT FOUND"; \
		fi; \
		rm -f .header_test.c
	@echo -n "PML headers (pml.h): "
	@echo '#include <pml.h>' > .header_test.c && \
		if $(CC) -E $(HEADER_TEST_FLAGS) .header_test.c >/dev/null 2>&1; then \
			echo "FOUND"; \
		else \
			echo "NOT FOUND"; \
		fi; \
		rm -f .header_test.c
	@echo -n "nmod_poly_mat_utils.h: "
	@echo '#include <nmod_poly_mat_utils.h>' > .header_test.c && \
		if $(CC) -E $(HEADER_TEST_FLAGS) .header_test.c >/dev/null 2>&1; then \
			echo "FOUND"; \
		else \
			echo "NOT FOUND"; \
		fi; \
		rm -f .header_test.c
	@echo -n "nmod_poly_mat_extra.h: "
	@echo '#include <nmod_poly_mat_extra.h>' > .header_test.c && \
		if $(CC) -E $(HEADER_TEST_FLAGS) .header_test.c >/dev/null 2>&1; then \
			echo "FOUND"; \
		else \
			echo "NOT FOUND"; \
		fi; \
		rm -f .header_test.c
	@echo "PML Available (all required headers + libraries): $(PML_AVAILABLE)"
	@echo ""
	@echo "=== Manual Path Search ==="
	@echo "Searching for FLINT headers in common locations..."
	@for path in /usr/include /usr/local/include ~/.local/include $(subst :, ,$(C_INCLUDE_PATH)); do \
		if [ -f "$$path/flint/flint.h" ]; then \
			echo "  FOUND: $$path/flint/flint.h"; \
		fi; \
	done
	@echo "Searching for PML headers in common locations..."
	@for path in /usr/include /usr/local/include ~/.local/include $(subst :, ,$(C_INCLUDE_PATH)); do \
		if [ -f "$$path/pml.h" ]; then \
			echo "  FOUND: $$path/pml.h"; \
		fi; \
	done
	@echo "Searching for nmod_poly_mat_utils.h in common locations..."
	@for path in /usr/include /usr/local/include ~/.local/include $(subst :, ,$(C_INCLUDE_PATH)); do \
		if [ -f "$$path/nmod_poly_mat_utils.h" ]; then \
			echo "  FOUND: $$path/nmod_poly_mat_utils.h"; \
		fi; \
	done
	@echo "Searching for nmod_poly_mat_extra.h in common locations..."
	@for path in /usr/include /usr/local/include ~/.local/include $(subst :, ,$(C_INCLUDE_PATH)); do \
		if [ -f "$$path/nmod_poly_mat_extra.h" ]; then \
			echo "  FOUND: $$path/nmod_poly_mat_extra.h"; \
		fi; \
	done

# Debug library detection with simple shell commands
debug-libs:
	@echo "=== Library Detection Debug ==="
	@echo ""
	@echo "=== System Library Paths ==="
	@echo "LD_LIBRARY_PATH: $$LD_LIBRARY_PATH"
	@echo "Detected paths: $(SYSTEM_LIB_PATHS)"
	@echo ""
	@echo "=== PML Library Search ==="
	@echo "Searching for PML libraries in all system paths..."
	@echo -n "Dynamic libraries (libpml.so*): "
	@found=no; for path in $(SYSTEM_LIB_PATHS); do \
		if [ -n "$$path" ] && [ -d "$$path" ]; then \
			if ls "$$path"/libpml.so* >/dev/null 2>&1; then \
				echo "FOUND"; \
				ls "$$path"/libpml.so* 2>/dev/null | sed 's/^/  /'; \
				found=yes; break; \
			fi; \
		fi; \
	done; if [ "$$found" = "no" ]; then echo "NOT FOUND"; fi
	@echo -n "Static libraries (libpml.a): "
	@found=no; for path in $(SYSTEM_LIB_PATHS); do \
		if [ -n "$$path" ] && [ -f "$$path/libpml.a" ]; then \
			echo "FOUND"; \
			echo "  $$path/libpml.a"; \
			found=yes; break; \
		fi; \
	done; if [ "$$found" = "no" ]; then echo "NOT FOUND"; fi
	@echo ""
	@echo "=== Detection Results ==="
	@echo "PML headers found: $(PML_HEADER_CHECK)"
	@echo "nmod_poly_mat_utils.h found: $(NMOD_POLY_MAT_UTILS_CHECK)"
	@echo "nmod_poly_mat_extra.h found: $(NMOD_POLY_MAT_EXTRA_CHECK)"
	@echo "PML dynamic library found: $(PML_DYNAMIC_LIB_CHECK)"
	@echo "PML static library found: $(PML_STATIC_LIB_CHECK)"
	@echo "PML available (all requirements met): $(PML_AVAILABLE)"
	@echo "Selected PML SO path: $(PML_SO_PATH)"
	@echo "Selected PML A path: $(PML_A_PATH)"

# Debug attack directory structure
debug-attack:
	@echo "=== Attack Directory Debug ==="
	@echo ""
	@echo "=== Attack Directory Analysis ==="
	@echo "Attack directory: $(ATTACK_DIR)"
	@echo -n "Directory exists: "
	@if [ -d "$(ATTACK_DIR)" ]; then \
		echo "YES"; \
		echo "Directory contents:"; \
		find $(ATTACK_DIR) -type f -name "*.c" 2>/dev/null | grep -v ".ipynb_checkpoints" | sed 's/^/  /' || echo "  No C files found"; \
		echo ""; \
		echo "Subdirectories:"; \
		find $(ATTACK_DIR) -type d 2>/dev/null | sed 's/^/  /' || echo "  No subdirectories"; \
		echo ""; \
		echo "All files (including checkpoints):"; \
		find $(ATTACK_DIR) -type f -name "*.c" 2>/dev/null | sed 's/^/  /' || echo "  No C files found"; \
	else \
		echo "NO"; \
	fi
	@echo ""
	@echo "=== Detected Variables ==="
	@echo "ATTACK_C_FILES: $(ATTACK_C_FILES)"
	@echo "ATTACK_EXECUTABLES: $(ATTACK_EXECUTABLES)"
	@echo ""
	@echo "=== Test Compilation Check ==="
	@echo "Current directory dixon library status:"
	@if [ -f "$(DIXON_STATIC_LIB)" ]; then \
		echo "  $(DIXON_STATIC_LIB): EXISTS"; \
	else \
		echo "  $(DIXON_STATIC_LIB): MISSING (run 'make static-lib' first)"; \
	fi
	@if [ -f "$(DIXON_SHARED_LIB)" ]; then \
		echo "  $(DIXON_SHARED_LIB): EXISTS"; \
	else \
		echo "  $(DIXON_SHARED_LIB): MISSING (run 'make dynamic-lib' first)"; \
	fi
	@echo ""
	@echo "=== Compiled Attack Programs Status ==="
	@if [ -n "$(ATTACK_EXECUTABLES)" ]; then \
		for exe in $(ATTACK_EXECUTABLES); do \
			if [ -f "$$exe" ]; then \
				echo "  $$exe: EXISTS"; \
				ls -la "$$exe" | sed 's/^/    /'; \
			else \
				echo "  $$exe: NOT COMPILED"; \
			fi; \
		done; \
	else \
		echo "  No attack executables expected"; \
	fi

debug-structure:
	@echo "=== Local Directory Structure Debug ==="
	@echo ""
	@echo "=== Current Directory ==="
	@echo "PWD: $(shell pwd)"
	@echo "Contents:"
	@ls -la . | sed 's/^/  /'
	@echo ""
	@echo "=== Source Directory ($(SRC_DIR)) ==="
	@echo -n "Directory exists: "
	@if [ -d "$(SRC_DIR)" ]; then \
		echo "YES"; \
		echo "Contents:"; \
		ls -la $(SRC_DIR) | sed 's/^/  /'; \
	else \
		echo "NO"; \
	fi
	@echo ""
	@echo "=== Include Directory ($(INCLUDE_DIR)) ==="
	@echo -n "Directory exists: "
	@if [ -d "$(INCLUDE_DIR)" ]; then \
		echo "YES"; \
		echo "Contents:"; \
		ls -la $(INCLUDE_DIR) | sed 's/^/  /'; \
	else \
		echo "NO"; \
	fi
	@echo ""
	@echo "=== Build Directory ($(BUILD_DIR)) ==="
	@echo -n "Directory exists: "
	@if [ -d "$(BUILD_DIR)" ]; then \
		echo "YES"; \
		echo "Contents:"; \
		ls -la $(BUILD_DIR) | sed 's/^/  /'; \
	else \
		echo "NO (will be created during build)"; \
	fi

# Help
help:
	@echo "Available targets:"
	@echo "  make (default)       - Build libraries first, then dixon with LTO (all sources compiled together)"
	@echo "  make all             - Same as default"
	@echo "  make lto             - Same as default - Build with Link Time Optimization"
	@echo "  make dynamic         - Build dixon with dynamic dixon library"
	@echo "  make static          - Build dixon with static dixon library (dynamic FLINT/PML)"
	@echo "  make static-pml      - Build dixon with static dixon+PML libraries (dynamic FLINT)"
	@echo "  make static-all      - Build dixon with all static libraries (fully static)"
	@echo "  make dynamic-lib     - Build dynamic dixon library only"
	@echo "  make static-lib      - Build static dixon library only"
	@echo "  make attack-programs - Build all C programs in ../Attack directory (using dynamic dixon library)"
	@echo "  make attack-programs-verbose - Build Attack programs with detailed output"
	@echo "  make attack-static   - Build all C programs in ../Attack directory (using static dixon library)"
	@echo "  make clean-attack    - Clean all compiled Attack programs"
	@echo "  make test-paths      - Test library path detection"
	@echo "  make test-attack     - Test Attack directory detection"
	@echo "  make info            - Show build configuration"
	@echo "  make debug-headers   - Debug header file detection (recommended)"
	@echo "  make debug-libs      - Debug external library detection"
	@echo "  make debug-structure - Debug local directory structure"
	@echo "  make debug-attack    - Debug Attack directory structure and C files"
	@echo "  make clean           - Clean all build artifacts (including Attack programs)"
	@echo "  make clean-build     - Clean only build directory"
	@echo "  make help            - Show this help"
	@echo ""
	@echo "Directory structure:"
	@echo "  $(SRC_DIR)/          - Source files (.c)"
	@echo "  $(INCLUDE_DIR)/      - Header files (.h)"
	@echo "  $(BUILD_DIR)/        - Object files (.o) [created during build]"
	@echo "  $(ATTACK_DIR)/       - Attack programs (.c files with main())"
	@echo "  ./               - Executables and libraries"
	@echo ""
	@echo "Attack programs workflow:"
	@echo "  1. Run 'make' to build dixon libraries AND all Attack programs automatically"
	@echo "  2. Or run 'make attack-programs' to build only Attack programs with dynamic library"
	@echo "  3. Or run 'make attack-programs-verbose' for detailed compilation output"
	@echo "  4. Or run 'make attack-static' to build with static dixon library"
	@echo "  5. Each .c file in ../Attack becomes an executable with the same name"
	@echo "  6. Use 'make debug-attack' to check compilation status"
	@echo "  7. .ipynb_checkpoints directories are automatically excluded"
	@echo ""
	@echo "Compilation strategy:"
	@echo "  default - Build libraries first, then compile all sources with LTO for maximum inlining"
	@echo "  dynamic - Traditional library-based compilation using pre-built library"
	@echo "  static  - Static dixon + dynamic FLINT/PML (needs rpath)"
	@echo "  static-pml  - Static dixon+PML + dynamic FLINT (needs rpath for FLINT)"
	@echo "  static-all  - Fully static (no runtime dependencies)"
	@echo ""
	@echo "Library structure:"
	@echo "  Dixon library: $(words $(MATH_SOURCES)) math source files"
	@echo "  Main program: dixon.c links against dixon library OR compiles with all sources"
	@echo "  Attack programs: Each .c in ../Attack links against dixon library + external libraries"
	@echo "  External deps: FLINT (required), PML (optional - auto-detected)"
	@echo ""
	@echo "PML Detection:"
	@echo "  PML support requires ALL of: pml.h, nmod_poly_mat_utils.h, nmod_poly_mat_extra.h"
	@echo "  PLUS at least one library file: libpml.so OR libpml.a"
	@echo "  If any requirement is missing, PML support is disabled automatically"
	@echo "  Use 'make test-paths' and 'make debug-libs' to check availability"

# Aliases for convenience
lto: $(DIXON_TARGET)-lto
dynamic: $(DIXON_TARGET)-dynamic

.PHONY: default all lto dynamic static static-pml static-all dynamic-lib static-lib attack-programs attack-programs-verbose attack-static clean-attack clean clean-build test-paths test-attack info debug-headers debug-libs debug-structure debug-attack help
