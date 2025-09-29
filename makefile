# 显式设置默认目标 - 必须在所有变量定义和 $(info) 之前
.DEFAULT_GOAL := default

# Compiler
CC = gcc

# Compiler flags
CFLAGS = -O3 -march=native -fopenmp -fPIC
# -fsanitize=address -static-libasan

# Linker flags
LDFLAGS = -fopenmp

# Library paths (corrected)
FLINT_LIB_PATH = /home/suohaohai02/mylinks
PML_LIB_PATH = /home/suohaohai02/mylinks

# Include paths - 使用环境变量 C_INCLUDE_PATH 中的路径
# 从环境变量中提取路径，如果未设置则使用默认值
FLINT_INCLUDE_PATH ?= /home/suohaohai02/apps/pml/flint-extras/include
PML_INCLUDE_PATH ?= /home/suohaohai02/apps/pml/flint-extras/include

# 检查是否设置了 C_INCLUDE_PATH 环境变量
ifdef C_INCLUDE_PATH
$(info Using C_INCLUDE_PATH: $(C_INCLUDE_PATH))
endif

# Local directories
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build

# 更智能的头文件检查 - 使用编译器来查找头文件
FLINT_HEADER_CHECK := $(shell echo '\#include <flint/flint.h>' | $(CC) -E $(INCLUDE_FLAGS) -x c - >/dev/null 2>&1 && echo yes || echo no)
PML_HEADER_CHECK := $(shell echo '\#include <pml.h>' | $(CC) -E $(INCLUDE_FLAGS) -x c - >/dev/null 2>&1 && echo yes || echo no)

# 旧的目录检查（作为后备）
FLINT_DIR_EXISTS := $(shell if [ -d "$(FLINT_INCLUDE_PATH)" ]; then echo yes; else echo no; fi)
PML_DIR_EXISTS := $(shell if [ -d "$(PML_INCLUDE_PATH)" ]; then echo yes; else echo no; fi)

# 验证 FLINT（必需）- 使用头文件检查而不是目录检查
ifneq ($(FLINT_HEADER_CHECK),yes)
$(error FLINT headers not found! Compiler cannot find flint/flint.h. Check C_INCLUDE_PATH or install FLINT development packages)
endif

# System libraries
SYSTEM_LIBS = -lmpfr -lgmp -lm -lpthread -lstdc++

# Include flags (including local include directory)
# 使用环境变量 C_INCLUDE_PATH，如果设置的话可以不显式指定路径
# 因为 GCC 会自动搜索 C_INCLUDE_PATH 中的目录
INCLUDE_FLAGS = -I./$(INCLUDE_DIR)
# 如果需要显式指定，可以取消注释下面的行
# INCLUDE_FLAGS += -I$(FLINT_INCLUDE_PATH)
# ifeq ($(PML_DIR_EXISTS),yes)
# INCLUDE_FLAGS += -I$(PML_INCLUDE_PATH)
# endif

# Library flags
FLINT_FLAGS = -DHAVE_FLINT
PML_FLAGS = 
# 使用头文件检查而不是目录检查
ifeq ($(PML_HEADER_CHECK),yes)
PML_FLAGS = -DHAVE_PML
endif

# Combined CFLAGS
ALL_CFLAGS = $(CFLAGS) $(INCLUDE_FLAGS) $(FLINT_FLAGS) $(PML_FLAGS)

# FLINT library linking
FLINT_LIBS = -L$(FLINT_LIB_PATH) -lflint
FLINT_STATIC_LIBS = $(FLINT_LIB_PATH)/libflint.a

# PML library linking (optional)
PML_LIBS = 
PML_STATIC_LIBS = 
# 使用头文件检查而不是目录检查
ifeq ($(PML_HEADER_CHECK),yes)
PML_LIBS = -L$(PML_LIB_PATH) -lpml
PML_STATIC_LIBS = $(PML_LIB_PATH)/libpml.a
endif

# Combined external libraries
EXTERNAL_LIBS = $(FLINT_LIBS) $(PML_LIBS) $(SYSTEM_LIBS)
EXTERNAL_STATIC_PML_LIBS = $(FLINT_LIBS) $(PML_STATIC_LIBS) $(SYSTEM_LIBS)
EXTERNAL_STATIC_ALL_LIBS = $(FLINT_STATIC_LIBS) $(PML_STATIC_LIBS) $(SYSTEM_LIBS)

# Runtime library path (rpath) flags
RPATH_FLAGS = -Wl,-rpath,$(FLINT_LIB_PATH) -Wl,-rpath,.
ifeq ($(PML_HEADER_CHECK),yes)
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

# Library names (in current directory)
DIXON_STATIC_LIB = libdixon.a
DIXON_SHARED_LIB = libdixon.so

# Output executable (in current directory)
DIXON_TARGET = dixon

# Create build directory
$(BUILD_DIR):
	@echo "Creating build directory..."
	mkdir -p $(BUILD_DIR)

# 默认目标 - 使用正确的依赖关系，完全模仿 static 的模式
default: $(DIXON_TARGET)-dynamic
	@echo "Built dixon with dynamic library (default)"

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

# Object file compilation (src/*.c -> build/*.o)
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	@echo "Compiling $<..."
	$(CC) $(ALL_CFLAGS) -c -o $@ $<

# Clean
clean:
	rm -f $(DIXON_TARGET) $(DIXON_STATIC_LIB) $(DIXON_SHARED_LIB)
	rm -rf $(BUILD_DIR)
	@echo "Cleaned all build artifacts"

# Clean only build directory (keep executables and libraries)
clean-build:
	rm -rf $(BUILD_DIR)
	@echo "Cleaned build directory"

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
	@echo "=== Library Status ==="
	@echo "FLINT headers found: $(FLINT_HEADER_CHECK)"
	@echo "FLINT directory exists: $(FLINT_DIR_EXISTS) at $(FLINT_INCLUDE_PATH)"
	@echo "PML headers found: $(PML_HEADER_CHECK)"
	@echo "PML directory exists: $(PML_DIR_EXISTS) at $(PML_INCLUDE_PATH)"
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
	@echo "=== Combined Settings ==="
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
	@if echo '\#include <flint/flint.h>' | $(CC) -E $(INCLUDE_FLAGS) -x c - >/dev/null 2>&1; then \
		echo "FOUND"; \
	else \
		echo "NOT FOUND"; \
	fi
	@echo -n "PML headers (pml.h): "
	@if echo '\#include <pml.h>' | $(CC) -E $(INCLUDE_FLAGS) -x c - >/dev/null 2>&1; then \
		echo "FOUND"; \
	else \
		echo "NOT FOUND"; \
	fi
	@echo ""
	@echo "=== Manual Path Search ==="
	@echo "Searching for FLINT headers in common locations..."
	@for path in /usr/include /usr/local/include ~/.local/include $(subst :, ,$(C_INCLUDE_PATH)); do \
		if [ -f "$path/flint/flint.h" ]; then \
			echo "  FOUND: $path/flint/flint.h"; \
		fi; \
	done
	@echo "Searching for PML headers in common locations..."
	@for path in /usr/include /usr/local/include ~/.local/include $(subst :, ,$(C_INCLUDE_PATH)); do \
		if [ -f "$path/pml.h" ]; then \
			echo "  FOUND: $path/pml.h"; \
		fi; \
	done

# Debug library detection
debug-libs:
	@echo "=== Library Detection Debug ==="
	@echo ""
	@echo "=== FLINT ==="
	@echo "Include path: $(FLINT_INCLUDE_PATH)"
	@echo -n "Directory exists: "
	@if [ -d "$(FLINT_INCLUDE_PATH)" ]; then \
		echo "YES"; \
		echo "Sample files:"; \
		ls "$(FLINT_INCLUDE_PATH)" | head -5 | sed 's/^/  /'; \
	else \
		echo "NO"; \
	fi
	@echo "Library path: $(FLINT_LIB_PATH)"
	@echo -n "libflint.so: "
	@if ls $(FLINT_LIB_PATH)/libflint.so* 2>/dev/null | head -1 >/dev/null; then \
		echo "FOUND"; \
		ls $(FLINT_LIB_PATH)/libflint.so* | sed 's/^/  /'; \
	else \
		echo "NOT FOUND"; \
	fi
	@echo -n "libflint.a: "
	@if [ -f "$(FLINT_LIB_PATH)/libflint.a" ]; then \
		echo "FOUND"; \
		ls -la $(FLINT_LIB_PATH)/libflint.a | sed 's/^/  /'; \
	else \
		echo "NOT FOUND"; \
	fi
	@echo ""
	@echo "=== PML ==="
	@echo "Include path: $(PML_INCLUDE_PATH)"
	@echo -n "Directory exists: "
	@if [ -d "$(PML_INCLUDE_PATH)" ]; then \
		echo "YES"; \
		echo "Sample files:"; \
		ls "$(PML_INCLUDE_PATH)" | head -5 | sed 's/^/  /'; \
	else \
		echo "NO"; \
	fi
	@echo "Library path: $(PML_LIB_PATH)"
	@echo -n "libpml.so: "
	@if ls $(PML_LIB_PATH)/libpml.so* 2>/dev/null | head -1 >/dev/null; then \
		echo "FOUND"; \
		ls $(PML_LIB_PATH)/libpml.so* | sed 's/^/  /'; \
	else \
		echo "NOT FOUND"; \
	fi
	@echo -n "libpml.a: "
	@if [ -f "$(PML_LIB_PATH)/libpml.a" ]; then \
		echo "FOUND"; \
		ls -la $(PML_LIB_PATH)/libpml.a | sed 's/^/  /'; \
	else \
		echo "NOT FOUND"; \
	fi

# Debug local directory structure
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
	@echo "  make (default)       - Build dixon with dynamic dixon library"
	@echo "  make static          - Build dixon with static dixon library (dynamic FLINT/PML)"
	@echo "  make static-pml      - Build dixon with static dixon+PML libraries (dynamic FLINT)"
	@echo "  make static-all      - Build dixon with all static libraries (fully static)"
	@echo "  make dynamic-lib     - Build dynamic dixon library only"
	@echo "  make static-lib      - Build static dixon library only"
	@echo "  make info            - Show build configuration"
	@echo "  make debug-headers   - Debug header file detection (recommended)"
	@echo "  make debug-libs      - Debug external library detection"
	@echo "  make debug-structure - Debug local directory structure"
	@echo "  make clean           - Clean all build artifacts"
	@echo "  make clean-build     - Clean only build directory"
	@echo "  make help            - Show this help"
	@echo ""
	@echo "Directory structure:"
	@echo "  $(SRC_DIR)/          - Source files (.c)"
	@echo "  $(INCLUDE_DIR)/      - Header files (.h)"
	@echo "  $(BUILD_DIR)/        - Object files (.o) [created during build]"
	@echo "  ./               - Executables and libraries"
	@echo ""
	@echo "Static compilation options:"
	@echo "  static      - Static dixon + dynamic FLINT/PML (needs rpath)"
	@echo "  static-pml  - Static dixon+PML + dynamic FLINT (needs rpath for FLINT)"
	@echo "  static-all  - Fully static (no runtime dependencies)"
	@echo ""
	@echo "Library structure:"
	@echo "  Dixon library: $(words $(MATH_SOURCES)) math source files"
	@echo "  Main program: dixon.c links against dixon library"
	@echo "  External deps: FLINT (required), PML (optional)"

.PHONY: default static static-pml static-all dynamic-lib static-lib clean clean-build info debug-headers debug-libs debug-structure help