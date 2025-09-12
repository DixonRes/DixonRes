# Compiler
CC = gcc

# Compiler flags
CFLAGS = -O3 -march=native -fopenmp
# -fsanitize=address -static-libasan

# Linker flags  
LDFLAGS = -fopenmp

# Base libraries (always required)
BASE_LIBS = -lflint -lmpfr -lgmp -lm -lpthread -lstdc++

# Custom library path
LIBPATH = -L/home/suohaohai02/mylinks

# Define library search paths
LIB_SEARCH_PATHS = /home/suohaohai02/mylinks /usr/local/lib /usr/lib /usr/lib64 /lib /lib64 . $(LD_LIBRARY_PATH)

# Function to check if PML exists in any of the paths
define check_pml
$(shell for path in $(LIB_SEARCH_PATHS); do \
    if [ -d "$$path" ]; then \
        if [ -f "$$path/libpml.so" ] || [ -f "$$path/libpml.a" ] || [ -f "$$path/libpml.dylib" ]; then \
            echo "$$path"; \
            exit 0; \
        fi; \
    fi; \
done; \
echo "")
endef

# Find PML library
PML_PATH := $(call check_pml)

# Check if PML can be linked
ifneq ($(PML_PATH),)
    # Found PML library file, just use it
    PML_LIBS = -L$(PML_PATH) -lpml
    PML_FLAGS = -DHAVE_PML
    $(info === PML library found in $(PML_PATH) - prime field optimizations ENABLED ===)
else
    PML_LIBS = 
    PML_FLAGS = 
    $(warning === PML library NOT found in any of the following locations: ===)
    $(warning === $(LIB_SEARCH_PATHS) ===)
    $(info === Continuing without PML optimizations... ===)
endif

# Combined libraries
LIBS = $(BASE_LIBS) $(PML_LIBS)

# Combined CFLAGS
ALL_CFLAGS = $(CFLAGS) $(PML_FLAGS)

# Source files
DIXON_SRC = dixon.c
RESCUE_SRC = rescue_attack_dixon.c

# Output executables
DIXON_TARGET = dixon
RESCUE_TARGET = rescue

# Default target (run with just 'make')
default: $(DIXON_TARGET)

# Build both targets
all: $(DIXON_TARGET) $(RESCUE_TARGET)

$(DIXON_TARGET): $(DIXON_SRC)
	@echo "Building $(DIXON_TARGET)..."
	$(CC) $(ALL_CFLAGS) -o $@ $^ $(LIBPATH) $(LIBS) $(LDFLAGS)
	@echo "Build complete: $(DIXON_TARGET)"

$(RESCUE_TARGET): $(RESCUE_SRC)
	@echo "Building $(RESCUE_TARGET)..."
	$(CC) $(ALL_CFLAGS) -o $@ $^ $(LIBPATH) $(LIBS) $(LDFLAGS)
	@echo "Build complete: $(RESCUE_TARGET)"

clean:
	rm -f $(DIXON_TARGET) $(RESCUE_TARGET)

# Print configuration
info:
	@echo "=== Build Configuration ==="
	@echo "CC: $(CC)"
	@echo "CFLAGS: $(ALL_CFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
	@echo "LIBPATH: $(LIBPATH)"
	@echo "LIBS: $(LIBS)"
	@echo ""
	@echo "=== PML Status ==="
ifneq ($(PML_PATH),)
	@echo "PML library: FOUND in $(PML_PATH)"
	@echo "Prime field optimizations: ENABLED"
else
	@echo "PML library: NOT FOUND"
	@echo "Searched paths: $(LIB_SEARCH_PATHS)"
	@echo "Prime field optimizations: DISABLED"
endif

# Debug target to help diagnose PML issues
debug-pml:
	@echo "=== PML Library Debug Information ==="
	@echo "Searching for PML in the following paths:"
	@for path in $(LIB_SEARCH_PATHS); do \
	    echo -n "  $path: "; \
	    if [ -d "$path" ]; then \
	        if [ -f "$path/libpml.so" ]; then \
	            echo "libpml.so FOUND"; \
	            ls -la "$path/libpml.so" 2>/dev/null || true; \
	        elif [ -f "$path/libpml.a" ]; then \
	            echo "libpml.a FOUND"; \
	            ls -la "$path/libpml.a" 2>/dev/null || true; \
	        elif [ -f "$path/libpml.dylib" ]; then \
	            echo "libpml.dylib FOUND"; \
	            ls -la "$path/libpml.dylib" 2>/dev/null || true; \
	        else \
	            echo "NOT FOUND"; \
	        fi; \
	    else \
	        echo "DIRECTORY DOES NOT EXIST"; \
	    fi; \
	done
	@echo ""
	@echo "PML_PATH detected: $(PML_PATH)"
	@echo ""
	@if [ -n "$(PML_PATH)" ]; then \
	    echo "Checking PML library dependencies:"; \
	    if [ -f "$(PML_PATH)/libpml.so" ]; then \
	        echo "Running ldd on $(PML_PATH)/libpml.so:"; \
	        ldd "$(PML_PATH)/libpml.so" 2>/dev/null || echo "ldd failed or not available"; \
	    fi; \
	    echo ""; \
	    echo "Testing compilation with PML:"; \
	    echo "#include <stdio.h>" > /tmp/test_pml.c; \
	    echo "int main() { printf(\"PML test\\n\"); return 0; }" >> /tmp/test_pml.c; \
	    echo "Command: $(CC) /tmp/test_pml.c -o /tmp/test_pml $(LIBPATH) -L$(PML_PATH) -lpml $(BASE_LIBS)"; \
	    if $(CC) /tmp/test_pml.c -o /tmp/test_pml $(LIBPATH) -L$(PML_PATH) -lpml $(BASE_LIBS) 2>&1; then \
	        echo "Compilation: SUCCESS"; \
	        rm -f /tmp/test_pml; \
	    else \
	        echo "Compilation: FAILED"; \
	    fi; \
	    rm -f /tmp/test_pml.c; \
	fi

.PHONY: default all clean info debug-pml