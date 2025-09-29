# Dixon Resultant & Polynomial System Solver

A high-performance C implementation for computing Dixon resultants and solving polynomial systems over finite fields using FLINT library.

## Features

- **Dixon Resultant Computation**: Eliminate variables from polynomial systems using the Dixon resultant method
- **Polynomial System Solver**: Solve n×n polynomial systems (n equations in n variables)
- **Dixon with Ideal Reduction**: Compute resultants with triangular ideal reduction
- **Finite Field Support**: 
  - Prime fields F_p (p < 2^63)
  - Extension fields F_{p^k}
  - Automatic field size parsing (e.g., `257`, `2^8`, `256`)
- **Flexible Input**: Command line arguments or file input
- **Silent Mode**: Background computation with file output only
- **Parallel Computation**: OpenMP optimization for performance

## Dependencies

### Required
- **FLINT** (Fast Library for Number Theory) - Headers and library must be accessible
- **GMP** (GNU Multiple Precision Arithmetic Library)
- **MPFR** (Multiple Precision Floating-Point Reliable Library)

### Optional
- **PML** (Polynomial Math Library) - Detected automatically if available

### System Requirements
- GCC compiler with OpenMP support
- Standard C libraries: pthread, math, stdc++

## Windows Version

For Windows users, a native GUI application is available at https://github.com/DixonRes/DixonRes-Windows

### Quick Start
```bash
make              # Build optimized executable + libraries
```

### Build Options
```bash
make                        # Default: Direct compilation (best performance)
make dixon-static-pml       # Static PML, dynamic FLINT
make dixon-static-all       # Fully static executable
make libs                   # Build both static and dynamic libraries
make static-lib             # Build static library only
make dynamic-lib            # Build dynamic library only
```

### Configuration
The Makefile automatically detects FLINT and PML headers using the compiler. If needed, set environment variables:

```bash
export C_INCLUDE_PATH="/path/to/flint/include:/path/to/pml/include"
export LIBRARY_PATH="/path/to/flint/lib:/path/to/pml/lib"
make
```

### Debug Information
```bash
make info              # Show build configuration
make debug-headers     # Debug header detection
make debug-libs        # Debug library detection
make debug-structure   # Debug directory structure
```

## Usage

### Polynomial System Solver (n×n systems)

**Command Line:**
```bash
# Solve 2×2 system
./dixon --solve "x + y - 3, 2*x - y" 257

# Solve 3×3 system
./dixon --solve "x^2 + y^2 + z^2 - 5, x + y + z - 3, x*y*z - 1" 257
```

**File Input:**
```bash
./dixon --solve input.dat
```

File format:
```
257
x + y - 3
2*x - y
```

### Dixon Resultant (Basic)

**Command Line:**
```bash
# 3 equations, eliminate 2 variables (x, y)
./dixon "x + y + z, x*y + y*z + z*x, x*y*z + 1" "x, y" 257
```

**File Input:**
```bash
./dixon input.dat
```

File format:
```
257
x + y + z
x*y + y*z + z*x
x*y*z + 1
x, y
```

**Important**: For Dixon mode, you must eliminate exactly (n-1) variables for n equations.

### Dixon with Ideal Reduction

```bash
./dixon "a1^2 + a2^2 + a3^2 + a4^2 - 10, a4^3 - a1 - a2*a3 - 5" \
        "a4" \
        "a2^3 = 2*a1 + 1, a3^3 = a1*a2 + 3, a4^3 = a1 + a2*a3 + 5" \
        "a1, a2, a3, a4" \
        257
```

### Extension Fields

```bash
# F_256 = F_2^8, 't' is the field generator
./dixon "x + y^2 + t, x*y + t*y + 1" "x" 2^8
```

### Silent Mode

```bash
# No console output except timing
./dixon --silent --solve "x^2 + y^2 - 1, x + y - 1" 7
./dixon --silent input.dat
```

### Field Size Formats

All of these are valid:
- `257` - Prime field
- `256` - Automatically detected as 2^8
- `2^8` - Extension field F_{2^8}
- `3^5` - Extension field F_{3^5}

## Examples

### Example 1: Simple Linear System
```bash
./dixon --solve "x + y - 5, x - y - 1" 7
```

### Example 2: Quadratic System
```bash
./dixon --solve "x^2 + y^2 - 5, x + y - 3" 257
```

### Example 3: Dixon Resultant
```bash
# Symmetric polynomials - eliminate x and y, keep z
./dixon "x + y + z - 6, x*y + y*z + z*x - 11, x*y*z - 6" "x, y" 257
```

### Example 4: Extension Field
```bash
# Binary extension field F_256
./dixon --solve "x^2 + x*t + t^2, x + t" 2^8
```

## Output

Results are saved to files:
- Command line input: `solution.dat`
- File input `example.dat`: `example_solution.dat`

Output includes:
- Field information
- Input polynomials
- Computation time
- Solutions or resultant

## Project Structure

```
.
├── dixon.c              # Main program
├── Makefile             # Build configuration
├── src/                 # Source files
│   ├── dixon_flint.c
│   ├── polynomial_system_solver.c
│   └── ...
├── include/             # Header files
│   ├── dixon_flint.h
│   ├── polynomial_system_solver.h
│   └── ...
└── build/               # Object files (created during build)
```

## Dixon vs Solver Mode

### Dixon Mode
- **Input**: n polynomials + (n-1) variables to eliminate
- **Output**: Resultant polynomial in remaining variables
- **Use case**: Variable elimination, implicitization

### Solver Mode
- **Input**: n polynomials in n variables
- **Output**: All solutions for all variables
- **Use case**: Finding roots, solving complete systems

## Optimization

The default build uses:
- `-O3` - Maximum optimization
- `-flto` - Link-Time Optimization
- `-march=native` - CPU-specific optimizations
- `-fopenmp` - Parallel computation
- Direct compilation for best cross-module inlining

## Troubleshooting

### FLINT not found
```bash
make debug-headers  # Check header detection
export C_INCLUDE_PATH="/path/to/flint/include:$C_INCLUDE_PATH"
```

### Library linking errors
```bash
export LD_LIBRARY_PATH="/path/to/flint/lib:$LD_LIBRARY_PATH"
```

### Compilation errors
```bash
make clean
make info  # Check configuration
make
```

## Testing

```bash
./dixon --test           # Run Dixon tests
./dixon --test-solver    # Run polynomial solver tests
```

## Cleaning

```bash
make clean         # Remove all build artifacts
make clean-build   # Keep executables, remove build directory
```

- Direct compilation (default) provides best performance through full inlining and LTO
- For large systems, parallel computation significantly reduces runtime
- Extension fields are slower than prime fields due to polynomial arithmetic
- Static linking increases startup time but eliminates runtime dependencies

