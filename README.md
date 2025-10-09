# Dixon Resultant & Polynomial System Solver

A high-performance C implementation for computing Dixon resultants and solving polynomial systems over finite fields using FLINT library.

## Implementation Notes for Version 1.0.0

This release (`v1.0.0`) is the version submitted to EUROCRYPT 2026. Please note that the `main` branch is under active development.

### Missing File `dixon_test.h`

**Due to our oversight, this version is missing the `dixon_test.h` file.** This does not affect the core functionality implementation. To compile and use the `./dixon` program normally, you can either download the latest dixon_test.h file from the main branch, or:

1. **Create an empty `dixon_test.h` file** in the `include/` directory containing only:
   ```c
   void dixon_test();
   ```

2. **Clear the `dixon_test.c` file** in the `src/` directory (or implement your own test functions as needed)

This will allow normal compilation and usage of the Dixon program.

### Automatic Algorithm Selection

This implementation automatically selects the optimal resultant computation method based on input configuration:

- **2 polynomials, 1 elimination variable** → FLINT's optimized bivariate resultant
- **3+ polynomials, n-1 elimination variables** → Dixon matrix construction

### Bivariate Resultant (FLINT Built-in)

For the bivariate case, we leverage FLINT's highly optimized resultant algorithms based on **optimizations of the subresultant algorithm**, which provide superior performance compared to Dixon matrix methods for two-polynomial systems.

#### Dixon Matrix Size Reporting (3+ Polynomials)

When using the Dixon method for systems with 3 or more polynomials, you'll observe two distinct matrix size reports:

#### Expected Matrix Size
```
Expected matrix size: 120 x 120
```
This theoretical maximum is computed as the product of `(degree + 1)` for each elimination variable, representing an upper bound before degree filtering. **This is not the Fuss-Catalan bound.**

#### Actual Dixon Matrix Size
```
Found 45 x-monomials and 45 ~x-monomials (after degree filtering)
```
This is the **actual Dixon matrix size** after applying total degree constraints. This size conforms to the Fuss-Catalan bound for multivariate resultants and may be further compressed to a maximal rank submatrix during computation.

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

## Windows GUI

For Windows users, a native GUI application is available at https://github.com/DixonRes/DixonRes-Windows


## Usage

### Polynomial System Solver (n×n systems)

**Command Line:**
```bash
# Solve 2×2 system
./dixon --solve "x + y - 3, 2*x - y" 257

# Solve 3×3 system
./dixon --solve "x^2 + y^2 + z^2 - 3, x + y + z - 3, x*y*z - 1" 257
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

