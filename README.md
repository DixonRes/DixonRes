# Dixon Resultant Implementation for Finite Fields

A high-performance implementation of the Dixon resultant method for polynomial system solving over finite fields (both prime fields and field extensions). This implementation includes optimizations for triangular ideal reduction and is designed for cryptanalysis applications, particularly for algebraic attacks on cryptographic systems like Rescue.

## Overview

The Dixon resultant is a powerful algebraic technique for eliminating variables from systems of polynomial equations. This implementation provides:

- Support for both prime fields F_p and extension fields F_{p^k}
- Triangular ideal reduction for improved efficiency
- Multiple determinant computation methods (recursive, Kronecker substitution, interpolation)
- Specialized attack framework for the Rescue hash function
- Command-line interface and file I/O support

## Features

### Core Functionality

1. **Basic Dixon Resultant**: Eliminates variables from polynomial systems
2. **Dixon with Ideal Reduction**: Uses triangular ideal structure to reduce polynomial degrees during computation
3. **Multiple Field Support**: 
   - Prime fields (e.g., F_257)
   - Binary extension fields (e.g., F_{2^8})
   - General extension fields (e.g., F_{p^k})

### Optimizations

- Parallel computation using OpenMP
- Batch polynomial operations for large systems
- Optimized matrix operations with early termination
- Memory-efficient sparse polynomial representations
- Adaptive algorithm selection based on problem size

## Building

### Prerequisites

- GCC compiler with C11 support
- FLINT library (2.9.0 or later)
- GMP library
- MPFR library
- OpenMP support
- pthreads

### Compilation

```bash
make          # Build the main Dixon tool
make all      # Build both Dixon tool and Rescue attack
make clean    # Clean build artifacts
```

Or compile manually:

```bash
gcc -O3 -march=native -fopenmp -o dixon dixon.c -lflint -lmpfr -lgmp -lpthread
gcc -O3 -march=native -fopenmp -o rescue rescue_attack_dixon.c -lflint -lmpfr -lgmp -lpthread
```

## Usage

### Command Line Interface

#### Basic Dixon Resultant

```bash
./dixon "polynomials" "eliminate_vars" field_size
```

Example:
```bash
./dixon "x + y + z, x*y + y*z + z*x, x*y*z + 1" "x, y" 257
```

This eliminates variables x and y, keeping z in the resultant over F_257.

#### Dixon with Ideal Reduction

```bash
./dixon "polynomials" "eliminate_vars" "ideal_generators" "all_variables" field_size
```

#### File Input

Create an input file with the following format:
```
field_size
polynomial_1
polynomial_2
...
polynomial_n
variables_to_eliminate
```

Then run:
```bash
./dixon input.dat
```

The result will be saved to `input_solution.dat`.

#### Silent Mode

For batch processing without console output:
```bash
./dixon --silent input.dat
```

### Field Specification

Fields can be specified in multiple formats:
- Prime field: `257` or `997`
- Prime power: `256` (automatically detected as 2^8)
- Explicit format: `2^8`, `3^5`, `7^2`

### Rescue Attack Example

The specialized Rescue attack tool demonstrates using Dixon resultants for cryptanalysis:

```bash
./rescue
```

This builds the Rescue polynomial system and performs iterative variable elimination.

## File Formats

### Input File Example (example.dat)

```
65537
31511*x0^9 + 5685*x0^8*x1 + 24547*x0^8 - 5355*x0^7*x1 - 3209*x0^7*x2^2 + ...
-2596*x0^9 - 24616*x0^8*x1^2 - 5730*x0^8*x2^2 - 691*x0^7*x1*x2^2 + ...
-28950*x0^10 + 4882*x0^9*x1 + 15386*x0^9 + 27295*x0^8*x1 + ...
x1, x2
```

### Output Format

The solution file contains:
- Field information
- Computation mode (basic or with ideal reduction)
- Input polynomials
- Variables eliminated
- Computation time
- Final resultant polynomial

## Implementation Details

### Key Components

1. **dixon_flint.h**: Core Dixon resultant implementation
2. **dixon_with_ideal_reduction.h**: Triangular ideal reduction algorithms
3. **fq_unified_interface.h**: Unified interface for different field types
4. **fq_multivariate_interpolation.h**: Multivariate polynomial interpolation
5. **fq_mpoly_mat_det.h**: Matrix determinant algorithms

### Algorithm Selection

The implementation automatically selects appropriate algorithms based on:
- Field characteristics (prime vs extension)
- Number of variables
- Number of parameters
- Polynomial density

### Performance Considerations

- For systems with many variables (>5), use ideal reduction when possible
- Binary extension fields (F_{2^n}) use optimized GF2 arithmetic
- Large prime fields benefit from Kronecker substitution
- Dense systems may require significant memory (>16GB for large problems)

## Testing

Run the built-in test suite:
```bash
./dixon --test
```

This executes various test cases including:
- Prime field examples
- Extension field examples
- Triangular ideal reduction tests
- Performance benchmarks

## Limitations

- Maximum practical system size depends on available memory
- Very high degree polynomials (>100) may require extended computation time
- The number of variables to eliminate must equal (number of equations - 1)

## Environment Variables

- `DIXON_DET_METHOD`: Override determinant computation method (0=recursive, 1=Kronecker, 2=interpolation)
- `OMP_NUM_THREADS`: Control parallel thread count

## Troubleshooting

### Common Issues

1. **Memory errors**: Reduce polynomial system size or increase system RAM
2. **Slow computation**: Enable ideal reduction or adjust OMP_NUM_THREADS
3. **Zero resultant**: Check that the system has a non-trivial solution

### Debug Mode

Compile with debug flags:
```bash
gcc -g -O0 -DDEBUG_OUTPUT_D=1 -o dixon_debug dixon.c -lflint -lmpfr -lgmp
```

## Academic References

This implementation is based on:
- Dixon, A.L. "The eliminant of three quantics in two independent variables" (1908)
- Kapur, D., Saxena, T., Yang, L. "Algebraic and geometric reasoning using Dixon resultants" (1994)

## License

This project is provided for academic and research purposes. Please cite appropriately if used in publications.

## Contributing

For bug reports or feature requests, please include:
- System specifications
- Input polynomials causing issues
- Expected vs actual output
- Compilation flags used
```

