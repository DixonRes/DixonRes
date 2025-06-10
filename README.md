# Dixon Resultant Computation in Magma

This program implements Dixon resultant calculation in Magma, supporting multiple computation modes for flexibility and performance optimization.

## Features

- Computes Dixon resultant for polynomial systems
- Multiple computation modes:
  - Magma's built-in determinant calculation
  - Interpolation-based approach
  - Matrix export for external computation
- Support for finite field computations
- Random polynomial generation for testing

## Installation

1. Ensure you have [Magma](http://magma.maths.usyd.edu.au/magma/) installed
2. For mode 2 (external computation), install:
   - [NTL](https://www.libntl.org/) or [FLINT](https://www.flintlib.org/)
   - [PML](https://github.com/vneiger/pml) (Polynomial Matrix Library)

## Usage

### Main Function

```magma
dixon(ps, vars, mode);
```

**Parameters:**
- `ps`: Polynomial system (sequence of polynomials)
- `vars`: Variables to eliminate (sequence of variables)
- `mode`: Computation mode:
  - `0`: Use Magma's built-in determinant calculation
  - `1`: Use interpolation method
  - `2`: Output Dixon matrix for external computation, Only supports well-posed systems (n equations with n variables)

### Example

```magma
load "dixon.m";

// Define polynomial ring
R<x0,x1,x2,y0,y1,y2> := PolynomialRing(GF(65537), 6);

// Define polynomial system
ps := [
    x0 + x1 + x2 - y0,
    x0*x1 + x0*x2 + x1*x2 - y1,
    x0*x1*x2 - y2
];

// Compute Dixon resultant eliminating x0, x1
resultant := dixon(ps, [x0,x1], 0);
```

### Testing

Refer to `test.txt` for examples using `random_polynomial_over_finite_field` to generate test cases or input your own polynomials.

## Computation Modes

1. **Mode 0**: Direct computation using Magma's determinant function
   - Fastest for small matrices
   - Limited by Magma's built-in capabilities

2. **Mode 1**: Interpolation-based determinant calculation
   - Better for larger matrices
   - More memory efficient

3. **Mode 2**: Export Dixon matrix for external computation
   - Use with NTL/FLINT+PML for best performance on large problems
   - Requires additional libraries

## Performance Notes

- For low-degree polynomial systems: Mode 0 recommended
- For high-degree polynomial systems: Mode 1 recommended
- For single-parameter systems: Mode 2 recommended (fastest option for n equations with n variables cases)



## Limitations

- Currently optimized for finite field computations
- Characteristic 0 support may require modifications

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
