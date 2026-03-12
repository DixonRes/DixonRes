# DixonRes: Dixon Resultant & Polynomial System Solver
A C implementation for computing Dixon resultants and solving polynomial systems over finite fields, based on the FLINT and PML library.

**Note:** This repository is under active development. The version submitted to Crypto 2026 is archived and available at:  
- <https://github.com/DixonRes/DixonRes/releases/tag/0.0.1>
  
Any updates after the submission deadline are not included in the above version.

## Features
- Dixon resultant computation for variable elimination
- Polynomial system solver for n×n systems
- Dixon with triangular ideal reduction
- Finite fields:
  - Prime fields F_p (p < 2^63): Implemented with FLINT modular arithmetic, optionally accelerated by PML.
  - Extension fields F_{p^k}: Further optimized for binary fields F_{2^n} with n in {8, 16, 32, 64, 128}.
- Complexity analysis — estimates Dixon matrix size, Bezout degree bound, and operation count before computing
- Command line input or file input. Automatic output to solution files

---

## Dependencies
- **FLINT** (recommended version: 3.4.0)  
  <https://github.com/flintlib/flint>

Optional:
- **PML** (used automatically if available)  
  <https://github.com/vneiger/pml>

---

## Build
```bash
./configure
make
make check                         # optional
make install                       # optional
```
For more options, run `./configure --help` or `make help`.

---

## Usage

### Dixon Resultant (Basic)
```bash
./dixon "polynomials" "eliminate_vars" field_size
```
Example:
```bash
./dixon "x+y+z, x*y+y*z+z*x, x*y*z+1" "x,y" 257
```

---

### Polynomial System Solver (n equations in n variables)
```bash
./dixon --solve "polynomials" field_size
```
Example:
```bash
./dixon --solve "x^2 + y^2 + z^2 - 6, x + y + z - 4, x*y*z - x - 1" 257
```

---

### Complexity Analysis
Estimates the difficulty of a Dixon resultant computation **without** performing it.
Reports equation count, variable count, degree sequence, Dixon matrix size
(via Hessenberg recurrence), Bezout degree bound, and complexity in bits.

```bash
./dixon --comp "polynomials" "eliminate_vars" field_size
./dixon -c     "polynomials" "eliminate_vars" field_size
```

Example:
```bash
./dixon --comp "x^3+y^3+z^3, x^2*y+y^2*z+z^2*x, x+y+z-1" "x,y" 257
```

**Custom omega** — set the matrix-multiplication exponent used in the complexity formula
(default: 2.3):
```bash
./dixon --comp --omega 2.373 "polynomials" "eliminate_vars" field_size
./dixon -c -w 2.0            "polynomials" "eliminate_vars" field_size
```

File input uses the same format as Dixon mode:
```bash
./dixon --comp example.dat          # output: example_comp.dat
```

---

### Extension Fields
```bash
./dixon "x + y^2 + t, x*y + t*y + 1" "x" 2^8
```
The default settings use `t` as the extension field generator and FLINT's built-in field polynomial.
```bash
./dixon --solve "x^2 + t*y, x*y + t^2" "2^8: t^8 + t^4 + t^3 + t + 1"
```
(with AES custom polynomial for F_256)

---

### Dixon with Ideal Reduction
```bash
./dixon "polynomials" "eliminate_vars" "ideal_generators" field_size
```
Example:
```bash
./dixon "a1^2 + a2^2 + a3^2 + a4^2 - 10, a4^3 - a1 - a2*a3 - 5" \
        "a4" \
        "a2^3 = 2*a1 + 1, a3^3 = a1*a2 + 3, a4^3 = a1 + a2*a3 + 5" \
        257
```

---

### Silent Mode
```bash
./dixon --silent [--solve|--comp|-c] <arguments>
```
No console output is produced; the solution/report file is still generated.

---

## Random Mode

Generate random polynomial systems with specified degrees for testing and benchmarking.

### Basic Usage
```bash
./dixon --random "[d1,d2,...,dn]" field_size
./dixon -r       "[d]*n"          field_size
```

- `[d1,d2,...,dn]`: degree list (comma-separated) for n polynomials
- `[d]*n`: all n polynomials have same degree d
- `field_size`: field size (prime or extension)

### Combine with Compute Flags
```bash
# Random + Dixon elimination
./dixon -r --solve "[d1,...,dn]" field_size

# Random + complexity analysis
./dixon -r --comp  "[d]*n" field_size
./dixon -r -c --omega 2.373 "[4]*5" 257   # custom omega

# Random + Dixon with ideal reduction
./dixon -r "[d1,d2,d3]" "ideal_generators" field_size
```

### Examples
```bash
# 3 polynomials (deg 3,3,2) in GF(257)
./dixon --random "[3,3,2]" 257

# Solve 3 quadratic system in GF(257)
./dixon -r --solve "[2]*3" 257

# Complexity analysis of 4 quartic polynomials
./dixon -r --comp --omega 2.373 "[4]*4" 257

# GF(2^8) with degrees 3 and 2
./dixon -r "[3,2]" 2^8
```

## File Input Format

### Dixon Mode / Complexity Mode (multiline)
```
Line 1 : field size
Line 2+: polynomials (comma-separated or multiline)
Last   : variables to ELIMINATE (comma-separated)
         (#eliminate = #equations - 1)
```
Example:
```bash
./dixon       example.dat
./dixon --comp example.dat
```

### Polynomial Solver Mode (multiline)
```
Line 1 : field size
Line 2+: polynomials
         (n equations in n variables)
```

---

## Output

| Mode | Command-line input | File input `example.dat` |
|---|---|---|
| Dixon / Solver | `solution_YYYYMMDD_HHMMSS.dat` | `example_solution.dat` |
| Complexity | `comp_YYYYMMDD_HHMMSS.dat` | `example_comp.dat` |

Each output file contains field information, input polynomials, computation time,
and the resultant, solutions, or complexity report.

### Complexity report contents
- Equation count, variable list, elimination variable list, remaining variables
- Degree sequence of input polynomials
- Bezout bound (product of degrees)
- Dixon matrix size (Hessenberg recurrence)
- Resultant degree estimate
- Complexity in log₂ bits (with the omega value used)

---

## Notes
- All computation modes generate a solution/report file by default
- Extension fields are slower than prime fields due to polynomial arithmetic
- The optional PML library only accelerates well-determined systems over prime fields
- Complexity analysis does not run any polynomial arithmetic; it parses only

---

## License
DixonRes is distributed under the GNU General Public License version 2.0 (GPL-2.0-or-later). See the file COPYING.
