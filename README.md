# DixonRes: Dixon Resultant & Polynomial System Solver
A C implementation for computing Dixon resultants and solving polynomial systems over finite fields and the rationals ℚ, based on the FLINT and PML libraries.

Website: <https://dixonres.github.io>

**Note:** This repository is under active development. The version submitted to Crypto 2026 is archived and available at:  
 - Release page: <https://github.com/DixonRes/DixonRes/releases/tag/0.0.1>
 - Code tree: <https://github.com/DixonRes/DixonRes/tree/0.0.1>

## Features
- Dixon resultant computation for variable elimination
- Polynomial system solver for n×n systems
- Dixon with triangular ideal reduction
- Finite fields:
  - Prime fields F_p (any size): Implemented with FLINT modular arithmetic, optionally accelerated by PML.
  - Extension fields F_{p^k}: Further optimized for binary fields F_{2^n} with n in {8, 16, 32, 64, 128}.
- Rational field ℚ: Rational reconstruction via multi-prime CRT. Set field_size = 0 to enable.
- Complexity analysis — estimates Dixon matrix size, Bezout degree bound, and operation count before computing
- Command line input or file input. Automatic output to solution files

---

## Dependencies
- **FLINT** (recommended version: 3.4.0)  
  <https://github.com/flintlib/flint>

```bash
git clone https://github.com/flintlib/flint.git && cd flint
./bootstrap.sh
./configure 
make
make install
```
  
Optional:
- **PML** (used automatically if available)  
  <https://github.com/vneiger/pml>
  
```bash
git clone https://github.com/vneiger/pml.git && cd pml/flint-extras
./bootstrap.sh
./configure
make
make install
```

---

## Build
```bash
git clone https://github.com/DixonRes/DixonRes.git && cd DixonRes
./configure
make
make check                         # optional
make install                       # optional
```
For more options, run `./configure --help` or `make help`.
We also provide a Windows GUI at [DixonRes-Windows](https://github.com/DixonRes/DixonRes-Windows) or [DixonRes-Cross](https://github.com/DixonRes/DixonRes-Cross).

---

## Usage

### Dixon Resultant (Basic)
```bash
./dixon "polynomials" "eliminate_vars" field_size
```
Examples:
```bash
./dixon "x+y+z, x*y+y*z+z*x, x*y*z+1" "x,y" 257
./dixon "x^2+y^2+z^2-1, x^2+y^2-2*z^2, x+y+z" "x,y" 0
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

## File Input Format

### Dixon Mode (multiline)
```
Line 1 : field size (prime or p^k; use 0 for Q; generator defaults to 't')
Line 2+: polynomials (comma-separated, may span multiple lines)
Last   : variables to ELIMINATE (comma-separated)
         (#eliminate = #equations - 1)
```
Example:
```bash
# example.dr
0
x0^3+x1^3+x2^3, x0*x1+x1*x2+x2*x1, x1*x2*x0+1
x0,x1
```
Run:
```bash
./dixon example.dr
```

### Polynomial Solver Mode (multiline)
```
Line 1 : field size
Line 2+: polynomials
         (n equations in n variables)
```
Example:
```bash
# example_solve.dr
257
x^2+y^2+z^2-6, x+y+z-4, x*y*z-x-1
```
Run:
```bash
./dixon --solve example_solve.dr
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

Examples:
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
./dixon --ideal "ideal_generators" "polynomials" "eliminate_vars" field_size
```
Example:
```bash
./dixon --ideal "a2^3=2*a1+1, a3^3=a1*a2+3" "a1^2+a2^2+a3^2-10, a3^3-a1*a2-3" "a3" 257
```

### Field-equation reduction mode 
After each multiplication, reduces x^q -> x for every variable.
```bash
./dixon --field-equation "polynomials" "eliminate_vars" field_size
./dixon --field-eqution -r "[d1,d2,...,dn]" field_size
```
Example:
```bash
./dixon --field-eqution "x0*x2+x1, x0*x1*x2+x2+1, x1*x2+x0+1" "x0,x1" 2
./dixon --field-eqution -r [3]*5 2
```

---

### Silent Mode
```bash
./dixon --silent [--solve|--comp|-c] <arguments>
```
No console output is produced; the solution/report file is still generated.

### Method Selection
```bash
./dixon --method <method> --threads <num> <args>
```
Available methods: 1.Recursive; 2.Kronecker; 3.Interpolation

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
- `field_size`: field size (prime or extension); use `0` for Q

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
```

---

## SageMath Interface

`dixon_sage_interface.sage` is a thin wrapper that lets you call DixonRes directly from a SageMath session.

### Quick start

```python
load("dixon_sage_interface.sage")
set_dixon_path("./dixon")   # set once per session

R.<x, y, z> = GF(257)[]
F = [x + y + z - 3, x*y + y*z + z*x - 3, x*y*z - 1]

res  = DixonResultant(F, [x, y])   # Dixon resultant, eliminating x and y
sols = DixonSolve(F)               # enumerate all solutions
info = DixonComplexity(F, [x, y])  # complexity estimate (no arithmetic)

# iterative elimination: output is a plain string, feed into the next call
res1 = DixonResultant([x+y+z, x*y+y*z+z*x+1], [x])
res2 = DixonResultant([res1, y*z-1], [y])
```

### API reference

| Function | Description | Returns |
|---|---|---|
| `DixonResultant(F, elim_vars, ...)` | Dixon resultant, eliminating the specified variables. | String or `None` |
| `DixonSolve(F, ...)` | Solve an n×n system, enumerate all solutions. | List of `{var: val}` dicts; `[]`; or `"infinite"` |
| `DixonComplexity(F, elim_vars, ...)` | Estimate complexity without any polynomial arithmetic. | Dict with `complexity_log2`, `bezout_bound`, `matrix_size`, … |
| `DixonIdeal(F, ideal_gens, elim_vars, ...)` | Dixon resultant with triangular ideal reduction. `ideal_gens`: list of strings like `"a^3=2*b+1"`. | String or `None` |
| `set_dixon_path(p)` / `get_dixon_path()` | Set / get the default path to the `dixon` binary. | — |
| `ToDixon(...)` / `ToDixonSolver(...)` | Write an input file without running the binary. | — |

`field_size` accepts an integer, `"p^k"` string, `GF(...)` object, or `0` for ℚ; inferred from the polynomial ring if omitted. All main functions also accept `debug=True` and `timeout` (seconds).

---

## Output

| Mode | Command-line input | File input `example.dr` |
|---|---|---|
| Dixon / Solver | `solution_YYYYMMDD_HHMMSS.dr` | `example_solution.dr` |
| Complexity | `comp_YYYYMMDD_HHMMSS.dr` | `example_comp.dr` |

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
- Over Q (`field_size=0`, `--ideal`, and `--field-equation` are not supported

---

## Feature Support by Field

| Feature | F_p (p<2^63) | F_p (p>2^63) | F_{p^k} (p<2^63) | Q |
|---|---|---|---|---|
| Dixon resultant | ✅ | ✅ | ✅ | ✅ |
| Complexity analysis (`--comp`) | ✅ | ✅ | ✅ | ✅ |
| Random mode (`-r`) | ✅ | ✅ | ✅ | ✅ |
| Polynomial solver (`--solve`) | ✅ | ❌ | ✅ | ✅ |
| Ideal reduction (`--ideal`) | ✅ | ❌ | ✅ | ❌ |
| Field-equation reduction | ✅ | ❌ | ✅ | ❌ |
| PML acceleration | ✅ | ✅ | ❌ | ✅ |

---

## License
DixonRes is distributed under the GNU General Public License version 2.0 (GPL-2.0-or-later). See the file COPYING.
