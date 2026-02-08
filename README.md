# Dixon Resultant & Polynomial System Solver

A C implementation for computing Dixon resultants and solving polynomial systems over finite fields, based on the FLINT library.

**Version:** v0.0.1

---

## Features

- Dixon resultant computation for variable elimination
- Polynomial system solver for n×n systems
- Dixon with triangular ideal reduction
- Finite fields:
  - Prime fields F_p (p < 2^63)
  - Extension fields F_{p^k} (e.g. 2^8, 3^5)
- Command line input or file input
- Silent mode for batch computation
- Automatic output to solution files

---

## Dependencies

- **FLINT**

Optional:
- **PML** (used automatically if available)

---

## Build

```bash
make
```

The Makefile automatically detects required libraries. For more options, run `make help`.

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

### Extension Fields

```bash
./dixon "x + y^2 + t, x*y + t*y + 1" "x" 2^8
```

`t` is used as the extension field generator by default.

```bash
./dixon --solve "x^2 + t*y, x*y + t^2" "2^8: t^8 + t^4 + t^3 + t + 1"
```

(with AES custom polynomial for F_256)

---

---

### Dixon with Ideal Reduction

```bash
./dixon "polynomials" "eliminate_vars" "ideal_generators" "all_variables" field_size
```

Example:

```bash
./dixon "a1^2 + a2^2 + a3^2 + a4^2 - 10, a4^3 - a1 - a2*a3 - 5" "a4" "a2^3 = 2*a1 + 1, a3^3 = a1*a2 + 3, a4^3 = a1 + a2*a3 + 5" 257"
```

```
### Silent Mode

```bash
./dixon --silent [--solve] <arguments>
```

No console output is produced; the solution file is still generated.

---

## File Input Format

### Dixon Mode (multiline)

```
Line 1 : field size
Line 2+: polynomials (comma-separated or multiline)
Last   : variables to ELIMINATE (comma-separated)
         (#eliminate = #equations - 1)
```

### Polynomial Solver Mode (multiline)

```
Line 1 : field size
Line 2+: polynomials
         (n equations in n variables)
```

---

## Output

* **Command line input**:
  `solution_YYYYMMDD_HHMMSS.dat`

* **File input `example.dat`**:
  `example_solution_YYYYMMDD_HHMMSS.dat`

Each output file contains:

* Field information
* Input polynomials
* Computation time
* Resultant or solutions

---

## Notes

* All computation modes generate a solution file by default
* Extension fields are slower than prime fields due to polynomial arithmetic
* The solver computes complete solutions for all variables

诉别人论文该怎么引你的代码。
```
