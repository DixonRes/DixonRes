"""
SageMath interface for DixonRes
================================

Usage
-----
    load("dixon_sage_interface.sage")

    R.<x, y, z> = GF(257)[]
    F = [x + y + z, x*y + y*z + z*x, x*y*z + 1]

    res  = DixonResultant(F, [x, y], dixon_path="./dixon")
    sols = DixonSolve(F, dixon_path="./dixon")
    info = DixonComplexity(F, [x, y], dixon_path="./dixon")

Pass debug=True to any function for verbose diagnostics.
"""

import os
import re
import glob
import subprocess
from itertools import product as _iproduct


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _field_size_str(field_size):
    """
    Convert various field-size representations to the string DixonRes expects.

    0 / "0"       -> "0"    (rational / Q mode)
    257           -> "257"
    (2, 8)        -> "2^8"
    GF(257)       -> "257"
    GF(2^8)       -> "2^8"
    "2^8"         -> "2^8"  (passed through)
    """
    try:
        from sage.rings.finite_rings.finite_field_base import FiniteField
        if isinstance(field_size, FiniteField):
            p = int(field_size.characteristic())
            k = int(field_size.degree())
            return "%d^%d" % (p, k) if k > 1 else str(p)
    except ImportError:
        pass

    if isinstance(field_size, tuple) and len(field_size) == 2:
        p, k = int(field_size[0]), int(field_size[1])
        return "%d^%d" % (p, k) if k > 1 else str(p)

    if isinstance(field_size, str):
        return field_size

    return str(int(field_size))


def _poly_to_str(f):
    """Sage polynomial -> string DixonRes can parse."""
    return str(f)


def _elim_vars_to_str(elim_vars):
    return ", ".join(str(v) for v in elim_vars)


def _run(cmd, timeout, debug):
    if debug:
        print("[debug] command: %s" % " ".join(cmd))
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if debug:
            print("[debug] return code: %d" % proc.returncode)
            if proc.stdout.strip():
                print("[debug] stdout:\n%s" % proc.stdout.rstrip())
            if proc.stderr.strip():
                print("[debug] stderr:\n%s" % proc.stderr.rstrip())
        return proc
    except FileNotFoundError:
        raise RuntimeError(
            "DixonRes binary not found: '%s'.\n"
            "Pass the correct path via dixon_path=..." % cmd[0]
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError("DixonRes timed out after %d s." % timeout)


def _find_output_file_tagged(finput, tag):
    """
    Reproduce generate_tagged_filename() from the C source:
    insert *tag* before the last dot-extension.

      /tmp/dixon_in.dat  +  "_solution"  ->  /tmp/dixon_in_solution.dat
      /tmp/dixon_in      +  "_solution"  ->  /tmp/dixon_in_solution
    """
    dot = finput.rfind(".")
    if dot != -1:
        return finput[:dot] + tag + finput[dot:]
    return finput + tag


def _locate_output_file(finput, tag, debug):
    """
    Try several candidate paths for the output file DixonRes wrote.

    Priority:
    1. The tagged path (file-input mode):
         generate_tagged_filename(finput, tag)
    2. The most recently modified  <tag[1:]>_*.dat  in the CWD
       (CLI mode uses generate_timestamped_filename which writes to CWD)
    3. Any <tag[1:]>*.dat in the same directory as finput

    Returns the path string if found, else None.
    """
    # 1. tagged path next to the input file
    tagged = _find_output_file_tagged(finput, tag)
    if debug:
        print("[debug] looking for tagged output file: %s" % tagged)
    if os.path.isfile(tagged):
        if debug:
            print("[debug] found tagged output file: %s" % tagged)
        return tagged

    # 2. timestamped file in CWD  (solution_*.dat / comp_*.dat)
    pattern_cwd = "%s_*.dat" % tag.lstrip("_")
    candidates  = sorted(glob.glob(pattern_cwd), key=os.path.getmtime)
    if debug:
        print("[debug] CWD glob '%s': %s" % (pattern_cwd, candidates))
    if candidates:
        if debug:
            print("[debug] using most recent CWD file: %s" % candidates[-1])
        return candidates[-1]

    # 3. same dir as finput
    d = os.path.dirname(finput) or "."
    pattern_dir = os.path.join(d, "%s_*.dat" % tag.lstrip("_"))
    candidates  = sorted(glob.glob(pattern_dir), key=os.path.getmtime)
    if debug:
        print("[debug] dir glob '%s': %s" % (pattern_dir, candidates))
    if candidates:
        return candidates[-1]

    if debug:
        print("[debug] output file NOT found for tag '%s'" % tag)
    return None


# ---------------------------------------------------------------------------
# File writers
# ---------------------------------------------------------------------------

def ToDixon(F, elim_vars, field_size=257, finput="/tmp/dixon_in.dat", debug=False):
    """
    Write a DixonRes input file (resultant / complexity mode).

    Format:
      Line 1      : field size
      Lines 2..n-1: one polynomial per line
      Line n      : comma-separated variables to ELIMINATE
    """
    A = F[0].parent()
    assert all(f.parent() == A for f in F), \
        "All polynomials must belong to the same ring."

    with open(finput, "w") as fd:
        fd.write(_field_size_str(field_size) + "\n")
        fd.write(", ".join(_poly_to_str(f) for f in F) + "\n")
        fd.write(_elim_vars_to_str(elim_vars) + "\n")

    if debug:
        print("[debug] wrote input file: %s" % finput)
        with open(finput) as fh:
            print("[debug] --- input file content ---")
            print(fh.read().rstrip())
            print("[debug] --- end ---")

    return finput


def ToDixonSolver(F, field_size=257, finput="/tmp/dixon_solve_in.dat", debug=False):
    """
    Write a DixonRes solver input file (no elimination-variable line).

    Format:
      Line 1     : field size
      Lines 2..n : one polynomial per line
    """
    A = F[0].parent()
    assert all(f.parent() == A for f in F), \
        "All polynomials must belong to the same ring."

    with open(finput, "w") as fd:
        fd.write(_field_size_str(field_size) + "\n")
        fd.write(", ".join(_poly_to_str(f) for f in F) + "\n")

    if debug:
        print("[debug] wrote solver input file: %s" % finput)
        with open(finput) as fh:
            print("[debug] --- input file content ---")
            print(fh.read().rstrip())
            print("[debug] --- end ---")

    return finput


# ---------------------------------------------------------------------------
# Output parsers
# ---------------------------------------------------------------------------

def _parse_resultant_file(foutput, debug):
    """
    Parse a DixonRes resultant output file.
    Returns the resultant as a raw string, or None if not found.
    """
    try:
        with open(foutput) as f:
            content = f.read()
    except FileNotFoundError:
        if debug:
            print("[debug] output file does not exist: %s" % foutput)
        return None

    if debug:
        print("[debug] --- raw output file content ---")
        print(content.rstrip())
        print("[debug] --- end ---")

    # DixonRes writes:  "\nResultant:\n<polynomial>\n"
    m = re.search(r"Resultant:\n(.*)", content, re.DOTALL)
    if not m:
        if debug:
            print("[debug] regex 'Resultant:\\n...' did NOT match")
        return None

    raw = m.group(1).strip()
    if debug:
        print("[debug] parsed resultant: %r" % raw[:120])
    return raw


def _parse_solutions_file(foutput, debug):
    """
    Parse a DixonRes solver output file.

    Returns
    -------
    list of dict  {var_name: value_str}   one dict per solution
    []                                    no solutions
    "infinite"                            positive-dimensional system
    None                                  parse failure
    """
    try:
        with open(foutput) as f:
            content = f.read()
    except FileNotFoundError:
        if debug:
            print("[debug] solver output file does not exist: %s" % foutput)
        return None

    if debug:
        print("[debug] --- raw solver output file ---")
        print(content.rstrip())
        print("[debug] --- end ---")

    if "has no solutions" in content:
        return []
    if "positive dimension" in content or "positive-dimensional" in content:
        return "infinite"

    solutions = []

    compat = re.search(
        r"=== Compatibility View ===(.*?)=== Solution Complete ===",
        content, re.DOTALL,
    )
    if compat:
        var_vals = {}
        for line in compat.group(1).strip().splitlines():
            m = re.match(r"(\w+)\s*=\s*\{(.*?)\}", line.strip())
            if m:
                var_name = m.group(1)
                raw = m.group(2).strip()
                var_vals[var_name] = [v.strip() for v in raw.split(",")] if raw else []
        if debug:
            print("[debug] compatibility view vars: %s" % list(var_vals.keys()))
        if var_vals:
            keys = list(var_vals.keys())
            for combo in _iproduct(*[var_vals[k] for k in keys]):
                solutions.append(dict(zip(keys, combo)))
            return solutions

    blocks = re.findall(
        r"Solution set \d+:\n(.*?)(?=\nSolution set|\n===|$)",
        content, re.DOTALL,
    )
    for block in blocks:
        sol = {}
        for line in block.strip().splitlines():
            m = re.match(r"\s*(\w+)\s*=\s*(.+)", line)
            if m:
                sol[m.group(1)] = m.group(2).strip()
        if sol:
            solutions.append(sol)

    return solutions if solutions else None


def _parse_complexity_file(foutput, debug):
    """Parse a DixonRes complexity output file."""
    try:
        with open(foutput) as f:
            content = f.read()
    except FileNotFoundError:
        if debug:
            print("[debug] complexity output file does not exist: %s" % foutput)
        return None

    if debug:
        print("[debug] --- raw complexity output file ---")
        print(content.rstrip())
        print("[debug] --- end ---")

    result = {}

    m = re.search(r"Complexity \(log2, omega=([\d.]+)\):\s*([\d.]+|inf)", content)
    if m:
        result["omega"]           = float(m.group(1))
        result["complexity_log2"] = float(m.group(2))

    m = re.search(r"Bezout bound.*?:\s*(\d+)", content)
    if m:
        result["bezout_bound"] = int(m.group(1))

    m = re.search(r"Dixon matrix size:\s*(\d+)", content)
    if m:
        result["matrix_size"] = int(m.group(1))

    m = re.search(r"Degree sequence:\s*\[(.*?)\]", content)
    if m:
        result["degrees"] = [int(d) for d in m.group(1).split(",")]

    return result or None


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def DixonResultant(
    F,
    elim_vars,
    field_size=0,
    dixon_path="./dixon",
    finput="/tmp/dixon_in.dat",
    debug=False,
    timeout=600,
):
    """
    Compute the Dixon resultant of a polynomial system.

    Parameters
    ----------
    F          : list of Sage polynomials
    elim_vars  : variables to eliminate (Sage vars or strings)
    field_size : prime, prime power, (p,k), GF(...), or 0 for Q
    dixon_path : path to the dixon executable
    finput     : temporary input file path
    debug      : print detailed diagnostics
    timeout    : seconds before aborting

    Returns
    -------
    str   raw resultant polynomial string
    None  on failure
    """
    ToDixon(F, elim_vars, field_size, finput, debug=debug)

    # Use file-input mode (not --silent) so stdout is captured but
    # the output file is still written next to the input file.
    cmd  = [dixon_path, finput]
    proc = _run(cmd, timeout, debug)

    if proc.returncode != 0:
        print("[DixonResultant] DixonRes exited with code %d" % proc.returncode)
        if proc.stderr:
            print(proc.stderr)
        return None

    foutput = _locate_output_file(finput, "_solution", debug)
    if foutput is None:
        print("[DixonResultant] could not locate output file")
        return None

    result = _parse_resultant_file(foutput, debug)

    if not debug:
        try:
            os.remove(foutput)
        except OSError:
            pass

    return result


def DixonSolve(
    F,
    field_size=None,
    dixon_path="./dixon",
    finput="/tmp/dixon_solve_in.dat",
    debug=False,
    timeout=600,
):
    """
    Solve a polynomial system with DixonRes.

    Parameters
    ----------
    F          : list of Sage polynomials
    field_size : prime or prime power; inferred from F[0].base_ring() if None
    dixon_path : path to the dixon executable
    finput     : temporary input file
    debug      : print detailed diagnostics
    timeout    : seconds before aborting

    Returns
    -------
    list of dict  {var_name: value_str}  one dict per solution
    []                                   no solutions
    "infinite"                           positive-dimensional system
    None                                 failure
    """
    if field_size is None:
        field_size = F[0].base_ring()

    ToDixonSolver(F, field_size, finput, debug=debug)

    cmd  = [dixon_path, "--solve", finput]
    proc = _run(cmd, timeout, debug)

    if proc.returncode != 0:
        print("[DixonSolve] DixonRes exited with code %d" % proc.returncode)
        if proc.stderr:
            print(proc.stderr)
        return None

    foutput   = _locate_output_file(finput, "_solution", debug)
    if foutput is None:
        print("[DixonSolve] could not locate output file")
        return None

    solutions = _parse_solutions_file(foutput, debug)

    if not debug:
        try:
            os.remove(foutput)
        except OSError:
            pass

    return solutions


def DixonComplexity(
    F,
    elim_vars,
    field_size=257,
    omega=None,
    dixon_path="./dixon",
    finput="/tmp/dixon_comp_in.dat",
    debug=False,
    timeout=120,
):
    """
    Run DixonRes complexity analysis.

    Parameters
    ----------
    F          : list of Sage polynomials
    elim_vars  : variables to eliminate
    field_size : prime or prime power
    omega      : matrix-multiplication exponent (default: DixonRes built-in)
    dixon_path : path to the dixon executable
    finput     : temporary input file
    debug      : print detailed diagnostics
    timeout    : seconds before aborting

    Returns
    -------
    dict with keys: complexity_log2, omega, bezout_bound, matrix_size, degrees
    None on failure
    """
    ToDixon(F, elim_vars, field_size, finput, debug=debug)

    cmd = [dixon_path, "--comp"]
    if omega is not None:
        cmd += ["--omega", str(omega)]
    cmd.append(finput)

    proc = _run(cmd, timeout, debug)

    if proc.returncode != 0:
        print("[DixonComplexity] DixonRes exited with code %d" % proc.returncode)
        return None

    foutput = _locate_output_file(finput, "_comp", debug)
    if foutput is None:
        print("[DixonComplexity] could not locate output file")
        return None

    result = _parse_complexity_file(foutput, debug)

    if not debug:
        try:
            os.remove(foutput)
        except OSError:
            pass

    return result


def DixonIdeal(
    F,
    ideal_gens,
    elim_vars,
    field_size=257,
    dixon_path="./dixon",
    debug=False,
    timeout=600,
):
    """
    Dixon resultant with ideal reduction.

    Parameters
    ----------
    F          : list of Sage polynomials
    ideal_gens : list of strings like "a^3=2*b+1", or Sage expressions
    elim_vars  : variables to eliminate
    field_size : prime or prime power
    dixon_path : path to the dixon executable
    debug      : print detailed diagnostics
    timeout    : seconds

    Returns
    -------
    str resultant string, or None on failure
    """
    ideal_str = ", ".join(str(g) for g in ideal_gens)
    polys_str = ", ".join(_poly_to_str(f) for f in F)
    elim_str  = _elim_vars_to_str(elim_vars)
    field_str = _field_size_str(field_size)

    if debug:
        print("[debug] ideal_str : %s" % ideal_str)
        print("[debug] polys_str : %s" % polys_str)
        print("[debug] elim_str  : %s" % elim_str)
        print("[debug] field_str : %s" % field_str)

    cmd = [dixon_path, "--ideal", ideal_str, polys_str, elim_str, field_str]
    proc = _run(cmd, timeout, debug)

    if proc.returncode != 0:
        print("[DixonIdeal] DixonRes exited with code %d" % proc.returncode)
        return None

    # CLI mode: timestamped solution_*.dat written to CWD
    candidates = sorted(glob.glob("solution_*.dat"), key=os.path.getmtime)
    if debug:
        print("[debug] solution_*.dat in CWD: %s" % candidates)
    if not candidates:
        print("[DixonIdeal] output file not found")
        return None

    foutput = candidates[-1]
    result  = _parse_resultant_file(foutput, debug)

    if not debug:
        try:
            os.remove(foutput)
        except OSError:
            pass

    return result
