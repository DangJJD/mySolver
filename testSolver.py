#!/usr/bin/env python3
"""
ECE 3340 - Programming Assignment 1 - Root Finding Test Suite

Runs ./mySolver against:
  1. Interface tests      - CLI argument handling, output format
  2. Known-answer tests   - hand-picked problems with exact roots, every solver
  3. General cases        - randomly generated problems (seeded, reproducible)
  4. Corner cases         - domain edges, steep/flat functions, endpoint roots
  5. Expected-break cases - inputs where a method SHOULD fail; the test proves it does

Every result is put in one of these buckets:
  PASS     - behaved correctly
  FAIL     - a bug you should fix (with a hint about where to look)
  WARN     - robustness / polish issue; worth fixing but not a wrong answer
  BROKE-OK - the method is *supposed* to fail here, and it did (shows why)
  NOTE     - expected to fail but succeeded anyway (not a bug, just informative)

Usage:
  python3 test_solver.py                     # everything
  python3 test_solver.py --solver ./mySolver.exe   # Windows
  python3 test_solver.py --random 50 --seed 7
  python3 test_solver.py --only Newton2      # one solver
  python3 test_solver.py --csv results.csv   # dump raw data for the report
  python3 test_solver.py -v                  # show every PASS too
"""

import argparse
import csv
import math
import os
import random
import re
import subprocess
import sys
import time

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

SOLVERS = [
    "BiSec",
    "FalsePos",
    "Newton1",
    "Newton2",
    "Bracket-Newton1",
    "Bracket-Newton2",
]
BRACKETING = {"BiSec", "FalsePos", "Bracket-Newton1", "Bracket-Newton2"}
NEWTON_PLAIN = {"Newton1", "Newton2"}
TIMEOUT_S = 5.0

OUTPUT_RE = re.compile(r"Solver output:\s*([-+0-9.eEinfnaINFNA]+)")
ITER_RE = re.compile(r"[Ii]terations?\s*[:=]\s*(\d+)")

PI = math.pi

# --------------------------------------------------------------------------
# Reference functions (what the C program is supposed to compute)
# --------------------------------------------------------------------------


def f_ref(problem, args, x):
    """f(x) for the problem, matching the assignment definitions."""
    if problem == "POLY":
        s = 0.0
        for c in args:  # Horner, highest power first
            s = s * x + c
        return s
    c = args[0]
    if problem == "COS":
        return math.cos(x) - c
    if problem == "SIN":
        return math.sin(x) - c
    if problem == "TAN":
        return math.tan(x) - c
    if problem == "EXP":
        return math.exp(x) - c
    if problem == "LOG":
        return math.log(x) - c
    raise ValueError(problem)


def df_ref(problem, args, x):
    if problem == "POLY":
        n = len(args) - 1
        d = 0.0
        for i, c in enumerate(args[:-1]):
            d = d * x + c * (n - i)
        return d
    if problem == "COS":
        return -math.sin(x)
    if problem == "SIN":
        return math.cos(x)
    if problem == "TAN":
        return 1.0 / math.cos(x) ** 2
    if problem == "EXP":
        return math.exp(x)
    if problem == "LOG":
        return 1.0 / x
    raise ValueError(problem)


def poly_from_roots(roots, lead=1.0):
    """Coefficients (highest first) of lead * prod(x - r)."""
    coeffs = [lead]
    for r in roots:
        new = coeffs + [0.0]
        for i in range(len(coeffs)):
            new[i + 1] -= r * coeffs[i]
        coeffs = new
    return coeffs


def x_tolerance(problem, args, root):
    """
    How close the returned x must be to the true root.

    The solver stops when |f(x)| < 1e-8, so the achievable x-error is about
    1e-8 / |f'(root)|.  We allow 100x that, with a floor tied to double
    precision and a ceiling of 1e-5 relative so a sloppy answer still fails.
    """
    try:
        slope = abs(df_ref(problem, args, root))
    except (ValueError, ZeroDivisionError, OverflowError):
        slope = 0.0
    scale = 1.0 + abs(root)
    floor = 1e-12 * scale
    ceil = 1e-5 * scale
    if slope == 0:
        return ceil
    return min(max(1e-6 / slope, floor), ceil)


# --------------------------------------------------------------------------
# Running the program
# --------------------------------------------------------------------------


class Run:
    def __init__(self, cmd, rc, out, err, secs, timed_out=False):
        self.cmd, self.rc, self.out, self.err = cmd, rc, out, err
        self.secs, self.timed_out = secs, timed_out
        self.value = None
        self.iters = None
        m = OUTPUT_RE.search(out)
        if m:
            try:
                self.value = float(m.group(1))
            except ValueError:
                pass
        m = ITER_RE.search(out + "\n" + err)
        if m:
            self.iters = int(m.group(1))

    @property
    def ok(self):
        return (not self.timed_out) and self.rc == 0 and self.value is not None

    @property
    def extra_stdout(self):
        """Lines on stdout other than the 'Solver output:' line."""
        return [
            l
            for l in self.out.splitlines()
            if l.strip() and not OUTPUT_RE.search(l) and not ITER_RE.search(l)
        ]

    def why(self):
        if self.timed_out:
            return f"timed out after {TIMEOUT_S}s (infinite loop?)"
        msg = self.err.strip() or self.out.strip() or "(no output)"
        msg = " | ".join(msg.splitlines()[-3:])
        return f"exit {self.rc}: {msg}"

    def cmdline(self):
        return " ".join(self.cmd)


def run_solver(exe, problem, solver, args):
    cmd = [exe, problem, solver] + [fmt(a) for a in args]
    t0 = time.perf_counter()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_S)
        dt = time.perf_counter() - t0
        return Run(cmd, p.returncode, p.stdout, p.stderr, dt)
    except subprocess.TimeoutExpired as e:
        dt = time.perf_counter() - t0
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
        err = e.stderr.decode() if isinstance(e.stderr, bytes) else (e.stderr or "")
        return Run(cmd, -1, out, err, dt, timed_out=True)


def fmt(a):
    if isinstance(a, str):
        return a
    return repr(float(a))  # full precision, e.g. 1.5707963267948966


# --------------------------------------------------------------------------
# Result bookkeeping
# --------------------------------------------------------------------------


class Results:
    COLORS = {
        "PASS": "\033[32m",
        "FAIL": "\033[31m",
        "WARN": "\033[33m",
        "BROKE-OK": "\033[36m",
        "NOTE": "\033[35m",
    }
    RESET = "\033[0m"

    def __init__(self, verbose, color):
        self.verbose, self.color = verbose, color
        self.rows = []  # for CSV / per-solver stats
        self.counts = {k: 0 for k in self.COLORS}
        self.problems = []  # FAIL/WARN details for the final summary

    def tag(self, status):
        if self.color:
            return f"{self.COLORS[status]}{status:<8}{self.RESET}"
        return f"{status:<8}"

    def record(
        self, section, name, solver, status, detail="", run=None, hint="", expected=None
    ):
        self.counts[status] += 1
        self.rows.append(
            {
                "section": section,
                "test": name,
                "solver": solver,
                "status": status,
                "expected": "" if expected is None else repr(expected),
                "got": "" if (run is None or run.value is None) else repr(run.value),
                "abs_error": (
                    ""
                    if (run is None or run.value is None or expected is None)
                    else abs(run.value - expected)
                ),
                "seconds": "" if run is None else f"{run.secs:.6f}",
                "iterations": "" if (run is None or run.iters is None) else run.iters,
                "detail": detail,
                "command": "" if run is None else run.cmdline(),
            }
        )
        show = status != "PASS" or self.verbose
        if show:
            print(f"  {self.tag(status)} [{solver}] {name}")
            if detail:
                print(f"           {detail}")
            if run is not None and status in ("FAIL", "WARN", "BROKE-OK"):
                print(f"           cmd: {run.cmdline()}")
            if hint and status in ("FAIL", "WARN"):
                print(f"           fix: {hint}")
        if status in ("FAIL", "WARN"):
            self.problems.append((status, section, solver, name, hint))


# --------------------------------------------------------------------------
# Test helpers
# --------------------------------------------------------------------------


def auto_hint(solver, run):
    """Best guess at the cause, from the solver's own error message."""
    err = (run.err + run.out).lower()
    if run.timed_out:
        return "a loop has no working exit - check max-iteration stop"
    if "max loops" in err and solver == "FalsePos":
        return (
            "plain false position keeps one endpoint fixed on curved "
            "functions and crawls; 100 iterations isn't enough. Use the "
            "Illinois modification (when the same endpoint survives twice, "
            "halve its stored f value) and/or add an x-based stop"
        )
    if "max loops" in err:
        return (
            "100 iterations weren't enough to reach |f| < 1e-8; also stop "
            "when the bracket width / step size gets tiny"
        )
    if "out of bounds" in err and solver.endswith("2"):
        return (
            "the numerical derivative probes x-h / x+h outside the "
            "problem's domain; shrink h near the edge or use a one-sided "
            "difference"
        )
    return ""


def check_root(
    R,
    exe,
    section,
    name,
    problem,
    solver,
    coeffs,
    a,
    b,
    root,
    fail_status="FAIL",
    hint="",
    tol=None,
    newton_may_diverge=True,
):
    """
    Expect the solver to find `root` in [a, b].
    Plain Newton is not bracketed, so if it errors out (diverges) that is
    reported as BROKE-OK unless newton_may_diverge=False (a real bug).
    A WRONG value is always a FAIL.
    """
    args = list(coeffs) + [a, b]
    run = run_solver(exe, problem, solver, args)
    if tol is None:
        tol = x_tolerance(problem, coeffs, root)
    if not run.ok:
        if solver in NEWTON_PLAIN and newton_may_diverge and not run.timed_out:
            R.record(
                section,
                name,
                solver,
                "BROKE-OK",
                f"plain Newton diverged from the midpoint start (no "
                f"bracket to save it) -> {run.why()}",
                run,
                expected=root,
            )
            return run
        R.record(
            section,
            name,
            solver,
            fail_status,
            f"expected root {root:.12g}, solver failed -> {run.why()}",
            run,
            hint or auto_hint(solver, run),
            expected=root,
        )
        return run
    err = abs(run.value - root)
    lo, hi = min(a, b), max(a, b)
    if not (lo - 1e-12 <= run.value <= hi + 1e-12) and solver in BRACKETING:
        R.record(
            section,
            name,
            solver,
            "FAIL",
            f"returned {run.value:.17g} which is OUTSIDE [{a}, {b}]",
            run,
            "a bracketing method must never leave its interval",
            expected=root,
        )
    elif err > tol:
        R.record(
            section,
            name,
            solver,
            fail_status,
            f"wrong root: got {run.value:.17g}, expected {root:.17g} "
            f"(err {err:.2e} > tol {tol:.1e})",
            run,
            hint or "check the update formula / stopping test",
            expected=root,
        )
    else:
        R.record(
            section,
            name,
            solver,
            "PASS",
            f"got {run.value:.17g}  err {err:.1e}  {run.secs * 1e3:.2f} ms",
            run,
            expected=root,
        )
    return run


def check_error(
    R,
    exe,
    section,
    name,
    problem,
    solver,
    cli_args,
    fail_status="FAIL",
    hint="",
    reason="",
):
    """
    The program SHOULD refuse / fail on this input.  Passing means it exited
    non-zero without hanging and without printing a bogus 'Solver output'.
    """
    run = run_solver(exe, problem, solver, cli_args)
    if run.timed_out:
        R.record(
            section,
            name,
            solver,
            "FAIL",
            f"should have reported an error but hung -> {run.why()}",
            run,
            "every loop needs a max-iteration exit",
        )
    elif run.rc != 0 and run.value is None:
        R.record(
            section,
            name,
            solver,
            "PASS",
            f"correctly rejected ({reason}) -> {run.why()}",
            run,
        )
    else:
        got = run.value if run.value is not None else "(no value)"
        R.record(
            section,
            name,
            solver,
            fail_status,
            f"should have failed ({reason}) but exited {run.rc} and printed {got}",
            run,
            hint,
        )
    return run


def check_breaks(
    R, exe, section, name, problem, solver, coeffs, a, b, reason, root=None
):
    """
    The METHOD is mathematically expected to fail here (not a coding bug).
    Prove it fails cleanly: non-zero exit, no hang, no fake answer.
    If it returns a value, that value must at least be correct.
    """
    args = list(coeffs) + [a, b]
    run = run_solver(exe, problem, solver, args)
    if run.timed_out:
        R.record(
            section,
            name,
            solver,
            "FAIL",
            f"expected a clean failure ({reason}) but it hung",
            run,
            "add/verify the max-iteration stop",
            expected=root,
        )
    elif run.rc != 0:
        R.record(
            section,
            name,
            solver,
            "BROKE-OK",
            f"{reason}. Program said -> {run.why()}",
            run,
            expected=root,
        )
    elif (
        run.value is not None
        and root is not None
        and abs(run.value - root) > x_tolerance(problem, coeffs, root) * 10
    ):
        R.record(
            section,
            name,
            solver,
            "FAIL",
            f"expected to fail ({reason}) but instead returned a WRONG "
            f"answer {run.value:.17g} (true root {root:.17g})",
            run,
            "a solver should exit with an error rather than return a non-root",
            expected=root,
        )
    else:
        R.record(
            section,
            name,
            solver,
            "NOTE",
            f"expected to break ({reason}) but converged to "
            f"{run.value!r} anyway - fine",
            run,
            expected=root,
        )
    return run


# --------------------------------------------------------------------------
# 1. Interface tests
# --------------------------------------------------------------------------


def test_interface(R, exe):
    sec = "interface"
    print("\n== 1. Command-line interface ==")

    check_error(
        R,
        exe,
        sec,
        "too few arguments",
        "POLY",
        "BiSec",
        ["1", "2"],
        reason="needs coefficients + a + b",
    )
    check_error(
        R,
        exe,
        sec,
        "unknown problem type",
        "QUAD",
        "BiSec",
        ["1", "0", "1"],
        reason="QUAD is not a problem type",
    )
    check_error(
        R,
        exe,
        sec,
        "unknown solver type",
        "COS",
        "Secant",
        ["0", "0", "3"],
        reason="Secant is not a solver",
    )
    check_error(
        R,
        exe,
        sec,
        "lowercase problem type",
        "cos",
        "BiSec",
        ["0", "0", "3"],
        fail_status="WARN",
        reason="names are case-sensitive in the spec",
        hint="either reject it (current spec) or document that "
        "case-insensitive names are accepted",
    )

    # Non-numeric arguments: atof("abc") silently returns 0.
    check_error(
        R,
        exe,
        sec,
        "non-numeric coefficient",
        "POLY",
        "BiSec",
        ["1", "abc", "-4", "0", "5"],
        fail_status="WARN",
        reason="'abc' is not a number",
        hint="atof() returns 0 for garbage. Use strtod() and check "
        "that endptr reached the end of the string",
    )
    check_error(
        R,
        exe,
        sec,
        "non-numeric bound",
        "COS",
        "BiSec",
        ["0.5", "zero", "3"],
        fail_status="WARN",
        reason="'zero' is not a number",
        hint="same strtod() validation for the interval bounds",
    )

    # COS/SIN/... take exactly one constant
    check_error(
        R,
        exe,
        sec,
        "COS given two constants",
        "COS",
        "BiSec",
        ["0.5", "0.7", "0", "3"],
        fail_status="WARN",
        reason="COS takes exactly: c a b",
        hint="for non-POLY problems, require length == 1 so a typo "
        "doesn't silently use the wrong value",
    )

    # Error message should name the right function
    run = run_solver(exe, "TAN", "BiSec", [1, 0, 2])
    if "SIN:" in run.err:
        R.record(
            sec,
            "TAN out-of-domain message",
            "BiSec",
            "WARN",
            f"TAN error is labelled SIN: {run.err.strip().splitlines()[-1]}",
            run,
            'copy-paste slip in findY(): change "SIN: X:" to '
            '"TAN: X:" in the TAN branch',
        )
    else:
        R.record(sec, "TAN out-of-domain message", "BiSec", "PASS", run=run)

    # Output format: exactly one parseable line on stdout
    for solver in SOLVERS:
        run = run_solver(exe, "POLY", solver, [1, -1, -2, 1, 3])
        if not run.ok:
            R.record(
                sec,
                "x^2-x-2 on [1,3] (assignment example)",
                solver,
                "FAIL",
                run.why(),
                run,
                "basic case must work",
                expected=2.0,
            )
            continue
        extra = run.extra_stdout
        if extra:
            R.record(
                sec,
                "stdout contains only the result",
                solver,
                "WARN",
                f"extra stdout lines: {extra[:4]}",
                run,
                "remove debug printf()s (or send them to stderr) - "
                "a grading script reading stdout will choke on them",
            )
        else:
            R.record(sec, "stdout contains only the result", solver, "PASS", run=run)

    # Max-args limit: argc > 20 means a degree >= 15 polynomial is refused
    coeffs = poly_from_roots([0.5] + [10 + i for i in range(15)])
    run = run_solver(exe, "POLY", "BiSec", coeffs + [0, 1])
    if run.ok:
        R.record(sec, "degree-16 polynomial accepted", "BiSec", "PASS", run=run)
    else:
        R.record(
            sec,
            "degree-16 polynomial accepted",
            "BiSec",
            "WARN",
            f"rejected: {run.why()}",
            run,
            "the 'max = 20' argc limit caps POLY at degree 14; the spec "
            "doesn't set a limit, so raise it or document it",
        )


# --------------------------------------------------------------------------
# 2. Known-answer tests (every solver, exact roots)
# --------------------------------------------------------------------------

KNOWN = [
    # name, problem, coeffs/c, a, b, exact root
    ("x^2 - x - 2 on [1,3]", "POLY", [1, -1, -2], 1, 3, 2.0),
    ("x^2 - 2 on [0,2]  (sqrt 2)", "POLY", [1, 0, -2], 0, 2, math.sqrt(2)),
    (
        "x^3 - 2x - 5 on [2,3] (Wallis)",
        "POLY",
        [1, 0, -2, -5],
        2,
        3,
        2.0945514815423265,
    ),
    ("2x - 3 linear on [0,5]", "POLY", [2, -3], 0, 5, 1.5),
    ("-x^2 + 4 on [0,3] (neg lead)", "POLY", [-1, 0, 4], 0, 3, 2.0),
    ("x^5 - x - 1 on [1,2]", "POLY", [1, 0, 0, 0, -1, -1], 1, 2, 1.1673039782614187),
    ("cos x = 0.5 on [0,pi]", "COS", [0.5], 0, PI, PI / 3),
    ("cos x = -0.3 on [0,pi]", "COS", [-0.3], 0, PI, math.acos(-0.3)),
    ("sin x = 0.5 on [-pi/2,pi/2]", "SIN", [0.5], -PI / 2, PI / 2, PI / 6),
    ("sin x = -0.8 on [-pi/2,pi/2]", "SIN", [-0.8], -PI / 2, PI / 2, math.asin(-0.8)),
    ("tan x = 1 on [-1.5,1.5]", "TAN", [1], -1.5, 1.5, PI / 4),
    ("tan x = -3 on [-1.5,1.5]", "TAN", [-3], -1.5, 1.5, math.atan(-3)),
    ("e^x = 2 on [0,1]", "EXP", [2], 0, 1, math.log(2)),
    ("e^x = 100 on [0,10]", "EXP", [100], 0, 10, math.log(100)),
    ("ln x = 1 on [1,4]", "LOG", [1], 1, 4, math.e),
    ("ln x = -2 on [0.01,1]", "LOG", [-2], 0.01, 1, math.exp(-2)),
]


def test_known(R, exe, solvers):
    sec = "known"
    print("\n== 2. Known-answer tests ==")
    for name, prob, coeffs, a, b, root in KNOWN:
        for s in solvers:
            check_root(R, exe, sec, name, prob, s, coeffs, a, b, root)

    # Reversed interval [b, a] - bracketing methods should still work
    for s in solvers:
        if s in BRACKETING:
            check_root(
                R,
                exe,
                sec,
                "reversed interval x^2-2 on [2,0]",
                "POLY",
                s,
                [1, 0, -2],
                2,
                0,
                math.sqrt(2),
                fail_status="WARN",
                hint="either swap a,b when a > b or reject the input "
                "with an error; don't silently mis-handle it",
            )


# --------------------------------------------------------------------------
# 3. General (random) cases
# --------------------------------------------------------------------------


def random_cases(rng, n):
    cases = []
    for i in range(n):
        kind = rng.choice(["POLY", "POLY", "COS", "SIN", "TAN", "EXP", "LOG"])
        if kind == "POLY":
            deg = rng.randint(1, 6)
            target = rng.uniform(-10, 10)
            # other roots are kept at least 1.5 away so [a,b] holds just one
            others = []
            while len(others) < deg - 1:
                r = rng.uniform(-15, 15)
                if abs(r - target) > 1.5 and all(abs(r - o) > 0.3 for o in others):
                    others.append(r)
            lead = rng.choice([-1, 1]) * rng.uniform(0.5, 3)
            coeffs = poly_from_roots([target] + others, lead)
            gap = min([abs(o - target) for o in others] + [1.5]) * 0.9
            a = target - rng.uniform(0.05, gap)
            b = target + rng.uniform(0.05, gap)
            cases.append((f"rand POLY deg {deg}", kind, coeffs, a, b, target))
        elif kind == "COS":
            c = rng.uniform(-0.98, 0.98)
            root = math.acos(c)
            cases.append(
                (
                    f"rand cos x = {c:.4f}",
                    kind,
                    [c],
                    rng.uniform(0, root * 0.9),
                    root + rng.uniform(0.01, (PI - root) * 0.95),
                    root,
                )
            )
        elif kind == "SIN":
            c = rng.uniform(-0.98, 0.98)
            root = math.asin(c)
            lo, hi = -PI / 2, PI / 2
            cases.append(
                (
                    f"rand sin x = {c:.4f}",
                    kind,
                    [c],
                    root - rng.uniform(0.01, (root - lo) * 0.95),
                    root + rng.uniform(0.01, (hi - root) * 0.95),
                    root,
                )
            )
        elif kind == "TAN":
            c = rng.uniform(-20, 20)
            root = math.atan(c)
            lo, hi = -PI / 2, PI / 2
            cases.append(
                (
                    f"rand tan x = {c:.4f}",
                    kind,
                    [c],
                    root - rng.uniform(0.01, (root - lo) * 0.9),
                    root + rng.uniform(0.01, (hi - root) * 0.9),
                    root,
                )
            )
        elif kind == "EXP":
            c = math.exp(rng.uniform(-5, 5))
            root = math.log(c)
            cases.append(
                (
                    f"rand e^x = {c:.4g}",
                    kind,
                    [c],
                    root - rng.uniform(0.1, 3),
                    root + rng.uniform(0.1, 3),
                    root,
                )
            )
        else:  # LOG
            c = rng.uniform(-3, 3)
            root = math.exp(c)
            cases.append(
                (
                    f"rand ln x = {c:.4f}",
                    kind,
                    [c],
                    root * rng.uniform(0.2, 0.95),
                    root * rng.uniform(1.05, 3),
                    root,
                )
            )
    return cases


def test_random(R, exe, solvers, n, seed):
    sec = "general"
    print(f"\n== 3. General cases ({n} random problems, seed={seed}) ==")
    rng = random.Random(seed)
    for name, prob, coeffs, a, b, root in random_cases(rng, n):
        # sanity: our generator must give a real bracket
        if f_ref(prob, coeffs, a) * f_ref(prob, coeffs, b) > 0:
            continue
        for s in solvers:
            if s in NEWTON_PLAIN:
                # Plain Newton is allowed to diverge from a mid-point start on
                # wide intervals; that is a property of the method, not a bug.
                run = run_solver(exe, prob, s, list(coeffs) + [a, b])
                tol = x_tolerance(prob, coeffs, root)
                if run.ok and abs(run.value - root) <= tol:
                    R.record(
                        sec,
                        name,
                        s,
                        "PASS",
                        f"err {abs(run.value - root):.1e}",
                        run,
                        expected=root,
                    )
                elif run.ok and not (min(a, b) <= run.value <= max(a, b)):
                    R.record(
                        sec,
                        name,
                        s,
                        "BROKE-OK",
                        f"converged to a different root {run.value:.10g} "
                        f"outside [{a:.4g},{b:.4g}] - plain Newton isn't "
                        f"bracketed, this is expected",
                        run,
                        expected=root,
                    )
                elif run.ok:
                    R.record(
                        sec,
                        name,
                        s,
                        "FAIL",
                        f"wrong answer {run.value:.17g} vs {root:.17g}",
                        run,
                        "check Newton update / derivative",
                        expected=root,
                    )
                elif run.timed_out:
                    R.record(
                        sec,
                        name,
                        s,
                        "FAIL",
                        run.why(),
                        run,
                        "max-iteration stop not working",
                        expected=root,
                    )
                else:
                    R.record(
                        sec,
                        name,
                        s,
                        "BROKE-OK",
                        f"Newton diverged from midpoint start -> {run.why()}",
                        run,
                        expected=root,
                    )
            else:
                check_root(R, exe, sec, name, prob, s, coeffs, a, b, root)


# --------------------------------------------------------------------------
# 4. Corner cases (should work - failures here are things to fix)
# --------------------------------------------------------------------------


def test_corners(R, exe, solvers):
    sec = "corner"
    print("\n== 4. Corner cases ==")

    for s in solvers:
        # Root exactly at an endpoint (f(a)*f(b) == 0 is allowed by the spec)
        check_root(
            R,
            exe,
            sec,
            "root at left endpoint x-1 on [1,3]",
            "POLY",
            s,
            [1, -1],
            1,
            3,
            1.0,
        )
        check_root(
            R,
            exe,
            sec,
            "root at right endpoint x-3 on [1,3]",
            "POLY",
            s,
            [1, -3],
            1,
            3,
            3.0,
        )
        # Root exactly at the midpoint
        check_root(
            R, exe, sec, "root at midpoint x-2 on [1,3]", "POLY", s, [1, -2], 1, 3, 2.0
        )
        # Domain endpoints of the trig problems
        check_root(
            R,
            exe,
            sec,
            "cos x = -1 -> root at pi (domain edge)",
            "COS",
            s,
            [-1],
            2.5,
            PI,
            PI,
            fail_status="WARN",
            tol=1e-6,
            hint="cos x = -1 is a double (flat) root: f ~ (x-pi)^2/2, "
            "so |f|<1e-8 only pins x to ~1e-4. Bracketing methods "
            "catch f(b)==0; Newton needs an x-based stop to do "
            "better",
        )
        check_root(
            R,
            exe,
            sec,
            "sin x = 0 on [-pi/2, pi/2]",
            "SIN",
            s,
            [0],
            -PI / 2,
            PI / 2,
            0.0,
        )
        check_root(R, exe, sec, "tan x = 0 on [-1, 1.5]", "TAN", s, [0], -1, 1.5, 0.0)
        # Tiny interval
        r = math.sqrt(2)
        check_root(
            R,
            exe,
            sec,
            "very narrow interval (width 2e-10)",
            "POLY",
            s,
            [1, 0, -2],
            r - 1e-10,
            r + 1e-10,
            r,
        )
        # Very wide interval / large root
        check_root(
            R,
            exe,
            sec,
            "x - 1e6 on [-1e7, 1e7] (wide)",
            "POLY",
            s,
            [1, -1e6],
            -1e7,
            1e7,
            1e6,
        )
        # Large coefficients (absolute |f| tolerance is hard to reach)
        check_root(
            R,
            exe,
            sec,
            "1e6 x^2 - 2e6 on [0, 3] (big coeffs)",
            "POLY",
            s,
            [1e6, 0, -2e6],
            0,
            3,
            math.sqrt(2),
        )
        # Tiny coefficients (absolute |f| tolerance is reached too early)
        check_root(
            R,
            exe,
            sec,
            "1e-9 (x - 0.7) on [0, 2] (tiny coeffs)",
            "POLY",
            s,
            [1e-9, -0.7e-9],
            0,
            2,
            0.7,
            hint="|f(x)| < 1e-8 is true EVERYWHERE for this function, "
            "so the solver stops at its first guess. Add an "
            "x-based stop (|b-a| or |dx| small) and don't rely on "
            "|f| alone",
        )
        # Steep function: |f| < 1e-8 is unreachable in double precision
        c = 1e8
        check_root(
            R,
            exe,
            sec,
            "tan x = 1e8 (root 1e-8 from pi/2, very steep)",
            "TAN",
            s,
            [c],
            1.5,
            PI / 2 - 1e-12,
            math.atan(c),
            hint="near the root f' ~ 1e16, so one ulp of x changes f "
            "by ~2; |f| < 1e-8 can never be met. Also stop when "
            "the bracket width (or Newton step) falls below "
            "~1e-15*|x|",
        )
        # Small root of LOG near the domain edge x > 0
        check_root(
            R,
            exe,
            sec,
            "ln x = -20 on [1e-10, 1e-8] (root 2.06e-9)",
            "LOG",
            s,
            [-20],
            1e-10,
            1e-8,
            math.exp(-20),
            newton_may_diverge=(s != "Newton2"),
            hint="numerical derivative step h = 1e-6*(1+|x|) is ~1e-6 "
            "here, far bigger than x itself, so x-h < 0 leaves the "
            "LOG domain. Use a relative step like "
            "h = sqrt(DBL_EPSILON)*max(|x|, tiny) and keep x-h, "
            "x+h inside the domain",
        )
        # Flat (triple) root: valid bracket, but |f|-only stopping is sloppy
        check_root(
            R,
            exe,
            sec,
            "(x-1)^3 on [0, 2.5]  (flat/triple root)",
            "POLY",
            s,
            [1, -3, 3, -1],
            0,
            2.5,
            1.0,
            fail_status="WARN",
            tol=1e-6,
            hint="|f|<1e-8 only guarantees |x-1| < ~2e-3 here. Stopping "
            "on interval width / step size gives much better x",
        )


# --------------------------------------------------------------------------
# 5. Inputs that SHOULD be rejected, and methods that SHOULD break
# --------------------------------------------------------------------------


def test_expected_failures(R, exe, solvers):
    sec = "expected-fail"
    print("\n== 5a. Bad inputs that must be rejected (all solvers) ==")

    # newton = how to judge PLAIN Newton (it never looks at f(a), f(b)):
    #   "error"  - no root exists at all, so Newton must still exit with error
    #   "domain" - interval is invalid for the problem; Newton may never touch
    #              the bad region, so only WARN that inputs aren't validated
    #   "breaks" - a root exists (maybe outside [a,b]); finding it is fine
    bad = [
        (
            "no sign change: x^2+1 on [-2,2]",
            "POLY",
            [1, 0, 1],
            -2,
            2,
            "f(a)*f(b) > 0, no real root",
            "error",
            None,
        ),
        (
            "double root (x-1)^2 on [0,3]",
            "POLY",
            [1, -2, 1],
            0,
            3,
            "f(a)*f(b) > 0; spec only guarantees simple roots",
            "breaks",
            1.0,
        ),
        (
            "cos x = 2 (c outside [-1,1])",
            "COS",
            [2],
            0,
            PI,
            "cos never reaches 2",
            "error",
            None,
        ),
        (
            "sin x = -1.5 (c outside [-1,1])",
            "SIN",
            [-1.5],
            -PI / 2,
            PI / 2,
            "sin never reaches -1.5",
            "error",
            None,
        ),
        (
            "e^x = -1 (c <= 0)",
            "EXP",
            [-1],
            -5,
            5,
            "e^x is always positive",
            "error",
            None,
        ),
        ("e^x = 0 (c <= 0)", "EXP", [0], -5, 5, "e^x is never 0", "error", None),
        (
            "cos interval leaves [0,pi]",
            "COS",
            [0.5],
            -1,
            4,
            "x outside COS domain [0, pi]",
            "domain",
            None,
        ),
        (
            "tan interval contains pi/2",
            "TAN",
            [1],
            0,
            2,
            "tan undefined at pi/2",
            "domain",
            None,
        ),
        (
            "log interval includes x <= 0",
            "LOG",
            [0],
            -1,
            5,
            "ln undefined for x <= 0",
            "domain",
            None,
        ),
        (
            "interval misses the root",
            "POLY",
            [1, -1, -2],
            3,
            5,
            "root 2 is not in [3,5]",
            "breaks",
            2.0,
        ),
    ]
    for name, prob, coeffs, a, b, reason, newton, root in bad:
        cli = [fmt(v) for v in list(coeffs) + [a, b]]
        for s in solvers:
            if s in BRACKETING:
                check_error(
                    R,
                    exe,
                    sec,
                    name,
                    prob,
                    s,
                    cli,
                    reason=reason,
                    hint="check f(a)*f(b) > 0 (and that a, b are in "
                    "the domain) up front, then exit with a "
                    "clear error",
                )
            elif newton == "error":
                hint = (
                    "Newton stopped because |f(x)| < 1e-8 at a point "
                    "that is NOT a root (e^-19 ~ 5.6e-9). An absolute "
                    "|f| test alone accepts fake roots; also require the "
                    "step |dx| to be small, or reject c <= 0 for EXP up "
                    "front"
                    if prob == "EXP"
                    else "Newton must end in max-iterations / tiny-derivative "
                    "/ non-finite error, never a value"
                )
                check_error(R, exe, sec, name, prob, s, cli, reason=reason, hint=hint)
            elif newton == "domain":
                check_error(
                    R,
                    exe,
                    sec,
                    name + " (input validation)",
                    prob,
                    s,
                    cli,
                    fail_status="WARN",
                    reason=reason,
                    hint="Newton starts at the midpoint so it never "
                    "evaluates the bad endpoint. Validate a and b "
                    "against the problem's domain in main() for "
                    "every solver",
                )
            else:
                check_breaks(
                    R,
                    exe,
                    sec,
                    name + " (plain Newton ignores brackets)",
                    prob,
                    s,
                    coeffs,
                    a,
                    b,
                    "Newton has no bracket, it may fail or find the root anyway",
                    root=root,
                )

    print("\n== 5b. Methods that are SUPPOSED to break here ==")
    # These show the weaknesses the assignment wants you to discuss.
    breaks = [
        # (name, problem, coeffs, a, b, true root, which solvers, reason)
        (
            "ln x = 0 on [0.1, 10] (Newton overshoots to x<0)",
            "LOG",
            [0],
            0.1,
            10,
            1.0,
            NEWTON_PLAIN,
            "from x0=5.05 the tangent crosses zero at x<0, outside ln's domain",
        ),
        (
            "e^x = 1 on [-50, 10] (tiny slope at start)",
            "EXP",
            [1],
            -50,
            10,
            0.0,
            NEWTON_PLAIN,
            "x0=-20 has slope 2e-9, so the step is ~5e8 and exp() overflows",
        ),
        (
            "x^3 - 3x - 3 on [-1, 3] (f'=0 at start)",
            "POLY",
            [1, 0, -3, -3],
            -1,
            3,
            2.1038034027355366,
            NEWTON_PLAIN,
            "the midpoint x0=1 is exactly where f'(x)=0",
        ),
        (
            "tan x = -3 on [-1.5, 1.5] (steep tangent overshoots)",
            "TAN",
            [-3],
            -1.5,
            1.5,
            math.atan(-3),
            NEWTON_PLAIN,
            "from x0=0 the slope is 1, so x1=-3, past -pi/2",
        ),
    ]
    for name, prob, coeffs, a, b, root, which, reason in breaks:
        for s in solvers:
            if s in which:
                check_breaks(
                    R, exe, sec, name, prob, s, coeffs, a, b, reason, root=root
                )
        # ...and show that the bracketed hybrids survive the same input
        for s in ("Bracket-Newton1", "Bracket-Newton2"):
            if s in solvers and which == NEWTON_PLAIN:
                check_root(
                    R,
                    exe,
                    sec,
                    name + "  [hybrid must survive]",
                    prob,
                    s,
                    coeffs,
                    a,
                    b,
                    root,
                    hint="Bracket-Newton must fall back to bisection "
                    "when the Newton step leaves [a,b] or doesn't "
                    "reduce |f|",
                )


# --------------------------------------------------------------------------
# Performance summary (for the Analysis section of the report)
# --------------------------------------------------------------------------


def perf_summary(R, exe, solvers, repeats):
    print("\n== 6. Performance on a common problem set ==")
    probs = [k for k in KNOWN]
    have_iters = False
    stats = {}
    for s in solvers:
        times, errs, ok, iters = [], [], 0, []
        for name, prob, coeffs, a, b, root in probs:
            cli = list(coeffs) + [a, b]
            best = run_solver(exe, prob, s, cli)
            for _ in range(repeats - 1):
                run = run_solver(exe, prob, s, cli)
                if run.secs < best.secs:
                    best = run
            if best.ok and best.value is not None:
                ok += 1
                times.append(best.secs)
                errs.append(abs(best.value - root))
                if best.iters is not None:
                    iters.append(best.iters)
                    have_iters = True
        stats[s] = (ok, len(probs), times, errs, iters)

    hdr = f"  {'solver':<16}{'solved':>8}{'median ms':>11}{'max |err|':>12}"
    if have_iters:
        hdr += f"{'avg iters':>11}"
    print(hdr)
    for s, (ok, n, times, errs, iters) in stats.items():
        med = sorted(times)[len(times) // 2] * 1e3 if times else float("nan")
        mx = max(errs) if errs else float("nan")
        line = f"  {s:<16}{ok:>4}/{n:<3}{med:>11.3f}{mx:>12.2e}"
        if have_iters:
            line += f"{(sum(iters) / len(iters) if iters else float('nan')):>11.1f}"
        print(line)
    print(
        "  (times include process start-up, which dominates; compare "
        "relative to each other)"
    )
    if not have_iters:
        print(
            "  tip: print 'Iterations: N' to stderr from each solver and "
            "this table will show iteration counts - that's the cleaner "
            "convergence-speed metric for your report."
        )


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--solver",
        default=None,
        help="path to the executable (default ./mySolver, or mySolver.exe on Windows)",
    )
    ap.add_argument(
        "--random",
        type=int,
        default=40,
        help="number of random general cases (default 40)",
    )
    ap.add_argument("--seed", type=int, default=3340)
    ap.add_argument(
        "--only",
        choices=SOLVERS,
        action="append",
        help="test only this solver (repeatable)",
    )
    ap.add_argument("--skip-perf", action="store_true")
    ap.add_argument(
        "--repeats", type=int, default=3, help="timing repeats per problem (default 3)"
    )
    ap.add_argument("--csv", help="write every result row to this CSV file")
    ap.add_argument(
        "-v", "--verbose", action="store_true", help="also print passing tests"
    )
    ap.add_argument("--no-color", action="store_true")
    opts = ap.parse_args()

    exe = opts.solver
    if exe is None:
        exe = "mySolver.exe" if os.name == "nt" else "./mySolver"
    if not os.path.exists(exe):
        sys.exit(
            f"Can't find solver executable '{exe}'. Run 'make' first or "
            f"pass --solver PATH."
        )

    solvers = opts.only or SOLVERS
    color = (not opts.no_color) and sys.stdout.isatty()
    R = Results(opts.verbose, color)

    print(f"Testing {exe}  solvers: {', '.join(solvers)}")
    test_interface(R, exe)
    test_known(R, exe, solvers)
    test_random(R, exe, solvers, opts.random, opts.seed)
    test_corners(R, exe, solvers)
    test_expected_failures(R, exe, solvers)
    if not opts.skip_perf:
        perf_summary(R, exe, solvers, opts.repeats)

    # ---- Summary --------------------------------------------------------
    print("\n" + "=" * 70)
    c = R.counts
    print(
        f"PASS {c['PASS']}   FAIL {c['FAIL']}   WARN {c['WARN']}   "
        f"BROKE-OK {c['BROKE-OK']}   NOTE {c['NOTE']}"
    )

    if R.problems:
        print("\nThings to fix, grouped by cause (most important first):")
        groups = {}
        for status, sec, solver, name, hint in R.problems:
            if sec == "interface" and "stdout" not in name:
                solver = "main()/all solvers"  # CLI handling, not a method
            key = (status, hint or name)
            g = groups.setdefault(key, {"solvers": set(), "tests": []})
            g["solvers"].add(solver)
            if name not in g["tests"]:
                g["tests"].append(name)
        order = sorted(
            groups.items(),
            key=lambda kv: (
                kv[0][0] != "FAIL",
                -len(kv[1]["tests"]) * len(kv[1]["solvers"]),
            ),
        )
        for n, ((status, cause), g) in enumerate(order, 1):
            who = sorted(
                g["solvers"], key=lambda w: SOLVERS.index(w) if w in SOLVERS else 99
            )
            tests = g["tests"]
            shown = "; ".join(tests[:3]) + (
                f"; +{len(tests) - 3} more" if len(tests) > 3 else ""
            )
            print(f"\n  {n}. [{status}] solvers: {', '.join(who)}")
            print(f"     seen in: {shown}")
            if cause != tests[0]:
                print(f"     fix: {cause}")
    else:
        print("\nNothing to fix. ")

    if opts.csv:
        with open(opts.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(R.rows[0].keys()))
            w.writeheader()
            w.writerows(R.rows)
        print(f"\nWrote {len(R.rows)} rows to {opts.csv}")

    sys.exit(1 if c["FAIL"] else 0)


if __name__ == "__main__":
    main()
