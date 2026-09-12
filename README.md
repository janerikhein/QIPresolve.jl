# QIPresolve

[![Build Status](https://github.com/janerikhein/QIPresolve.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/janerikhein/QIPresolve.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/janerikhein/QIPresolve.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/janerikhein/QIPresolve.jl)

## Presolve strategy comparison

Run all instances listed in the main benchmark directory:

```sh
julia --project=. scripts/presolve_main_benchmark_comparison.jl
```

The script compares parity alone, residue alone, the alternating combination,
SCIP alone, and the alternating combination followed by SCIP. Each independent
strategy starts from the same imported model. SCIP runs `SCIPpresolve` with its
default settings and no search; `--scip-config FILE.set` applies the same parameter
file to both SCIP strategies. Use `--help` for reduction settings and directory
options. A small run can be written separately with:

```sh
julia --project=. scripts/presolve_main_benchmark_comparison.jl --limit 2 --output-dir /tmp/qipresolve-comparison
```

Results go to `results/main_benchmark_presolve_comparison/` by default:

- `per_instance.csv`: eight reduction metrics and five strategy statuses, one row
  per instance.
- `aggregated.csv`: arithmetic means of the metrics by instance type, plus the
  instance count; status columns are omitted. This file is updated after each
  completed instance, so it also describes a partially completed run.
- `diagnostics.log`: SCIP termination statuses, search-node checks, and any
  surviving constraint rewrites whose bounds could not be compared reliably.
- `run_config.txt`: effective experiment options, Julia/SCIP versions, and the
  contents of an optional SCIP parameter file.

Instance types come from `random_instances.csv` and `embedding_instances.csv` in
the input directory. Files are processed in filename order. Re-running into the
same output directory replaces the previous results.

Domain reduction is `(L_before - L_after) / L_before`, where
`L = sum(log(ub - lb + 1))` over active variables, including helper variables.
Every strategy uses the original model as its baseline. A zero baseline gives
zero reduction; detected infeasibility uses `L_after = 0`. This is a measure of
the transformed domain box, not the number of feasible assignments. Introducing
helper variables or changing the representation can enlarge that box, so the
reported reduction is not clamped and the full strategy need not dominate the
combined strategy numerically.

Bound tightening averages contributions over the original logical constraints.
Complementary LP inequalities are rejoined into ranged constraints before taking
the baseline. A removed original constraint contributes **1.0**, including an
equality or one-sided constraint. A surviving constraint contributes
`(original_width - final_width) / original_width`; surviving equalities and
one-sided constraints contribute zero. Empty constraint collections give zero.
Values above one are preserved when tightened bounds cross. Coefficient scaling
is undone before comparing widths, and newly added helper constraints are
excluded. A replacement constraint is followed rather than counted as removed;
aggregation of a redundant original constraint counts as removal. Constraints
removed in either stage of the full strategy are counted once.

For surviving SCIP constraints, bounds are compared through affine variable
substitutions and proportional expression matching. If a rewrite cannot be
verified, the last verified interval is retained and the case is logged; these
cases can underestimate tightening. The residue experiment's range formula is
used directly, rather than the core statistics average over tightened rows.

Statuses are `infeasible` when presolving proves infeasibility, `feasible` when
all constraints disappear with consistent variable domains or SCIP reports the
problem solved, and `reduced` otherwise (including unchanged instances). Finding
an incumbent alone does not qualify. A solved status does not automatically
assign a reduction of one to surviving domains or constraints.

The public entry point also supports isolated reductions:

```julia
presolve!(model; enable_parity = true, enable_residue = false)
presolve!(model; enable_parity = false, enable_residue = true)
```

Both switches default to `true`. Shared normalization and postsolve tracking
remain enabled; setting both switches to `false` runs normalization only.
