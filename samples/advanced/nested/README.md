# Nested Block Samples

This directory contains AMG4PSBLAS examples for PSBLAS nested matrices.  The
samples read block MatrixMarket files, assemble a PSBLAS `psb_d_nest_matrix`, and
solve the resulting saddle-point systems through the usual PSBLAS Krylov
interface.

The Stokes sample uses the AMG4PSBLAS-side block preconditioner by default
(`STOKES_PREC=AMG_BLOCK`), with AMG on the velocity block and a Schur
preconditioner for pressure.  `STOKES_PREC=MASS_BLOCK` instead uses the
block-diagonal preconditioner `diag(A,Mp)`, where `Mp` is the exported pressure
mass matrix; BJAC/ILU is applied to both diagonal blocks.  The optimal-control
sample uses the PSBLAS built-in diagonal nested composition by default (`DIAG`),
which is the robust default for the current KKT MatrixMarket data.

The additional `amg_d_nested_block_prec_mod.F90` file provides an
AMG4PSBLAS-side block preconditioner for comparison runs.  This keeps PSBLAS
independent of AMG4PSBLAS while still allowing nested/block tests to use AMG on
selected field blocks.

## Added Preconditioner

`amg_d_nested_block_prec_type` extends `psb_dprec_type` and is selected by the
sample drivers with:

- `STOKES_PREC=AMG_BLOCK`, `AMG_SCHUR`, or `SCHUR_AMG`
- `OPTIMAL_PREC=AMG_BLOCK`, `AMG_KKT`, or `KKT_AMG`

The preconditioner supports two explicit block structures:

- Stokes 2x2 Schur form

  ```text
  [ A   Bt ]
  [ B    0 ]
  ```

  AMG is built for the velocity block `A`.  The pressure Schur action is applied
  matrix-free as `-B A^{-1} Bt`, with AMG used for the `A^{-1}` action.  The
  Schur solve uses a small inner BiCGSTAB iteration with a probed diagonal Schur
  approximation as preconditioner.

- Optimal-control 3x3 KKT form

  ```text
  [ A11   0   A13 ]
  [  0   A22  A23 ]
  [ A31  A32  A33 ]
  ```

  AMG is built for the first two diagonal field blocks, `A11` and `A22`.  The
  constraint Schur action is applied matrix-free as
  `-(A31 A11^{-1} A13 + A32 A22^{-1} A23)`, again using AMG for the field solves.

A monolithic AMG path is also available by setting `STOKES_PREC=AMG` or
`OPTIMAL_PREC=AMG`, but this is generally not appropriate for these indefinite
saddle-point/KKT systems.

## Sample Problems

### Stokes File Test

Executable:

```text
runs/amg_d_nest_stokes_file_test
```

Source:

```text
amg_d_nest_stokes_file_test.F90
```

The Stokes sample reads a 2x2 block system exported by the deal.II Stokes test.
Expected files in `STOKES_DIR` are:

```text
A_ref<ref>.mtx
Bt_ref<ref>.mtx
B_ref<ref>.mtx
rhs_u_ref<ref>.mtx
rhs_p_ref<ref>.mtx
Mp_ref<ref>.mtx (required by STOKES_PREC=MASS_BLOCK)
```


Useful environment variables:

```text
STOKES_DIR=/path/to/dealii-test/build
STOKES_REF=0
STOKES_METHOD=BICGSTAB
STOKES_PREC=AMG_BLOCK|MASS_BLOCK|NEST|AMG
STOKES_COMPOSITION=SCHUR_FULL
STOKES_SCHUR_SOLVE=MATRIX_FREE|SELF|SELFP
STOKES_SCHUR_MAXIT=10
STOKES_SCHUR_TOL=0
STOKES_ITMAX=2000
STOKES_EPS=1e-8
STOKES_RESIDUAL_TOL=1e-6
```

Example:

```sh
cd samples/advanced/nested
make all
STOKES_DIR=$HOME/dealii-test/build \
./runs/amg_d_nest_stokes_file_test

# Compare with the pressure-mass block-diagonal route:
STOKES_DIR=$HOME/dealii-test/build STOKES_PREC=MASS_BLOCK \
./runs/amg_d_nest_stokes_file_test
```

For `STOKES_PREC=AMG_BLOCK`, `SELF` and `SELFP` use the lumped pressure
Schur diagonal built from `B diag(A)^{-1} Bt`; `MATRIX_FREE` keeps the
existing inner Schur iteration using AMG applications on the velocity block.
For `STOKES_PREC=NEST`, `SELF` and `SELFP` are aliases for the PSBLAS-built
approximate Schur matrix and its internal block preconditioner.


### Optimal-Control File Test

Executable:

```text
runs/amg_d_nest_optimal_control_test
```

Source:

```text
amg_d_nest_optimal_control_test.F90
```

The optimal-control sample reads a 3x3 KKT system exported by the deal.II
optimal-control test.  Expected files in `OPTIMAL_DIR` are named with the suffix
`_refinements<ref>_degree<degree>.mtx`, for example:

```text
A11_refinements5_degree1.mtx
A13_refinements5_degree1.mtx
A22_refinements5_degree1.mtx
A23_refinements5_degree1.mtx
A31_refinements5_degree1.mtx
A32_refinements5_degree1.mtx
A33_refinements5_degree1.mtx
b1_refinements5_degree1.mtx
b2_refinements5_degree1.mtx
b3_refinements5_degree1.mtx
```

Useful environment variables:

```text
OPTIMAL_DIR=/path/to/dealii-test/build_optimal
OPTIMAL_REF=5
OPTIMAL_DEGREE=1
OPTIMAL_METHOD=BICGSTAB
OPTIMAL_PREC=NEST|AMG_BLOCK|AMG
OPTIMAL_COMPOSITION=DIAG|SCHUR_PDE_CONTROL
OPTIMAL_SCHUR_SOLVE=MATRIX_FREE|A33
OPTIMAL_SCHUR_MAXIT=12
OPTIMAL_SCHUR_TOL=0
OPTIMAL_FIELD3_BLOCK_SOLVE=NONE|DIAG|BJAC
OPTIMAL_ITMAX=4000
OPTIMAL_EPS=1e-3
OPTIMAL_CHECK_TOL=1e-2
OPTIMAL_PROFILE=MASS_DIAG
```

Example:

```sh
cd samples/advanced/nested
make all
OPTIMAL_DIR=$HOME/dealii-test/build_optimal \
OPTIMAL_PREC=NEST \
./runs/amg_d_nest_optimal_control_test
```

The `OPTIMAL_PREC=NEST` default uses PSBLAS nested diagonal field solves.
Setting `OPTIMAL_COMPOSITION=SCHUR_PDE_CONTROL` switches to the 3-field
PDE-control Schur complement on field 3:

```text
[ A11   0   A13 ]       grouped diagonal block: diag(A11,A22)
[  0   A22  A23 ]       Schur field: field 3
[ A31  A32  A33 ]
```

The run uses `BICGSTAB`, total size 3267 (`n_u=n_y=n_p=1089`),
`OPTIMAL_EPS=1e-3`, and `OPTIMAL_CHECK_TOL=1e-2`.

The `OPTIMAL_PREC=AMG_BLOCK` path is implemented and exercises AMG field solves
and a matrix-free KKT Schur action, but it does not yet converge robustly for the
current optimal-control matrices.  The KKT Schur term is very strongly scaled
because `A13/A31` are order-one blocks while `A11` is a small mass-like block, so
a diagonal/probed Schur preconditioner is too weak.  A stronger Schur
preconditioner is still needed for this case.

## Benchmark Notes

### Purpose and test data

These benchmarks are functional and algorithmic comparisons of the nested
preconditioner routes.  They are intended to expose convergence behavior and
the relative cost of the current implementations, not to establish
machine-independent performance rankings.

The Stokes benchmark uses the deal.II export at `STOKES_REF=0`.  Its block
operator is

```text
[ A  Bt ] [u] = [f]
[ B   0 ] [p]   [g]
```

with 594 velocity unknowns, 85 pressure unknowns, and 679 unknowns in total.
The MatrixMarket blocks have the following dimensions and stored nonzeros:

| Block | Dimensions | Nonzeros |
|---|---:|---:|
| `A` | 594 x 594 | 14,514 |
| `B` | 85 x 594 | 3,078 |
| `Bt` | 594 x 85 | 3,078 |
| `Mp` | 85 x 85 | 637 |

`Mp` is the pressure mass matrix exported from the deal.II block
preconditioner matrix.  It is not inserted into the physical `(2,2)` block;
for `MASS_BLOCK` it is used only in the block-diagonal preconditioner
`diag(A,Mp)`.

The optimal-control benchmark uses `OPTIMAL_REF=5` and `OPTIMAL_DEGREE=1`.
It has three equally sized fields, `n_u=n_y=n_p=1089`, for 3,267 unknowns in
total.

### Algorithmic choices

All reported runs use preconditioned BiCGSTAB through the PSBLAS Krylov
interface.  The stopping test is based on the relative residual
`||b-Ax||_2/||b||_2`.

The Stokes comparison uses `STOKES_EPS=1e-8`, a maximum of 2,000 outer
iterations, and an independent post-solve acceptance threshold of `1e-6`.
The tested preconditioners are:

- `AMG_BLOCK` with `MATRIX_FREE`: AMG approximates the velocity inverse.  The
  pressure Schur action `-B A^{-1} Bt` is evaluated matrix-free and solved by
  an inner BiCGSTAB iteration, with at most 10 inner iterations.
- `AMG_BLOCK` with `SELF`: AMG is again used for velocity, while the pressure
  solve uses the lumped approximation formed from `B diag(A)^{-1} Bt`.
- `MASS_BLOCK`: the separate preconditioning operator is `diag(A,Mp)`.
  PSBLAS nested diagonal composition applies BJAC/ILU(0) to both blocks.
- `NEST`: the built-in PSBLAS nested Schur route uses BJAC/ILU(0) for its
  block solves and the same matrix-free Schur controls.
- `NEST` with `SELF` or `SELFP`: PSBLAS explicitly assembles the approximate
  Schur matrix from the nested blocks and builds a separate BJAC/ILU(0)
  preconditioner for it.  `SELF` and `SELFP` select the same implementation.

The optimal-control comparison uses `NEST`, `OPTIMAL_EPS=1e-3`, a maximum of
4,000 outer iterations, and an independent acceptance threshold of `1e-2`.
The two configurations are:

- `DIAG`: independent diagonal field solves.  With the `MASS_DIAG` profile,
  fields 1 and 2 use diagonal approximations to their mass-like blocks and
  field 3 uses the `NONE` fallback.
- `SCHUR_PDE_CONTROL`: fields 1 and 2 again use the `MASS_DIAG` profile, while
  field 3 is treated through the matrix-free PDE-control Schur action.  The
  Schur solve uses at most 200 inner CG iterations.  No separate field-3 block
  or inner solver is enabled.

### Execution setup

The main results below were collected on one Apple arm64 process (`np=1`) on
macOS 26.5.1, using GNU Fortran 15.2.0, Open MPI 5.0.9, PSBLAS, AMG4PSBLAS, and
OpenBLAS.  Each row is one representative run.  `Prec time` is preconditioner
construction time and `Solve time` is the Krylov solve time reported by the
driver; file reading and nested-matrix assembly are excluded.  Because the
problems are small and each configuration was run only once, timing differences
should be treated cautiously.  Iteration counts and final residuals are the
more useful comparison here.

### Stokes benchmark report

| Route | Outer iterations | Relative residual | Prec time (s) | Solve time (s) |
|---|---:|---:|---:|---:|
| `AMG_BLOCK`, `MATRIX_FREE` | 22 | 6.6543e-09 | 0.002422 | 0.074329 |
| `AMG_BLOCK`, `SELF` | 40 | 4.2237e-09 | 0.004471 | 0.008340 |
| `MASS_BLOCK` | 56 | 8.3966e-09 | 0.000901 | 0.003952 |
| `NEST`, `MATRIX_FREE` | 137 | 5.2265e-09 | 0.001289 | 0.087068 |
| `NEST`, `SELF/SELFP` | 40 | 8.3981e-09 | 0.001154 | 0.004906 |

On this Stokes case, the matrix-free AMG/Schur route requires the fewest outer
iterations.  The pressure-mass route converges in 56 iterations and is a useful
spectrally motivated block-diagonal baseline.  Its lower time in this single
small run should not be interpreted as a scalability result: the matrix-free
route performs nested iterative work per outer iteration, while `MASS_BLOCK`
uses local ILU(0) block applications.

The direct PSBLAS `NEST + SELF/SELFP` route was also checked after moving the
internally assembled Schur matrix and block preconditioner into each nested
preconditioner instance and replaying the configured Schur options before its
build.

To reproduce the principal Stokes benchmark routes:

```sh
STOKES_PREC=AMG_BLOCK STOKES_SCHUR_SOLVE=MATRIX_FREE \
./runs/amg_d_nest_stokes_file_test

STOKES_PREC=MASS_BLOCK \
./runs/amg_d_nest_stokes_file_test

STOKES_PREC=NEST STOKES_SCHUR_SOLVE=SELFP \
./runs/amg_d_nest_stokes_file_test
```

### Optimal-control benchmark report

All rows use one process, the refinement-5 degree-1 matrices
(`n_u=n_y=n_p=1089`), `OPTIMAL_PROFILE=MASS_DIAG`, and the stopping controls
described above.  The first two rows use `BICGSTAB`; the third uses `MINRES`
with the symmetric diagonal nested preconditioner.

| Krylov method | Route | Schur inner limit | Outer iterations | Solver error | Relative residual | Prec time (s) | Solve time (s) | Result |
|---|---|---:|---:|---:|---:|---:|---:|---|
| `BICGSTAB` | `NEST`, `DIAG` | -- | 1,134 | 9.8249e-04 | 9.8249e-04 | 0.000143 | 0.18755 | PASS |
| `BICGSTAB` | `NEST`, `SCHUR_PDE_CONTROL`, `MATRIX_FREE` | 200 | 260 | 7.7616e-04 | 7.7616e-04 | 0.000046 | 4.3951 | PASS |
| `MINRES` | `NEST`, `DIAG` | -- | 705 | 9.9964e-04 | 9.0928e-04 | 0.000400 | 0.060146 | PASS |

The PDE-control Schur configuration reduces the outer iteration count from
1,134 to 260.  Its solve time is higher for this small case because every outer
preconditioner application can perform up to 200 matrix-free inner Schur
iterations; the iteration reduction therefore does not imply a reduction in
wall-clock time.

With the diagonal nested preconditioner, `MINRES` reduces the iteration count
from 1,134 to 705 and the representative solve time from 0.18755 s to
0.060146 s.  The reported MINRES solver estimate (`9.9964e-04`) differs from
the independently computed true relative residual (`9.0928e-04`); both satisfy
their configured thresholds.

To reproduce the three optimal-control rows:

```sh
OPTIMAL_PREC=NEST \
OPTIMAL_COMPOSITION=DIAG \
OPTIMAL_PROFILE=MASS_DIAG \
OPTIMAL_EPS=1e-3 \
OPTIMAL_CHECK_TOL=1e-2 \
./runs/amg_d_nest_optimal_control_test

OPTIMAL_PREC=NEST \
OPTIMAL_COMPOSITION=SCHUR_PDE_CONTROL \
OPTIMAL_PROFILE=MASS_DIAG \
OPTIMAL_MASS_INNER_SOLVE=NONE \
OPTIMAL_FIELD3_BLOCK_SOLVE=NONE \
OPTIMAL_FIELD3_INNER_SOLVE=NONE \
OPTIMAL_SCHUR_SOLVE=MATRIX_FREE \
OPTIMAL_SCHUR_MAXIT=200 \
OPTIMAL_SCHUR_TOL=0 \
OPTIMAL_EPS=1e-3 \
OPTIMAL_CHECK_TOL=1e-2 \
./runs/amg_d_nest_optimal_control_test

OPTIMAL_METHOD=MINRES \
OPTIMAL_PREC=NEST \
OPTIMAL_COMPOSITION=DIAG \
OPTIMAL_PROFILE=MASS_DIAG \
OPTIMAL_MASS_INNER_SOLVE=NONE \
OPTIMAL_FIELD3_BLOCK_SOLVE=NONE \
OPTIMAL_FIELD3_INNER_SOLVE=NONE \
OPTIMAL_EPS=1e-3 \
OPTIMAL_CHECK_TOL=1e-2 \
./runs/amg_d_nest_optimal_control_test
```

## Build

From this directory:

```sh
make all
```

The executables are written under `runs/`.

To clean generated objects, modules, and executables:

```sh
make clean
```

A standalone CMake file is also provided for builds against installed PSBLAS and
AMG4PSBLAS packages.
