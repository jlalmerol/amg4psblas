# Nested Block Samples

This directory contains AMG4PSBLAS examples for PSBLAS nested matrices.  The
samples read block MatrixMarket files, assemble a PSBLAS `psb_d_nest_matrix`, and
solve the resulting saddle-point systems through the usual PSBLAS Krylov
interface.

The Stokes sample uses the AMG4PSBLAS-side block preconditioner by default
(`STOKES_PREC=AMG_BLOCK`), with AMG on the velocity block and a Schur
preconditioner for pressure.  `STOKES_PREC=MASS_BLOCK` instead uses the
block-diagonal preconditioner `diag(A,Mp)`, where `Mp` is the exported pressure
mass matrix; BJAC/ILU is applied to both diagonal blocks.
`STOKES_PREC=AMG_MASS_BLOCK` implements equation (21) with a fixed symmetric
AMG cycle for `A` and the lumped diagonal inverse of `Mp`. The optimal-control
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
- `STOKES_PREC=AMG_MASS_BLOCK` for the equation-(21) block-diagonal route
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
Mp_ref<ref>.mtx (required by STOKES_PREC=MASS_BLOCK or AMG_MASS_BLOCK)
```


Useful environment variables:

```text
STOKES_DIR=/path/to/dealii-test/build
STOKES_REF=0
STOKES_METHOD=BICGSTAB
STOKES_PREC=AMG_BLOCK|AMG_MASS_BLOCK|MASS_BLOCK|NEST|AMG
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

```sh
# Equation (21): fixed AMG velocity solve and lumped pressure mass with MINRES:
STOKES_DIR=$HOME/dealii-test/build STOKES_METHOD=MINRES \
STOKES_PREC=AMG_MASS_BLOCK ./runs/amg_d_nest_stokes_file_test
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
OPTIMAL_PREC=NEST|AMG_BLOCK|AMG_KKT_DIAG|MGW_EXACT|AMG
OPTIMAL_COMPOSITION=DIAG|SCHUR_PDE_CONTROL|SCHUR_PDE_CONTROL_DIAG
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

`OPTIMAL_COMPOSITION=SCHUR_PDE_CONTROL_DIAG` applies the structurally
symmetric block-diagonal form `diag(B1,B2,BS)`.  It omits the lower and upper
coupling corrections used by `SCHUR_PDE_CONTROL`.  For MINRES, all three
solves must additionally be fixed symmetric positive-definite operators.
`OPTIMAL_SCHUR_SOLVE=A33` satisfies that requirement with the current
`A33=0`/`NONE` identity fallback.  The current `MATRIX_FREE` Schur option uses
an inner CG iteration and is not a fixed linear preconditioner, so it should
not be combined with strict MINRES.

`OPTIMAL_PREC=AMG_KKT_DIAG` selects the sample-level prototype of the shifted
SPD approximation from equation (17).  It applies

```text
diag(M^{-1}, (alpha M)^{-1}, H^{-1} M H^{-1}),
H = K + alpha^{-1/2} M.
```

The exported matrices retain constrained boundary unknowns, making `A13`
nonsymmetric by itself (`A31=A13^T`).  The prototype projects the stiffness
action onto the unconstrained degrees of freedom before adding the mass shift.
The two mass inverses use fixed diagonal actions. The prototype assembles
`H` as an ordinary distributed square PSBLAS matrix and builds one reusable
AMG hierarchy for it; the same fixed AMG action is used for both `H^{-1}`
factors. Jacobi smoothing and a fixed block-Jacobi coarse solve keep the
complete block preconditioner suitable for MINRES.

`OPTIMAL_PREC=MGW_EXACT` selects a serial reference implementation of the
ideal Murphy--Golub--Wathen block-diagonal preconditioner:

```text
P_MGW = diag(A11, A22, S),
S = A31 A11^{-1} A13 + A32 A22^{-1} A23 - A33.
```

For the distributed-control matrices, `A11=M`, `A22=alpha*M`, and
`A33=0`, so this is the exact Schur complement
`S=K M^{-1} K + alpha^{-1} M`. The route converts the blocks to dense
matrices, factors `A11`, `A22`, and `S` with LAPACK Cholesky, and applies
three fixed exact solves. It is intended to verify the ideal three-eigenvalue
MINRES result; it is deliberately rejected when more than one MPI process is
used and is not a scalable production preconditioner.

Run the exact reference case with:

```sh
OPTIMAL_METHOD=MINRES \
OPTIMAL_PREC=MGW_EXACT \
OPTIMAL_EPS=1e-3 \
OPTIMAL_CHECK_TOL=1e-2 \
./runs/amg_d_nest_optimal_control_test
```

At tighter tolerances, floating-point roundoff can require more than three
iterations even though the exact-arithmetic preconditioned operator has only
three distinct eigenvalues.

For example, the currently available fixed-SPD configuration is:

```sh
OPTIMAL_METHOD=MINRES \
OPTIMAL_PREC=NEST \
OPTIMAL_COMPOSITION=SCHUR_PDE_CONTROL_DIAG \
OPTIMAL_SCHUR_SOLVE=A33 \
OPTIMAL_PROFILE=MASS_DIAG \
OPTIMAL_MASS_INNER_SOLVE=NONE \
OPTIMAL_FIELD3_BLOCK_SOLVE=NONE \
OPTIMAL_FIELD3_INNER_SOLVE=NONE \
./runs/amg_d_nest_optimal_control_test
```

The run uses `MINRES`, total size 3267 (`n_u=n_y=n_p=1089`),
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

The original benchmark rows use preconditioned BiCGSTAB unless stated
otherwise. The stopping test is based on the relative residual
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
- `AMG_MASS_BLOCK`: the equation-(21) block-diagonal preconditioner uses a
  fixed AMG V-cycle with matching Jacobi pre/post smoothing for velocity and
  the lumped diagonal inverse of `Mp` for pressure. It is currently a serial
  reference route and is run with MINRES.
- `NEST`: the built-in PSBLAS nested Schur route uses BJAC/ILU(0) for its
  block solves and the same matrix-free Schur controls.
- `NEST` with `SELF` or `SELFP`: PSBLAS explicitly assembles the approximate
  Schur matrix from the nested blocks and builds a separate BJAC/ILU(0)
  preconditioner for it.  `SELF` and `SELFP` select the same implementation.

The optimal-control comparison uses `NEST`, `OPTIMAL_EPS=1e-3`, a maximum of
4,000 outer iterations, and an independent acceptance threshold of `1e-2`.
Relevant configurations include:

- `DIAG`: independent diagonal field solves.  With the `MASS_DIAG` profile,
  fields 1 and 2 use diagonal approximations to their mass-like blocks and
  field 3 uses the `NONE` fallback.
- `SCHUR_PDE_CONTROL`: fields 1 and 2 again use the `MASS_DIAG` profile, while
  field 3 is treated through the matrix-free PDE-control Schur action.  The
  Schur solve uses at most 200 inner CG iterations.  No separate field-3 block
  or inner solver is enabled.
- `MGW_EXACT`: exact dense solves are used for `A11`, `A22`, and the exact
  Schur complement. This serial-only route is a correctness reference for the
  three-iteration result, not a scalable alternative.


### Execution setup

All results below were rerun with one MPI process (`np=1`) on the current
machine: Ubuntu 26.04 LTS under WSL2 (`x86_64`), on an AMD Ryzen 7 7735HS.
The build uses GNU Fortran 15.2.0, Open MPI 5.0.10, CMake 4.2.3, and the Ubuntu
reference BLAS/LAPACK 3.12.1 packages. The tested source revisions are
AMG4PSBLAS `97c7fae3` and PSBLAS `157d5ea6`, including the uncommitted sample
changes in this working tree.

The current input directories are:

```sh
export STOKES_DIR=/home/jenny/dealii-test/build
export OPTIMAL_DIR=/home/jenny/dealii-test/build_optimal
```

Each row is one representative run. `Prec time` is preconditioner construction
time and `Solve time` is the Krylov solve time reported by the driver; file
reading and nested-matrix assembly are excluded. Because the problems are small
and each configuration was run only once, timing differences should be treated
cautiously. Iteration counts and final residuals are the more useful comparison.

### Stokes benchmark report

| Route | Outer iterations | Relative residual | Prec time (s) | Solve time (s) |
|---|---:|---:|---:|---:|
| `AMG_BLOCK`, `MATRIX_FREE` | 23 | 1.0235e-09 | 0.004373 | 0.084733 |
| `AMG_BLOCK`, `SELF` | 41 | 4.8287e-09 | 0.004254 | 0.015657 |
| `MASS_BLOCK` | 59 | 7.8751e-09 | 0.002866 | 0.006058 |
| `NEST`, `MATRIX_FREE` | 135 | 7.1175e-09 | 0.001978 | 0.099826 |
| `NEST`, `SELF/SELFP` | 41 | 5.6304e-09 | 0.003033 | 0.006594 |
| `MINRES`, `AMG_MASS_BLOCK` | 112 | 7.3804e-09 | 0.002385 | 0.019207 |

On this Stokes case, the matrix-free AMG/Schur route requires the fewest outer
iterations.  The pressure-mass route converges in 59 iterations and is a useful
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

STOKES_METHOD=MINRES STOKES_PREC=AMG_MASS_BLOCK \
./runs/amg_d_nest_stokes_file_test
```

### Optimal-control benchmark report

All rows use one process, the refinement-5 degree-1 matrices
(`n_u=n_y=n_p=1089`), `OPTIMAL_PROFILE=MASS_DIAG`, and the stopping controls
described above. The first two rows use `BICGSTAB`; the remaining rows use
`MINRES` with structurally symmetric block-diagonal preconditioners.

| Krylov method | Route | Schur inner limit | Outer iterations | Solver error | Relative residual | Prec time (s) | Solve time (s) | Result |
|---|---|---:|---:|---:|---:|---:|---:|---|
| `BICGSTAB` | `NEST`, `DIAG` | -- | 847 | 9.4861e-04 | 9.4861e-04 | 0.000149 | 0.18018 | PASS |
| `BICGSTAB` | `NEST`, `SCHUR_PDE_CONTROL`, `MATRIX_FREE` | 200 | 4,000 | 8.7289e+02 | 8.7289e+02 | 0.000237 | 94.735 | FAIL |
| `MINRES` | `NEST`, `DIAG` | -- | 701 | 9.9986e-04 | 9.0965e-04 | 0.000129 | 0.079291 | PASS |
| `MINRES` | `NEST`, `SCHUR_PDE_CONTROL_DIAG`, `A33` | -- | 701 | 9.9986e-04 | 9.0965e-04 | 0.000154 | 0.081658 | PASS |
| `MINRES` | `AMG_KKT_DIAG`, shifted SPD/AMG | -- | 9,488 | 1.9968e-05 | 3.2601e-03 | 0.002796 | 2.5530 | PASS |
| `MINRES` | `MGW_EXACT`, exact dense Schur | -- | 3 | 4.2410e-10 | 7.7501e-10 | 3.3199 | 0.017473 | PASS |

On the current machine, the experimental PDE-control Schur configuration did
not reproduce the earlier convergent result. It reached the 4,000-iteration
limit with relative residual `8.7289e+02`; its row is therefore reported as a
failure and should be treated as a regression target. Every outer
preconditioner application can perform up to 200 matrix-free inner Schur
iterations, which accounts for the long failed-run solve time.

With the diagonal nested preconditioner, `MINRES` reduces the iteration count
from 847 to 701 and the representative solve time from 0.18018 s to
0.079291 s. The reported MINRES solver estimate (`9.9986e-04`) differs from
the independently computed true relative residual (`9.0965e-04`); both satisfy
their configured thresholds.

The new `SCHUR_PDE_CONTROL_DIAG + A33` route also converges in 701 MINRES
iterations.  Here `A33=0`, so the configured `NONE` field-3 block provides the
fixed identity fallback.  The measured read and assembly times for this run
were 0.051312 s and 0.008227 s, respectively. Its identical iteration count
and residual confirm that, with this fallback, it has the same mathematical
action as the existing diagonal baseline while exercising the new
block-diagonal Schur composition.

The shifted SPD/AMG prototype needs a tighter MINRES estimate (`OPTIMAL_EPS=2e-5`)
to meet the independently checked residual threshold.  It validates the new
operator, the reusable AMG hierarchy, and MINRES compatibility.  It is faster
than the earlier Jacobi-factor prototype but is not yet competitive with the
simpler diagonal baseline in outer iterations or total solve time.

The serial `MGW_EXACT` route terminates in three MINRES iterations at
`OPTIMAL_EPS=1e-3`, reproducing the ideal Murphy--Golub--Wathen result. Its
3.3199-second setup time is dominated by forming and factoring the dense exact
Schur complement. At `OPTIMAL_EPS=1e-12`, the same run required five
iterations and reached a true relative residual of `1.3814e-13` because of
finite-precision loss of the exact three-eigenvalue property.

To reproduce the six optimal-control rows:

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

OPTIMAL_METHOD=MINRES \
OPTIMAL_PREC=NEST \
OPTIMAL_COMPOSITION=SCHUR_PDE_CONTROL_DIAG \
OPTIMAL_SCHUR_SOLVE=A33 \
OPTIMAL_PROFILE=MASS_DIAG \
OPTIMAL_MASS_INNER_SOLVE=NONE \
OPTIMAL_FIELD3_BLOCK_SOLVE=NONE \
OPTIMAL_FIELD3_INNER_SOLVE=NONE \
./runs/amg_d_nest_optimal_control_test

OPTIMAL_METHOD=MINRES \
OPTIMAL_PREC=AMG_KKT_DIAG \
OPTIMAL_EPS=2e-5 \
OPTIMAL_CHECK_TOL=1e-2 \
OPTIMAL_ITMAX=10000 \
./runs/amg_d_nest_optimal_control_test
```

The exact reference row is reproduced separately with one process:

```sh
OPTIMAL_METHOD=MINRES \
OPTIMAL_PREC=MGW_EXACT \
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
AMG4PSBLAS packages. On the current WSL machine, configure and build it with:

```sh
cmake -S . -B build \
  -DAMG4PSBLAS_INSTALL_DIR=/home/jenny/amg4psblas-install \
  -DPSBLAS_INSTALL_DIR=/home/jenny/psblas3-install
cmake --build build -j
```

CMake writes the executables under `build/runs/`. Run them from this source
directory so that the documented relative input defaults resolve correctly, or
set the absolute input directories shown in the execution setup:

```sh
STOKES_DIR=/home/jenny/dealii-test/build ./build/runs/amg_d_nest_stokes_file_test
OPTIMAL_DIR=/home/jenny/dealii-test/build_optimal \
  ./build/runs/amg_d_nest_optimal_control_test
```

## Automated regression checks

The standalone CMake build registers analytic Stokes/KKT regressions on 1, 2, 3,
and 4 MPI ranks, with SELF and MATRIX_FREE Schur solves and unit/tiny RHS.
It also checks independently assembled shifted operators. No external files
are required for these checks.

Set `STOKES_TEST_DATA_DIR` and `OPTIMAL_TEST_DATA_DIR` during configuration to
register external examples: Stokes AMG routes, optimal-control NEST/DIAG with
MINRES, and the serial exact reference. Run:

```sh
ctest --test-dir build --output-on-failure
```

`NEST_TEST_RANKS` controls MPI sizes. CTest locks file runs because they write
solutions beside the inputs. Readers support real/general coordinate matrices
and real/general single-column array vectors; other storage headers are rejected.

Schur diagonal probing visits global columns collectively: it is an MPI
reference implementation and is expensive for large systems. Analytic tests
do not establish convergence of every inexact Schur/outer-solver combination.
Experimental optimal-control matrix-free factorizations and shifted-AMG must
still pass the independent true-residual check on each problem and MPI size.
