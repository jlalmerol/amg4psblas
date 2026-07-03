# Nested Block Samples

This directory contains AMG4PSBLAS examples for PSBLAS nested matrices.  The
samples read block MatrixMarket files, assemble a PSBLAS `psb_d_nest_matrix`, and
solve the resulting saddle-point systems through the usual PSBLAS Krylov
interface.

The Stokes sample uses the AMG4PSBLAS-side block preconditioner by default
(`STOKES_PREC=AMG_BLOCK`), with AMG on the velocity block and a Schur
preconditioner for pressure.  The optimal-control sample uses the PSBLAS built-in diagonal
nested composition by default (`DIAG`), which is the robust default for the
current KKT MatrixMarket data.

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
```


Useful environment variables:

```text
STOKES_DIR=/path/to/dealii-test/build
STOKES_REF=0
STOKES_METHOD=BICGSTAB
STOKES_PREC=AMG_BLOCK|NEST|AMG
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
```

Current status on the available `STOKES_REF=0` test data:

```text
default: STOKES_PREC=AMG_BLOCK STOKES_SCHUR_SOLVE=MATRIX_FREE iterations = 23    relative residual = 1.0235E-09
STOKES_PREC=NEST                                             iterations = 135   relative residual = 7.1175E-09
STOKES_PREC=AMG_BLOCK STOKES_SCHUR_SOLVE=SELF        iterations = 41    relative residual = 4.8287E-09
```

For `STOKES_PREC=AMG_BLOCK`, `SELF` and `SELFP` use the lumped pressure
Schur diagonal built from `B diag(A)^{-1} Bt`; `MATRIX_FREE` keeps the
existing inner Schur iteration using AMG applications on the velocity block.


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

Current status on the available default optimal-control data (`OPTIMAL_REF=5`,
`OPTIMAL_DEGREE=1`):

```text
OPTIMAL_PREC=NEST OPTIMAL_COMPOSITION=DIAG iterations = 847   relative residual = 9.4861E-04
```

The run uses `BICGSTAB`, total size 3267 (`n_u=n_y=n_p=1089`),
`OPTIMAL_EPS=1e-3`, and `OPTIMAL_CHECK_TOL=1e-2`.

The `OPTIMAL_PREC=AMG_BLOCK` path is implemented and exercises AMG field solves
and a matrix-free KKT Schur action, but it does not yet converge robustly for the
current optimal-control matrices.  The KKT Schur term is very strongly scaled
because `A13/A31` are order-one blocks while `A11` is a small mass-like block, so
a diagonal/probed Schur preconditioner is too weak.  A stronger Schur
preconditioner is still needed for this case.

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
