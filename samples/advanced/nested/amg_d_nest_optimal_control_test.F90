!
! File: psb_d_nest_optimal_control_test.F90
!
! Read the 3x3 block MatrixMarket files exported by
! dealii-test/optimal_control.cpp, assemble the KKT operator as a PSBLAS
! nested matrix, and solve it through the standard PSBLAS Krylov interface.
!
! Expected input files in OPTIMAL_DIR, for OPTIMAL_REF/OPTIMAL_DEGREE:
!   A11_*.mtx, A13_*.mtx, A22_*.mtx, A23_*.mtx,
!   A31_*.mtx, A32_*.mtx, A33_*.mtx, b1_*.mtx, b2_*.mtx, b3_*.mtx
!
program amg_d_nest_optimal_control_test
  use psb_base_mod
  use psb_prec_mod
  use psb_linsolve_mod
  use psb_d_nest_mod
  use amg_d_prec_mod, only : amg_dprec_type
  use amg_d_nested_block_prec_mod, only : amg_d_nested_block_prec_type
  implicit none

  type mm_matrix
    integer(psb_lpk_) :: nrow = 0_psb_lpk_
    integer(psb_lpk_) :: ncol = 0_psb_lpk_
    integer(psb_lpk_) :: nnz  = 0_psb_lpk_
    integer(psb_lpk_), allocatable :: row(:), col(:)
    real(psb_dpk_),    allocatable :: val(:)
  end type mm_matrix

  type mm_vector
    integer(psb_lpk_) :: n = 0_psb_lpk_
    real(psb_dpk_), allocatable :: val(:)
  end type mm_vector

  type(psb_ctxt_type)     :: context
  type(psb_d_nest_matrix), target :: nested_matrix
  class(psb_dprec_type), allocatable :: preconditioner
  type(psb_d_vect_type)   :: rhs, x_solution, residual

  type(mm_matrix) :: a11, a13, a22, a23, a31, a32, a33
  type(mm_vector) :: b1, b2, b3

  integer(psb_ipk_) :: my_rank, num_procs, info
  integer(psb_ipk_) :: max_iter, trace_level, stop_criterion, n_iter, n_insert
  integer(psb_lpk_) :: n_u, n_y, n_p, n_total
  integer(psb_lpk_), allocatable :: rows(:), cols(:)
  real(psb_dpk_),    allocatable :: vals(:), x_global(:)
  real(psb_dpk_) :: eps, check_tol, final_residual, residual_norm, rhs_norm, inner_tol, schur_tol
  real(psb_dpk_) :: t0, t_read, t_assemble, t_prec, t_solve
  character(len=256) :: input_dir, method, ptype, composition, block_solve, sub_solve
  character(len=256) :: optimal_profile, mass_inner_solve, field3_inner_solve, field3_block_solve
  character(len=256) :: schur_solve
  character(len=64)  :: suffix
  integer :: refinement, degree, sub_fillin, mass_inner_maxit, field3_inner_maxit, schur_maxit

  call psb_init(context)
  call psb_info(context, my_rank, num_procs)

  call get_optimal_dir(input_dir)
  call get_string_env('OPTIMAL_METHOD', method, 'BICGSTAB')
  call get_string_env('OPTIMAL_PREC', ptype, 'NEST')
  call get_string_env('OPTIMAL_COMPOSITION', composition, 'DIAG')
  call get_string_env('OPTIMAL_BLOCK_SOLVE', block_solve, 'DIAG')
  call get_string_env('OPTIMAL_SUB_SOLVE', sub_solve, 'ILU')
  call get_string_env('OPTIMAL_PROFILE', optimal_profile, 'MASS_DIAG')
  call get_string_env('OPTIMAL_MASS_INNER_SOLVE', mass_inner_solve, 'NONE')
  call get_string_env('OPTIMAL_FIELD3_INNER_SOLVE', field3_inner_solve, 'NONE')
  call get_string_env('OPTIMAL_FIELD3_BLOCK_SOLVE', field3_block_solve, 'NONE')
  call get_string_env('OPTIMAL_SCHUR_SOLVE', schur_solve, 'MATRIX_FREE')
  call get_int_env('OPTIMAL_REF', refinement, 5)
  call get_int_env('OPTIMAL_DEGREE', degree, 1)
  call get_int_env('OPTIMAL_ITMAX', max_iter, 4000)
  call get_int_env('OPTIMAL_TRACE', trace_level, 0)
  call get_int_env('OPTIMAL_SUB_FILLIN', sub_fillin, 0)
  call get_int_env('OPTIMAL_MASS_INNER_MAXIT', mass_inner_maxit, 8)
  call get_int_env('OPTIMAL_FIELD3_INNER_MAXIT', field3_inner_maxit, 12)
  call get_int_env('OPTIMAL_SCHUR_MAXIT', schur_maxit, 12)
  call get_real_env('OPTIMAL_EPS', eps, 1.0e-3_psb_dpk_)
  call get_real_env('OPTIMAL_CHECK_TOL', check_tol, 1.0e-2_psb_dpk_)
  call get_real_env('OPTIMAL_INNER_TOL', inner_tol, 1.0e-2_psb_dpk_)
  call get_real_env('OPTIMAL_SCHUR_TOL', schur_tol, 0.0_psb_dpk_)
  stop_criterion = 2
  write(suffix, '(a,i0,a,i0,a)') '_refinements', refinement, '_degree', degree, '.mtx'

  if (my_rank == 0) then
    write(*,'(a)') 'Optimal-control nested MatrixMarket test'
    write(*,'(a,a)') '  input dir  : ', trim(input_dir)
    write(*,'(a,a)') '  suffix     : ', trim(suffix)
    write(*,'(a,a)') '  method     : ', trim(method)
    write(*,'(a,a)') '  prec       : ', trim(ptype)
    write(*,'(a,a)') '  composition: ', trim(composition)
    write(*,'(a,a)') '  block solve: ', trim(block_solve)
    write(*,'(a,a)') '  opt profile : ', trim(optimal_profile)
    write(*,'(a,a)') '  sub solve  : ', trim(sub_solve)
    write(*,'(a,i0)') '  sub fillin : ', sub_fillin
    write(*,'(a,a,a,i0)') '  mass inner : ', trim(mass_inner_solve), ' maxit=', mass_inner_maxit
    write(*,'(a,a)') '  fld3 block : ', trim(field3_block_solve)
    write(*,'(a,a,a,i0)') '  fld3 inner : ', trim(field3_inner_solve), ' maxit=', field3_inner_maxit
    write(*,'(a,es12.4)') '  inner tol  : ', inner_tol
    write(*,'(a,a)') '  schur solve: ', trim(schur_solve)
    write(*,'(a,i0)') '  schur maxit: ', schur_maxit
    write(*,'(a,es12.4)') '  schur tol  : ', schur_tol
    write(*,'(a,es12.4)') '  tolerance  : ', eps
    write(*,'(a,es12.4)') '  check tol  : ', check_tol
  end if

  t0 = psb_wtime()
  call read_mm_matrix(path_join(input_dir, 'A11' // trim(suffix)), a11)
  call read_mm_matrix(path_join(input_dir, 'A13' // trim(suffix)), a13)
  call read_mm_matrix(path_join(input_dir, 'A22' // trim(suffix)), a22)
  call read_mm_matrix(path_join(input_dir, 'A23' // trim(suffix)), a23)
  call read_mm_matrix(path_join(input_dir, 'A31' // trim(suffix)), a31)
  call read_mm_matrix(path_join(input_dir, 'A32' // trim(suffix)), a32)
  call read_mm_matrix(path_join(input_dir, 'A33' // trim(suffix)), a33)
  call read_mm_vector(path_join(input_dir, 'b1' // trim(suffix)), b1)
  call read_mm_vector(path_join(input_dir, 'b2' // trim(suffix)), b2)
  call read_mm_vector(path_join(input_dir, 'b3' // trim(suffix)), b3)
  t_read = psb_wtime() - t0

  n_u = a11%nrow
  n_y = a22%nrow
  n_p = a33%nrow
  n_total = n_u + n_y + n_p

  if ((a11%ncol /= n_u) .or. (a13%nrow /= n_u) .or. (a13%ncol /= n_p) .or. &
      (a22%ncol /= n_y) .or. (a23%nrow /= n_y) .or. (a23%ncol /= n_p) .or. &
      (a31%nrow /= n_p) .or. (a31%ncol /= n_u) .or. &
      (a32%nrow /= n_p) .or. (a32%ncol /= n_y) .or. &
      (a33%ncol /= n_p) .or. (b1%n /= n_u) .or. (b2%n /= n_y) .or. (b3%n /= n_p)) then
    if (my_rank == 0) write(*,*) 'FAIL: inconsistent optimal-control block dimensions'
    call psb_abort(context)
  end if

  t0 = psb_wtime()
  call nested_matrix%init(context, [n_u, n_y, n_p], info)
  call check_info(info, 'nested_matrix%init')

  call insert_block(1, 1, a11, 'A11')
  call insert_block(1, 3, a13, 'A13')
  call insert_block(2, 2, a22, 'A22')
  call insert_block(2, 3, a23, 'A23')
  call insert_block(3, 1, a31, 'A31')
  call insert_block(3, 2, a32, 'A32')
  call insert_block(3, 3, a33, 'A33')

  call nested_matrix%asb(info)
  call check_info(info, 'nested_matrix%asb')
  t_assemble = psb_wtime() - t0

  call psb_geall(rhs, nested_matrix%desc_glob, info)
  call check_info(info, 'psb_geall(rhs)')
  call insert_owned_rhs(b1, b2, b3, n_u, n_y, nested_matrix%get_owned_rows(1), &
       & nested_matrix%get_owned_rows(2), nested_matrix%get_owned_rows(3), rhs, &
       & nested_matrix%desc_glob, info)
  call check_info(info, 'insert rhs')
  call psb_geasb(rhs, nested_matrix%desc_glob, info, dupl=psb_dupl_add_)
  call check_info(info, 'psb_geasb(rhs)')

  call psb_geall(x_solution, nested_matrix%desc_glob, info)
  call check_info(info, 'psb_geall(x)')
  call x_solution%zero()
  call psb_geasb(x_solution, nested_matrix%desc_glob, info)
  call check_info(info, 'psb_geasb(x)')

  t0 = psb_wtime()
  call allocate_preconditioner()
  call check_info(info, 'prec%init')
  if (psb_toupper(trim(ptype)) == 'NEST') then
    call preconditioner%set('COMPOSITION', trim(composition), info)
    call check_info(info, 'prec%set(COMPOSITION)')
    call preconditioner%set('SCHUR_SOLVE', trim(schur_solve), info)
    call check_info(info, 'prec%set(SCHUR_SOLVE)')
    call preconditioner%set('SCHUR_MAXIT', schur_maxit, info)
    call check_info(info, 'prec%set(SCHUR_MAXIT)')
    call preconditioner%set('SCHUR_TOL', schur_tol, info)
    call check_info(info, 'prec%set(SCHUR_TOL)')
    select case (psb_toupper(trim(optimal_profile)))
    case ('MASS_DIAG','POISSON_CONTROL','OPTIMAL_CONTROL','OPTIMAL')
      call configure_optimal_profile()
    case default
      call preconditioner%set('BLOCK_SOLVE', trim(block_solve), info)
      call check_info(info, 'prec%set(BLOCK_SOLVE)')
      call preconditioner%set('SUB_SOLVE', trim(sub_solve), info)
      call check_info(info, 'prec%set(SUB_SOLVE)')
      call preconditioner%set('SUB_FILLIN', sub_fillin, info)
      call check_info(info, 'prec%set(SUB_FILLIN)')
    end select
  end if
  if (my_rank == 0) then
    write(*,'(a)') '  resolved operator composition:'
    write(*,'(a)') '    row 1: [A11,   0, A13]'
    write(*,'(a)') '    row 2: [  0, A22, A23]'
    write(*,'(a)') '    row 3: [A31, A32, A33]'
    write(*,'(a,a)') '  resolved preconditioner: ', trim(ptype)
    select case (psb_toupper(trim(ptype)))
    case ('NEST')
      write(*,'(a,a)') '    composition: ', trim(composition)
      select case (psb_toupper(trim(optimal_profile)))
      case ('MASS_DIAG','POISSON_CONTROL','OPTIMAL_CONTROL','OPTIMAL')
        write(*,'(a,a)') '    block 1 (A11): DIAG; inner=', trim(mass_inner_solve)
        if (psb_toupper(trim(mass_inner_solve)) /= 'NONE') &
             & write(*,'(a,i0,a,es12.4)') '      maxit=', mass_inner_maxit, ', tol=', inner_tol
        write(*,'(a,a)') '    block 2 (A22): DIAG; inner=', trim(mass_inner_solve)
        if (psb_toupper(trim(mass_inner_solve)) /= 'NONE') &
             & write(*,'(a,i0,a,es12.4)') '      maxit=', mass_inner_maxit, ', tol=', inner_tol
        write(*,'(a,a,a,a)') '    block 3 (A33): ', trim(field3_block_solve), &
             & '; inner=', trim(field3_inner_solve)
        if (psb_toupper(trim(field3_block_solve)) == 'BJAC') &
             & write(*,'(a,a,a,i0)') '      subsolve=', trim(sub_solve), ', fillin=', sub_fillin
        if (psb_toupper(trim(field3_inner_solve)) /= 'NONE') &
             & write(*,'(a,i0,a,es12.4)') '      maxit=', field3_inner_maxit, ', tol=', inner_tol
      case default
        write(*,'(a,a)') '    blocks 1-3: ', trim(block_solve)
        if (psb_toupper(trim(block_solve)) == 'BJAC') &
             & write(*,'(a,a,a,i0)') '      subsolve=', trim(sub_solve), ', fillin=', sub_fillin
      end select
      if (index(psb_toupper(trim(composition)), 'SCHUR') == 1) &
           & write(*,'(a,a,a,i0,a,es12.4)') '    Schur: ', trim(schur_solve), &
           & ', maxit=', schur_maxit, ', tol=', schur_tol
    case ('AMG_KKT')
      write(*,'(a)') '    composition: KKT block factorization'
      write(*,'(a)') '    blocks 1 and 2: AMG (ML)'
      write(*,'(a,a)') '    block 3: KKT Schur solve ', trim(schur_solve)
    case ('AMG_KKT_DIAG')
      write(*,'(a)') '    composition: KKT block diagonal'
      write(*,'(a)') '    blocks 1 and 2: inverse mass diagonals'
      write(*,'(a)') '    block 3: AMG (ML) on dominant Schur approximation'
    case ('MGW_EXACT')
      write(*,'(a)') '    composition: exact Murphy-Golub-Wathen block diagonal'
      write(*,'(a)') '    block 1: exact Cholesky solve with A11=M'
      write(*,'(a)') '    block 2: exact Cholesky solve with A22=alpha*M'
      write(*,'(a)') '    block 3: exact Cholesky solve with the dense Schur complement'
    case ('ML')
      write(*,'(a)') '    all blocks: global AMG (ML)'
    end select
  end if
  call preconditioner%build(nested_matrix%a_glob, nested_matrix%desc_glob, info)
  call check_info(info, 'prec%build')
  t_prec = psb_wtime() - t0

  t0 = psb_wtime()
  call psb_krylov(trim(method), nested_matrix%a_glob, preconditioner, rhs, x_solution, eps, &
       & nested_matrix%desc_glob, info, itmax=max_iter, iter=n_iter, err=final_residual, &
       & itrace=trace_level, istop=stop_criterion)
  call check_info(info, 'psb_krylov')
  t_solve = psb_wtime() - t0

  call psb_geall(residual, nested_matrix%desc_glob, info)
  call check_info(info, 'psb_geall')
  call psb_geasb(residual, nested_matrix%desc_glob, info)
  call check_info(info, 'psb_geasb')
  call psb_geaxpby(done, rhs, dzero, residual, nested_matrix%desc_glob, info)
  call check_info(info, 'psb_geaxpby')
  call psb_spmm(-done, nested_matrix%a_glob, x_solution, done, residual, nested_matrix%desc_glob, info)
  call check_info(info, 'residual spmm')
  residual_norm = psb_genrm2(residual, nested_matrix%desc_glob, info)
  call check_info(info, 'norm')
  rhs_norm      = psb_genrm2(rhs,      nested_matrix%desc_glob, info)
  call check_info(info, 'norm')

  allocate(x_global(n_total))
  call psb_gather(x_global, x_solution, nested_matrix%desc_glob, info, root=psb_root_)
  call check_info(info, 'psb_gather(solution)')

  if (my_rank == 0) then
    call write_mm_vector(path_join(input_dir, 'psblas_solution_u' // trim(suffix)), x_global(1:n_u))
    call write_mm_vector(path_join(input_dir, 'psblas_solution_y' // trim(suffix)), x_global(n_u+1:n_u+n_y))
    call write_mm_vector(path_join(input_dir, 'psblas_solution_p' // trim(suffix)), x_global(n_u+n_y+1:n_total))

    write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') '  block sizes: n_u=', n_u, ' n_y=', n_y, &
         & ' n_p=', n_p, ' total=', n_total, ' np=', num_procs
    write(*,'(a,i8,a,es12.4)') '  iterations=', n_iter, '  solver err=', final_residual
    write(*,'(a,es12.4)') '  ||b-Ax||_2 / ||b||_2=', residual_norm / max(rhs_norm, tiny(done))
    write(*,'(a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
         & '  times read=', t_read, ' assemble=', t_assemble, ' prec=', t_prec, ' solve=', t_solve
    if (residual_norm / max(rhs_norm, tiny(done)) > check_tol) then
      write(*,'(a)') '[FAIL] amg_d_nest_optimal_control_test: residual above OPTIMAL_CHECK_TOL'
      call psb_abort(context)
    end if
    write(*,'(a)') '[PASS] amg_d_nest_optimal_control_test'
  end if

9999 continue
  if (allocated(preconditioner)) call preconditioner%free(info)
  call psb_gefree(residual, nested_matrix%desc_glob, info)
  call check_info(info, 'psb_gefree')
  call psb_gefree(x_solution, nested_matrix%desc_glob, info)
  call check_info(info, 'psb_gefree')
  call psb_gefree(rhs, nested_matrix%desc_glob, info)
  call check_info(info, 'psb_gefree')
  call nested_matrix%free(info)
  call check_info(info, 'nested_matrix%free')
  call psb_exit(context)

contains

  subroutine allocate_preconditioner()
    select case (psb_toupper(trim(ptype)))
    case ('AMG_KKT_DIAG','KKT_DIAG','KKT_BLOCK_DIAG','OPTIMAL_DIAG')
      allocate(amg_d_nested_block_prec_type :: preconditioner, stat=info)
      call check_info(info, 'allocate AMG KKT block-diagonal preconditioner')
      ptype = 'AMG_KKT_DIAG'
    case ('MGW_EXACT')
      allocate(amg_d_nested_block_prec_type :: preconditioner, stat=info)
      call check_info(info, 'allocate exact MGW preconditioner')
      ptype = 'MGW_EXACT'
    case ('AMG_BLOCK','AMG_KKT','KKT_AMG')
      allocate(amg_d_nested_block_prec_type :: preconditioner, stat=info)
      call check_info(info, 'allocate AMG KKT preconditioner')
      ptype = 'AMG_KKT'
    case ('AMG','ML','MULTILEVEL')
      allocate(amg_dprec_type :: preconditioner, stat=info)
      call check_info(info, 'allocate AMG preconditioner')
      ptype = 'ML'
    case default
      allocate(psb_dprec_type :: preconditioner, stat=info)
      call check_info(info, 'allocate PSBLAS preconditioner')
    end select
    call preconditioner%init(context, trim(ptype), info)
  end subroutine allocate_preconditioner

  subroutine insert_block(i_block, j_block, mat, label)
    integer, intent(in) :: i_block, j_block
    type(mm_matrix), intent(in) :: mat
    character(len=*), intent(in) :: label

    call select_owned_entries(mat, nested_matrix%get_owned_rows(i_block), rows, cols, vals, n_insert)
    call nested_matrix%ins(i_block, j_block, n_insert, rows, cols, vals, info)
    call check_info(info, 'insert ' // trim(label))
    call clear_triplets(rows, cols, vals)
  end subroutine insert_block

  subroutine configure_optimal_profile()
    call preconditioner%set('BLOCK_SOLVE', 'DIAG', info, idx=1)
    call check_info(info, 'prec%set(BLOCK_SOLVE,1)')
    call preconditioner%set('BLOCK_SOLVE', 'DIAG', info, idx=2)
    call check_info(info, 'prec%set(BLOCK_SOLVE,2)')

    if (psb_toupper(trim(mass_inner_solve)) /= 'NONE') then
      call preconditioner%set('INNER_SOLVE', trim(mass_inner_solve), info, idx=1)
      call check_info(info, 'prec%set(INNER_SOLVE,1)')
      call preconditioner%set('INNER_SOLVE', trim(mass_inner_solve), info, idx=2)
      call check_info(info, 'prec%set(INNER_SOLVE,2)')
      call preconditioner%set('INNER_MAXIT', mass_inner_maxit, info, idx=1)
      call check_info(info, 'prec%set(INNER_MAXIT,1)')
      call preconditioner%set('INNER_MAXIT', mass_inner_maxit, info, idx=2)
      call check_info(info, 'prec%set(INNER_MAXIT,2)')
      call preconditioner%set('INNER_TOL', inner_tol, info, idx=1)
      call check_info(info, 'prec%set(INNER_TOL,1)')
      call preconditioner%set('INNER_TOL', inner_tol, info, idx=2)
      call check_info(info, 'prec%set(INNER_TOL,2)')
    end if

    call preconditioner%set('BLOCK_SOLVE', trim(field3_block_solve), info, idx=3)
    call check_info(info, 'prec%set(BLOCK_SOLVE,3)')
    if (psb_toupper(trim(field3_block_solve)) == 'BJAC') then
      call preconditioner%set('SUB_SOLVE', trim(sub_solve), info, idx=3)
      call check_info(info, 'prec%set(SUB_SOLVE,3)')
      call preconditioner%set('SUB_FILLIN', sub_fillin, info, idx=3)
      call check_info(info, 'prec%set(SUB_FILLIN,3)')
    end if

    if (psb_toupper(trim(field3_inner_solve)) /= 'NONE') then
      call preconditioner%set('INNER_SOLVE', trim(field3_inner_solve), info, idx=3)
      call check_info(info, 'prec%set(INNER_SOLVE,3)')
      call preconditioner%set('INNER_MAXIT', field3_inner_maxit, info, idx=3)
      call check_info(info, 'prec%set(INNER_MAXIT,3)')
      call preconditioner%set('INNER_TOL', inner_tol, info, idx=3)
      call check_info(info, 'prec%set(INNER_TOL,3)')
    end if
  end subroutine configure_optimal_profile
  subroutine check_info(info_value, label)
    integer(psb_ipk_), intent(in) :: info_value
    character(len=*), intent(in) :: label
    if (info_value /= psb_success_) then
      if (my_rank == 0) write(*,*) 'FAIL: ', trim(label), ' info=', info_value
      call psb_abort(context)
    end if
  end subroutine check_info

  function path_join(dir, name) result(path)
    character(len=*), intent(in) :: dir, name
    character(len=512) :: path
    integer :: n
    n = len_trim(dir)
    if (n > 0 .and. dir(n:n) == '/') then
      path = trim(dir) // trim(name)
    else
      path = trim(dir) // '/' // trim(name)
    end if
  end function path_join

  subroutine get_optimal_dir(value)
    character(len=*), intent(out) :: value
    integer :: status

    call get_environment_variable('OPTIMAL_DIR', value, status=status)
    if (status == 0 .and. len_trim(value) > 0) return

    value = '../../../../dealii-test/build_optimal'
  end subroutine get_optimal_dir

  subroutine get_string_env(name, value, default_value)
    character(len=*), intent(in) :: name, default_value
    character(len=*), intent(out) :: value
    integer :: status
    call get_environment_variable(name, value, status=status)
    if (status /= 0 .or. len_trim(value) == 0) value = default_value
  end subroutine get_string_env

  subroutine get_int_env(name, value, default_value)
    character(len=*), intent(in) :: name
    integer, intent(out) :: value
    integer, intent(in) :: default_value
    character(len=64) :: text
    integer :: status
    call get_environment_variable(name, text, status=status)
    if (status == 0 .and. len_trim(text) > 0) then
      read(text,*) value
    else
      value = default_value
    end if
  end subroutine get_int_env

  subroutine get_real_env(name, value, default_value)
    character(len=*), intent(in) :: name
    real(psb_dpk_), intent(out) :: value
    real(psb_dpk_), intent(in) :: default_value
    character(len=64) :: text
    integer :: status
    call get_environment_variable(name, text, status=status)
    if (status == 0 .and. len_trim(text) > 0) then
      read(text,*) value
    else
      value = default_value
    end if
  end subroutine get_real_env

  subroutine read_mm_matrix(filename, mat)
    character(len=*), intent(in) :: filename
    type(mm_matrix), intent(out) :: mat
    character(len=256) :: header, line
    integer :: unit, ios
    integer(psb_lpk_) :: i

    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Could not open matrix file: ', trim(filename)
      stop 2
    end if
    read(unit,'(a)') header
    if (trim(psb_toupper(header)) /= '%%MATRIXMARKET MATRIX COORDINATE REAL GENERAL') then
      write(*,*) '[FAIL] unsupported MatrixMarket matrix format: ', trim(header)
      call psb_abort(context)
    end if
    do
      read(unit,'(a)', iostat=ios) line
      if (ios /= 0) error stop 'Unexpected end of MatrixMarket header'
      if (line(1:1) /= '%') exit
    end do
    read(line,*) mat%nrow, mat%ncol, mat%nnz
    allocate(mat%row(mat%nnz), mat%col(mat%nnz), mat%val(mat%nnz))
    do i = 1_psb_lpk_, mat%nnz
      read(unit,*,iostat=ios) mat%row(i), mat%col(i), mat%val(i)
      if (ios /= 0) error stop 'Bad MatrixMarket coordinate row'
    end do
    close(unit)
  end subroutine read_mm_matrix

  subroutine read_mm_vector(filename, vec)
    character(len=*), intent(in) :: filename
    type(mm_vector), intent(out) :: vec
    character(len=256) :: header, line
    integer :: unit, ios
    integer(psb_lpk_) :: ncol, i

    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Could not open vector file: ', trim(filename)
      stop 2
    end if
    read(unit,'(a)') header
    if (trim(psb_toupper(header)) /= '%%MATRIXMARKET MATRIX ARRAY REAL GENERAL') then
      write(*,*) '[FAIL] unsupported MatrixMarket vector format: ', trim(header)
      call psb_abort(context)
    end if
    do
      read(unit,'(a)', iostat=ios) line
      if (ios /= 0) error stop 'Unexpected end of MatrixMarket vector header'
      if (line(1:1) /= '%') exit
    end do
    read(line,*) vec%n, ncol
    if (ncol /= 1) call psb_abort(context)
    allocate(vec%val(vec%n))
    do i = 1_psb_lpk_, vec%n
      read(unit,*,iostat=ios) vec%val(i)
      if (ios /= 0) error stop 'Bad MatrixMarket vector row'
    end do
    close(unit)
  end subroutine read_mm_vector

  subroutine write_mm_vector(filename, values)
    character(len=*), intent(in) :: filename
    real(psb_dpk_), intent(in) :: values(:)
    integer :: unit, i

    open(newunit=unit, file=trim(filename), status='replace', action='write')
    write(unit,'(a)') '%%MatrixMarket matrix array real general'
    write(unit,*) size(values), 1
    do i = 1, size(values)
      write(unit,'(es26.18)') values(i)
    end do
    close(unit)
  end subroutine write_mm_vector

  subroutine select_owned_entries(mat, owned_rows, out_rows, out_cols, out_vals, count)
    type(mm_matrix), intent(in) :: mat
    integer(psb_lpk_), intent(in) :: owned_rows(:)
    integer(psb_lpk_), allocatable, intent(out) :: out_rows(:), out_cols(:)
    real(psb_dpk_), allocatable, intent(out) :: out_vals(:)
    integer(psb_ipk_), intent(out) :: count
    logical, allocatable :: is_owned(:)
    integer(psb_lpk_) :: i
    integer(psb_ipk_) :: k

    count = 0
    if (size(owned_rows) == 0) then
      allocate(out_rows(0), out_cols(0), out_vals(0))
      return
    end if

    allocate(is_owned(mat%nrow))
    is_owned = .false.
    do i = 1_psb_lpk_, int(size(owned_rows), psb_lpk_)
      is_owned(owned_rows(i)) = .true.
    end do

    do i = 1_psb_lpk_, mat%nnz
      if (is_owned(mat%row(i))) count = count + 1
    end do

    allocate(out_rows(count), out_cols(count), out_vals(count))
    k = 0
    do i = 1_psb_lpk_, mat%nnz
      if (is_owned(mat%row(i))) then
        k = k + 1
        out_rows(k) = mat%row(i)
        out_cols(k) = mat%col(i)
        out_vals(k) = mat%val(i)
      end if
    end do
    deallocate(is_owned)
  end subroutine select_owned_entries

  subroutine insert_owned_rhs(b1_in, b2_in, b3_in, n1, n2, owned_u, owned_y, owned_p, rhs_vec, desc, info)
    type(mm_vector), intent(in) :: b1_in, b2_in, b3_in
    integer(psb_lpk_), intent(in) :: n1, n2
    integer(psb_lpk_), intent(in) :: owned_u(:), owned_y(:), owned_p(:)
    type(psb_d_vect_type), intent(inout) :: rhs_vec
    type(psb_desc_type), intent(in) :: desc
    integer(psb_ipk_), intent(out) :: info
    integer(psb_lpk_), allocatable :: idx(:)
    real(psb_dpk_), allocatable :: values(:)
    integer :: i, n

    n = size(owned_u) + size(owned_y) + size(owned_p)
    allocate(idx(n), values(n))
    do i = 1, size(owned_u)
      idx(i) = owned_u(i)
      values(i) = b1_in%val(owned_u(i))
    end do
    do i = 1, size(owned_y)
      idx(size(owned_u) + i) = n1 + owned_y(i)
      values(size(owned_u) + i) = b2_in%val(owned_y(i))
    end do
    do i = 1, size(owned_p)
      idx(size(owned_u) + size(owned_y) + i) = n1 + n2 + owned_p(i)
      values(size(owned_u) + size(owned_y) + i) = b3_in%val(owned_p(i))
    end do
    call psb_geins(n, idx, values, rhs_vec, desc, info)
    deallocate(idx, values)
  end subroutine insert_owned_rhs

  subroutine clear_triplets(in_rows, in_cols, in_vals)
    integer(psb_lpk_), allocatable, intent(inout) :: in_rows(:), in_cols(:)
    real(psb_dpk_), allocatable, intent(inout) :: in_vals(:)
    if (allocated(in_rows)) deallocate(in_rows)
    if (allocated(in_cols)) deallocate(in_cols)
    if (allocated(in_vals)) deallocate(in_vals)
  end subroutine clear_triplets

end program amg_d_nest_optimal_control_test
