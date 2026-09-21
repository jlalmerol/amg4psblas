!
! File: psb_d_nest_stokes_file_test.F90
!
! Read the block MatrixMarket files exported by dealii-test/stokes.cpp,
! assemble the Stokes operator as a PSBLAS nested matrix, and solve it with
! the standard PSBLAS Krylov/preconditioner interface.
!
! Expected input files in STOKES_DIR, for refinement STOKES_REF:
!   A_ref<ref>.mtx, Bt_ref<ref>.mtx, B_ref<ref>.mtx
!   rhs_u_ref<ref>.mtx, rhs_p_ref<ref>.mtx
! For STOKES_PREC=MASS_BLOCK, also Mp_ref<ref>.mtx.
!
program amg_d_nest_stokes_file_test
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
  type(psb_d_nest_matrix), target :: mass_prec_matrix
  class(psb_dprec_type), allocatable :: preconditioner
  type(psb_d_vect_type)   :: rhs, x_solution, residual

  type(mm_matrix) :: a_block, bt_block, b_block, mp_block
  type(mm_vector) :: rhs_u, rhs_p

  integer(psb_ipk_) :: my_rank, num_procs, info
  integer(psb_ipk_) :: max_iter, trace_level, stop_criterion, n_iter
  integer(psb_ipk_) :: n_insert
  integer(psb_lpk_) :: n_u, n_p, n_total
  integer(psb_lpk_), allocatable :: rows(:), cols(:)
  real(psb_dpk_),    allocatable :: vals(:)
  real(psb_dpk_),    allocatable :: x_global(:)
  real(psb_dpk_) :: eps, final_residual, residual_norm, rhs_norm, residual_tol, schur_tol
  real(psb_dpk_) :: t0, t_read, t_assemble, t_prec, t_solve
  character(len=256) :: input_dir, method, ptype, ref_text, composition, schur_solve
  character(len=512) :: filename
  integer :: refinement, schur_maxit, k
  logical :: use_pressure_mass

  call psb_init(context)
  call psb_info(context, my_rank, num_procs)

  call get_string_env('STOKES_DIR', input_dir, '../../../../dealii-test/build')
  call get_string_env('STOKES_METHOD', method, 'BICGSTAB')
  call get_string_env('STOKES_PREC', ptype, 'AMG_BLOCK')
  use_pressure_mass = (psb_toupper(trim(ptype)) == 'MASS_BLOCK') .or. &
       & (psb_toupper(trim(ptype)) == 'PRESSURE_MASS') .or. &
       & (psb_toupper(trim(ptype)) == 'AMG_MASS_BLOCK')
  if ((psb_toupper(trim(ptype)) == 'AMG_MASS_BLOCK') .and. (num_procs /= 1)) then
    if (my_rank == 0) write(*,'(a)') 'FAIL: AMG_MASS_BLOCK currently requires one MPI process'
    call psb_abort(context)
  end if

  call get_string_env('STOKES_COMPOSITION', composition, 'SCHUR_FULL')
  call get_string_env('STOKES_SCHUR_SOLVE', schur_solve, 'MATRIX_FREE')
  call get_int_env('STOKES_REF', refinement, 0)
  call get_int_env('STOKES_ITMAX', max_iter, 2000)
  call get_int_env('STOKES_SCHUR_MAXIT', schur_maxit, 10)
  call get_real_env('STOKES_EPS', eps, 1.0e-8_psb_dpk_)
  call get_real_env('STOKES_SCHUR_TOL', schur_tol, 0.0_psb_dpk_)
  call get_real_env('STOKES_RESIDUAL_TOL', residual_tol, 1.0e-6_psb_dpk_)
  trace_level    = 0
  stop_criterion = 2
  write(ref_text, '(i0)') refinement

  if (my_rank == 0) then
    write(*,'(a)') 'Stokes nested MatrixMarket test'
    write(*,'(a,a)') '  input dir : ', trim(input_dir)
    write(*,'(a,a)') '  refinement: ', trim(ref_text)
    write(*,'(a,a)') '  method    : ', trim(method)
    write(*,'(a,a)') '  prec      : ', trim(ptype)
    if (use_pressure_mass) write(*,'(a)') '  pressure block: Mp (pressure mass matrix)'
    write(*,'(a,a)') '  composition: ', trim(composition)
    write(*,'(a,a)') '  schur solve: ', trim(schur_solve)
    write(*,'(a,i0)') '  schur maxit: ', schur_maxit
    write(*,'(a,es12.4)') '  schur tol  : ', schur_tol
    write(*,'(a,es12.4)') '  residual tol: ', residual_tol
  end if

  t0 = psb_wtime()
  call read_mm_matrix(path_join(input_dir, 'A_ref'  // trim(ref_text) // '.mtx'), a_block)
  call read_mm_matrix(path_join(input_dir, 'Bt_ref' // trim(ref_text) // '.mtx'), bt_block)
  call read_mm_matrix(path_join(input_dir, 'B_ref'  // trim(ref_text) // '.mtx'), b_block)
  if (use_pressure_mass) then
    call read_mm_matrix(path_join(input_dir, 'Mp_ref' // trim(ref_text) // '.mtx'), mp_block)
  end if
  call read_mm_vector(path_join(input_dir, 'rhs_u_ref' // trim(ref_text) // '.mtx'), rhs_u)
  call read_mm_vector(path_join(input_dir, 'rhs_p_ref' // trim(ref_text) // '.mtx'), rhs_p)
  t_read = psb_wtime() - t0

  n_u = a_block%nrow
  n_p = b_block%nrow
  n_total = n_u + n_p

  if ((a_block%ncol /= n_u) .or. (bt_block%nrow /= n_u) .or. &
      (bt_block%ncol /= n_p) .or. (b_block%ncol /= n_u) .or. &
      (rhs_u%n /= n_u) .or. (rhs_p%n /= n_p)) then
    if (my_rank == 0) write(*,*) 'FAIL: inconsistent Stokes block dimensions'
    call psb_abort(context)
  end if
  if (use_pressure_mass) then
    if ((mp_block%nrow /= n_p) .or. (mp_block%ncol /= n_p)) then
      if (my_rank == 0) write(*,*) 'FAIL: pressure mass matrix has inconsistent dimensions'
      call psb_abort(context)
    end if
  end if

  t0 = psb_wtime()
  call nested_matrix%init(context, [n_u, n_p], info)
  call check_info(info, 'nested_matrix%init')

  call select_owned_entries(a_block,  nested_matrix%get_owned_rows(1), rows, cols, vals, n_insert)
  call nested_matrix%ins(1, 1, n_insert, rows, cols, vals, info)
  call check_info(info, 'insert A')
  call clear_triplets(rows, cols, vals)

  call select_owned_entries(bt_block, nested_matrix%get_owned_rows(1), rows, cols, vals, n_insert)
  call nested_matrix%ins(1, 2, n_insert, rows, cols, vals, info)
  call check_info(info, 'insert Bt')
  call clear_triplets(rows, cols, vals)

  call select_owned_entries(b_block,  nested_matrix%get_owned_rows(2), rows, cols, vals, n_insert)
  call nested_matrix%ins(2, 1, n_insert, rows, cols, vals, info)
  call check_info(info, 'insert B')
  call clear_triplets(rows, cols, vals)

  call nested_matrix%asb(info)
  call check_info(info, 'nested_matrix%asb')

  if (use_pressure_mass) then
    call mass_prec_matrix%init(context, [n_u, n_p], info)
    call check_info(info, 'mass_prec_matrix%init')
    call select_owned_entries(a_block, mass_prec_matrix%get_owned_rows(1), rows, cols, vals, n_insert)
    call mass_prec_matrix%ins(1, 1, n_insert, rows, cols, vals, info)
    call check_info(info, 'insert preconditioner A')
    call clear_triplets(rows, cols, vals)
    call select_owned_entries(bt_block, mass_prec_matrix%get_owned_rows(1), rows, cols, vals, n_insert)
    call mass_prec_matrix%ins(1, 2, n_insert, rows, cols, vals, info)
    call check_info(info, 'insert preconditioner Bt')
    call clear_triplets(rows, cols, vals)
    call select_owned_entries(b_block, mass_prec_matrix%get_owned_rows(2), rows, cols, vals, n_insert)
    call mass_prec_matrix%ins(2, 1, n_insert, rows, cols, vals, info)
    call check_info(info, 'insert preconditioner B')
    call clear_triplets(rows, cols, vals)
    call select_owned_entries(mp_block, mass_prec_matrix%get_owned_rows(2), rows, cols, vals, n_insert)
    call mass_prec_matrix%ins(2, 2, n_insert, rows, cols, vals, info)
    call check_info(info, 'insert pressure mass matrix')
    call clear_triplets(rows, cols, vals)
    call mass_prec_matrix%asb(info)
    call check_info(info, 'mass_prec_matrix%asb')
  end if
  t_assemble = psb_wtime() - t0

  call psb_geall(rhs, nested_matrix%desc_glob, info)
  call check_info(info, 'psb_geall(rhs)')
  call insert_rhs(rhs_u, rhs_p, n_u, rhs, nested_matrix%desc_glob, info)
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
  select case (psb_toupper(trim(ptype)))
  case ('NEST')
    if (use_pressure_mass) then
      call preconditioner%set('COMPOSITION', 'DIAG', info)
      call check_info(info, 'preconditioner%set')
    else
      call preconditioner%set('COMPOSITION', trim(composition), info)
      call check_info(info, 'preconditioner%set')
    end if
    call check_info(info, 'prec%set(COMPOSITION)')
    if (.not. use_pressure_mass) then
      call preconditioner%set('SCHUR_SOLVE', trim(schur_solve), info)
      call check_info(info, 'prec%set(SCHUR_SOLVE)')
      call preconditioner%set('SCHUR_MAXIT', schur_maxit, info)
      call check_info(info, 'prec%set(SCHUR_MAXIT)')
      call preconditioner%set('SCHUR_TOL', schur_tol, info)
      call check_info(info, 'prec%set(SCHUR_TOL)')
    end if
    call preconditioner%set('BLOCK_SOLVE', 'BJAC', info)
    call check_info(info, 'prec%set(BLOCK_SOLVE)')
    call preconditioner%set('SUB_SOLVE', 'ILU', info)
    call check_info(info, 'prec%set(SUB_SOLVE)')
    call preconditioner%set('SUB_FILLIN', 0, info)
    call check_info(info, 'prec%set(SUB_FILLIN)')
  case ('AMG_STOKES')
    call preconditioner%set('SCHUR_SOLVE', trim(schur_solve), info)
    call check_info(info, 'prec%set(SCHUR_SOLVE)')
    call preconditioner%set('SCHUR_MAXIT', schur_maxit, info)
    call check_info(info, 'prec%set(SCHUR_MAXIT)')
    call preconditioner%set('SCHUR_TOL', schur_tol, info)
    call check_info(info, 'prec%set(SCHUR_TOL)')
  end select
  if (my_rank == 0) then
    write(*,'(a)') '  resolved operator composition:'
    write(*,'(a)') '    row 1: [A, Bt]'
    write(*,'(a)') '    row 2: [B,  0]'
    write(*,'(a,a)') '  resolved preconditioner: ', trim(ptype)
    select case (psb_toupper(trim(ptype)))
    case ('NEST')
      if (use_pressure_mass) then
        write(*,'(a)') '    composition: DIAG'
      else
        write(*,'(a,a)') '    composition: ', trim(composition)
      end if
      write(*,'(a)') '    block 1 (velocity, A): BJAC / ILU(0)'
      if (use_pressure_mass) then
        write(*,'(a)') '    block 2 (pressure, Mp): BJAC / ILU(0)'
      else
        write(*,'(a)') '    block 2 (pressure, A22=0): NONE; handled by Schur solve'
        write(*,'(a,a,a,i0,a,es12.4)') '      Schur: ', trim(schur_solve), &
             & ', maxit=', schur_maxit, ', tol=', schur_tol
      end if
    case ('AMG_STOKES')
      write(*,'(a)') '    composition: Stokes block factorization'
      write(*,'(a)') '    block 1 (velocity, A): AMG (ML)'
      write(*,'(a,a,a,i0,a,es12.4)') '    block 2 (pressure): Schur ', trim(schur_solve), &
           & ', maxit=', schur_maxit, ', tol=', schur_tol
    case ('AMG_STOKES_MASS')
      write(*,'(a)') '    composition: DIAG (equation 21)'
      write(*,'(a)') '    block 1 (velocity, A): fixed AMG (ML)'
      write(*,'(a)') '    block 2 (pressure, Mp): lumped diagonal inverse'
    case ('ML')
      write(*,'(a)') '    all blocks: global AMG (ML)'
    end select
  end if
  if (use_pressure_mass .and. (psb_toupper(trim(ptype)) /= 'AMG_STOKES_MASS')) then
    call preconditioner%build(mass_prec_matrix%a_glob, mass_prec_matrix%desc_glob, info)
    call check_info(info, 'preconditioner%build')
  else
    call preconditioner%build(nested_matrix%a_glob, nested_matrix%desc_glob, info)
    call check_info(info, 'preconditioner%build')
  end if
  call check_info(info, 'prec%build')
  if (psb_toupper(trim(ptype)) == 'AMG_STOKES_MASS') then
    select type (amg_prec => preconditioner)
    type is (amg_d_nested_block_prec_type)
      if (size(amg_prec%schur_diag) /= n_p) then
        if (my_rank == 0) write(*,'(a)') 'FAIL: pressure diagonal size mismatch'
        call psb_abort(context)
      end if
      amg_prec%schur_diag(:) = dzero
      do k = 1, mp_block%nnz
        if (mp_block%row(k) == mp_block%col(k)) &
             & amg_prec%schur_diag(mp_block%row(k)) = mp_block%val(k)
      end do
      if (any(amg_prec%schur_diag <= dzero)) then
        if (my_rank == 0) write(*,'(a)') 'FAIL: pressure mass diagonal must be positive'
        call psb_abort(context)
      end if
    end select
  end if
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
  call check_info(info, 'psb_spmm')
  residual_norm = psb_genrm2(residual, nested_matrix%desc_glob, info)
  call check_info(info, 'norm')
  rhs_norm      = psb_genrm2(rhs,      nested_matrix%desc_glob, info)
  call check_info(info, 'norm')
  if (residual_norm / max(rhs_norm, tiny(done)) > residual_tol) then
    if (my_rank == 0) write(*,'(a)') '[FAIL] amg_d_nest_stokes_file_test: residual above STOKES_RESIDUAL_TOL'
    call psb_abort(context)
  end if

  allocate(x_global(n_total))
  call psb_gather(x_global, x_solution, nested_matrix%desc_glob, info, root=psb_root_)
  call check_info(info, 'psb_gather(solution)')

  if (my_rank == 0) then
    call write_mm_vector(path_join(input_dir, 'psblas_solution_u_ref' // trim(ref_text) // '.mtx'), &
         & x_global(1:n_u))
    call write_mm_vector(path_join(input_dir, 'psblas_solution_p_ref' // trim(ref_text) // '.mtx'), &
         & x_global(n_u+1:n_total))


    write(*,'(a,i0,a,i0,a,i0)') '  block sizes: n_u=', n_u, ' n_p=', n_p, ' total=', n_total
    write(*,'(a,i8,a,es12.4)') '  iterations=', n_iter, '  solver err=', final_residual
    write(*,'(a,es12.4)') '  ||b-Ax||_2 / ||b||_2=', residual_norm / max(rhs_norm, tiny(done))
    write(*,'(a,es12.4,a,es12.4,a,es12.4,a,es12.4)') &
         & '  times read=', t_read, ' assemble=', t_assemble, ' prec=', t_prec, ' solve=', t_solve
    write(*,'(a)') '[PASS] amg_d_nest_stokes_file_test'
  end if

9999 continue
  if (allocated(preconditioner)) call preconditioner%free(info)
  call psb_gefree(residual, nested_matrix%desc_glob, info)
  call check_info(info, 'psb_gefree')
  call psb_gefree(x_solution, nested_matrix%desc_glob, info)
  call check_info(info, 'psb_gefree')
  call psb_gefree(rhs, nested_matrix%desc_glob, info)
  call check_info(info, 'psb_gefree')
  if (use_pressure_mass) call mass_prec_matrix%free(info)
  call nested_matrix%free(info)
  call check_info(info, 'nested_matrix%free')
  call psb_exit(context)

contains

  subroutine allocate_preconditioner()
    select case (psb_toupper(trim(ptype)))
    case ('MASS_BLOCK','PRESSURE_MASS')
      allocate(psb_dprec_type :: preconditioner, stat=info)
      call check_info(info, 'allocate pressure-mass preconditioner')
      ptype = 'NEST'
    case ('AMG_BLOCK','AMG_SCHUR','SCHUR_AMG')
      allocate(amg_d_nested_block_prec_type :: preconditioner, stat=info)
      call check_info(info, 'allocate AMG block preconditioner')
      ptype = 'AMG_STOKES'
    case ('AMG','ML','MULTILEVEL')
      allocate(amg_dprec_type :: preconditioner, stat=info)
      call check_info(info, 'allocate AMG preconditioner')
      ptype = 'ML'
    case ('AMG_MASS_BLOCK')
      allocate(amg_d_nested_block_prec_type :: preconditioner, stat=info)
      call check_info(info, 'allocate AMG pressure-mass block preconditioner')
      ptype = 'AMG_STOKES_MASS'
    case default
      allocate(psb_dprec_type :: preconditioner, stat=info)
      call check_info(info, 'allocate PSBLAS preconditioner')
    end select
    call preconditioner%init(context, trim(ptype), info)
  end subroutine allocate_preconditioner

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
    integer(psb_lpk_) :: i, first_owned, last_owned
    integer(psb_ipk_) :: k

    count = 0
    if (size(owned_rows) == 0) then
      allocate(out_rows(0), out_cols(0), out_vals(0))
      return
    end if

    first_owned = minval(owned_rows)
    last_owned  = maxval(owned_rows)
    do i = 1_psb_lpk_, mat%nnz
      if (mat%row(i) >= first_owned .and. mat%row(i) <= last_owned) count = count + 1
    end do

    allocate(out_rows(count), out_cols(count), out_vals(count))
    k = 0
    do i = 1_psb_lpk_, mat%nnz
      if (mat%row(i) >= first_owned .and. mat%row(i) <= last_owned) then
        k = k + 1
        out_rows(k) = mat%row(i)
        out_cols(k) = mat%col(i)
        out_vals(k) = mat%val(i)
      end if
    end do
  end subroutine select_owned_entries

  subroutine insert_rhs(rhs_u_in, rhs_p_in, velocity_size, rhs_vec, desc, info)
    type(mm_vector), intent(in) :: rhs_u_in, rhs_p_in
    integer(psb_lpk_), intent(in) :: velocity_size
    type(psb_d_vect_type), intent(inout) :: rhs_vec
    type(psb_desc_type), intent(in) :: desc
    integer(psb_ipk_), intent(out) :: info
    integer(psb_lpk_), allocatable :: idx(:)
    real(psb_dpk_), allocatable :: values(:)
    integer(psb_lpk_) :: i

    allocate(idx(rhs_u_in%n), values(rhs_u_in%n))
    do i = 1_psb_lpk_, rhs_u_in%n
      idx(i) = i
      values(i) = rhs_u_in%val(i)
    end do
    call psb_geins(int(rhs_u_in%n, psb_ipk_), idx, values, rhs_vec, desc, info)
    if (info /= psb_success_) return
    deallocate(idx, values)

    allocate(idx(rhs_p_in%n), values(rhs_p_in%n))
    do i = 1_psb_lpk_, rhs_p_in%n
      idx(i) = velocity_size + i
      values(i) = rhs_p_in%val(i)
    end do
    call psb_geins(int(rhs_p_in%n, psb_ipk_), idx, values, rhs_vec, desc, info)
    deallocate(idx, values)
  end subroutine insert_rhs

  subroutine clear_triplets(in_rows, in_cols, in_vals)
    integer(psb_lpk_), allocatable, intent(inout) :: in_rows(:), in_cols(:)
    real(psb_dpk_), allocatable, intent(inout) :: in_vals(:)
    if (allocated(in_rows)) deallocate(in_rows)
    if (allocated(in_cols)) deallocate(in_cols)
    if (allocated(in_vals)) deallocate(in_vals)
  end subroutine clear_triplets

end program amg_d_nest_stokes_file_test
