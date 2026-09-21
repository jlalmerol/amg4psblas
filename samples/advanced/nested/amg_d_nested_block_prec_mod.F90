module amg_d_nested_block_prec_mod
  use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_epk_, psb_dpk_, psb_success_, &
       & psb_err_invalid_input_, psb_err_invalid_mat_state_, psb_err_alloc_dealloc_, &
       & psb_errpush, psb_toupper, psb_geall, psb_geasb, psb_gefree, &
       & psb_spall, psb_spins, psb_spasb, psb_dupl_add_, &
       & psb_halo, psb_gedot, psb_sum, psb_cdall, psb_cdins, psb_cdasb, &
       & done, dzero, psb_ctxt_type, psb_info
  use psb_prec_mod, only : psb_dprec_type
  use psb_d_mat_mod, only : psb_dspmat_type, psb_d_get_diag, psb_d_csgetrow
  use psb_d_vect_mod, only : psb_d_vect_type
  use psb_desc_mod, only : psb_desc_type
  use psb_d_base_mat_mod, only : psb_d_base_sparse_mat
  use psb_d_base_vect_mod, only : psb_d_base_vect_type
  use psb_i_base_vect_mod, only : psb_i_base_vect_type
  use psb_d_nest_mod, only : psb_d_nest_base_mat, psb_d_nest_get_n_fields, &
       & psb_d_nest_get_block, psb_d_nest_get_field_desc, &
       & psb_d_nest_restrict_field, psb_d_nest_restrict_field_local, &
       & psb_d_nest_prolong_field, psb_d_nest_apply_block
  use amg_d_prec_mod, only : amg_dprec_type
  implicit none

  private
  public :: amg_d_nested_block_prec_type

  type, extends(psb_dprec_type) :: amg_d_nested_block_prec_type
    character(len=16) :: mode = 'AUTO'
    integer(psb_ipk_) :: nfields = 0
    type(psb_d_nest_base_mat), pointer :: nest_op => null()
    type(amg_dprec_type), allocatable :: field_amg(:)
    type(amg_dprec_type) :: schur_amg
    type(psb_dspmat_type) :: shifted_mat
    type(psb_desc_type) :: shifted_desc
    logical, allocatable :: use_amg(:)
    logical :: schur_amg_built = .false.
    logical :: mgw_exact_built = .false.
    real(psb_dpk_), allocatable :: mgw_fact11(:,:), mgw_fact22(:,:), mgw_fact_s(:,:)

    real(psb_dpk_), allocatable :: schur_diag(:)
    real(psb_dpk_), allocatable :: mass_diag1(:), mass_diag2(:)
    real(psb_dpk_), allocatable :: schur_active(:)
    character(len=16) :: schur_solve = 'MATRIX_FREE'
    integer(psb_ipk_) :: schur_maxit = 8
    real(psb_dpk_) :: schur_tol = dzero
    logical :: wrk_allocated = .false.
  contains
    procedure, pass(prec) :: init   => amg_d_nested_block_init
    procedure, pass(prec) :: build  => amg_d_nested_block_build
    procedure, pass(prec) :: free   => amg_d_nested_block_free
    procedure, pass(prec) :: csetc  => amg_d_nested_block_csetc
    procedure, pass(prec) :: cseti  => amg_d_nested_block_cseti
    procedure, pass(prec) :: csetr  => amg_d_nested_block_csetr
    procedure, pass(prec) :: allocate_wrk => amg_d_nested_block_allocate_wrk
    procedure, pass(prec) :: free_wrk => amg_d_nested_block_free_wrk
    procedure, pass(prec) :: deallocate_wrk => amg_d_nested_block_free_wrk
    procedure, pass(prec) :: is_allocated_wrk => amg_d_nested_block_is_allocated_wrk
    procedure, pass(prec) :: psb_d_apply2v => amg_d_nested_block_apply2v
    procedure, pass(prec) :: psb_d_apply1v => amg_d_nested_block_apply1v

    procedure, pass(prec) :: psb_d_apply2_vect => amg_d_nested_block_apply2_vect
    procedure, pass(prec) :: psb_d_apply1_vect => amg_d_nested_block_apply1_vect
    procedure, pass(prec) :: sizeof => amg_d_nested_block_sizeof
  end type amg_d_nested_block_prec_type

  interface
    subroutine dpotrf(uplo, n, a, lda, lapack_info)
      import :: psb_dpk_
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, lda
      real(psb_dpk_), intent(inout) :: a(lda,*)
      integer, intent(out) :: lapack_info
    end subroutine dpotrf
    subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, lapack_info)
      import :: psb_dpk_
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, nrhs, lda, ldb
      real(psb_dpk_), intent(in) :: a(lda,*)
      real(psb_dpk_), intent(inout) :: b(ldb,*)
      integer, intent(out) :: lapack_info
    end subroutine dpotrs
  end interface

contains

  ! Distributed block product. Only owned input entries are authoritative;
  ! each column field may have a different union halo and local ordering.
  subroutine amg_d_apply_block(nest, i, j, alpha, x, beta, y, info)
    type(psb_d_nest_base_mat), intent(in) :: nest
    integer(psb_ipk_), intent(in) :: i, j
    real(psb_dpk_), intent(in) :: alpha, beta, x(:)
    real(psb_dpk_), intent(inout) :: y(:)
    integer(psb_ipk_), intent(out) :: info
    type(psb_desc_type), pointer :: desc
    real(psb_dpk_), allocatable :: xh(:)
    integer(psb_ipk_) :: nr
    desc => psb_d_nest_get_field_desc(nest, j)
    info = psb_err_invalid_mat_state_
    if (.not. associated(desc)) return
    nr = desc%get_local_rows()
    if (size(x) < nr) return
    allocate(xh(desc%get_local_cols()), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; return
    end if
    xh = dzero
    xh(1:nr) = x(1:nr)
    call psb_halo(xh, desc, info)
    if (info == psb_success_) &
         & call psb_d_nest_apply_block(nest, i, j, alpha, xh, beta, y, info)
  end subroutine amg_d_apply_block

  subroutine amg_d_nested_block_init(ctxt, prec, ptype, info)
    type(psb_ctxt_type), intent(in) :: ctxt
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    character(len=*), intent(in) :: ptype
    integer(psb_ipk_), intent(out) :: info

    call prec%free(info)
    if (info /= psb_success_) return
    prec%ctxt = ctxt
    select case (psb_toupper(trim(ptype)))
    case ('AMG_STOKES','STOKES')
      prec%mode = 'STOKES'
    case ('AMG_STOKES_MASS','STOKES_MASS','AMG_MASS_BLOCK')
      prec%mode = 'STOKES_MASS'
    case ('AMG_KKT','KKT','OPTIMAL','OPTIMAL_CONTROL')
      prec%mode = 'KKT'
    case ('AMG_KKT_DIAG','KKT_DIAG','KKT_BLOCK_DIAG','OPTIMAL_DIAG')
      prec%mode = 'KKT_DIAG'
    case ('MGW_EXACT')
      prec%mode = 'MGW_EXACT'
    case default
      prec%mode = 'AUTO'
    end select
    prec%schur_solve = 'MATRIX_FREE'
    prec%schur_maxit = 8
    prec%schur_tol = dzero
    info = psb_success_
  end subroutine amg_d_nested_block_init

  subroutine amg_d_nested_block_csetc(prec, what, string, info, ilev, ilmax, pos, idx)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    character(len=*), intent(in) :: what, string
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), optional, intent(in) :: ilev, ilmax, idx
    character(len=*), optional, intent(in) :: pos

    info = psb_success_
    select case (psb_toupper(trim(what)))
    case ('SCHUR_SOLVE','NEST_SCHUR_SOLVE')
      select case (psb_toupper(trim(string)))
      case ('MATRIX_FREE','MATFREE','MF')
        prec%schur_solve = 'MATRIX_FREE'
      case ('SELF','SELFP')
        prec%schur_solve = 'SELF'
      case default
        info = psb_err_invalid_input_
        call psb_errpush(info, 'amg_nested_block_csetc', a_err='SCHUR_SOLVE')
      end select
    case default
      info = psb_err_invalid_input_
      call psb_errpush(info, 'amg_nested_block_csetc', a_err=trim(what))
    end select
  end subroutine amg_d_nested_block_csetc

  subroutine amg_d_nested_block_cseti(prec, what, val, info, ilev, ilmax, pos, idx)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    character(len=*), intent(in) :: what
    integer(psb_ipk_), intent(in) :: val
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), optional, intent(in) :: ilev, ilmax, idx
    character(len=*), optional, intent(in) :: pos

    info = psb_success_
    select case (psb_toupper(trim(what)))
    case ('SCHUR_MAXIT','NEST_SCHUR_MAXIT')
      prec%schur_maxit = max(1, val)
    case default
      info = psb_err_invalid_input_
      call psb_errpush(info, 'amg_nested_block_cseti', a_err=trim(what))
    end select
  end subroutine amg_d_nested_block_cseti

  subroutine amg_d_nested_block_csetr(prec, what, val, info, ilev, ilmax, pos, idx)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    character(len=*), intent(in) :: what
    real(psb_dpk_), intent(in) :: val
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), optional, intent(in) :: ilev, ilmax, idx
    character(len=*), optional, intent(in) :: pos

    info = psb_success_
    select case (psb_toupper(trim(what)))
    case ('SCHUR_TOL','NEST_SCHUR_TOL')
      prec%schur_tol = max(dzero, val)
    case default
      info = psb_err_invalid_input_
      call psb_errpush(info, 'amg_nested_block_csetr', a_err=trim(what))
    end select
  end subroutine amg_d_nested_block_csetr

  subroutine amg_d_nested_block_build(a, desc_a, prec, info, amold, vmold, imold)
    type(psb_dspmat_type), intent(inout), target :: a
    type(psb_desc_type), intent(inout), target :: desc_a
    class(amg_d_nested_block_prec_type), intent(inout), target :: prec
    integer(psb_ipk_), intent(out) :: info
    class(psb_d_base_sparse_mat), intent(in), optional :: amold
    class(psb_d_base_vect_type), intent(in), optional :: vmold
    class(psb_i_base_vect_type), intent(in), optional :: imold

    integer(psb_ipk_) :: i
    type(psb_dspmat_type), pointer :: blk
    type(psb_desc_type), pointer :: fdesc

    info = psb_success_
    select type (ap => a%a)
    type is (psb_d_nest_base_mat)
      prec%nest_op => ap
    class default
      info = psb_err_invalid_input_
      call psb_errpush(info, 'amg_nested_block_build', a_err='requires NEST matrix')
      return
    end select

    prec%nfields = psb_d_nest_get_n_fields(prec%nest_op)
    if (prec%nfields /= 2 .and. prec%nfields /= 3) then
      info = psb_err_invalid_input_
      call psb_errpush(info, 'amg_nested_block_build', a_err='expected 2 or 3 fields')
      return
    end if
    if (prec%mode == 'AUTO') then
      if (prec%nfields == 2) prec%mode = 'STOKES'
      if (prec%nfields == 3) prec%mode = 'KKT'
    end if

    allocate(prec%field_amg(prec%nfields), prec%use_amg(prec%nfields), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info, 'amg_nested_block_build', a_err='field preconditioners')
      return
    end if
    prec%use_amg(:) = .false.

    do i = 1, prec%nfields
      blk => psb_d_nest_get_block(prec%nest_op, i, i)
      if (.not. associated(blk)) cycle
      fdesc => psb_d_nest_get_field_desc(prec%nest_op, i)
      if (.not. associated(fdesc)) then
        info = psb_err_invalid_mat_state_
        call psb_errpush(info, 'amg_nested_block_build', a_err='missing field descriptor')
        return
      end if
      if ((i < prec%nfields) .and. (trim(prec%mode) /= 'KKT_DIAG') .and. &
           & (trim(prec%mode) /= 'MGW_EXACT')) then
        call prec%field_amg(i)%init(fdesc%get_ctxt(), 'ML', info)
        if (info /= psb_success_) return
        if (trim(prec%mode) == 'STOKES_MASS') then
          call prec%field_amg(i)%set('SMOOTHER_TYPE', 'JACOBI', info)
          if (info /= psb_success_) return
          call prec%field_amg(i)%set('SMOOTHER_SWEEPS', 4, info)
          if (info /= psb_success_) return
          call prec%field_amg(i)%set('SMOOTHER_TYPE', 'JACOBI', info, pos='post')
          if (info /= psb_success_) return
          call prec%field_amg(i)%set('SMOOTHER_SWEEPS', 4, info, pos='post')
          if (info /= psb_success_) return
          call prec%field_amg(i)%set('COARSE_SOLVE', 'JACOBI', info)
          if (info /= psb_success_) return
          call prec%field_amg(i)%set('COARSE_SWEEPS', 8, info)
          if (info /= psb_success_) return
        end if
        call prec%field_amg(i)%build(blk, fdesc, info, amold=amold, vmold=vmold, imold=imold)
        if (info /= psb_success_) return
        prec%use_amg(i) = .true.
      end if
    end do

    select case (trim(prec%mode))
    case ('STOKES')
      call amg_d_build_stokes_schur_diag(prec, info)
    case ('STOKES_MASS')
      call amg_d_build_stokes_mass_diag(prec, info)
    case ('KKT')
      call amg_d_build_kkt_schur_diag(prec, info)
    case ('KKT_DIAG')
      call amg_d_build_kkt_dominant_schur(prec, info, amold, vmold, imold)
    case ('MGW_EXACT')
      call amg_d_build_mgw_exact(prec, info)
    case default
      info = psb_err_invalid_input_
      call psb_errpush(info, 'amg_nested_block_build', a_err='unknown mode')
    end select
  end subroutine amg_d_nested_block_build

  subroutine amg_d_nested_block_allocate_wrk(prec, info, vmold, desc)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out) :: info
    class(psb_d_base_vect_type), intent(in), optional :: vmold
    type(psb_desc_type), intent(in), optional :: desc

    info = psb_success_
    prec%wrk_allocated = .true.
  end subroutine amg_d_nested_block_allocate_wrk

  subroutine amg_d_nested_block_free_wrk(prec, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out) :: info

    integer(psb_ipk_) :: i, linfo

    info = psb_success_
    if (allocated(prec%field_amg) .and. allocated(prec%use_amg)) then
      do i = 1, min(size(prec%field_amg), size(prec%use_amg))
        if (.not. prec%use_amg(i)) cycle
        call prec%field_amg(i)%free_wrk(linfo)
        if (linfo /= psb_success_ .and. info == psb_success_) info = linfo
      end do
    end if
    prec%wrk_allocated = .false.
  end subroutine amg_d_nested_block_free_wrk

  function amg_d_nested_block_is_allocated_wrk(prec) result(res)
    class(amg_d_nested_block_prec_type), intent(in) :: prec
    logical :: res

    res = prec%wrk_allocated
  end function amg_d_nested_block_is_allocated_wrk

  subroutine amg_d_apply_field(prec, ifld, rhs, sol, desc, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in) :: ifld
    real(psb_dpk_), intent(in) :: rhs(:)
    real(psb_dpk_), intent(out) :: sol(:)
    type(psb_desc_type), intent(in) :: desc
    integer(psb_ipk_), intent(out) :: info

    type(psb_d_vect_type) :: xv, yv
    integer(psb_ipk_) :: nloc, linfo

    info = psb_success_
    nloc = desc%get_local_rows()
    sol(:) = dzero
    if (.not. prec%use_amg(ifld)) then
      sol(1:nloc) = rhs(1:nloc)
      return
    end if
    call psb_geall(xv, desc, info)
    if (info == psb_success_) call psb_geall(yv, desc, info)
    if (info == psb_success_) call psb_geasb(xv, desc, info)
    if (info == psb_success_) call psb_geasb(yv, desc, info)
    if (info == psb_success_) call xv%zero()
    if (info == psb_success_) call yv%zero()
    if (info == psb_success_) call xv%set(rhs(1:min(size(rhs), nloc)))
    if (info == psb_success_) call prec%field_amg(ifld)%apply(xv, yv, desc, info, trans='N')
    if (info == psb_success_) sol(1:nloc) = yv%get_vect(nloc)
    call psb_gefree(xv, desc, linfo)
    call psb_gefree(yv, desc, linfo)
  end subroutine amg_d_apply_field

  subroutine amg_d_schur_action(prec, x, y, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(in) :: x(:)
    real(psb_dpk_), intent(out) :: y(:)
    integer(psb_ipk_), intent(out) :: info

    type(psb_desc_type), pointer :: d1, d2
    real(psb_dpk_), allocatable :: t1(:), w1(:), t2(:), w2(:)
    integer(psb_ipk_) :: n1, n2

    info = psb_success_
    y(:) = dzero

    call amg_d_apply_block(prec%nest_op, prec%nfields, prec%nfields, done, x, dzero, y, info)
    if (info /= psb_success_) return
    select case (trim(prec%mode))
    case ('STOKES')
      d1 => psb_d_nest_get_field_desc(prec%nest_op, 1)
      if (.not. associated(d1)) then
        info = psb_err_invalid_mat_state_; return
      end if
      n1 = d1%get_local_cols()
      allocate(t1(n1), w1(n1), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_; return
      end if
      t1(:) = dzero; w1(:) = dzero
      call amg_d_apply_block(prec%nest_op, 1, 2, done, x, dzero, t1, info)
      if (info == psb_success_) call amg_d_apply_field(prec, 1, t1, w1, d1, info)
      if (info == psb_success_) call amg_d_apply_block(prec%nest_op, 2, 1, -done, w1, done, y, info)
      deallocate(t1, w1)

    case ('KKT')
      d1 => psb_d_nest_get_field_desc(prec%nest_op, 1)
      d2 => psb_d_nest_get_field_desc(prec%nest_op, 2)
      if (.not. associated(d1) .or. .not. associated(d2)) then
        info = psb_err_invalid_mat_state_; return
      end if
      n1 = d1%get_local_cols(); n2 = d2%get_local_cols()
      allocate(t1(n1), w1(n1), t2(n2), w2(n2), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_; return
      end if
      t1(:) = dzero; w1(:) = dzero; t2(:) = dzero; w2(:) = dzero
      call amg_d_apply_block(prec%nest_op, 1, 3, done, x, dzero, t1, info)
      if (info == psb_success_) call amg_d_apply_field(prec, 1, t1, w1, d1, info)
      if (info == psb_success_) call amg_d_apply_block(prec%nest_op, 3, 1, -done, w1, done, y, info)
      if (info == psb_success_) call amg_d_apply_block(prec%nest_op, 2, 3, done, x, dzero, t2, info)
      if (info == psb_success_) call amg_d_apply_field(prec, 2, t2, w2, d2, info)
      if (info == psb_success_) call amg_d_apply_block(prec%nest_op, 3, 2, -done, w2, done, y, info)
      deallocate(t1, w1, t2, w2)

    case default
      info = psb_err_invalid_input_
    end select
  end subroutine amg_d_schur_action

  subroutine amg_d_schur_solve(prec, rhs, sol, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(in) :: rhs(:)
    real(psb_dpk_), intent(out) :: sol(:)
    integer(psb_ipk_), intent(out) :: info

    real(psb_dpk_), allocatable :: r(:), rr(:), p(:), v(:), s(:), t(:), ph(:), sh(:)
    type(psb_desc_type), pointer :: desc
    integer(psb_ipk_) :: k, n, nr
    real(psb_dpk_) :: alpha, beta, omega, rho, rho_old, denom, rhsn, rn, floorv

    info = psb_success_
    if (psb_toupper(trim(prec%schur_solve)) == 'SELF') then
      call amg_d_schur_diag_apply(prec, rhs, sol)
      sol(:) = -sol(:)
      return
    end if

    desc => psb_d_nest_get_field_desc(prec%nest_op, prec%nfields)
    nr = desc%get_local_rows()
    n = size(rhs)
    allocate(r(n), rr(n), p(n), v(n), s(n), t(n), ph(n), sh(n), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; return
    end if

    floorv = sqrt(tiny(done))
    sol(:) = dzero
    r(:) = dzero
    r(1:nr) = rhs(1:nr)
    rr(:) = r(:)
    p(:) = dzero
    v(:) = dzero
    alpha = done
    omega = done
    rho_old = done
    rhsn = sqrt(max(dzero, psb_gedot(rhs, rhs, desc, info)))
    if (info /= psb_success_) goto 100
    if (rhsn == dzero) goto 100
    r(1:nr) = rhs(1:nr) / rhsn
    rr = r

    do k = 1, max(1, prec%schur_maxit)
      rho = psb_gedot(rr, r, desc, info)
      if (info /= psb_success_) exit
      if (abs(rho) <= floorv) exit
      if (k == 1) then
        p(:) = r(:)
      else
        if (abs(rho_old) <= floorv .or. abs(omega) <= floorv) exit
        beta = (rho / rho_old) * (alpha / omega)
        p(:) = r(:) + beta * (p(:) - omega * v(:))
      end if

      call amg_d_schur_diag_apply(prec, p, ph)
      call amg_d_schur_action(prec, ph, v, info)
      if (info /= psb_success_) exit
      denom = psb_gedot(rr, v, desc, info)
      if (info /= psb_success_) exit
      if (abs(denom) <= floorv) exit
      alpha = rho / denom
      s(:) = r(:) - alpha * v(:)
      rn = sqrt(psb_gedot(s, s, desc, info))
      if (info /= psb_success_) exit
      if (rn <= max(prec%schur_tol, floorv)) then
        sol(:) = sol(:) + alpha * ph(:)
        exit
      end if

      call amg_d_schur_diag_apply(prec, s, sh)
      call amg_d_schur_action(prec, sh, t, info)
      if (info /= psb_success_) exit
      denom = psb_gedot(t, t, desc, info)
      if (info /= psb_success_) exit
      if (abs(denom) <= floorv) exit
      omega = psb_gedot(t, s, desc, info) / denom
      if (info /= psb_success_) exit
      sol(:) = sol(:) + alpha * ph(:) + omega * sh(:)
      r(:) = s(:) - omega * t(:)
      rn = sqrt(psb_gedot(r, r, desc, info))
      if (info /= psb_success_) exit
      if (rn <= max(prec%schur_tol, floorv)) exit
      rho_old = rho
    end do

100 continue
    sol(:) = rhsn * sol(:)
    deallocate(r, rr, p, v, s, t, ph, sh)
  end subroutine amg_d_schur_solve

  subroutine amg_d_schur_diag_apply(prec, rhs, sol)
    class(amg_d_nested_block_prec_type), intent(in) :: prec
    real(psb_dpk_), intent(in) :: rhs(:)
    real(psb_dpk_), intent(out) :: sol(:)
    real(psb_dpk_) :: floorv

    floorv = sqrt(tiny(done))
    if (allocated(prec%schur_diag) .and. size(prec%schur_diag) == size(rhs)) then
      sol(:) = rhs(:) / max(abs(prec%schur_diag(:)), floorv)
    else
      sol(:) = rhs(:)
    end if
  end subroutine amg_d_schur_diag_apply
  subroutine amg_d_build_stokes_schur_diag(prec, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out) :: info
    call amg_d_build_schur_diag(prec, info)
  end subroutine amg_d_build_stokes_schur_diag

  subroutine amg_d_build_stokes_mass_diag(prec, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out) :: info
    type(psb_dspmat_type), pointer :: pressure_mass
    type(psb_desc_type), pointer :: pressure_desc
    real(psb_dpk_), allocatable :: diag(:)

    info = psb_success_
    pressure_mass => psb_d_nest_get_block(prec%nest_op, 2, 2)
    pressure_desc => psb_d_nest_get_field_desc(prec%nest_op, 2)
    if (.not. associated(pressure_desc)) then
      info = psb_err_invalid_mat_state_
      return
    end if
    if (associated(pressure_mass)) then
      diag = psb_d_get_diag(pressure_mass, info)
      if (info /= psb_success_) return
      allocate(prec%schur_diag(pressure_desc%get_local_cols()), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_; return
      end if
      prec%schur_diag = done
      prec%schur_diag(1:size(diag)) = diag
      call psb_halo(prec%schur_diag, pressure_desc, info)
      if (info /= psb_success_) return
    else
      allocate(prec%schur_diag(pressure_desc%get_local_cols()), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_
        return
      end if
      prec%schur_diag(:) = done
    end if
  end subroutine amg_d_build_stokes_mass_diag

  subroutine amg_d_build_kkt_schur_diag(prec, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out) :: info
    call amg_d_build_schur_diag(prec, info)
  end subroutine amg_d_build_kkt_schur_diag


  subroutine amg_d_build_schur_diag(prec, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out) :: info
    type(psb_desc_type), pointer :: ds, di
    type(psb_dspmat_type), pointer :: aii
    integer(psb_lpk_), allocatable :: owned(:)
    real(psb_dpk_), allocatable :: e(:), v(:), w(:), diag(:)
    integer(psb_lpk_) :: g
    integer(psb_ipk_) :: ns, nr, i, k, ni
    ds => psb_d_nest_get_field_desc(prec%nest_op, prec%nfields)
    ns = ds%get_local_cols(); nr = ds%get_local_rows()
    owned = ds%get_global_indices(owned=.true.)
    allocate(prec%schur_diag(ns), e(ns), v(ns), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; return
    end if
    prec%schur_diag = dzero
    aii => psb_d_nest_get_block(prec%nest_op, prec%nfields, prec%nfields)
    if (associated(aii)) then
      diag = psb_d_get_diag(aii, info)
      if (info /= psb_success_) return
      prec%schur_diag(1:nr) = diag(1:nr)
      deallocate(diag)
    end if
    do i = 1, prec%nfields-1
      di => psb_d_nest_get_field_desc(prec%nest_op, i)
      ni = di%get_local_rows()
      allocate(w(di%get_local_cols()), diag(ni), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_; return
      end if
      diag = done
      aii => psb_d_nest_get_block(prec%nest_op, i, i)
      if (associated(aii)) then
        diag = psb_d_get_diag(aii, info)
        if (info /= psb_success_) return
      end if
      do g = 1, ds%get_global_rows()
        e = dzero; w = dzero; v = dzero
        do k = 1, nr
          if (owned(k) == g) e(k) = done
        end do
        call amg_d_apply_block(prec%nest_op, i, prec%nfields, done, e, dzero, w, info)
        if (info /= psb_success_) return
        w(1:ni) = w(1:ni) / max(abs(diag), sqrt(tiny(done)))
        call amg_d_apply_block(prec%nest_op, prec%nfields, i, done, w, dzero, v, info)
        if (info /= psb_success_) return
        do k = 1, nr
          if (owned(k) == g) prec%schur_diag(k) = prec%schur_diag(k) - v(k)
        end do
      end do
      deallocate(w, diag)
    end do
    prec%schur_diag(1:nr) = max(abs(prec%schur_diag(1:nr)), sqrt(tiny(done)))
    call psb_halo(prec%schur_diag, ds, info)
  end subroutine amg_d_build_schur_diag

  ! Build the AMG approximation used by the shifted SPD Schur block
  !
  !   S_alpha = H M^{-1} H,       S_alpha^{-1} = H^{-1} M H^{-1},
  !   H = K + alpha^{-1/2} M.
  !
  ! For the PDE-control ordering used by this sample, K=A13, M=A11, and
  ! A22=alpha*M.  The shift also regularizes retained Dirichlet rows for which
  ! the exported stiffness block has a zero diagonal.
  subroutine amg_d_build_kkt_dominant_schur(prec, info, amold, vmold, imold)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out) :: info
    class(psb_d_base_sparse_mat), intent(in), optional :: amold
    class(psb_d_base_vect_type), intent(in), optional :: vmold
    class(psb_i_base_vect_type), intent(in), optional :: imold
    type(psb_dspmat_type), pointer :: stiffness, mass, alpha_mass
    type(psb_desc_type), pointer :: desc1, desc3
    real(psb_dpk_), allocatable :: dmass(:), dalpha_mass(:), dstiffness(:)
    real(psb_dpk_) :: alpha_control, shift, totals(2)
    integer(psb_lpk_), allocatable :: owned1(:), owned3(:)
    integer(psb_ipk_) :: mismatch

    info = psb_success_
    desc1 => psb_d_nest_get_field_desc(prec%nest_op, 1)
    desc3 => psb_d_nest_get_field_desc(prec%nest_op, 3)
    stiffness => psb_d_nest_get_block(prec%nest_op, 1, 3)
    mass => psb_d_nest_get_block(prec%nest_op, 1, 1)
    alpha_mass => psb_d_nest_get_block(prec%nest_op, 2, 2)
    if ((.not. associated(desc1)) .or. (.not. associated(desc3)) .or. &
         & (.not. associated(stiffness)) .or. (.not. associated(mass)) .or. &
         & (.not. associated(alpha_mass))) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'amg_nested_kkt_diag_build', &
           & a_err='requires A11, A13, A22 and field descriptors')
      return
    end if
    owned1 = desc1%get_global_indices(owned=.true.)
    owned3 = desc3%get_global_indices(owned=.true.)
    mismatch = 0
    if (size(owned1) /= size(owned3)) then
      mismatch = 1
    else if (any(owned1 /= owned3)) then
      mismatch = 1
    end if
    call psb_sum(prec%ctxt, mismatch)
    if (mismatch /= 0) then
      info = psb_err_invalid_input_
      call psb_errpush(info, 'amg_nested_kkt_diag_build', &
           & a_err='dominant Schur prototype requires matching fields')
      return
    end if

    dmass = psb_d_get_diag(mass, info)
    if (info /= psb_success_) return
    dalpha_mass = psb_d_get_diag(alpha_mass, info)
    if (info /= psb_success_) return
    totals = [sum(abs(dmass)), sum(abs(dalpha_mass))]
    call psb_sum(prec%ctxt, totals)
    if (any(totals <= tiny(done))) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'amg_nested_kkt_diag_build', &
           & a_err='cannot infer positive control regularization')
      return
    end if
    alpha_control = totals(2) / totals(1)
    shift = done / sqrt(alpha_control)

    dstiffness = psb_d_get_diag(stiffness, info)
    if (info /= psb_success_) return
    if (allocated(prec%schur_active)) deallocate(prec%schur_active)
    allocate(prec%schur_active(desc3%get_local_cols()), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      return
    end if
    prec%schur_active(:) = dzero
    prec%schur_active(1:size(dstiffness)) = merge(done, dzero, &
         & abs(dstiffness(:)) > sqrt(tiny(done)))
    call psb_halo(prec%schur_active, desc3, info)
    if (info /= psb_success_) return
    dstiffness(:) = dstiffness(:) + shift * dmass(:)
    if (minval(abs(dstiffness)) <= tiny(done)) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'amg_nested_kkt_diag_build', &
           & a_err='shifted stiffness has a zero diagonal')
      return
    end if

    ! Assemble the shifted operator as an ordinary distributed square PSBLAS
    ! matrix before passing it to AMG.  AMG cannot build directly from the
    ! local storage used for a nested off-diagonal block.
    call amg_d_assemble_shifted_matrix(prec, stiffness, mass, desc1, shift, info)
    if (info /= psb_success_) return

    call prec%schur_amg%init(desc1%get_ctxt(), 'ML', info)
    if (info /= psb_success_) return
    call prec%schur_amg%set('SMOOTHER_TYPE', 'JACOBI', info)
    if (info /= psb_success_) return
    call prec%schur_amg%set('COARSE_SOLVE', 'BJAC', info)
    if (info /= psb_success_) return
    call prec%schur_amg%build(prec%shifted_mat, prec%shifted_desc, info, amold=amold, &
         & vmold=vmold, imold=imold)
    if (info /= psb_success_) return
    prec%schur_amg_built = .true.

    ! Retain the diagonal for diagnostics and a possible Jacobi fallback.
    if (allocated(prec%schur_diag)) deallocate(prec%schur_diag)
    allocate(prec%schur_diag(size(dstiffness)), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      return
    end if
    prec%schur_diag(:) = abs(dstiffness(:))
    if (allocated(prec%mass_diag1)) deallocate(prec%mass_diag1)
    if (allocated(prec%mass_diag2)) deallocate(prec%mass_diag2)
    allocate(prec%mass_diag1(size(dmass)), &
         & prec%mass_diag2(size(dalpha_mass)), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      return
    end if
    prec%mass_diag1(:) = abs(dmass(:))
    prec%mass_diag2(:) = abs(dalpha_mass(:))
  end subroutine amg_d_build_kkt_dominant_schur
  ! Build the ideal Murphy-Golub-Wathen block diagonal preconditioner.
  ! This deliberately dense route is a serial reference implementation used
  ! to verify the three-eigenvalue MINRES result, not a scalable solver.
  subroutine amg_d_build_mgw_exact(prec, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out) :: info
    type(psb_dspmat_type), pointer :: a11, a22, a13, a23, a31, a32, a33
    type(psb_desc_type), pointer :: d1, d2, d3
    real(psb_dpk_), allocatable :: b13(:,:), b23(:,:), b31(:,:), b32(:,:), b33(:,:)
    real(psb_dpk_), allocatable :: x13(:,:), x23(:,:)
    integer(psb_ipk_) :: iam, np
    integer :: n1, n2, n3, lapack_info

    info = psb_success_
    call psb_info(prec%ctxt, iam, np)
    if (np /= 1) then
      info = psb_err_invalid_input_
      call psb_errpush(info, 'amg_nested_mgw_exact_build', &
           & a_err='MGW_EXACT is a serial reference preconditioner')
      return
    end if

    d1 => psb_d_nest_get_field_desc(prec%nest_op, 1)
    d2 => psb_d_nest_get_field_desc(prec%nest_op, 2)
    d3 => psb_d_nest_get_field_desc(prec%nest_op, 3)
    a11 => psb_d_nest_get_block(prec%nest_op, 1, 1)
    a22 => psb_d_nest_get_block(prec%nest_op, 2, 2)
    a13 => psb_d_nest_get_block(prec%nest_op, 1, 3)
    a23 => psb_d_nest_get_block(prec%nest_op, 2, 3)
    a31 => psb_d_nest_get_block(prec%nest_op, 3, 1)
    a32 => psb_d_nest_get_block(prec%nest_op, 3, 2)
    a33 => psb_d_nest_get_block(prec%nest_op, 3, 3)
    if ((.not. associated(d1)) .or. (.not. associated(d2)) .or. &
         & (.not. associated(d3)) .or. (.not. associated(a11)) .or. &
         & (.not. associated(a22)) .or. (.not. associated(a13)) .or. &
         & (.not. associated(a23)) .or. (.not. associated(a31)) .or. &
         & (.not. associated(a32))) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'amg_nested_mgw_exact_build', &
           & a_err='requires A11,A22,A13,A23,A31,A32')
      return
    end if

    n1 = d1%get_local_cols()
    n2 = d2%get_local_cols()
    n3 = d3%get_local_cols()
    if ((n1 <= 0) .or. (n2 <= 0) .or. (n3 <= 0)) then
      info = psb_err_invalid_input_
      return
    end if

    call amg_d_block_to_dense(a11, n1, n1, prec%mgw_fact11, info)
    if (info /= psb_success_) return
    call amg_d_block_to_dense(a22, n2, n2, prec%mgw_fact22, info)
    if (info /= psb_success_) return
    call amg_d_block_to_dense(a13, n1, n3, b13, info)
    if (info /= psb_success_) return
    call amg_d_block_to_dense(a23, n2, n3, b23, info)
    if (info /= psb_success_) return
    call amg_d_block_to_dense(a31, n3, n1, b31, info)
    if (info /= psb_success_) return
    call amg_d_block_to_dense(a32, n3, n2, b32, info)
    if (info /= psb_success_) return
    if (associated(a33)) then
      call amg_d_block_to_dense(a33, n3, n3, b33, info)
      if (info /= psb_success_) return
    else
      allocate(b33(n3,n3), stat=info)
      if (info /= 0) then
        info = psb_err_alloc_dealloc_
        return
      end if
      b33(:,:) = dzero
    end if

    call dpotrf('L', n1, prec%mgw_fact11, n1, lapack_info)
    if (lapack_info /= 0) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'amg_nested_mgw_exact_build', &
           & a_err='A11 is not SPD')
      return
    end if
    call dpotrf('L', n2, prec%mgw_fact22, n2, lapack_info)
    if (lapack_info /= 0) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'amg_nested_mgw_exact_build', &
           & a_err='A22 is not SPD')
      return
    end if

    x13 = b13
    x23 = b23
    call dpotrs('L', n1, n3, prec%mgw_fact11, n1, x13, n1, lapack_info)
    if (lapack_info /= 0) then
      info = psb_err_invalid_mat_state_
      return
    end if
    call dpotrs('L', n2, n3, prec%mgw_fact22, n2, x23, n2, lapack_info)
    if (lapack_info /= 0) then
      info = psb_err_invalid_mat_state_
      return
    end if

    allocate(prec%mgw_fact_s(n3,n3), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      return
    end if
    prec%mgw_fact_s = matmul(b31, x13) + matmul(b32, x23) - b33
    ! Remove roundoff-level asymmetry before Cholesky factorization.
    prec%mgw_fact_s = 0.5_psb_dpk_ * &
         & (prec%mgw_fact_s + transpose(prec%mgw_fact_s))
    call dpotrf('L', n3, prec%mgw_fact_s, n3, lapack_info)
    if (lapack_info /= 0) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'amg_nested_mgw_exact_build', &
           & a_err='exact MGW Schur complement is not SPD')
      return
    end if
    prec%mgw_exact_built = .true.
  end subroutine amg_d_build_mgw_exact

  subroutine amg_d_block_to_dense(block, nrow, ncol, dense, info)
    type(psb_dspmat_type), intent(in) :: block
    integer, intent(in) :: nrow, ncol
    real(psb_dpk_), allocatable, intent(out) :: dense(:,:)
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), allocatable :: ia(:), ja(:)
    real(psb_dpk_), allocatable :: val(:)
    integer(psb_ipk_) :: nz
    integer :: k

    info = psb_success_
    allocate(dense(nrow,ncol), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      return
    end if
    dense(:,:) = dzero
    call psb_d_csgetrow(1, nrow, block, nz, ia, ja, val, info)
    if (info /= psb_success_) return
    do k = 1, nz
      if ((ia(k) < 1) .or. (ia(k) > nrow) .or. &
           & (ja(k) < 1) .or. (ja(k) > ncol)) then
        info = psb_err_invalid_mat_state_
        call psb_errpush(info, 'amg_nested_block_to_dense', &
             & a_err='nonlocal index in serial block')
        return
      end if
      dense(ia(k),ja(k)) = dense(ia(k),ja(k)) + val(k)
    end do
  end subroutine amg_d_block_to_dense

  subroutine amg_d_assemble_shifted_matrix(prec, stiffness, mass, desc, shift, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    type(psb_dspmat_type), intent(in) :: stiffness, mass
    type(psb_desc_type), intent(inout) :: desc
    real(psb_dpk_), intent(in) :: shift
    integer(psb_ipk_), intent(out) :: info
    type(psb_desc_type), pointer :: desc3
    integer(psb_ipk_), allocatable :: ia(:), ja(:)
    integer(psb_lpk_), allocatable :: gia(:), gja(:), owned(:)
    real(psb_dpk_), allocatable :: val(:), vals(:)
    integer(psb_ipk_) :: nz, nk, k, nrows
    info = psb_success_
    nrows = desc%get_local_rows()
    desc3 => psb_d_nest_get_field_desc(prec%nest_op, 3)
    owned = desc%get_global_indices(owned=.true.)
    call psb_cdall(prec%ctxt, prec%shifted_desc, info, vl=owned)
    if (info /= psb_success_) return
    call psb_d_csgetrow(1, nrows, stiffness, nz, ia, ja, val, info)
    if (info /= psb_success_) return
    nk = nz
    allocate(gia(nk+mass%get_nzeros()), gja(nk+mass%get_nzeros()), &
         & vals(nk+mass%get_nzeros()), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; return
    end if
    do k = 1, nk
      call desc%l2g(ia(k), gia(k), info)
      if (info /= psb_success_) return
      call desc3%l2g(ja(k), gja(k), info)
      if (info /= psb_success_) return
      vals(k) = val(k) * prec%schur_active(ia(k)) * prec%schur_active(ja(k))
    end do
    call psb_d_csgetrow(1, nrows, mass, nz, ia, ja, val, info)
    if (info /= psb_success_) return
    do k = 1, nz
      call desc%l2g(ia(k), gia(nk+k), info)
      if (info /= psb_success_) return
      call desc%l2g(ja(k), gja(nk+k), info)
      if (info /= psb_success_) return
      vals(nk+k) = shift * val(k)
    end do
    nz = nk + nz
    call psb_cdins(nz, gja(1:nz), prec%shifted_desc, info)
    if (info /= psb_success_) return
    call psb_cdasb(prec%shifted_desc, info)
    if (info /= psb_success_) return
    call psb_spall(prec%shifted_mat, prec%shifted_desc, info, nnz=nz)
    if (info /= psb_success_) return
    call psb_spins(nz, gia, gja, vals, prec%shifted_mat, prec%shifted_desc, info)
    if (info /= psb_success_) return
    call psb_spasb(prec%shifted_mat, prec%shifted_desc, info, dupl=psb_dupl_add_)
  end subroutine amg_d_assemble_shifted_matrix

  subroutine amg_d_apply_shifted_amg(prec, rhs, sol, desc, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(in) :: rhs(:)
    real(psb_dpk_), intent(out) :: sol(:)
    type(psb_desc_type), intent(in) :: desc
    integer(psb_ipk_), intent(out) :: info
    type(psb_d_vect_type) :: xv, yv
    integer(psb_ipk_) :: nloc, linfo

    info = psb_success_
    nloc = prec%shifted_desc%get_local_rows()
    sol(:) = dzero
    call psb_geall(xv, prec%shifted_desc, info)
    if (info == psb_success_) call psb_geall(yv, prec%shifted_desc, info)
    if (info == psb_success_) call psb_geasb(xv, prec%shifted_desc, info)
    if (info == psb_success_) call psb_geasb(yv, prec%shifted_desc, info)
    if (info == psb_success_) call xv%zero()
    if (info == psb_success_) call yv%zero()
    if (info == psb_success_) call xv%set(rhs(1:min(size(rhs), nloc)))
    if (info == psb_success_) &
         & call prec%schur_amg%apply(xv, yv, prec%shifted_desc, info, trans='N')
    if (info == psb_success_) sol(1:nloc) = yv%get_vect(nloc)
    call psb_gefree(xv, prec%shifted_desc, linfo)
    call psb_gefree(yv, prec%shifted_desc, linfo)
  end subroutine amg_d_apply_shifted_amg

  subroutine amg_d_nested_block_apply2v(prec, x, y, desc_data, info, trans, work)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    type(psb_desc_type), intent(in) :: desc_data
    real(psb_dpk_), intent(inout) :: x(:)
    real(psb_dpk_), intent(inout) :: y(:)
    integer(psb_ipk_), intent(out) :: info
    character(len=1), optional :: trans
    real(psb_dpk_), intent(inout), optional, target :: work(:)

    select case (trim(prec%mode))
    case ('STOKES')
      call amg_d_apply_stokes(prec, x, y, desc_data, info)
    case ('STOKES_MASS')
      call amg_d_apply_stokes_mass(prec, x, y, desc_data, info)
    case ('KKT')
      call amg_d_apply_kkt(prec, x, y, desc_data, info)
    case ('KKT_DIAG')
      call amg_d_apply_kkt_block_diag(prec, x, y, desc_data, info)
    case ('MGW_EXACT')
      call amg_d_apply_mgw_exact(prec, x, y, desc_data, info)
    case default
      info = psb_err_invalid_input_
      call psb_errpush(info, 'amg_nested_block_apply', a_err='unknown mode')
    end select
  end subroutine amg_d_nested_block_apply2v

  subroutine amg_d_apply_stokes(prec, x, y, desc_data, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(inout) :: x(:), y(:)
    type(psb_desc_type), intent(in) :: desc_data
    integer(psb_ipk_), intent(out) :: info
    type(psb_desc_type), pointer :: d1, d2
    real(psb_dpk_), allocatable :: r1(:), r2(:), z1(:), z2(:), t2(:), c1(:), yh(:)
    integer(psb_ipk_) :: n1, n2

    info = psb_success_
    d1 => psb_d_nest_get_field_desc(prec%nest_op, 1)
    d2 => psb_d_nest_get_field_desc(prec%nest_op, 2)
    if (.not. associated(d1) .or. .not. associated(d2)) then
      info = psb_err_invalid_mat_state_; return
    end if
    n1 = d1%get_local_cols(); n2 = d2%get_local_cols()
    allocate(r1(n1), r2(n2), z1(n1), z2(n2), t2(n2), c1(n1), yh(size(y)), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; return
    end if
    y(:) = dzero; yh(:) = dzero
    call psb_d_nest_restrict_field(prec%nest_op, 1, x, r1, info); if (info /= psb_success_) goto 100
    call psb_d_nest_restrict_field(prec%nest_op, 2, x, r2, info); if (info /= psb_success_) goto 100
    z1(:) = dzero
    call amg_d_apply_field(prec, 1, r1, z1, d1, info)
    if (info /= psb_success_) goto 100
    t2(:) = r2
    call amg_d_apply_block(prec%nest_op, 2, 1, -done, z1, done, t2, info)
    if (info /= psb_success_) goto 100
    call amg_d_schur_solve(prec, t2, z2, info); if (info /= psb_success_) goto 100
    call psb_d_nest_prolong_field(prec%nest_op, 2, z2, yh, info); if (info /= psb_success_) goto 100
    call amg_d_apply_block(prec%nest_op, 1, 2, done, z2, dzero, c1, info)
    if (info /= psb_success_) goto 100
    call amg_d_apply_field(prec, 1, c1, r1, d1, info)
    if (info /= psb_success_) goto 100
    z1(:) = z1(:) - r1(:)
    call psb_d_nest_prolong_field(prec%nest_op, 1, z1, yh, info)
    if (info == psb_success_) y(:) = yh(:)
100 continue


    deallocate(r1, r2, z1, z2, t2, c1, yh)
  end subroutine amg_d_apply_stokes

  subroutine amg_d_apply_stokes_mass(prec, x, y, desc_data, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(inout) :: x(:), y(:)
    type(psb_desc_type), intent(in) :: desc_data
    integer(psb_ipk_), intent(out) :: info
    type(psb_desc_type), pointer :: d1, d2
    real(psb_dpk_), allocatable :: r1(:), r2(:), z1(:), z2(:)
    integer(psb_ipk_) :: n1, n2

    info = psb_success_
    d1 => psb_d_nest_get_field_desc(prec%nest_op, 1)
    d2 => psb_d_nest_get_field_desc(prec%nest_op, 2)
    if ((.not. associated(d1)) .or. (.not. associated(d2)) .or. &
         & (.not. allocated(prec%schur_diag))) then
      info = psb_err_invalid_mat_state_
      return
    end if
    n1 = d1%get_local_cols()
    n2 = d2%get_local_cols()
    if (size(prec%schur_diag) /= n2) then
      info = psb_err_invalid_mat_state_
      return
    end if
    allocate(r1(n1), r2(n2), z1(n1), z2(n2), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      return
    end if
    y(:) = dzero
    call psb_d_nest_restrict_field(prec%nest_op, 1, x, r1, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_restrict_field(prec%nest_op, 2, x, r2, info)
    if (info /= psb_success_) goto 100
    call amg_d_apply_field(prec, 1, r1, z1, d1, info)
    if (info /= psb_success_) goto 100
    z2(:) = r2(:) / max(prec%schur_diag(:), sqrt(tiny(done)))
    call psb_d_nest_prolong_field(prec%nest_op, 1, z1, y, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_prolong_field(prec%nest_op, 2, z2, y, info)
100 continue
    deallocate(r1, r2, z1, z2)
  end subroutine amg_d_apply_stokes_mass

  subroutine amg_d_apply_kkt(prec, x, y, desc_data, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(inout) :: x(:), y(:)
    type(psb_desc_type), intent(in) :: desc_data
    integer(psb_ipk_), intent(out) :: info
    type(psb_desc_type), pointer :: d1, d2, d3
    real(psb_dpk_), allocatable :: r1(:), r2(:), r3(:), z1(:), z2(:), z3(:), t3(:), c1(:), c2(:)
    integer(psb_ipk_) :: n1, n2, n3

    info = psb_success_
    d1 => psb_d_nest_get_field_desc(prec%nest_op, 1)
    d2 => psb_d_nest_get_field_desc(prec%nest_op, 2)
    d3 => psb_d_nest_get_field_desc(prec%nest_op, 3)
    if (.not. associated(d1) .or. .not. associated(d2) .or. .not. associated(d3)) then
      info = psb_err_invalid_mat_state_; return
    end if
    n1 = d1%get_local_cols(); n2 = d2%get_local_cols(); n3 = d3%get_local_cols()
    allocate(r1(n1), r2(n2), r3(n3), z1(n1), z2(n2), z3(n3), t3(n3), c1(n1), c2(n2), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; return
    end if
    y(:) = dzero
    call psb_d_nest_restrict_field(prec%nest_op, 1, x, r1, info); if (info /= psb_success_) goto 100
    call psb_d_nest_restrict_field(prec%nest_op, 2, x, r2, info); if (info /= psb_success_) goto 100
    call psb_d_nest_restrict_field(prec%nest_op, 3, x, r3, info); if (info /= psb_success_) goto 100
    z1(:) = dzero; z2(:) = dzero
    call amg_d_apply_field(prec, 1, r1, z1, d1, info); if (info /= psb_success_) goto 100
    call amg_d_apply_field(prec, 2, r2, z2, d2, info); if (info /= psb_success_) goto 100
    t3(:) = r3
    call amg_d_apply_block(prec%nest_op, 3, 1, -done, z1, done, t3, info); if (info /= psb_success_) goto 100
    call amg_d_apply_block(prec%nest_op, 3, 2, -done, z2, done, t3, info); if (info /= psb_success_) goto 100
    call amg_d_schur_solve(prec, t3, z3, info); if (info /= psb_success_) goto 100
    call amg_d_apply_block(prec%nest_op, 1, 3, done, z3, dzero, c1, info); if (info /= psb_success_) goto 100
    call amg_d_apply_field(prec, 1, c1, r1, d1, info); if (info /= psb_success_) goto 100
    z1(:) = z1(:) - r1(:)
    call amg_d_apply_block(prec%nest_op, 2, 3, done, z3, dzero, c2, info); if (info /= psb_success_) goto 100
    call amg_d_apply_field(prec, 2, c2, r2, d2, info); if (info /= psb_success_) goto 100
    z2(:) = z2(:) - r2(:)
    call psb_d_nest_prolong_field(prec%nest_op, 1, z1, y, info); if (info /= psb_success_) goto 100
    call psb_d_nest_prolong_field(prec%nest_op, 2, z2, y, info); if (info /= psb_success_) goto 100
    call psb_d_nest_prolong_field(prec%nest_op, 3, z3, y, info)
100 continue


    deallocate(r1, r2, r3, z1, z2, z3, t3, c1, c2)
  end subroutine amg_d_apply_kkt

  ! Apply diag(M^{-1}, (alpha M)^{-1}, S_alpha^{-1}) with
  ! S_alpha^{-1}=H^{-1} M H^{-1}.  There are no triangular coupling corrections,
  ! so this is the symmetric block-diagonal form intended for MINRES.
  subroutine amg_d_apply_kkt_block_diag(prec, x, y, desc_data, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(inout) :: x(:), y(:)
    type(psb_desc_type), intent(in) :: desc_data
    integer(psb_ipk_), intent(out) :: info
    type(psb_desc_type), pointer :: d1, d2, d3
    real(psb_dpk_), allocatable :: r1(:), r2(:), r3(:), z1(:), z2(:), &
         & z3(:), ktmp(:), mtmp(:)
    integer(psb_ipk_) :: n1, n2, n3

    info = psb_success_
    d1 => psb_d_nest_get_field_desc(prec%nest_op, 1)
    d2 => psb_d_nest_get_field_desc(prec%nest_op, 2)
    d3 => psb_d_nest_get_field_desc(prec%nest_op, 3)
    if ((.not. associated(d1)) .or. (.not. associated(d2)) .or. &
         & (.not. associated(d3)) .or. (.not. allocated(prec%schur_diag)) .or. &
         & (.not. allocated(prec%mass_diag1)) .or. &
         & (.not. allocated(prec%mass_diag2)) .or. &
         & (.not. allocated(prec%schur_active)) .or. &
         & (.not. prec%schur_amg_built)) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'amg_nested_kkt_diag_apply', &
           & a_err='block-diagonal Schur preconditioner not built')
      return
    end if
    n1 = d1%get_local_cols()
    n2 = d2%get_local_cols()
    n3 = d3%get_local_cols()
    if ((d1%get_local_rows() /= d3%get_local_rows()) .or. &
         & (d2%get_local_rows() /= d3%get_local_rows())) then
      info = psb_err_invalid_input_
      return
    end if
    allocate(r1(n1), r2(n2), r3(n3), z1(n1), z2(n2), z3(n3), &
         & ktmp(n1), mtmp(n1), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      return
    end if

    y(:) = dzero
    r1(:) = dzero; r2(:) = dzero; r3(:) = dzero
    z1(:) = dzero; z2(:) = dzero; z3(:) = dzero
    ktmp(:) = dzero; mtmp(:) = dzero
    call psb_d_nest_restrict_field(prec%nest_op, 1, x, r1, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_restrict_field(prec%nest_op, 2, x, r2, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_restrict_field(prec%nest_op, 3, x, r3, info)
    if (info /= psb_success_) goto 100

    z1(1:size(prec%mass_diag1)) = r1(1:size(prec%mass_diag1)) / max(prec%mass_diag1, sqrt(tiny(done)))
    z2(1:size(prec%mass_diag2)) = r2(1:size(prec%mass_diag2)) / max(prec%mass_diag2, sqrt(tiny(done)))
    call amg_d_apply_shifted_amg(prec, r3, ktmp, d1, info)
    if (info /= psb_success_) goto 100
    call amg_d_apply_block(prec%nest_op, 1, 1, done, ktmp, &
         & dzero, mtmp, info)
    if (info /= psb_success_) goto 100
    call amg_d_apply_shifted_amg(prec, mtmp, z3, d1, info)
    if (info /= psb_success_) goto 100

    call psb_d_nest_prolong_field(prec%nest_op, 1, z1, y, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_prolong_field(prec%nest_op, 2, z2, y, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_prolong_field(prec%nest_op, 3, z3, y, info)

100 continue
    deallocate(r1, r2, r3, z1, z2, z3, ktmp, mtmp)
  end subroutine amg_d_apply_kkt_block_diag

  subroutine amg_d_apply_mgw_exact(prec, x, y, desc_data, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    real(psb_dpk_), intent(inout) :: x(:), y(:)
    type(psb_desc_type), intent(in) :: desc_data
    integer(psb_ipk_), intent(out) :: info
    type(psb_desc_type), pointer :: d1, d2, d3
    real(psb_dpk_), allocatable :: r1(:), r2(:), r3(:)
    real(psb_dpk_), allocatable :: z1(:,:), z2(:,:), z3(:,:)
    integer :: n1, n2, n3, lapack_info

    info = psb_success_
    if (.not. prec%mgw_exact_built) then
      info = psb_err_invalid_mat_state_
      call psb_errpush(info, 'amg_nested_mgw_exact_apply', &
           & a_err='MGW_EXACT factors are not built')
      return
    end if
    d1 => psb_d_nest_get_field_desc(prec%nest_op, 1)
    d2 => psb_d_nest_get_field_desc(prec%nest_op, 2)
    d3 => psb_d_nest_get_field_desc(prec%nest_op, 3)
    if ((.not. associated(d1)) .or. (.not. associated(d2)) .or. &
         & (.not. associated(d3))) then
      info = psb_err_invalid_mat_state_
      return
    end if
    n1 = d1%get_local_cols()
    n2 = d2%get_local_cols()
    n3 = d3%get_local_cols()
    allocate(r1(n1), r2(n2), r3(n3), z1(n1,1), z2(n2,1), z3(n3,1), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      return
    end if

    call psb_d_nest_restrict_field(prec%nest_op, 1, x, r1, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_restrict_field(prec%nest_op, 2, x, r2, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_restrict_field(prec%nest_op, 3, x, r3, info)
    if (info /= psb_success_) goto 100
    z1(:,1) = r1
    z2(:,1) = r2
    z3(:,1) = r3
    call dpotrs('L', n1, 1, prec%mgw_fact11, n1, z1, n1, lapack_info)
    if (lapack_info /= 0) then
      info = psb_err_invalid_mat_state_
      goto 100
    end if
    call dpotrs('L', n2, 1, prec%mgw_fact22, n2, z2, n2, lapack_info)
    if (lapack_info /= 0) then
      info = psb_err_invalid_mat_state_
      goto 100
    end if
    call dpotrs('L', n3, 1, prec%mgw_fact_s, n3, z3, n3, lapack_info)
    if (lapack_info /= 0) then
      info = psb_err_invalid_mat_state_
      goto 100
    end if

    y(:) = dzero
    call psb_d_nest_prolong_field(prec%nest_op, 1, z1(:,1), y, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_prolong_field(prec%nest_op, 2, z2(:,1), y, info)
    if (info /= psb_success_) goto 100
    call psb_d_nest_prolong_field(prec%nest_op, 3, z3(:,1), y, info)

100 continue
    deallocate(r1, r2, r3, z1, z2, z3)
  end subroutine amg_d_apply_mgw_exact

  subroutine amg_d_nested_block_apply1v(prec, x, desc_data, info, trans)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    type(psb_desc_type), intent(in) :: desc_data
    real(psb_dpk_), intent(inout) :: x(:)
    integer(psb_ipk_), intent(out) :: info
    character(len=1), optional :: trans
    real(psb_dpk_), allocatable :: y(:)
    allocate(y(size(x)), stat=info)
    if (info /= 0) return
    call prec%psb_d_apply2v(x, y, desc_data, info, trans)
    if (info == psb_success_) x(:) = y(:)
    deallocate(y)
  end subroutine amg_d_nested_block_apply1v

  subroutine amg_d_nested_block_apply2_vect(prec, x, y, desc_data, info, trans)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    type(psb_desc_type), intent(in) :: desc_data
    type(psb_d_vect_type), intent(inout) :: x, y
    integer(psb_ipk_), intent(out) :: info
    character(len=1), optional :: trans
    real(psb_dpk_), allocatable :: xb(:), yb(:)
    integer(psb_ipk_) :: ncol
    ncol = desc_data%get_local_cols()
    xb = x%get_vect(ncol)
    allocate(yb(max(y%get_nrows(), ncol)), stat=info)
    if (info /= 0) return
    yb(:) = dzero
    if (present(trans)) then
      call amg_d_nested_block_apply2v(prec, xb, yb, desc_data, info, trans)
    else
      call amg_d_nested_block_apply2v(prec, xb, yb, desc_data, info)
    end if
    if (info == psb_success_) call y%set(yb(1:y%get_nrows()))
    deallocate(yb)
  end subroutine amg_d_nested_block_apply2_vect

  subroutine amg_d_nested_block_apply1_vect(prec, x, desc_data, info, trans)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    type(psb_desc_type), intent(in) :: desc_data
    type(psb_d_vect_type), intent(inout) :: x
    integer(psb_ipk_), intent(out) :: info
    character(len=1), optional :: trans
    type(psb_d_vect_type) :: y
    real(psb_dpk_), allocatable :: yb(:)
    integer(psb_ipk_) :: ncol
    ncol = desc_data%get_local_cols()
    allocate(yb(ncol), stat=info)
    if (info /= 0) return
    yb(:) = dzero
    call psb_geall(y, desc_data, info)
    if (info == psb_success_) call psb_geasb(y, desc_data, info)
    if (info == psb_success_) call prec%psb_d_apply2_vect(x, y, desc_data, info, trans)
    if (info == psb_success_) then
      yb = y%get_vect(ncol)
      call x%set(yb(1:x%get_nrows()))
    end if
    call psb_gefree(y, desc_data, info)
    deallocate(yb)
  end subroutine amg_d_nested_block_apply1_vect

  subroutine amg_d_nested_block_free(prec, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_) :: i, linfo
    info = psb_success_
    if (allocated(prec%field_amg)) then
      do i = 1, size(prec%field_amg)
        call prec%field_amg(i)%free(linfo)
        if (linfo /= psb_success_ .and. info == psb_success_) info = linfo
      end do
      deallocate(prec%field_amg)
    end if
    if (allocated(prec%use_amg)) deallocate(prec%use_amg)
    if (prec%schur_amg_built) then
      call prec%schur_amg%free(linfo)
      if (linfo /= psb_success_ .and. info == psb_success_) info = linfo
    end if
    prec%schur_amg_built = .false.
    call prec%shifted_mat%free()
    call prec%shifted_desc%free(linfo)
    if (info == psb_success_) info = linfo
    if (allocated(prec%schur_diag)) deallocate(prec%schur_diag)
    if (allocated(prec%mass_diag1)) deallocate(prec%mass_diag1)
    if (allocated(prec%mass_diag2)) deallocate(prec%mass_diag2)
    if (allocated(prec%schur_active)) deallocate(prec%schur_active)
    prec%schur_solve = 'MATRIX_FREE'
    prec%schur_maxit = 8
    if (allocated(prec%mgw_fact11)) deallocate(prec%mgw_fact11)
    if (allocated(prec%mgw_fact22)) deallocate(prec%mgw_fact22)
    if (allocated(prec%mgw_fact_s)) deallocate(prec%mgw_fact_s)
    prec%mgw_exact_built = .false.
    prec%schur_tol = dzero
    prec%nest_op => null()
    prec%nfields = 0
    prec%wrk_allocated = .false.
  end subroutine amg_d_nested_block_free

  function amg_d_nested_block_sizeof(prec, global) result(val)
    class(amg_d_nested_block_prec_type), intent(in) :: prec
    logical, intent(in), optional :: global
    integer(psb_epk_) :: val
    integer(psb_ipk_) :: i
    val = 0_psb_epk_
    if (allocated(prec%field_amg)) then
      do i = 1, size(prec%field_amg)
        val = val + prec%field_amg(i)%sizeof()
      end do
    end if
  end function amg_d_nested_block_sizeof

end module amg_d_nested_block_prec_mod
