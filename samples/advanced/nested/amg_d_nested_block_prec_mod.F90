module amg_d_nested_block_prec_mod
  use psb_base_mod, only : psb_ipk_, psb_epk_, psb_dpk_, psb_success_, &
       & psb_err_invalid_input_, psb_err_invalid_mat_state_, psb_err_alloc_dealloc_, &
       & psb_errpush, psb_toupper, psb_geall, psb_geasb, psb_gefree, done, dzero, psb_ctxt_type
  use psb_prec_mod, only : psb_dprec_type
  use psb_d_mat_mod, only : psb_dspmat_type, psb_d_get_diag
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
    logical, allocatable :: use_amg(:)
    real(psb_dpk_), allocatable :: schur_diag(:)
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

contains

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
    case ('AMG_KKT','KKT','OPTIMAL','OPTIMAL_CONTROL')
      prec%mode = 'KKT'
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
      if (i < prec%nfields) then
        call prec%field_amg(i)%init(fdesc%get_ctxt(), 'ML', info)
        if (info /= psb_success_) return
        call prec%field_amg(i)%build(blk, fdesc, info, amold=amold, vmold=vmold, imold=imold)
        if (info /= psb_success_) return
        prec%use_amg(i) = .true.
      end if
    end do

    select case (trim(prec%mode))
    case ('STOKES')
      call amg_d_build_stokes_schur_diag(prec, info)
    case ('KKT')
      call amg_d_build_kkt_schur_diag(prec, info)
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
    nloc = desc%get_local_cols()
    call psb_geall(xv, desc, info)
    if (info == psb_success_) call psb_geall(yv, desc, info)
    if (info == psb_success_) call psb_geasb(xv, desc, info)
    if (info == psb_success_) call psb_geasb(yv, desc, info)
    if (info == psb_success_) call xv%zero()
    if (info == psb_success_) call yv%zero()
    if (info == psb_success_) call xv%set(rhs(1:min(size(rhs), nloc)))
    if (info == psb_success_) call prec%field_amg(ifld)%apply(xv, yv, desc, info, trans='N')
    if (info == psb_success_) sol(:) = yv%get_vect(size(sol))
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
      call psb_d_nest_apply_block(prec%nest_op, 1, 2, done, x, dzero, t1, info)
      if (info == psb_success_) call amg_d_apply_field(prec, 1, t1, w1, d1, info)
      if (info == psb_success_) call psb_d_nest_apply_block(prec%nest_op, 2, 1, -done, w1, dzero, y, info)
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
      call psb_d_nest_apply_block(prec%nest_op, 1, 3, done, x, dzero, t1, info)
      if (info == psb_success_) call amg_d_apply_field(prec, 1, t1, w1, d1, info)
      if (info == psb_success_) call psb_d_nest_apply_block(prec%nest_op, 3, 1, -done, w1, dzero, y, info)
      if (info == psb_success_) call psb_d_nest_apply_block(prec%nest_op, 2, 3, done, x, dzero, t2, info)
      if (info == psb_success_) call amg_d_apply_field(prec, 2, t2, w2, d2, info)
      if (info == psb_success_) call psb_d_nest_apply_block(prec%nest_op, 3, 2, -done, w2, done, y, info)
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
    integer(psb_ipk_) :: k, n
    real(psb_dpk_) :: alpha, beta, omega, rho, rho_old, denom, rhsn, rn, floorv

    info = psb_success_
    if (psb_toupper(trim(prec%schur_solve)) == 'SELF') then
      call amg_d_schur_diag_apply(prec, rhs, sol)
      if (trim(prec%mode) == 'STOKES') sol(:) = -sol(:)
      return
    end if

    n = size(rhs)
    allocate(r(n), rr(n), p(n), v(n), s(n), t(n), ph(n), sh(n), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; return
    end if

    floorv = sqrt(tiny(done))
    sol(:) = dzero
    r(:) = rhs(:)
    rr(:) = r(:)
    p(:) = dzero
    v(:) = dzero
    alpha = done
    omega = done
    rho_old = done
    rhsn = sqrt(sum(rhs(:) * rhs(:)))
    if (rhsn <= floorv) goto 100

    do k = 1, max(1, prec%schur_maxit)
      rho = sum(rr(:) * r(:))
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
      denom = sum(rr(:) * v(:))
      if (abs(denom) <= floorv) exit
      alpha = rho / denom
      s(:) = r(:) - alpha * v(:)
      rn = sqrt(sum(s(:) * s(:)))
      if (rn <= max(prec%schur_tol * rhsn, floorv)) then
        sol(:) = sol(:) + alpha * ph(:)
        exit
      end if

      call amg_d_schur_diag_apply(prec, s, sh)
      call amg_d_schur_action(prec, sh, t, info)
      if (info /= psb_success_) exit
      denom = sum(t(:) * t(:))
      if (abs(denom) <= floorv) exit
      omega = sum(t(:) * s(:)) / denom
      sol(:) = sol(:) + alpha * ph(:) + omega * sh(:)
      r(:) = s(:) - omega * t(:)
      rn = sqrt(sum(r(:) * r(:)))
      if (rn <= max(prec%schur_tol * rhsn, floorv)) exit
      rho_old = rho
    end do

100 continue
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
    type(psb_dspmat_type), pointer :: a11, a12, a21
    type(psb_desc_type), pointer :: descu, descp
    real(psb_dpk_), allocatable :: d11(:), ep(:), wu(:), vp(:)
    integer(psb_ipk_) :: n, nu, k
    real(psb_dpk_) :: floorv

    info = psb_success_
    descu => psb_d_nest_get_field_desc(prec%nest_op, 1)
    descp => psb_d_nest_get_field_desc(prec%nest_op, 2)
    if (.not. associated(descu) .or. .not. associated(descp)) then
      info = psb_err_invalid_mat_state_; return
    end if
    n = descp%get_local_cols(); nu = descu%get_local_cols()
    allocate(prec%schur_diag(n), ep(n), wu(nu), vp(n), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; return
    end if
    prec%schur_diag(:) = done

    a11 => psb_d_nest_get_block(prec%nest_op, 1, 1)
    a12 => psb_d_nest_get_block(prec%nest_op, 1, 2)
    a21 => psb_d_nest_get_block(prec%nest_op, 2, 1)
    if (associated(a11) .and. associated(a12) .and. associated(a21)) then
      d11 = psb_d_get_diag(a11, info); if (info /= psb_success_) return
      floorv = sqrt(tiny(done))
      do k = 1, n
        ep(:) = dzero; wu(:) = dzero; vp(:) = dzero
        ep(k) = done
        call psb_d_nest_apply_block(prec%nest_op, 1, 2, done, ep, dzero, wu, info)
        if (info /= psb_success_) return
        wu(:) = wu(:) / max(abs(d11(1:min(size(d11), size(wu)))), floorv)
        call psb_d_nest_apply_block(prec%nest_op, 2, 1, done, wu, dzero, vp, info)
        if (info /= psb_success_) return
        prec%schur_diag(k) = max(abs(vp(k)), floorv)
      end do
    end if
  end subroutine amg_d_build_stokes_schur_diag

  subroutine amg_d_build_kkt_schur_diag(prec, info)
    class(amg_d_nested_block_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out) :: info
    type(psb_dspmat_type), pointer :: a11, a22, a13, a23, a31, a32
    type(psb_desc_type), pointer :: desc1, desc2, desc3
    real(psb_dpk_), allocatable :: d11(:), d22(:), e3(:), w1(:), w2(:), v3(:)
    integer(psb_ipk_) :: n1, n2, n3, k
    real(psb_dpk_) :: floorv

    info = psb_success_
    desc1 => psb_d_nest_get_field_desc(prec%nest_op, 1)
    desc2 => psb_d_nest_get_field_desc(prec%nest_op, 2)
    desc3 => psb_d_nest_get_field_desc(prec%nest_op, 3)
    if (.not. associated(desc1) .or. .not. associated(desc2) .or. .not. associated(desc3)) then
      info = psb_err_invalid_mat_state_; return
    end if
    n1 = desc1%get_local_cols(); n2 = desc2%get_local_cols(); n3 = desc3%get_local_cols()
    allocate(prec%schur_diag(n3), e3(n3), w1(n1), w2(n2), v3(n3), stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_; return
    end if
    prec%schur_diag(:) = done

    a11 => psb_d_nest_get_block(prec%nest_op, 1, 1)
    a22 => psb_d_nest_get_block(prec%nest_op, 2, 2)
    a13 => psb_d_nest_get_block(prec%nest_op, 1, 3)
    a23 => psb_d_nest_get_block(prec%nest_op, 2, 3)
    a31 => psb_d_nest_get_block(prec%nest_op, 3, 1)
    a32 => psb_d_nest_get_block(prec%nest_op, 3, 2)
    if (associated(a11) .and. associated(a22) .and. associated(a13) .and. &
        & associated(a23) .and. associated(a31) .and. associated(a32)) then
      d11 = psb_d_get_diag(a11, info); if (info /= psb_success_) return
      d22 = psb_d_get_diag(a22, info); if (info /= psb_success_) return
      floorv = sqrt(tiny(done))
      do k = 1, n3
        e3(:) = dzero; w1(:) = dzero; w2(:) = dzero; v3(:) = dzero
        e3(k) = done
        call psb_d_nest_apply_block(prec%nest_op, 1, 3, done, e3, dzero, w1, info)
        if (info /= psb_success_) return
        w1(:) = w1(:) / max(abs(d11(1:min(size(d11), size(w1)))), floorv)
        call psb_d_nest_apply_block(prec%nest_op, 3, 1, done, w1, done, v3, info)
        if (info /= psb_success_) return
        call psb_d_nest_apply_block(prec%nest_op, 2, 3, done, e3, dzero, w2, info)
        if (info /= psb_success_) return
        w2(:) = w2(:) / max(abs(d22(1:min(size(d22), size(w2)))), floorv)
        call psb_d_nest_apply_block(prec%nest_op, 3, 2, done, w2, done, v3, info)
        if (info /= psb_success_) return
        prec%schur_diag(k) = max(abs(v3(k)), floorv)
      end do
    end if
  end subroutine amg_d_build_kkt_schur_diag
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
    case ('KKT')
      call amg_d_apply_kkt(prec, x, y, desc_data, info)
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
    call psb_d_nest_apply_block(prec%nest_op, 2, 1, -done, z1, done, t2, info)
    if (info /= psb_success_) goto 100
    call amg_d_schur_solve(prec, t2, z2, info); if (info /= psb_success_) goto 100
    call psb_d_nest_prolong_field(prec%nest_op, 2, z2, yh, info); if (info /= psb_success_) goto 100
    call psb_d_nest_apply_block(prec%nest_op, 1, 2, done, z2, dzero, c1, info)
    if (info /= psb_success_) goto 100
    call amg_d_apply_field(prec, 1, c1, c1, d1, info)
    if (info /= psb_success_) goto 100
    z1(:) = z1(:) - c1(:)
    call psb_d_nest_prolong_field(prec%nest_op, 1, z1, yh, info)
    if (info == psb_success_) y(:) = yh(:)
100 continue


    deallocate(r1, r2, z1, z2, t2, c1, yh)
  end subroutine amg_d_apply_stokes

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
    call psb_d_nest_apply_block(prec%nest_op, 3, 1, -done, z1, done, t3, info); if (info /= psb_success_) goto 100
    call psb_d_nest_apply_block(prec%nest_op, 3, 2, -done, z2, done, t3, info); if (info /= psb_success_) goto 100
    call amg_d_schur_solve(prec, t3, z3, info); if (info /= psb_success_) goto 100
    call psb_d_nest_apply_block(prec%nest_op, 1, 3, done, z3, dzero, c1, info); if (info /= psb_success_) goto 100
    call amg_d_apply_field(prec, 1, c1, c1, d1, info); if (info /= psb_success_) goto 100
    z1(:) = z1(:) - c1(:)
    call psb_d_nest_apply_block(prec%nest_op, 2, 3, done, z3, dzero, c2, info); if (info /= psb_success_) goto 100
    call amg_d_apply_field(prec, 2, c2, c2, d2, info); if (info /= psb_success_) goto 100
    z2(:) = z2(:) - c2(:)
    call psb_d_nest_prolong_field(prec%nest_op, 1, z1, y, info); if (info /= psb_success_) goto 100
    call psb_d_nest_prolong_field(prec%nest_op, 2, z2, y, info); if (info /= psb_success_) goto 100
    call psb_d_nest_prolong_field(prec%nest_op, 3, z3, y, info)
100 continue


    deallocate(r1, r2, r3, z1, z2, z3, t3, c1, c2)
  end subroutine amg_d_apply_kkt

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
    call prec%psb_d_apply2v(xb, yb, desc_data, info, trans)
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
    if (allocated(prec%schur_diag)) deallocate(prec%schur_diag)
    prec%schur_solve = 'MATRIX_FREE'
    prec%schur_maxit = 8
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
