module weno_tensor_product_mod

  use poly_utils_mod
  use math_mod

  implicit none

  private

  public weno_tensor_product_type

  type weno_tensor_product_type
    logical :: initialized = .false.
    integer :: id     = 0                  ! Sub-stencil ID
    integer :: nd     = 0                  ! Dimension number
    integer :: sw     = 0                  ! Stencil width
    integer :: sub_sw = 0                  ! Sub-stencil width
    integer :: nc     = 0                  ! Number of cells
    integer :: ns     = 0                  ! Number of sub-stencils
    integer :: npt    = 0                  ! Number of evaluation points
    integer :: di0    = 0                  ! Origin shift of sub-stencil relative to stencil
    integer :: dj0    = 0                  ! Origin shift of sub-stencil relative to stencil
    integer :: di     = 0                  ! Index shift of sub-stencil relative to stencil
    integer :: dj     = 0                  ! Index shift of sub-stencil relative to stencil
    integer, allocatable :: cell_mask(:,:) ! Mask unavailable cells by 0
    integer, allocatable :: poly_mask(:,:) ! Mask unavailable polynomial terms by 0
    real(8), allocatable :: x        (:)   ! X coordinate of evaluation point
    real(8), allocatable :: y        (:)   ! Y coordinate of evaluation point
    real(8), allocatable :: iA       (:,:) ! A * a = f, iA = inverse(A)
    real(8), allocatable :: iAp      (:,:) ! iA * p
    real(8), allocatable :: ic       (:,:) ! Ideal coefficients for combining sub-stencils (only on stencil)
    real(16), allocatable :: iAp_r16 (:,:) ! Working arrays
    type(weno_tensor_product_type), allocatable :: subs(:) ! Sub-stencils
  contains
    procedure :: init                  => weno_tensor_product_init
    procedure :: add_point             => weno_tensor_product_add_point
    procedure :: calc_recon_matrix     => weno_tensor_product_calc_recon_matrix
    procedure :: calc_ideal_coefs      => weno_tensor_product_calc_ideal_coefs
    procedure :: release_unused_memory => weno_tensor_product_release_unused_memory
    procedure :: clear                 => weno_tensor_product_clear
    final :: weno_tensor_product_final
  end type weno_tensor_product_type

contains

  subroutine weno_tensor_product_init(this, nd, sw, di0, dj0, di, dj, mask, id, has_subs)

    class(weno_tensor_product_type), intent(inout) :: this
    integer, intent(in) :: nd
    integer, intent(in) :: sw
    integer, intent(in), optional :: di0
    integer, intent(in), optional :: dj0
    integer, intent(in), optional :: di
    integer, intent(in), optional :: dj
    integer, intent(in), optional :: mask(:,:) ! Cell mask
    integer, intent(in), optional :: id
    logical, intent(in), optional :: has_subs

    integer i, j, k

    call this%clear()

    this%sw = sw
    this%nd = nd
    this%nc = sw**nd

    allocate(this%cell_mask(sw,sw**(nd-1)))
    allocate(this%poly_mask(sw,sw**(nd-1)))

    if (present(id )) this%id  = id
    if (present(di0)) this%di0 = di0
    if (present(dj0)) this%dj0 = dj0
    if (present(di )) this%di  = di
    if (present(dj )) this%dj  = dj

    if (present(mask)) then
      this%cell_mask = mask
    else
      this%cell_mask = 1
    end if

    ! Set polynomial mask from cell mask.
    if (this%cell_mask(1,1) == 0) then
      this%poly_mask = this%cell_mask(this%sw:1:-1,this%sw**(this%nd-1):1:-1)
    else if (this%cell_mask(this%sw,1) == 0) then
      this%poly_mask = this%cell_mask(:,this%sw**(this%nd-1):1:-1)
    else if (this%cell_mask(1,this%sw**(this%nd-1)) == 0) then
      this%poly_mask = this%cell_mask(this%sw:1:-1,:)
    else
      this%poly_mask = this%cell_mask
    end if

    ! Initialize sub-stencils.
    if (merge(has_subs, .true., present(has_subs))) then
      ! Set sub-stencil size.
      select case (sw)
      case (5)
        this%sub_sw = 3
      end select
      this%ns = this%sub_sw**nd
      allocate(this%subs(this%ns))
      k = 1
      do j = 1, this%sub_sw**(nd - 1)
        do i = 1, this%sub_sw
          call this%subs(k)%init(nd=nd, sw=this%sub_sw, di0=i-int(this%sw/2), dj0=j-int(this%sw/2), &
                                 di=i-1, dj=j-1, id=k, has_subs=.false., &
                                 mask=this%cell_mask(i:i+this%sub_sw-1,j:j+this%sub_sw**(nd-1)-1))
          k = k + 1
        end do
      end do
    end if

    this%initialized = .true.

  end subroutine weno_tensor_product_init

  subroutine weno_tensor_product_add_point(this, x, y)

    class(weno_tensor_product_type), intent(inout) :: this
    real(8), intent(in) :: x
    real(8), intent(in), optional :: y

    real(8), allocatable :: tmp(:)
    integer i, k

    this%npt = this%npt + 1

    allocate(tmp(this%npt))
    do i = 1, this%npt - 1
      tmp(i) = this%x(i)
    end do
    tmp(this%npt) = x - this%di0 ! Convert to local coordinate.
    if (allocated(this%x)) deallocate(this%x)
    allocate(this%x(this%npt))
    this%x = tmp
    if (present(y)) then ! 1D case does not has y.
      do i = 1, this%npt - 1
        tmp(i) = this%y(i)
      end do
      tmp(this%npt) = y - this%dj0 ! Convert to local coordinate.
      if (allocated(this%y)) deallocate(this%y)
      allocate(this%y(this%npt))
      this%y = tmp
    end if
    deallocate(tmp)

    if (allocated(this%subs)) then
      ! Add point to sub-stencils.
      do k = 1, this%ns
        call this%subs(k)%add_point(x, y)
      end do
    end if

  end subroutine weno_tensor_product_add_point

  subroutine weno_tensor_product_calc_recon_matrix(this, ierr)

    class(weno_tensor_product_type), intent(inout) :: this
    integer, intent(out) :: ierr

    ! Local double double arrays for preserving precision.
    real(16), allocatable, dimension(:,:) :: poly, A, iA
    integer , allocatable :: idx_map(:)
    integer i, j, k, p, n

    ierr = 0

    if (allocated(this%iA)) deallocate(this%iA)
    allocate(this%iA(this%nc,this%nc))
    if (allocated(this%iAp)) deallocate(this%iAp)
    allocate(this%iAp(this%nc,this%npt))
    if (allocated(this%iAp_r16)) deallocate(this%iAp_r16)
    allocate(this%iAp_r16(this%nc,this%npt))

    allocate(poly   (this%nc,this%npt)); poly = 0
    allocate( A     (this%nc,this%nc ))
    allocate(iA     (this%nc,this%nc ))
    allocate(idx_map(this%nc         )); idx_map = 0

    ! Set index maps.
    k = 1; n = 1
    do j = 1, this%sw**(this%nd-1)
      do i = 1, this%sw
        if (this%cell_mask(i,j) == 1) then
          idx_map(n) = k; k = k + 1
        end if
        n = n + 1
      end do
    end do

    ! Set the p for each evaluation point.
    ! Select monomials according to mask.
    select case (this%nd)
    case (1)
      do p = 1, this%npt
        k = 1
        do i = 1, this%sw
          if (this%poly_mask(i,1) == 1) then
            call calc_poly_tensor_product_monomial(this%x(p), i - 1, poly(k,p))
            k = k + 1
          end if
        end do
      end do
    case (2)
      do p = 1, this%npt
        k = 1
        do j = 1, this%sw**(this%nd-1)
          do i = 1, this%sw
            if (this%poly_mask(i,j) == 1) then
              call calc_poly_tensor_product_monomial(this%x(p), this%y(p), i - 1, j - 1, poly(k,p))
              k = k + 1
            end if
          end do
        end do
      end do
    end select

    ! Calculate inverse of integral coefficient matrix.
    select case (this%nd)
    case (1)
      call calc_poly_tensor_product_integral_coef_matrix(this%sw, A, this%cell_mask(:,1), this%poly_mask(:,1))
    case (2)
      call calc_poly_tensor_product_integral_coef_matrix(this%sw, this%sw, A, this%cell_mask, this%poly_mask)
    end select
    n = count(this%cell_mask == 1)
    call inverse_matrix(A(1:n,1:n), iA(1:n,1:n), ierr)
    if (ierr /= 0) then
      deallocate(A, iA, idx_map)
      return
    end if

    this%iAp_r16(1:n,:) = matmul(iA(1:n,1:n), poly(1:n,:))

    ! Copy double double iA into double iA.
    this%iA = iA

    ! Rearrange iA and iAp.
    do n = this%nc, 1, -1
      if (idx_map(n) /= 0) then
        this%iA     (n,:) = this%iA     (idx_map(n),:)
        this%iAp_r16(n,:) = this%iAp_r16(idx_map(n),:)
      else
        this%iA     (n,:) = 0
        this%iAp_r16(n,:) = 0
      end if
    end do

    this%iAp = this%iAp_r16

    deallocate(poly, A, iA, idx_map)

  end subroutine weno_tensor_product_calc_recon_matrix

  subroutine weno_tensor_product_calc_ideal_coefs(this, ierr)

    class(weno_tensor_product_type), intent(inout) :: this
    integer, intent(out) :: ierr

    ! Local double double arrays for preserving precision.
    real(16), allocatable, dimension(:,:) :: A, AtA, iAtA, ic
    integer i, j, k, p, m, n

    ierr = 0

    if (.not. allocated(this%subs)) then
      ierr = 1
      return
    end if

    do k = 1, this%ns
      call this%subs(k)%calc_recon_matrix(ierr)
      if (ierr /= 0) return
    end do
    call this%calc_recon_matrix(ierr)

    if (allocated(this%ic)) deallocate(this%ic)
    allocate(this%ic(this%ns,this%npt))

    allocate(   A(this%nc,this%ns )); A = 0
    allocate( AtA(this%ns,this%ns ))
    allocate(iAtA(this%ns,this%ns ))
    allocate(  ic(this%ns,this%npt)); ic = 0
    do p = 1, this%npt
      do k = 1, this%ns
        do j = 1, this%sub_sw**(this%nd-1)
          do i = 1, this%sub_sw
            if (this%subs(k)%cell_mask(i,j) == 1) then
              m = (j + this%subs(k)%dj - 1) * this%sw     + i + (this%subs(k)%di)
              n = (j                   - 1) * this%sub_sw + i
              A(m,k) = this%subs(k)%iAp_r16(n,p)
            end if
          end do
        end do
      end do
      AtA = matmul(transpose(A), A)
      call inverse_matrix(AtA, iAtA, ierr)
      if (ierr /= 0) return
      ic(:,p) = matmul(matmul(iAtA, transpose(A)), this%iAp_r16(:,p))
    end do
    if (abs(sum(matmul(A, ic(:,1)) - this%iAp_r16(:,1))) > 1.0e-30) then
      ierr = 3
    end if
    this%ic = ic
    deallocate(A, AtA, iAtA, ic)

  end subroutine weno_tensor_product_calc_ideal_coefs

  subroutine weno_tensor_product_release_unused_memory(this)

    class(weno_tensor_product_type), intent(inout) :: this

    if (allocated(this%poly_mask)) deallocate(this%poly_mask)
    if (allocated(this%x        )) deallocate(this%x        )
    if (allocated(this%y        )) deallocate(this%y        )
    if (allocated(this%iAp_r16  )) deallocate(this%iAp_r16  )

  end subroutine weno_tensor_product_release_unused_memory

  subroutine weno_tensor_product_clear(this)

    class(weno_tensor_product_type), intent(inout) :: this

    this%id     = 0
    this%nd     = 0
    this%sw     = 0
    this%sub_sw = 0
    this%nc     = 0
    this%ns     = 0
    this%npt    = 0
    this%di0    = 0
    this%dj0    = 0
    this%di     = 0
    this%dj     = 0

    if (allocated(this%cell_mask)) deallocate(this%cell_mask)
    if (allocated(this%poly_mask)) deallocate(this%poly_mask)
    if (allocated(this%x        )) deallocate(this%x        )
    if (allocated(this%y        )) deallocate(this%y        )
    if (allocated(this%iAp      )) deallocate(this%iAp      )
    if (allocated(this%iAp_r16  )) deallocate(this%iAp_r16  )
    if (allocated(this%ic       )) deallocate(this%ic       )
    if (allocated(this%subs     )) deallocate(this%subs     )

    this%initialized = .false.

  end subroutine weno_tensor_product_clear

  subroutine weno_tensor_product_final(this)

    type(weno_tensor_product_type), intent(inout) :: this

    call this%clear()

  end subroutine weno_tensor_product_final

end module weno_tensor_product_mod
