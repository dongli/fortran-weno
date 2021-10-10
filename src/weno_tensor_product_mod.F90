module weno_tensor_product_mod

  use poly_utils_mod
  use math_mod

  implicit none

  private

  public weno_tensor_product_type

  type weno_tensor_product_type
    logical :: initialized = .false.
    integer :: id     = 0             ! Sub-stencil ID
    integer :: nd     = 0             ! Dimension number
    integer :: sw     = 0             ! Stencil width
    integer :: sub_sw = 0             ! Sub-stencil width
    integer :: nc     = 0             ! Number of cells
    integer :: ns     = 0             ! Number of sub-stencils
    integer :: npt    = 0             ! Number of evaluation points
    integer :: di0    = 0
    integer :: dj0    = 0
    integer :: di     = 0
    integer :: dj     = 0
    integer , allocatable :: mask(:,:) ! Mask unavailable nodes by 0
    real(8 ), allocatable :: x(:)      ! X coordinate of evaluation point
    real(8 ), allocatable :: y(:)      ! Y coordinate of evaluation point
    real(16), allocatable :: p(:,:)    ! Monomial terms on each evaluation point
    real(16), allocatable :: iA_p(:,:) ! p * iA
    real(16), allocatable :: ic(:,:)   ! Ideal coefficients for combining sub-stencils (only on stencil)
    type(weno_tensor_product_type), allocatable :: subs(:) ! Sub-stencils
  contains
    procedure :: init              => weno_tensor_product_init
    procedure :: add_point         => weno_tensor_product_add_point
    procedure :: calc_recon_matrix => weno_tensor_product_calc_recon_matrix
    procedure :: calc_ideal_coefs  => weno_tensor_product_calc_ideal_coefs
    procedure :: clear             => weno_tensor_product_clear
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
    integer, intent(in), optional :: mask(:,:)
    integer, intent(in), optional :: id
    logical, intent(in), optional :: has_subs

    integer i, j, k

    call this%clear()

    this%sw = sw
    this%nd = nd
    this%nc = sw**nd

    allocate(this%mask(sw,sw**(nd-1)))

    if (present(id )) this%id  = id
    if (present(di0)) this%di0 = di0
    if (present(dj0)) this%dj0 = dj0
    if (present(di )) this%di  = di
    if (present(dj )) this%dj  = dj

    if (present(mask)) then
      this%mask = mask
    else
      this%mask = 1
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
                                 di=i-1, dj=j-1, id=k, has_subs=.false.)
          this%subs(k)%mask = this%mask(i:i+this%sub_sw-1,j:j+this%sub_sw-1)
          k = k + 1
        end do
      end do
    end if

    this%initialized = .true.

  end subroutine weno_tensor_product_init

  subroutine weno_tensor_product_add_point(this, x, y)

    class(weno_tensor_product_type), intent(inout) :: this
    real(8), intent(in) :: x
    real(8), intent(in) :: y

    real(8), allocatable :: tmp(:)
    integer i, j, k

    this%npt = this%npt + 1

    allocate(tmp(this%npt))
    do i = 1, this%npt - 1
      tmp(i) = this%x(i)
    end do
    tmp(this%npt) = x - this%di0 ! Convert to local coordinate.
    if (allocated(this%x)) deallocate(this%x)
    allocate(this%x(this%npt))
    this%x = tmp
    do i = 1, this%npt - 1
      tmp(i) = this%y(i)
    end do
    tmp(this%npt) = y - this%dj0 ! Convert to local coordinate.
    if (allocated(this%y)) deallocate(this%y)
    allocate(this%y(this%npt))
    this%y = tmp
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

    real(16), allocatable :: A(:,:), iA(:,:)
    integer , allocatable :: kmap(:)
    logical mirror_i, mirror_j
    integer i, j, k, p, n

    if (allocated(this%p   )) deallocate(this%p   )
    if (allocated(this%iA_p)) deallocate(this%iA_p)
    allocate(this%p   (this%nc,this%npt)); this%p    = 0
    allocate(this%iA_p(this%nc,this%npt)); this%iA_p = 0

    allocate(kmap(this%nc)); kmap = 0

    if (this%mask(1,1) == 0) then
      mirror_i = .true. ; mirror_j = .true.
    else if (this%mask(this%sw,1) == 0) then
      mirror_i = .false.; mirror_j = .true.
    else if (this%mask(1,this%sw) == 0) then
      mirror_i = .true. ; mirror_j = .false.
    else
      mirror_i = .false.; mirror_j = .false.
    end if

    ! Set the p for each evaluation point.
    ! Select monomials according to mask.
    do p = 1, this%npt
      k = 1; n = 1
      do j = 1, this%sw
        do i = 1, this%sw
          if (this%mask(merge(this%sw-i+1, i, mirror_i),merge(this%sw-j+1, j, mirror_j)) == 1) then
            call calc_poly_tensor_product_monomial(this%x(p), this%y(p), i - 1, j - 1, this%p(k,p))
            kmap(n) = k; k = k + 1
          end if
          n = n + 1
        end do
      end do
    end do

    ! Calculate inverse of integral coefficient matrix.
    allocate( A(this%nc,this%nc))
    allocate(iA(this%nc,this%nc))
    call calc_poly_tensor_product_integral_coef_matrix(this%sw, this%sw, A, this%mask)
    n = count(this%mask == 1)
    call inverse_matrix(A(1:n,1:n), iA(1:n,1:n), ierr)
    ! Rearrange iA and p.
    do n = this%nc, 1, -1
      if (kmap(n) /= 0) then
        iA(:,n) = iA(:,kmap(n))
      else
        iA(:,n) = 0
      end if
    end do
    do n = this%nc, 1, -1
      if (kmap(n) /= 0) then
        iA(n,:) = iA(kmap(n),:)
      else
        iA(n,:) = 0
      end if
    end do
    do n = this%nc, 1, -1
      if (kmap(n) /= 0) then
        this%p(n,:) = this%p(kmap(n),:)
      else
        this%p(n,:) = 0
      end if
    end do
    if (ierr /= 0) then
      deallocate(A, iA, kmap)
      return
    end if

    this%iA_p = matmul(iA, this%p)
    ! Rearrange iA_p.

    deallocate(A, iA, kmap)

  end subroutine weno_tensor_product_calc_recon_matrix

  subroutine weno_tensor_product_calc_ideal_coefs(this, ierr)

    class(weno_tensor_product_type), intent(inout) :: this
    integer, intent(out) :: ierr

    real(16), allocatable, dimension(:,:) :: A, AtA, iAtA
    integer i, j, k, p, m, n

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
    allocate(this%ic(this%ns,this%npt)); this%ic = 0

    allocate(   A(this%nc,this%ns)); A = 0
    allocate( AtA(this%ns,this%ns))
    allocate(iAtA(this%ns,this%ns))
    do p = 1, this%npt
      do k = 1, this%ns
        do j = 1, this%sub_sw
          do i = 1, this%sub_sw
            if (this%subs(k)%mask(i,j) == 1) then
              m = (j + this%subs(k)%dj - 1) * this%sw     + i + (this%subs(k)%di)
              n = (j                   - 1) * this%sub_sw + i
              A(m,k) = this%subs(k)%iA_p(n,p)
            end if
          end do
        end do
      end do
      AtA = matmul(transpose(A), A)
      call inverse_matrix(AtA, iAtA, ierr)
      if (ierr /= 0) return
      this%ic(:,p) = matmul(matmul(iAtA, transpose(A)), this%iA_p(:,p))
    end do
    do i = 1, this%ns
      write(*, '(F10.5)') this%ic(i,1)
    end do
    if (abs(sum(matmul(A, this%ic(:,1)) - this%iA_p(:,1))) > 1.0e-30) then
      ierr = 3
    end if
    deallocate(A, AtA, iAtA)

  end subroutine weno_tensor_product_calc_ideal_coefs

  subroutine weno_tensor_product_clear(this)

    class(weno_tensor_product_type), intent(inout) :: this

    this%sw  = 0
    this%nd  = 0
    this%nc  = 0
    this%ns  = 0
    this%npt = 0
    this%di0 = 0
    this%dj0 = 0

    if (allocated(this%mask)) deallocate(this%mask)
    if (allocated(this%x   )) deallocate(this%x   )
    if (allocated(this%y   )) deallocate(this%y   )
    if (allocated(this%p   )) deallocate(this%p   )
    if (allocated(this%iA_p)) deallocate(this%iA_p)
    if (allocated(this%ic  )) deallocate(this%ic  )
    if (allocated(this%subs)) deallocate(this%subs)

    this%initialized = .false.

  end subroutine weno_tensor_product_clear

  subroutine weno_tensor_product_final(this)

    type(weno_tensor_product_type), intent(inout) :: this

    call this%clear()

  end subroutine weno_tensor_product_final

end module weno_tensor_product_mod
