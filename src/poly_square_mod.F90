module poly_square_mod

  use poly_utils_mod
  use math_mod

  implicit none

  private

  public poly_square_type

  type poly_square_type
    logical :: initialized = .false.
    integer :: nd     = 0                ! Dimension number
    integer :: sw     = 0                ! Stencil width
    integer :: nc     = 0                ! Number of cells
    integer :: ns     = 0                ! Number of sub-stencils
    integer :: npt    = 0                ! Number of evaluation points
    integer :: is     = 0                ! Start index of subarray
    integer :: ie     = 0                ! End index of subarray
    integer :: js     = 0                ! Start index of subarray
    integer :: je     = 0                ! End index of subarray
    real(8), allocatable :: xc     (:)   ! X coordinate of cell centroids
    real(8), allocatable :: yc     (:)   ! Y coordinate of cell centroids
    real(8), allocatable :: x      (:)   ! X coordinate of evaluation point
    real(8), allocatable :: y      (:)   ! Y coordinate of evaluation point
    real(8), allocatable :: poly   (:,:) ! Polynomial terms on each evaluation point
    real(8), allocatable :: poly_iA(:,:) ! poly * iA^T
  contains
    procedure :: init                  => poly_square_init
    procedure :: add_point             => poly_square_add_point
    procedure :: calc_recon_matrix     => poly_square_calc_recon_matrix
    procedure :: release_unused_memory => poly_square_release_unused_memory
    procedure :: clear                 => poly_square_clear
    procedure :: poly_square_reconstruct_1d
    procedure :: poly_square_reconstruct_2d
    generic :: reconstruct => poly_square_reconstruct_1d, &
                              poly_square_reconstruct_2d
    final :: poly_square_final
  end type poly_square_type

contains

  subroutine poly_square_init(this, nd, sw, xc, yc, is, ie, js, je)

    class(poly_square_type), intent(inout) :: this
    integer, intent(in) :: nd
    integer, intent(in) :: sw
    real(8), intent(in), optional :: xc(sw)
    real(8), intent(in), optional :: yc(sw)
    integer, intent(in), optional :: is
    integer, intent(in), optional :: ie
    integer, intent(in), optional :: js
    integer, intent(in), optional :: je

    integer i, j, k

    call this%clear()

    this%sw = sw
    this%nd = nd
    this%nc = sw**nd

    if (present(is) .and. present(ie) .and. present(js) .and. present(je)) then
      this%is = is; this%ie = ie; this%js = js; this%je = je
    else if (present(is) .and. present(ie)) then
      this%is = is; this%ie = ie; this%js =  1; this%je = sw**(nd-1)
    else
      this%is =  1; this%ie = sw; this%js =  1; this%je = sw**(nd-1)
    end if

    allocate(this%xc(this%is:this%ie))
    allocate(this%yc(this%js:this%je))
    if (present(xc) .and. present(yc)) then
      this%xc = xc
      this%yc = yc
    else
      ! Set coordinates of cells on the large stencil with origin at center.
      do i = this%is, this%ie
        this%xc(i) = -int(sw / 2) + i - 1
      end do
      do j = this%js, this%je
        this%yc(j) = -int(sw / 2) + j - 1
      end do
    end if

    this%initialized = .true.

  end subroutine poly_square_init

  subroutine poly_square_add_point(this, x, y)

    class(poly_square_type), intent(inout) :: this
    real(8), intent(in) :: x
    real(8), intent(in), optional :: y

    real(8), allocatable :: tmp(:)
    integer i, k

    this%npt = this%npt + 1

    allocate(tmp(this%npt))
    do i = 1, this%npt - 1
      tmp(i) = this%x(i)
    end do
    tmp(this%npt) = x
    if (allocated(this%x)) deallocate(this%x)
    allocate(this%x(this%npt))
    this%x = tmp
    if (present(y)) then ! 1D case does not has y.
      do i = 1, this%npt - 1
        tmp(i) = this%y(i)
      end do
      tmp(this%npt) = y
      if (allocated(this%y)) deallocate(this%y)
      allocate(this%y(this%npt))
      this%y = tmp
    end if
    deallocate(tmp)

  end subroutine poly_square_add_point

  subroutine poly_square_calc_recon_matrix(this, ierr)

    class(poly_square_type), intent(inout) :: this
    integer, intent(out) :: ierr

    ! Local double double arrays for preserving precision.
    real(16), allocatable, dimension(:,:) :: A, iA
    integer i, j, k, ipt, n

    ierr = 0

    if (allocated(this%poly   )) deallocate(this%poly   )
    if (allocated(this%poly_iA)) deallocate(this%poly_iA)

    allocate(this%poly   (this%nc ,this%npt))
    allocate(this%poly_iA(this%npt,this%nc ))

    allocate( A(this%nc,this%nc))
    allocate(iA(this%nc,this%nc))

    ! Set the poly for each evaluation point.
    ! Select monomials according to mask.
    select case (this%nd)
    case (1)
      do ipt = 1, this%npt
        do k = 1, this%sw
          call calc_monomial(this%x(ipt), k - 1, this%poly(k,ipt))
        end do
      end do
    case (2)
      do ipt = 1, this%npt
        k = 1
        do j = 1, this%sw
          do i = 1, this%sw
            call calc_monomial(this%x(ipt), this%y(ipt), i - 1, j - 1, this%poly(k,ipt))
            k = k + 1
          end do
        end do
      end do
    end select

    ! Calculate inverse of integral coefficient matrix.
    select case (this%nd)
    case (1)
      call calc_poly_tensor_product_integral_matrix(this%sw, this%xc, A)
    case (2)
      call calc_poly_tensor_product_integral_matrix(this%sw, this%sw, this%xc, this%yc, A)
    end select
    call inverse_matrix(A, iA, ierr)
    if (ierr /= 0) then
      deallocate(A, iA)
      return
    end if

    this%poly_iA = transpose(matmul(iA, this%poly))

    deallocate(A, iA)

  end subroutine poly_square_calc_recon_matrix

  subroutine poly_square_reconstruct_1d(this, fi, fo, ierr)

    class(poly_square_type), intent(inout) :: this
    real(8), intent(in ) :: fi(:)   ! Cell averaged function values
    real(8), intent(out) :: fo(:)   ! Reconstructed function values on evaluation points
    integer, intent(out) :: ierr

    ierr = 0

#ifndef NDEBUG
    if (size(fi) /= this%nc) then
      ierr = 1
      return
    end if
    if (size(fo) /= this%npt) then
      ierr = 1
      return
    end if
#endif

    fo = matmul(this%poly_iA, fi)

  end subroutine poly_square_reconstruct_1d

  subroutine poly_square_reconstruct_2d(this, fi, fo, ierr)

    class(poly_square_type), intent(inout) :: this
    real(8), intent(in ) :: fi(:,:) ! Cell averaged function values
    real(8), intent(out) :: fo(:)   ! Reconstructed function values on evaluation points
    integer, intent(out) :: ierr

    ierr = 0

#ifndef NDEBUG
    if (size(fi) /= this%nc) then
      ierr = 1
      return
    end if
    if (size(fo) /= this%npt) then
      ierr = 1
      return
    end if
#endif

    fo = matmul(this%poly_iA, pack(fi, .true.))

  end subroutine poly_square_reconstruct_2d

  subroutine poly_square_release_unused_memory(this)

    class(poly_square_type), intent(inout) :: this

    type(poly_square_type), allocatable :: subs(:)
    real(8), allocatable :: ic(:,:)
    integer i, k, ns

    if (allocated(this%xc  )) deallocate(this%xc  )
    if (allocated(this%yc  )) deallocate(this%yc  )
    if (allocated(this%x   )) deallocate(this%x   )
    if (allocated(this%y   )) deallocate(this%y   )
    if (allocated(this%poly)) deallocate(this%poly)

  end subroutine poly_square_release_unused_memory

  subroutine poly_square_clear(this)

    class(poly_square_type), intent(inout) :: this

    this%nd     = 0
    this%sw     = 0
    this%nc     = 0
    this%ns     = 0
    this%npt    = 0
    this%is     = 0
    this%ie     = 0
    this%js     = 0
    this%je     = 0

    if (allocated(this%xc     )) deallocate(this%xc     )
    if (allocated(this%yc     )) deallocate(this%yc     )
    if (allocated(this%x      )) deallocate(this%x      )
    if (allocated(this%y      )) deallocate(this%y      )
    if (allocated(this%poly   )) deallocate(this%poly   )
    if (allocated(this%poly_iA)) deallocate(this%poly_iA)

    this%initialized = .false.

  end subroutine poly_square_clear

  subroutine poly_square_final(this)

    type(poly_square_type), intent(inout) :: this

    call this%clear()

  end subroutine poly_square_final

end module poly_square_mod
