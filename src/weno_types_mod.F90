module weno_types_mod

  implicit none

  private

  public weno_tensor_product_type

  type weno_tensor_product_type
    logical :: initialized = .false.
    integer :: nd  = 0                      ! Dimension number
    integer :: sw  = 0                      ! Stencil width
    integer :: npt = 0                      ! Number of evaluation points
    integer, allocatable :: mask(:,:)       ! Mask unavailable nodes by 0
    real(8), allocatable :: iA(:,:)         ! Polynomial integral coefficient matrix inverse
    real(8), allocatable :: xp(:)           ! X coordinate of evaluation point
    real(8), allocatable :: yp(:)           ! Y coordinate of evaluation point
    real(8), allocatable :: coefs(:,:)      ! Reconstruction coefficients for each point (only on sub-stencils)
    type(weno_tensor_product_type), allocatable :: subs(:) ! Sub-stencils
  contains
    procedure :: init      => weno_tensor_product_init
    procedure :: add_point => weno_tensor_product_add_point
    procedure :: clear     => weno_tensor_product_clear
    final :: weno_tensor_product_final
  end type weno_tensor_product_type

contains

  subroutine weno_tensor_product_init(this, nd, sw, no_subs)

    class(weno_tensor_product_type), intent(inout) :: this
    integer, intent(in) :: nd
    integer, intent(in) :: sw
    logical, intent(in), optional :: no_subs

    integer i

    call this%clear()

    this%sw = sw
    this%nd = nd

    select case (nd)
    case (1)
      allocate(this%mask(sw,1 ))
    case (2)
      allocate(this%mask(sw,sw))
    end select
    allocate(this%iA  (sw**nd,sw**nd))

    ! Initialize sub-stencils.
    if (merge(.not. no_subs, .true., present(no_subs))) then
      select case (nd)
      case (1)
        select case (sw)
        case (5)
          allocate(this%subs(3))
        end select
      case (2)
        select case (sw)
        case (5)
          allocate(this%subs(9))
          do i = 1, 9
            call this%subs(i)%init(nd, 3, no_subs=.true.)
          end do
        end select
      end select
    end if

    this%initialized = .true.

  end subroutine weno_tensor_product_init

  subroutine weno_tensor_product_add_point(this, x, y)

    class(weno_tensor_product_type), intent(inout) :: this
    real(8), intent(in) :: x
    real(8), intent(in) :: y

    real(8), allocatable :: tmp(:)
    integer i

    this%npt = this%npt + 1

    allocate(tmp(this%npt))
    do i = 1, this%npt - 1
      tmp(i) = this%xp(i)
    end do
    tmp(this%npt) = x
    if (allocated(this%xp)) deallocate(this%xp)
    allocate(this%xp(this%npt))
    this%xp = tmp
    do i = 1, this%npt - 1
      tmp(i) = this%yp(i)
    end do
    tmp(this%npt) = y
    if (allocated(this%yp)) deallocate(this%yp)
    allocate(this%yp(this%npt))
    this%yp = tmp

    if (allocated(this%coefs)) deallocate(this%coefs)
    allocate(this%coefs(this%sw**this%nd,this%npt))

  end subroutine weno_tensor_product_add_point

  subroutine weno_tensor_product_clear(this)

    class(weno_tensor_product_type), intent(inout) :: this

    this%sw  = 0
    this%nd  = 0
    this%npt = 0

    if (allocated(this%mask )) deallocate(this%mask )
    if (allocated(this%iA   )) deallocate(this%iA   )
    if (allocated(this%xp   )) deallocate(this%xp   )
    if (allocated(this%yp   )) deallocate(this%yp   )
    if (allocated(this%coefs)) deallocate(this%coefs)
    if (allocated(this%subs )) deallocate(this%subs )

    this%initialized = .false.

  end subroutine weno_tensor_product_clear

  subroutine weno_tensor_product_final(this)

    type(weno_tensor_product_type), intent(inout) :: this

    call this%clear()

  end subroutine weno_tensor_product_final

end module weno_types_mod
