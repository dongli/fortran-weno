module weno_rhombus_mod

  use poly_utils_mod
  use math_mod

  implicit none

  private

  public weno_rhombus_type

  ! 1: i,j,x^?,y^?
  ! 2: cell
  ! 3: sub-stencil
  !                            1 2
  integer, target :: stencil_3(4,5) = [ &
               [2,3,0,2]           , & !       y^2
    [1,2,1,0], [2,2,0,0], [3,2,2,0], & !   x    1    y
               [2,1,0,1]             & !       x^2
  ]
  !                               1 2 3
  integer, target :: substencil_3(4,3,4) = [ &
    [[2,2,0,0], [1,2,1,0], [2,1,0,1]], & ! sub-stencil 1
    [[2,2,0,0], [3,2,1,0], [2,1,0,1]], & ! sub-stencil 2
    [[2,2,0,0], [3,2,1,0], [2,3,0,1]], & ! sub-stencil 3
    [[2,2,0,0], [1,2,1,0], [2,3,0,1]]  & ! sub-stencil 4
  ]
  !                            1  2
  integer, target :: stencil_5(4,13) = [ &
                          [3,5,0,4]                      , & !                     y^4
               [2,4,1,2], [3,4,0,3], [4,4,2,2]           , & !          x y^2      y^3      x^2 y^2
    [1,3,2,0], [2,3,1,0], [3,3,0,0], [4,3,3,0], [5,3,4,0], & !  x^2     x           1       x^3          x^4
               [2,2,1,1], [3,2,0,1], [4,2,2,1]           , & !          x y         y       x^2 y
                          [3,1,0,2]                        & !                     y^2
  ]
  !                               1 2 3
  integer, target :: substencil_5(4,5,9) = [ &
    [[3,3,0,0], [2,3,1,0], [4,3,0,1], [3,2,2,0], [3,4,0,2]], &
    [[2,3,0,0], [1,3,1,0], [3,3,0,1], [2,2,2,0], [2,4,0,2]], &
    [[3,2,0,0], [2,2,1,0], [4,2,0,1], [3,3,2,0], [3,1,0,2]], &
    [[4,3,0,0], [3,3,1,0], [5,3,0,1], [4,2,2,0], [4,4,0,2]], &
    [[3,4,0,0], [2,4,1,0], [4,4,0,1], [3,3,2,0], [3,5,0,2]], &
    [[3,5,0,0], [3,4,1,0], [3,3,0,1], [2,3,2,0], [4,3,0,2]], &
    [[1,3,0,0], [2,3,1,0], [3,3,0,1], [3,2,2,0], [3,4,0,2]], &
    [[3,1,0,0], [3,2,1,0], [3,3,0,1], [2,3,2,0], [4,3,0,2]], &
    [[5,3,0,0], [4,3,1,0], [3,3,0,1], [3,4,2,0], [3,2,0,2]]  &
  ]

  integer, target :: stencil_7(4,25) = [ &
                                     [4,7,0,6]                                 , &
                          [3,6,1,4], [4,6,0,4], [5,6,2,4]                      , &
               [2,5,3,2], [3,5,1,2], [4,5,0,2], [5,5,2,2], [6,5,4,2]           , &
    [1,4,5,0], [2,4,3,0], [3,4,1,0], [4,4,0,0], [5,4,2,0], [6,4,4,0], [7,4,6,0], &
               [2,3,3,1], [3,3,1,1], [4,3,0,1], [5,3,2,1], [6,3,4,1]           , &
                          [3,2,1,3], [4,2,0,3], [5,2,2,3]                      , &
                                     [4,1,0,5]                                   &
  ]
  !integer, target :: substencil_7(4,9,25) = [ &
  !  [[1,4,0,0], [2,4,1,0], [3,4,2,0], [4,4,3,0], [5,4,4,0], [3,6,0,1], [3,5,0,2], [3,3,0,3], [3,2,0,4]], &
  !  [[2,3,0,0], [3,3,1,0], [4,3,2,0], [5,3,3,0], [6,3,4,0], [4,5,0,1], [4,4,0,2], [4,2,0,3], [4,1,0,4]], &
  !  [[3,4,0,0], [4,4,1,0], [5,4,2,0], [6,4,3,0], [7,4,4,0], [5,6,0,1], [5,5,0,2], [5,3,0,3], [5,2,0,4]], &
  !  [[2,5,0,0], [3,5,1,0], [4,5,2,0], [5,5,3,0], [6,5,4,0], [4,7,0,1], [4,6,0,2], [4,4,0,3], [4,3,0,4]], &
  !  [[2,4,0,0], [3,4,1,0], [4,4,2,0], [5,4,3,0], [6,4,4,0], [3,6,0,1], [3,5,0,2], [3,3,0,3], [3,2,0,4]], &
  !  [[2,3,0,0], [3,3,1,0], [4,3,2,0], [5,3,3,0], [6,3,4,0], [4,6,0,1], [4,5,0,2], [4,4,0,3], [4,2,0,4]], &
  !  [[2,4,0,0], [3,4,1,0], [4,4,2,0], [5,4,3,0], [6,4,4,0], [5,6,0,1], [5,5,0,2], [5,3,0,3], [5,2,0,4]], &
  !  [[2,5,0,0], [3,5,1,0], [4,5,2,0], [5,5,3,0], [6,5,4,0], [4,6,0,1], [4,4,0,2], [4,3,0,3], [4,2,0,4]], &
  !  [[3,4,0,0], [4,4,1,0], [5,4,2,0], [6,4,3,0], [7,4,4,0], [3,6,0,1], [3,5,0,2], [3,3,0,3], [3,2,0,4]], &
  !  [[2,3,0,0], [3,3,1,0], [4,3,2,0], [5,3,3,0], [6,3,4,0], [4,7,0,1], [4,6,0,2], [4,5,0,3], [4,4,0,4]], &
  !  [[1,4,0,0], [2,4,1,0], [3,4,2,0], [4,4,3,0], [5,4,4,0], [5,6,0,1], [5,5,0,2], [5,3,0,3], [5,2,0,4]], &
  !  [[2,5,0,0], [3,5,1,0], [4,5,2,0], [5,5,3,0], [6,5,4,0], [4,4,0,1], [4,3,0,2], [4,2,0,3], [4,1,0,4]], &
  !  [[2,4,0,0], [3,4,1,0], [4,4,2,0], [5,4,3,0], [6,4,4,0], [4,6,0,1], [4,5,0,2], [4,3,0,3], [4,2,0,4]], &
  !  [[2,5,0,0], [3,5,1,0], [4,6,2,0], [4,5,3,0], [4,4,4,0], [4,3,0,1], [4,2,0,2], [5,3,0,3], [6,3,0,4]], &
  !  [[2,3,0,0], [3,3,1,0], [4,6,2,0], [4,5,3,0], [4,4,4,0], [4,3,0,1], [4,2,0,2], [5,5,0,3], [6,5,0,4]], &
  !  [[2,4,0,0], [3,4,1,0], [4,4,2,0], [5,4,3,0], [6,4,4,0], [3,2,0,1], [3,3,0,2], [5,6,0,3], [5,5,0,4]], &
  !  [[2,4,0,0], [3,4,1,0], [4,4,2,0], [5,4,3,0], [6,4,4,0], [3,6,0,1], [3,5,0,2], [5,3,0,3], [5,2,0,4]], &
  !  [[1,4,0,0], [2,4,1,0], [3,4,2,0], [4,4,3,0], [5,4,4,0], [4,1,0,1], [4,2,0,2], [4,3,0,3], [4,5,0,4]], &
  !  [[3,4,0,0], [4,4,1,0], [5,4,2,0], [6,4,3,0], [7,4,4,0], [4,1,0,1], [4,2,0,2], [4,3,0,3], [4,5,0,4]], &
  !  [[3,4,0,0], [4,4,1,0], [5,4,2,0], [6,4,3,0], [7,4,4,0], [4,3,0,1], [4,5,0,2], [4,6,0,3], [4,7,0,4]], &
  !  [[1,4,0,0], [2,4,1,0], [3,4,2,0], [4,4,3,0], [5,4,4,0], [4,3,0,1], [4,5,0,2], [4,6,0,3], [4,7,0,4]], &
  !  [[1,4,0,0], [2,4,1,0], [3,4,2,0], [4,4,3,0], [5,4,4,0], [3,5,0,1], [3,6,0,2], [5,2,0,3], [5,3,0,4]], &
  !  [[4,1,0,0], [4,2,1,0], [4,3,2,0], [4,4,3,0], [4,5,4,0], [2,3,0,1], [3,3,0,2], [5,5,0,3], [6,5,0,4]], &
  !  [[3,4,0,0], [4,4,1,0], [5,4,2,0], [6,4,3,0], [7,4,4,0], [3,5,0,1], [3,6,0,2], [5,2,0,3], [5,3,0,4]], &
  !  [[4,3,0,0], [4,4,1,0], [4,5,2,0], [4,6,3,0], [4,7,4,0], [2,5,0,1], [3,5,0,2], [5,3,0,3], [6,3,0,4]]  &
  !]

  integer, target :: substencil_7(4,7,20) = [ &
    [[1,4,0,0], [2,4,1,0], [3,4,2,0], [4,4,3,0], [4,1,0,1], [4,2,0,2], [4,3,0,3]], &
    [[4,4,0,0], [5,4,1,0], [6,4,2,0], [7,4,3,0], [4,1,0,1], [4,2,0,2], [4,3,0,3]], &
    [[4,4,0,0], [5,4,1,0], [6,4,2,0], [7,4,3,0], [4,5,0,1], [4,6,0,2], [4,7,0,3]], &
    [[1,4,0,0], [2,4,1,0], [3,4,2,0], [4,4,3,0], [4,5,0,1], [4,6,0,2], [4,7,0,3]], &
    
    [[3,3,0,0], [3,4,1,0], [3,5,2,0], [4,3,0,1], [4,4,0,2], [4,5,1,1], [5,4,1,2]], &
    [[3,5,0,0], [4,5,1,0], [5,5,2,0], [3,4,0,1], [4,4,0,2], [5,4,1,1], [4,3,2,1]], &
    [[5,3,0,0], [5,4,1,0], [5,5,2,0], [4,3,0,1], [4,4,0,2], [4,5,1,1], [3,4,1,2]], &
    [[3,3,0,0], [4,3,1,0], [5,3,2,0], [3,4,0,1], [4,4,0,2], [5,4,1,1], [4,5,2,1]], &
    
    [[1,4,0,0], [2,4,1,0], [3,4,2,0], [4,4,3,0], [5,4,4,0], [2,3,0,1], [2,5,0,2]], &
    [[4,1,0,0], [4,2,1,0], [4,3,2,0], [4,4,0,1], [4,5,0,2], [3,2,0,3], [5,2,0,4]], &
    [[3,4,0,0], [4,4,1,0], [5,4,2,0], [6,4,3,0], [7,4,4,0], [6,3,0,1], [6,5,0,2]], &
    [[4,3,0,0], [4,4,1,0], [4,5,2,0], [4,6,0,1], [4,7,0,2], [3,6,0,3], [5,6,0,4]], &
    
    [[3,2,0,0], [3,3,1,0], [3,4,2,0], [3,5,0,1], [3,6,0,2], [4,4,0,3], [5,4,0,4]], &
    [[2,3,0,0], [3,3,1,0], [4,3,2,0], [5,3,3,0], [6,3,4,0], [4,4,0,1], [4,5,0,2]], &
    [[3,4,0,0], [4,4,1,0], [5,4,2,0], [5,2,0,1], [5,3,0,2], [5,5,0,3], [5,6,0,4]], &
    [[2,5,0,0], [3,5,1,0], [4,5,2,0], [5,5,3,0], [6,5,4,0], [4,3,0,1], [4,4,0,2]], &
    
    [[3,4,0,0], [4,4,1,0], [5,4,2,0], [4,5,0,1], [3,6,0,2], [4,6,1,1], [5,6,2,1]], &
    [[3,4,0,0], [4,4,1,0], [5,4,2,0], [4,3,0,1], [3,2,0,2], [4,2,1,1], [5,2,2,1]], &
    [[2,3,0,0], [2,4,1,0], [2,5,2,0], [3,4,0,1], [4,3,0,2], [4,4,1,1], [4,5,1,2]], &
    [[6,3,0,0], [6,4,1,0], [6,5,2,0], [5,4,0,1], [4,3,0,2], [4,4,1,1], [4,5,1,2]]  &
  ]
    
  type weno_rhombus_type
    logical :: initialized = .false.
    integer :: id     = 0                  ! Sub-stencil ID
    integer :: sw     = 0                  ! Stencil width
    integer :: swx    = 0
    integer :: swy    = 0
    integer :: nc     = 0                  ! Number of cells
    integer :: ns     = 0                  ! Number of sub-stencils
    integer :: npt    = 0                  ! Number of evaluation points
    integer :: is     = 0                  ! Start index of subarray
    integer :: ie     = 0                  ! End index of subarray
    integer :: js     = 0                  ! Start index of subarray
    integer :: je     = 0                  ! End index of subarray
    integer, allocatable :: cell_mask(:,:)
    integer, allocatable :: ij_to_1d (:,:) ! Map i,j index to 1D cell index
    real(8), allocatable :: xc       (:)   ! X coordinate of cell centroids
    real(8), allocatable :: yc       (:)   ! Y coordinate of cell centroids
    real(8), allocatable :: x        (:)   ! X coordinate of evaluation point
    real(8), allocatable :: y        (:)   ! Y coordinate of evaluation point
    real(8), allocatable :: poly     (:,:) ! Polynomial terms on each evaluation point
    real(8), allocatable :: iA       (:,:) ! A * a = f, iA = inverse(A)
    real(8), allocatable :: poly_iA (:,:) ! poly * iA
    real(8), allocatable :: gamma    (:,:) ! Ideal coefficients for combining sub-stencils (only on stencil)
    integer, pointer :: ijxy(:,:)
    type(weno_rhombus_type), allocatable :: subs(:)
  contains
    procedure :: init              => weno_rhombus_init
    procedure :: add_point         => weno_rhombus_add_point
    procedure :: calc_recon_matrix => weno_rhombus_calc_recon_matrix
    procedure :: calc_ideal_coefs  => weno_rhombus_calc_ideal_coefs
    procedure :: clear             => weno_rhombus_clear
    final :: weno_rhombus_final
  end type weno_rhombus_type

contains

  recursive subroutine weno_rhombus_init(this, sw, swx, swy, xc, yc, is, ie, js, je, id)

    class(weno_rhombus_type), intent(inout) :: this
    integer, intent(in), optional :: sw
    integer, intent(in), optional :: swx
    integer, intent(in), optional :: swy
    real(8), intent(in), optional :: xc(:)
    real(8), intent(in), optional :: yc(:)
    integer, intent(in), optional :: is
    integer, intent(in), optional :: ie
    integer, intent(in), optional :: js
    integer, intent(in), optional :: je
    integer, intent(in), optional :: id

    real(8), allocatable :: x(:), y(:)
    integer i, j, k, sub_swx, sub_swy, sub_is, sub_ie, sub_js, sub_je

    call this%clear()

    if (present(sw)) then
      this%swx = sw ; this%swy = sw
    else if (present(swx) .and. present(swy)) then
      this%swx = swx; this%swy = swy
    end if

    if (present(id)) this%id = id

    if (present(is) .and. present(ie) .and. present(js) .and. present(je)) then
      this%is = is; this%ie =       ie; this%js = js; this%je =       je
    else if (present(is) .and. present(ie)) then
      this%is = is; this%ie =       ie; this%js =  1; this%je = this%swy
    else
      this%is =  1; this%ie = this%swx; this%js =  1; this%je = this%swy
    end if

    allocate(this%cell_mask(this%is:this%ie,this%js:this%je)); this%cell_mask = 0

    if (this%id == 0) then
      select case (sw)
      case (3)
        this%ns = size(substencil_3, 3)
        this%ijxy => stencil_3
      case (5)
        this%ns = size(substencil_5, 3)
        this%ijxy => stencil_5
      case (7)
        this%ns = size(substencil_7, 3)
        this%ijxy => stencil_7
      end select
    end if
    this%nc = size(this%ijxy, 2)

    ! Set cell coordinates.
    allocate(this%xc(this%nc))
    allocate(this%yc(this%nc))
    allocate(x(this%is:this%ie))
    allocate(y(this%js:this%je))
    if (present(xc) .and. present(yc)) then
      x = xc
      y = yc
    else
      ! Set coordinates of cells on the large stencil with origin at center.
      do i = this%is, this%ie
        x(i) = -int(this%swx / 2) + i - 1
      end do
      do j = this%js, this%je
        y(j) = -int(this%swy / 2) + j - 1
      end do
    end if
    do k = 1, this%nc
      this%xc(k) = x(this%ijxy(1,k))
      this%yc(k) = y(this%ijxy(2,k))
    end do

    ! Initialize sub-stencils.
    if (this%id == 0) then
      allocate(this%subs(this%ns))
      do k = 1, this%ns
        select case (sw)
        case (3)
          sub_is  = minval(substencil_3(1,:,k)); sub_ie  = maxval(substencil_3(1,:,k))
          sub_js  = minval(substencil_3(2,:,k)); sub_je  = maxval(substencil_3(2,:,k))
          sub_swx = sub_ie - sub_is + 1; sub_swy = sub_je - sub_js + 1
          this%subs(k)%ijxy => substencil_3(:,:,k)
        case (5)
          sub_is  = minval(substencil_5(1,:,k)); sub_ie  = maxval(substencil_5(1,:,k))
          sub_js  = minval(substencil_5(2,:,k)); sub_je  = maxval(substencil_5(2,:,k))
          sub_swx = sub_ie - sub_is + 1; sub_swy = sub_je - sub_js + 1
          this%subs(k)%ijxy => substencil_5(:,:,k)
        case (7)
          sub_is  = minval(substencil_7(1,:,k)); sub_ie  = maxval(substencil_7(1,:,k))
          sub_js  = minval(substencil_7(2,:,k)); sub_je  = maxval(substencil_7(2,:,k))
          sub_swx = sub_ie - sub_is + 1; sub_swy = sub_je - sub_js + 1
          this%subs(k)%ijxy => substencil_7(:,:,k)
        end select
        call this%subs(k)%init(swx=sub_swx, swy=sub_swy, xc=x(sub_is:sub_ie), yc=y(sub_js:sub_je), &
                               is=sub_is, ie=sub_ie, js=sub_js, je=sub_je, id=k)
      end do
      ! Set index map.
      allocate(this%ij_to_1d(sw,sw)); this%ij_to_1d = 0
      do k = 1, this%nc
        this%ij_to_1d(this%ijxy(1,k),this%ijxy(2,k)) = k
      end do
    end if

    this%initialized = .true.

  end subroutine weno_rhombus_init

  recursive subroutine weno_rhombus_add_point(this, x, y)

    class(weno_rhombus_type), intent(inout) :: this
    real(8), intent(in) :: x
    real(8), intent(in) :: y

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
    do i = 1, this%npt - 1
      tmp(i) = this%y(i)
    end do
    tmp(this%npt) = y
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

  end subroutine weno_rhombus_add_point

  subroutine weno_rhombus_calc_recon_matrix(this, ierr)

    class(weno_rhombus_type), intent(inout) :: this
    integer, intent(out) :: ierr

    real(16), allocatable, dimension(:,:) :: A, iA
    integer ic, ipt

    ierr = 0

    if (allocated(this%poly   )) deallocate(this%poly   )
    if (allocated(this%iA     )) deallocate(this%iA     )
    if (allocated(this%poly_iA)) deallocate(this%poly_iA)

    allocate(this%poly   (this%npt,this%nc))
    allocate(this%iA     (this%nc ,this%nc))
    allocate(this%poly_iA(this%npt,this%nc))

    allocate( A(this%nc,this%nc))
    allocate(iA(this%nc,this%nc))

    do ipt = 1, this%npt
      do ic = 1, this%nc
        call calc_monomial(this%x(ipt), this%y(ipt), this%ijxy(3,ic), this%ijxy(4,ic), this%poly(ipt,ic))
      end do
    end do

    call calc_poly_integral_matrix(this%nc, this%nc, this%ijxy(3:4,:), this%xc, this%yc, A)
    A = transpose(A) ! nc x np
    call inverse_matrix(A, iA, ierr)
    if (ierr /= 0) then
      deallocate(A, iA)
      return
    end if

    this%poly_iA = matmul(this%poly, iA)

    deallocate(A, iA)

  end subroutine weno_rhombus_calc_recon_matrix

  subroutine weno_rhombus_calc_ideal_coefs(this, ierr)

    class(weno_rhombus_type), intent(inout) :: this
    integer, intent(out) :: ierr

    ! Local double double arrays for preserving precision.
    real(16), allocatable, dimension(:,:) :: A, AtA, iAtA
    integer ic, k, ipt, i, j

    real(16), parameter :: eps = 1.0e-17

    ierr = 0

    if (.not. this%initialized) then
      ierr = 1
      return
    end if

    do k = 1, this%ns
      call this%subs(k)%calc_recon_matrix(ierr)
      if (ierr /= 0) return
    end do
    call this%calc_recon_matrix(ierr)
    if (ierr /= 0) return

    if (allocated(this%gamma)) deallocate(this%gamma)
    allocate(this%gamma(this%ns,this%npt))

    allocate(    A(this%nc,this%ns )); A = 0
    allocate(  AtA(this%ns,this%ns ))
    allocate( iAtA(this%ns,this%ns )); iAtA = 0
    do ipt = 1, this%npt
      substencil: do k = 1, this%ns
        monomial: do ic = 1, this%subs(k)%nc
          i = this%subs(k)%ijxy(1,ic)
          j = this%subs(k)%ijxy(2,ic)
          A(this%ij_to_1d(i,j),k) = this%subs(k)%poly_iA(ipt,ic)
        end do monomial
      end do substencil

      !print *, 'A'
      !do i = 1, this%nc
      !  do j = 1, this%ns
      !    write(*, '(F20.15)', advance='no') A(i,j)
      !  end do
      !  write(*, *)
      !end do

      AtA = matmul(transpose(A), A)
      do i = 1, this%ns
        AtA(i,i) = AtA(i,i) + eps
      end do

      !print *, 'AtA'
      !do i = 1, this%ns
      !  do j = 1, this%ns
      !    write(*, '(F20.15)', advance='no') AtA(i,j)
      !  end do
      !  write(*, *)
      !end do

      call inverse_matrix(AtA, iAtA, ierr)
      if (ierr /= 0) return

      !print *, 'iAtA'
      !do i = 1, this%ns
      !  do j = 1, this%ns
      !    write(*, '(F20.1)', advance='no') iAtA(i,j)
      !  end do
      !  write(*, *)
      !end do
      this%gamma(:,ipt) = matmul(matmul(iAtA, transpose(A)), this%poly_iA(ipt,:))
      ! Set near-zero values to zero exactly.
      do k = 1, this%ns
        if (abs(this%gamma(k,ipt)) < 1.0e-15) this%gamma(k,ipt) = 0
      end do
    end do
    if (abs(sum(matmul(A, this%gamma(:,1)) - this%poly_iA(1,:))) > 1.0e-15) then
      print *, abs(sum(matmul(A, this%gamma(:,1)) - this%poly_iA(1,:)))
      print *, this%gamma(:,1)
      ierr = 3
      this%gamma = 0
    end if
    deallocate(A, AtA, iAtA)

  end subroutine weno_rhombus_calc_ideal_coefs

  subroutine weno_rhombus_clear(this)

    class(weno_rhombus_type), intent(inout) :: this

    if (allocated(this%cell_mask)) deallocate(this%cell_mask)
    if (allocated(this%ij_to_1d )) deallocate(this%ij_to_1d )
    if (allocated(this%xc       )) deallocate(this%xc       )
    if (allocated(this%yc       )) deallocate(this%yc       )
    if (allocated(this%x        )) deallocate(this%x        )
    if (allocated(this%y        )) deallocate(this%y        )
    if (allocated(this%poly     )) deallocate(this%poly     )
    if (allocated(this%iA       )) deallocate(this%iA       )
    if (allocated(this%poly_iA  )) deallocate(this%poly_iA  )
    if (allocated(this%gamma    )) deallocate(this%gamma    )
    if (allocated(this%subs     )) deallocate(this%subs     )

    this%initialized = .false.

  end subroutine weno_rhombus_clear

  subroutine weno_rhombus_final(this)

    type(weno_rhombus_type), intent(inout) :: this

    call this%clear()

  end subroutine weno_rhombus_final

end module weno_rhombus_mod
