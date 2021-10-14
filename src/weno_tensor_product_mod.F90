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
    integer :: is     = 0                  ! Start index of subarray
    integer :: ie     = 0                  ! End index of subarray
    integer :: js     = 0                  ! Start index of subarray
    integer :: je     = 0                  ! End index of subarray
    real(8) :: si     = 0                  ! Smoothness indicator
    integer, allocatable :: cell_mask(:,:) ! Mask unavailable cells by 0
    integer, allocatable :: poly_mask(:,:) ! Mask unavailable polynomial terms by 0
    real(8), allocatable :: xc       (:)   ! X coordinate of cell centroids
    real(8), allocatable :: yc       (:)   ! Y coordinate of cell centroids
    real(8), allocatable :: x        (:)   ! X coordinate of evaluation point
    real(8), allocatable :: y        (:)   ! Y coordinate of evaluation point
    real(8), allocatable :: a        (:)   ! Reconstruction coefficients
    real(8), allocatable :: p        (:,:) ! Polynomial terms on each evaluation point
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
    procedure :: calc_smooth_indicator => weno_tensor_product_calc_smooth_indicator
    procedure :: reconstruct           => weno_tensor_product_reconstruct
    procedure :: release_unused_memory => weno_tensor_product_release_unused_memory
    procedure :: clear                 => weno_tensor_product_clear
    procedure, private :: weno_tensor_product_assign
    generic :: assignment(=) => weno_tensor_product_assign
    final :: weno_tensor_product_final
  end type weno_tensor_product_type

contains

  subroutine weno_tensor_product_init(this, nd, sw, xc, yc, is, ie, js, je, mask, id, has_subs)

    class(weno_tensor_product_type), intent(inout) :: this
    integer, intent(in) :: nd
    integer, intent(in) :: sw
    real(8), intent(in), optional :: xc(sw)
    real(8), intent(in), optional :: yc(sw)
    integer, intent(in), optional :: is
    integer, intent(in), optional :: ie
    integer, intent(in), optional :: js
    integer, intent(in), optional :: je
    integer, intent(in), optional :: mask(:,:) ! Cell mask
    integer, intent(in), optional :: id
    logical, intent(in), optional :: has_subs

    integer i, j, k, sub_ie, sub_je

    call this%clear()

    this%sw = sw
    this%nd = nd
    this%nc = sw**nd

    if (present(id)) this%id = id

    if (present(is) .and. present(ie) .and. present(js) .and. present(je)) then
      this%is = is; this%ie = ie; this%js = js; this%je = je
    else if (present(is) .and. present(ie)) then
      this%is = is; this%ie = ie; this%js =  1; this%je = sw**(nd-1)
    else
      this%is =  1; this%ie = sw; this%js =  1; this%je = sw**(nd-1)
    end if

    allocate(this%cell_mask(this%is:this%ie,this%js:this%je))
    allocate(this%poly_mask(sw,sw**(nd-1)))
    allocate(this%a(this%nc))

    if (present(mask)) then
      this%cell_mask = mask
    else
      this%cell_mask = 1
    end if

    ! Set polynomial mask from cell mask.
    if (this%cell_mask(this%is,this%js) == 0) then
      this%poly_mask = this%cell_mask(this%ie:this%is:-1,this%je:this%js:-1)
    else if (this%cell_mask(this%ie,this%js) == 0) then
      this%poly_mask = this%cell_mask(:,this%je:this%js:-1)
    else if (this%cell_mask(this%is,this%je) == 0) then
      this%poly_mask = this%cell_mask(this%ie:this%is:-1,:)
    else
      this%poly_mask = this%cell_mask
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

    ! Initialize sub-stencils.
    if (merge(has_subs, .true., present(has_subs))) then
      this%sub_sw = int((sw + 1) / 2) ! Set sub-stencil size.
      this%ns = this%sub_sw**nd
      allocate(this%subs(this%ns))
      k = 1
      do j = this%js, this%js + this%sub_sw**(nd - 1) - 1
        do i = this%is, this%is + this%sub_sw - 1
          sub_ie = i + this%sub_sw - 1
          sub_je = j + this%sub_sw**(nd-1) - 1
          call this%subs(k)%init(nd=nd, sw=this%sub_sw, xc=this%xc(i:sub_ie), yc=this%yc(j:sub_je), &
                                 is=i, ie=sub_ie, js=j, je=sub_je, id=k, has_subs=.false., &
                                 mask=this%cell_mask(i:sub_ie,j:sub_je))
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
    real(16), allocatable, dimension(:,:) :: A, iA
    integer , allocatable :: idx_map(:)
    integer i, j, k, p, n

    ierr = 0

    if (allocated(this%p      )) deallocate(this%p      )
    if (allocated(this%iA     )) deallocate(this%iA     )
    if (allocated(this%iAp    )) deallocate(this%iAp    )
    if (allocated(this%iAp_r16)) deallocate(this%iAp_r16)

    allocate(this%p      (this%nc,this%npt))
    allocate(this%iA     (this%nc,this%nc))
    allocate(this%iAp    (this%nc,this%npt))
    allocate(this%iAp_r16(this%nc,this%npt))

    allocate( A     (this%nc,this%nc ))
    allocate(iA     (this%nc,this%nc ))
    allocate(idx_map(this%nc         )); idx_map = 0

    ! Set index maps.
    k = 1; n = 1
    do j = this%js, this%je
      do i = this%is, this%ie
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
            call calc_poly_tensor_product_monomial(this%x(p), i - 1, this%p(k,p))
            k = k + 1
          end if
        end do
      end do
    case (2)
      do p = 1, this%npt
        k = 1
        do j = 1, this%sw
          do i = 1, this%sw
            if (this%poly_mask(i,j) == 1) then
              call calc_poly_tensor_product_monomial(this%x(p), this%y(p), i - 1, j - 1, this%p(k,p))
              k = k + 1
            end if
          end do
        end do
      end do
    end select

    ! Calculate inverse of integral coefficient matrix.
    select case (this%nd)
    case (1)
      call calc_poly_tensor_product_integral_coef_matrix(this%sw, this%xc, A, this%cell_mask(:,1), this%poly_mask(:,1))
    case (2)
      call calc_poly_tensor_product_integral_coef_matrix(this%sw, this%sw, this%xc, this%yc, A, this%cell_mask, this%poly_mask)
    end select
    n = count(this%cell_mask == 1)
    call inverse_matrix(A(1:n,1:n), iA(1:n,1:n), ierr)

    if (ierr /= 0) then
      deallocate(A, iA, idx_map)
      return
    end if

    this%iAp_r16(1:n,:) = matmul(iA(1:n,1:n), this%p(1:n,:))

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

    deallocate(A, iA, idx_map)

  end subroutine weno_tensor_product_calc_recon_matrix

  subroutine weno_tensor_product_calc_ideal_coefs(this, ierr)

    class(weno_tensor_product_type), intent(inout) :: this
    integer, intent(out) :: ierr

    ! Local double double arrays for preserving precision.
    real(16), allocatable, dimension(:,:) :: A, AtA, iAtA, ic
    integer i, j, k, p, m, n, di, dj

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
        di = this%subs(k)%is - this%is
        dj = this%subs(k)%js - this%js
        do j = 1, this%sub_sw**(this%nd-1)
          do i = 1, this%sub_sw
            if (this%subs(k)%cell_mask(this%subs(k)%is+i-1,this%subs(k)%js+j-1) == 1) then
              m = (j + dj - 1) * this%sw     + i + di
              n = (j      - 1) * this%sub_sw + i
              A(m,k) = this%subs(k)%iAp_r16(n,p)
            end if
          end do
        end do
      end do
      AtA = matmul(transpose(A), A)
      call inverse_matrix(AtA, iAtA, ierr)
      if (ierr /= 0) return
      ic(:,p) = matmul(matmul(iAtA, transpose(A)), this%iAp_r16(:,p))
      ! Set near-zero values to zero exactly.
      do k = 1, this%ns
        if (abs(ic(k,p)) < 1.0e-15) ic(k,p) = 0
      end do
    end do
    if (abs(sum(matmul(A, ic(:,1)) - this%iAp_r16(:,1))) > 1.0e-15) then
      ierr = 3
      this%ic = 0
    else
      this%ic = ic
    end if
    deallocate(A, AtA, iAtA, ic)

  end subroutine weno_tensor_product_calc_ideal_coefs

  subroutine weno_tensor_product_calc_smooth_indicator(this)

    class(weno_tensor_product_type), intent(inout) :: this

    associate (a => this%a)
    select case (this%nd)
    case (1)
    case (2)
      select case (this%sw)
      case (3)
        this%si = (  720 * a(2) * a(2) + 3120 * a(3) * a(3) + 720  * a(4) * a(4) + 840   * a(5) * a(5) &
                  +  120 * a(4) * a(6) + 3389 * a(6) * a(6) + 3120 * a(7) * a(7) + 120   * a(2) * a(8) &
                  + 3389 * a(8) * a(8) + 520  * a(3) * a(9) + 520  * a(7) * a(9) + 13598 * a(9) * a(9) ) / 720
      case (4)
      end select
    end select
    end associate

  end subroutine weno_tensor_product_calc_smooth_indicator

  subroutine weno_tensor_product_reconstruct(this, fi, fo, ierr)

    class(weno_tensor_product_type), intent(inout) :: this
    real(8), intent(in ) :: fi(:,:) ! Cell averaged function values
    real(8), intent(out) :: fo(:)   ! Reconstructed function values on evaluation points
    integer, intent(out) :: ierr

    integer k, p
    real(8), parameter :: theta = 3
    real(8), parameter :: eps = 1.0d-6
    real(8) fs(this%ns,this%npt)
    real(8) ic_p(this%ns), sc_p, w_p(this%ns)
    real(8) ic_n(this%ns), sc_n, w_n(this%ns)

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

    ! Calculate point values from each available sub-stencils and smoothness
    ! indicators for those sub-stencils.
    do k = 1, this%ns
      this%subs(k)%a = matmul(this%subs(k)%iA, pack(fi(this%is:this%ie,this%js:this%je), .true.))
      fs(k,:) = matmul(this%subs(k)%a, this%subs(k)%p)
      call this%subs(k)%calc_smooth_indicator()
    end do

    do p = 1, this%npt
      ! Handle negative ideal coefficients by splitting method (Shi et al., 2002).
      ic_p = (this%ic(:,p) + theta * abs(this%ic(:,p))) / 2.0d0 ! Positive part
      ic_n = ic_p - this%ic(:,p)                                ! Negative part
      sc_p = sum(ic_p); ic_p = ic_p / sc_p
      sc_n = sum(ic_n); ic_n = ic_n / sc_n
      w_p = ic_p / (eps + this%subs(:)%si)**2
      w_n = ic_n / (eps + this%subs(:)%si)**2
      fo(p) = sc_p * sum(w_p * fs(:,p)) + sc_n * sum(w_n * fs(:,p))
    end do

  end subroutine weno_tensor_product_reconstruct

  subroutine weno_tensor_product_release_unused_memory(this)

    class(weno_tensor_product_type), intent(inout) :: this

    type(weno_tensor_product_type), allocatable :: subs(:)
    real(8), allocatable :: ic(:,:)
    integer i, k, ns

    if (allocated(this%poly_mask)) deallocate(this%poly_mask)
    if (allocated(this%xc       )) deallocate(this%xc       )
    if (allocated(this%yc       )) deallocate(this%yc       )
    if (allocated(this%x        )) deallocate(this%x        )
    if (allocated(this%y        )) deallocate(this%y        )
    if (allocated(this%iAp_r16  )) deallocate(this%iAp_r16  )

    ! Shrink subs and ic arrays to only contain unmasked sub-stencils.
    ns = this%ns
    do k = 1, this%ns
      if (any(this%subs(k)%cell_mask == 0)) ns = ns - 1
    end do
    if (ns < this%ns) then
      allocate(subs(ns), ic(ns,this%npt))
      i = 1
      do k = 1, this%ns
        if (any(this%subs(k)%cell_mask == 0)) cycle
        subs(i) = this%subs(k)
        ic(i,:) = this%ic(k,:)
        !print *, k, i, ic(i,:)
        i = i + 1
      end do
      deallocate(this%subs, this%ic)
      this%ns = ns
      allocate(this%subs(ns), this%ic(ns,this%npt))
      do k = 1, this%ns
        this%subs(k) = subs(k)
        this%ic(k,:) = ic(k,:)
        !print *, k, ic(k,:)
      end do
      deallocate(subs, ic)
      !stop
    end if

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
    this%is     = 0
    this%ie     = 0
    this%js     = 0
    this%je     = 0

    if (allocated(this%cell_mask)) deallocate(this%cell_mask)
    if (allocated(this%poly_mask)) deallocate(this%poly_mask)
    if (allocated(this%xc       )) deallocate(this%xc       )
    if (allocated(this%yc       )) deallocate(this%yc       )
    if (allocated(this%a        )) deallocate(this%a        )
    if (allocated(this%x        )) deallocate(this%x        )
    if (allocated(this%y        )) deallocate(this%y        )
    if (allocated(this%p        )) deallocate(this%p        )
    if (allocated(this%iA       )) deallocate(this%iA       )
    if (allocated(this%iAp      )) deallocate(this%iAp      )
    if (allocated(this%iAp_r16  )) deallocate(this%iAp_r16  )
    if (allocated(this%ic       )) deallocate(this%ic       )
    if (allocated(this%subs     )) deallocate(this%subs     )

    this%initialized = .false.

  end subroutine weno_tensor_product_clear

  subroutine weno_tensor_product_assign(this, other)

    class(weno_tensor_product_type), intent(inout) :: this
    class(weno_tensor_product_type), intent(in) :: other

    call this%clear()

    this%id     = other%id
    this%nd     = other%nd
    this%sw     = other%sw
    this%sub_sw = other%sub_sw
    this%nc     = other%nc
    this%ns     = other%ns
    this%npt    = other%npt
    this%is     = other%is
    this%ie     = other%ie
    this%js     = other%js
    this%je     = other%je

    allocate(this%cell_mask(this%sw,this%sw**(this%nd-1)))
    allocate(this%a  (this%nc))
    allocate(this%p  (this%nc,this%npt))
    allocate(this%iA (this%nc,this%nc))
    allocate(this%iAp(this%nc,this%npt))

    this%cell_mask = other%cell_mask
    this%p         = other%p
    this%iA        = other%iA
    this%iAp       = other%iAp

    this%initialized = .true.

  end subroutine weno_tensor_product_assign

  subroutine weno_tensor_product_final(this)

    type(weno_tensor_product_type), intent(inout) :: this

    call this%clear()

  end subroutine weno_tensor_product_final

end module weno_tensor_product_mod
