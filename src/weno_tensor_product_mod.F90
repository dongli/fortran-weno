module weno_tensor_product_mod

  use poly_utils_mod
  use math_mod
  use smooth_indicators_mod

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
    real(8) :: beta   = 0                  ! Smoothness indicator
    integer, allocatable :: cell_mask(:,:) ! Mask unavailable cells by 0
    integer, allocatable :: poly_mask(:,:) ! Mask unavailable polynomial terms by 0
    real(8), allocatable :: xc       (:)   ! X coordinate of cell centroids
    real(8), allocatable :: yc       (:)   ! Y coordinate of cell centroids
    real(8), allocatable :: x        (:)   ! X coordinate of evaluation point
    real(8), allocatable :: y        (:)   ! Y coordinate of evaluation point
    real(8), allocatable :: a        (:)   ! Reconstruction coefficients
    real(8), allocatable :: poly     (:,:) ! Polynomial terms on each evaluation point
    real(8), allocatable :: iA       (:,:) ! A * a = f, iA = inverse(A)
    real(8), allocatable :: poly_iA  (:,:) ! poly * iA
    real(8), allocatable :: gamma    (:,:) ! Ideal coefficients for combining sub-stencils (only on stencil)
    real(16), allocatable :: iA_poly (:,:) ! Working arrays
    type(weno_tensor_product_type), allocatable :: subs(:) ! Sub-stencils
    procedure(smooth_indicator_interface), nopass, pointer :: smooth_indicator => null()
  contains
    procedure :: init                  => weno_tensor_product_init
    procedure :: add_point             => weno_tensor_product_add_point
    procedure :: calc_recon_matrix     => weno_tensor_product_calc_recon_matrix
    procedure :: calc_ideal_coefs      => weno_tensor_product_calc_ideal_coefs
    procedure :: release_unused_memory => weno_tensor_product_release_unused_memory
    procedure :: clear                 => weno_tensor_product_clear
    procedure, private :: weno_tensor_product_reconstruct_1d
    procedure, private :: weno_tensor_product_reconstruct_2d
    generic :: reconstruct => weno_tensor_product_reconstruct_1d, &
                              weno_tensor_product_reconstruct_2d
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
          select case (nd)
          case (1)
            select case (this%sub_sw)
            case (3)
              this%subs(k)%smooth_indicator => smooth_indicator_3
            end select
          case (2)
            select case (this%sub_sw)
            case (3)
              this%subs(k)%smooth_indicator => smooth_indicator_3x3
            end select
          end select
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
    integer i, j, k, ipt, n

    ierr = 0

    if (allocated(this%poly   )) deallocate(this%poly   )
    if (allocated(this%iA     )) deallocate(this%iA     )
    if (allocated(this%iA_poly)) deallocate(this%iA_poly)
    if (allocated(this%poly_iA)) deallocate(this%poly_iA)

    allocate(this%poly   (this%nc ,this%npt))
    allocate(this%iA     (this%nc ,this%nc ))
    allocate(this%iA_poly(this%nc ,this%npt))
    allocate(this%poly_iA(this%npt,this%nc ))

    allocate( A     (this%nc,this%nc))
    allocate(iA     (this%nc,this%nc))
    allocate(idx_map(this%nc        )); idx_map = 0

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

    ! Set the poly for each evaluation point.
    ! Select monomials according to mask.
    select case (this%nd)
    case (1)
      do ipt = 1, this%npt
        k = 1
        do i = 1, this%sw
          if (this%poly_mask(i,1) == 1) then
            call calc_monomial(this%x(ipt), i - 1, this%poly(k,ipt))
            k = k + 1
          end if
        end do
      end do
    case (2)
      do ipt = 1, this%npt
        k = 1
        do j = 1, this%sw
          do i = 1, this%sw
            if (this%poly_mask(i,j) == 1) then
              call calc_monomial(this%x(ipt), this%y(ipt), i - 1, j - 1, this%poly(k,ipt))
              k = k + 1
            end if
          end do
        end do
      end do
    end select

    ! Calculate inverse of integral coefficient matrix.
    select case (this%nd)
    case (1)
      call calc_poly_tensor_product_integral_matrix(this%sw, this%xc, A, this%cell_mask(:,1), this%poly_mask(:,1))
    case (2)
      call calc_poly_tensor_product_integral_matrix(this%sw, this%sw, this%xc, this%yc, A, this%cell_mask, this%poly_mask)
    end select
    n = count(this%cell_mask == 1)
    call inverse_matrix(A(1:n,1:n), iA(1:n,1:n), ierr)

    if (ierr /= 0) then
      deallocate(A, iA, idx_map)
      return
    end if

    this%iA_poly(1:n,:) = matmul(iA(1:n,1:n), this%poly(1:n,:))

    ! Copy double double iA into double iA.
    this%iA = transpose(iA) ! Transpose iA for later convenience.

    ! Rearrange iA and iA_poly.
    do n = this%nc, 1, -1
      if (idx_map(n) /= 0) then
        this%iA     (:,n) = this%iA     (:,idx_map(n))
        this%iA_poly(n,:) = this%iA_poly(idx_map(n),:)
      else
        this%iA     (:,n) = 0
        this%iA_poly(n,:) = 0
      end if
    end do

    this%poly_iA = transpose(this%iA_poly)

    deallocate(A, iA, idx_map)

  end subroutine weno_tensor_product_calc_recon_matrix

  subroutine weno_tensor_product_calc_ideal_coefs(this, ierr)

    class(weno_tensor_product_type), intent(inout) :: this
    integer, intent(out) :: ierr

    ! Local double double arrays for preserving precision.
    real(16), allocatable, dimension(:,:) :: A, AtA, iAtA, gamma
    integer i, j, k, ipt, m, n, di, dj

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
    allocate( iAtA(this%ns,this%ns ))
    allocate(gamma(this%ns,this%npt)); gamma = 0
    do ipt = 1, this%npt
      do k = 1, this%ns
        di = this%subs(k)%is - this%is
        dj = this%subs(k)%js - this%js
        do j = 1, this%sub_sw**(this%nd-1)
          do i = 1, this%sub_sw
            if (this%subs(k)%cell_mask(this%subs(k)%is+i-1,this%subs(k)%js+j-1) == 1) then
              m = (j + dj - 1) * this%sw     + i + di
              n = (j      - 1) * this%sub_sw + i
              A(m,k) = this%subs(k)%iA_poly(n,ipt)
            end if
          end do
        end do
      end do
      AtA = matmul(transpose(A), A)
      call inverse_matrix(AtA, iAtA, ierr)
      if (ierr /= 0) return
      gamma(:,ipt) = matmul(matmul(iAtA, transpose(A)), this%iA_poly(:,ipt))
      ! Set near-zero values to zero exactly.
      do k = 1, this%ns
        if (abs(gamma(k,ipt)) < 1.0e-15) gamma(k,ipt) = 0
      end do
    end do
    if (abs(sum(matmul(A, gamma(:,1)) - this%iA_poly(:,1))) > 1.0e-15) then
      ierr = 3
      this%gamma = 0
    else
      this%gamma = gamma
    end if
    deallocate(A, AtA, iAtA, gamma)

  end subroutine weno_tensor_product_calc_ideal_coefs

  subroutine weno_tensor_product_reconstruct_1d(this, fi, fo, ierr)

    class(weno_tensor_product_type), intent(inout) :: this
    real(8), intent(in ) :: fi(:)   ! Cell averaged function values
    real(8), intent(out) :: fo(:)   ! Reconstructed function values on evaluation points
    integer, intent(out) :: ierr

    integer k, ipt
    real(8), parameter :: theta = 3
    real(8), parameter :: eps = 1.0d-40
    real(8) fs(this%ns,this%npt), tau, tmp(this%ns), alpha(this%ns)
    real(8) gamma_p(this%ns), sigma_p, alpha_p(this%ns)
    real(8) gamma_n(this%ns), sigma_n, alpha_n(this%ns)

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
      this%subs(k)%a = matmul(this%subs(k)%iA, fi(this%subs(k)%is:this%subs(k)%ie))
      fs(k,:) = matmul(this%subs(k)%a, this%subs(k)%poly)
      this%subs(k)%beta = this%subs(k)%smooth_indicator(this%subs(k)%a)
    end do

    ! WENO-Z
    tau = abs(this%subs(1)%beta - this%subs(this%ns-1)%beta)

    do ipt = 1, this%npt
      ! Handle negative ideal coefficients by splitting method (Shi et al., 2002).
      gamma_p = (this%gamma(:,ipt) + theta * abs(this%gamma(:,ipt))) / 2.0d0 ! Positive part
      gamma_n = gamma_p - this%gamma(:,ipt)                                  ! Negative part
      sigma_p = sum(gamma_p); gamma_p = gamma_p / sigma_p
      sigma_n = sum(gamma_n); gamma_n = gamma_n / sigma_n
      tmp     = 1 + (tau / (eps + this%subs(:)%beta))**2
      alpha_p = gamma_p * tmp
      alpha_n = gamma_n * tmp
      fo(ipt) = sigma_p * sum(alpha_p / sum(alpha_p) * fs(:,ipt)) - sigma_n * sum(alpha_n / sum(alpha_n) * fs(:,ipt))
    end do

  end subroutine weno_tensor_product_reconstruct_1d

  subroutine weno_tensor_product_reconstruct_2d(this, fi, fo, ierr)

    class(weno_tensor_product_type), intent(inout) :: this
    real(8), intent(in ) :: fi(:,:) ! Cell averaged function values
    real(8), intent(out) :: fo(:)   ! Reconstructed function values on evaluation points
    integer, intent(out) :: ierr

    integer i, j, k, ipt, n
    real(8), parameter :: theta = 3
    real(8), parameter :: eps = 1.0d-40
    real(8) fs(this%ns,this%npt), tau, tmp(this%ns)
    real(8) gamma_p(this%ns), sigma_p, alpha_p(this%ns)
    real(8) gamma_n(this%ns), sigma_n, alpha_n(this%ns)

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
      fs(k,:) = matmul(this%subs(k)%a, this%subs(k)%poly)
      this%subs(k)%beta = this%subs(k)%smooth_indicator(this%subs(k)%a)
    end do

    ! WENO-Z
    tau = 0; n = 0
    do j = 1, this%ns - 1
      do i = j + 1, this%ns
        tau = tau + abs(this%subs(i)%beta - this%subs(j)%beta)
        n = n + 1
      end do
    end do
    tau = tau / n

    do ipt = 1, this%npt
      ! Handle negative ideal coefficients by splitting method (Shi et al., 2002).
      gamma_p = (this%gamma(:,ipt) + theta * abs(this%gamma(:,ipt))) / 2.0d0 ! Positive part
      gamma_n = gamma_p - this%gamma(:,ipt)                                  ! Negative part
      sigma_p = sum(gamma_p); gamma_p = gamma_p / sigma_p
      sigma_n = sum(gamma_n); gamma_n = gamma_n / sigma_n
      tmp     = 1 + (tau / (eps + this%subs(:)%beta))**2
      alpha_p = gamma_p * tmp
      alpha_n = gamma_n * tmp
      fo(ipt) = sigma_p * sum(alpha_p / sum(alpha_p) * fs(:,ipt)) - sigma_n * sum(alpha_n / sum(alpha_n) * fs(:,ipt))
    end do

  end subroutine weno_tensor_product_reconstruct_2d

  subroutine weno_tensor_product_release_unused_memory(this)

    class(weno_tensor_product_type), intent(inout) :: this

    type(weno_tensor_product_type), allocatable :: subs(:)
    real(8), allocatable :: gamma(:,:)
    integer i, k, ns

    if (allocated(this%poly_mask)) deallocate(this%poly_mask)
    if (allocated(this%xc       )) deallocate(this%xc       )
    if (allocated(this%yc       )) deallocate(this%yc       )
    if (allocated(this%x        )) deallocate(this%x        )
    if (allocated(this%y        )) deallocate(this%y        )
    if (allocated(this%iA_poly  )) deallocate(this%iA_poly  )

    ! Shrink subs and gamma arrays to only contain unmasked sub-stencils.
    ns = this%ns
    do k = 1, this%ns
      if (any(this%subs(k)%cell_mask == 0)) ns = ns - 1
    end do
    if (ns < this%ns) then
      allocate(subs(ns), gamma(ns,this%npt))
      i = 1
      do k = 1, this%ns
        if (any(this%subs(k)%cell_mask == 0)) cycle
        subs(i) = this%subs(k)
        gamma(i,:) = this%gamma(k,:)
        i = i + 1
      end do
      deallocate(this%subs, this%gamma)
      this%ns = ns
      allocate(this%subs(ns), this%gamma(ns,this%npt))
      do k = 1, this%ns
        this%subs(k) = subs(k)
        this%gamma(k,:) = gamma(k,:)
      end do
      deallocate(subs, gamma)
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

    this%smooth_indicator => null()

    if (allocated(this%cell_mask)) deallocate(this%cell_mask)
    if (allocated(this%poly_mask)) deallocate(this%poly_mask)
    if (allocated(this%xc       )) deallocate(this%xc       )
    if (allocated(this%yc       )) deallocate(this%yc       )
    if (allocated(this%a        )) deallocate(this%a        )
    if (allocated(this%x        )) deallocate(this%x        )
    if (allocated(this%y        )) deallocate(this%y        )
    if (allocated(this%poly     )) deallocate(this%poly     )
    if (allocated(this%iA       )) deallocate(this%iA       )
    if (allocated(this%poly_iA  )) deallocate(this%poly_iA  )
    if (allocated(this%iA_poly  )) deallocate(this%iA_poly  )
    if (allocated(this%gamma    )) deallocate(this%gamma    )
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
    allocate(this%a        (this%nc))
    allocate(this%poly     (this%nc,this%npt))
    allocate(this%iA       (this%nc,this%nc))
    allocate(this%poly_iA  (this%nc,this%npt))

    this%cell_mask = other%cell_mask
    this%poly      = other%poly
    this%iA        = other%iA
    this%poly_iA   = other%poly_iA

    this%initialized = .true.

  end subroutine weno_tensor_product_assign

  subroutine weno_tensor_product_final(this)

    type(weno_tensor_product_type), intent(inout) :: this

    call this%clear()

  end subroutine weno_tensor_product_final

end module weno_tensor_product_mod
