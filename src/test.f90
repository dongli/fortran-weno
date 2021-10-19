program test

  use weno

  implicit none

  real(8), parameter :: pi = 4.0d0 * atan(1.0d0)

  integer mask_5(5,1)
  integer mask_3x3(3,3)
  integer mask_5x5(5,5)
  integer mask_7x7(7,7)

  call test_poly_1d(5)

  mask_7x7 = 1
  call test_weno_2d(mask_7x7, 1)

  mask_3x3 = 1
  call test_weno_2d(mask_3x3, 1)

  mask_5x5 = 1
  call test_weno_2d(mask_5x5, 1)

  mask_5x5 = 1
  mask_5x5(5,5) = 0
  call test_weno_2d(mask_5x5, 1)

  mask_5x5 = 1
  mask_5x5(1,1) = 0
  call test_weno_2d(mask_5x5, 2)

  mask_5x5 = 1
  mask_5x5(5,5) = 0; mask_5x5(4,5) = 0
  call test_weno_2d(mask_5x5, 1)

  mask_5x5 = 1
  mask_5x5(1,1) = 0; mask_5x5(2,1) = 0
  call test_weno_2d(mask_5x5, 2)

  mask_5x5 = 1
  mask_5x5(5,5) = 0; mask_5x5(5,4) = 0
  call test_weno_2d(mask_5x5, 1)

  mask_5x5 = 1
  mask_5x5(1,1) = 0; mask_5x5(1,2) = 0
  call test_weno_2d(mask_5x5, 2)

  mask_5x5 = 1
  mask_5x5(5,5) = 0; mask_5x5(4,5) = 0; mask_5x5(5,4) = 0; mask_5x5(4,4) = 0
  call test_weno_2d(mask_5x5, 1)

  mask_5x5 = 1
  mask_5x5(1,1) = 0; mask_5x5(1,2) = 0; mask_5x5(2,1) = 0; mask_5x5(2,2) = 0
  call test_weno_2d(mask_5x5, 2)

  mask_5x5 = 1
  mask_5x5(5,1) = 0; mask_5x5(5,2) = 0; mask_5x5(4,1) = 0; mask_5x5(4,2) = 0
  call test_weno_2d(mask_5x5, 1)

  mask_5x5 = 1
  mask_5x5(1,5) = 0; mask_5x5(2,5) = 0; mask_5x5(1,4) = 0; mask_5x5(2,4) = 0
  call test_weno_2d(mask_5x5, 2)

  mask_5 = 1
  call test_weno_1d(mask_5, 1)

  mask_5 = 1
  call test_weno_1d(mask_5, 2)

  mask_5 = 1
  mask_5(5,1) = 0
  call test_weno_1d(mask_5, 2)

  mask_5 = 1
  mask_5(1,1) = 0
  call test_weno_1d(mask_5, 1)

  mask_5 = 1
  mask_5(5,1) = 0
  call test_weno_1d(mask_5, 1)

  mask_5 = 1
  mask_5(1,1) = 0
  call test_weno_1d(mask_5, 2)

  mask_5 = 1
  mask_5(1,1) = 0; mask_5(2,1) = 0
  call test_weno_1d(mask_5, 2)

  mask_5 = 1
  mask_5(5,1) = 0; mask_5(4,1) = 0
  call test_weno_1d(mask_5, 1)

  call test_acker2016()

  call test_sine()

contains

  subroutine test_weno_1d(mask, check_point_idx)

    integer, intent(in) :: mask(:,:)
    integer, intent(in) :: check_point_idx

    type(weno_tensor_product_type), allocatable :: weno1d
    integer sw, i, k, ierr

    sw = size(mask, 1)

    write(*, *) 'Mask:'
    do i = 1, sw
      write(*, '(I5)', advance='no') mask(i,1)
    end do
    write(*, *)

    allocate(weno1d)

    !  _____________________________
    ! |     |     |     |     |     |
    ! |     |     |  0  x     |     |
    ! |_____|_____|_____|_____|_____|
    !
    call weno1d%init(nd=1, sw=sw, mask=mask)
    call weno1d%add_point(x=-0.5d0)
    call weno1d%add_point(x= 0.5d0)
    call weno1d%calc_ideal_coefs(ierr)
    if (ierr /= 0) then
      write(*, *) 'Failed to calculate WENO ideal coefficients!'
      stop 1
    end if

    call weno1d%release_unused_memory()

    write(*, *) 'Ideal sub-stencil coefficients:'
    do i = 1, weno1d%sub_sw
      do k = 1, weno1d%ns
        if (weno1d%subs(k)%id == i) then
          write(*, '(F10.5)', advance='no') weno1d%gamma(k,check_point_idx)
          exit
        end if
      end do
      if (k == weno1d%ns + 1) write(*, '(A10)', advance='no') ' XXXXXXXXX'
    end do
    write(*, *)
    write(*, *)

    deallocate(weno1d)

  end subroutine test_weno_1d

  subroutine test_weno_2d(mask, check_point_idx)

    integer, intent(in) :: mask(:,:)
    integer, intent(in) :: check_point_idx

    type(weno_tensor_product_type), allocatable :: weno2d
    integer sw, i, j, k, ierr

    sw = size(mask, 1)

    write(*, *) 'Mask:'
    do j = sw, 1, -1
      do i = 1, sw
        write(*, '(I5)', advance='no') mask(i,j)
      end do
      write(*, *)
    end do
    write(*, *)

    allocate(weno2d)

    !  _____________________________
    ! |     |     |     |     |     |
    ! |     |     |     |     |     |
    ! |_____|_____|_____|_____|_____|
    ! |     |     |     |     |     |
    ! |     |     |     |     |     |
    ! |_____|_____|__x__|_____|_____|
    ! |     |     |     |     |     |
    ! |     |     xi0,j0x     |     |
    ! |_____|_____|__x__|_____|_____|
    ! |     |     |     |     |     |
    ! |     |     |     |     |     |
    ! |_____|_____|_____|_____|_____|
    ! |     |     |     |     |     |
    ! |     |     |     |     |     |
    ! |_____|_____|_____|_____|_____|

    call weno2d%init(nd=2, sw=sw, is=-int(sw/2), ie=int(sw/2), js=-int(sw/2), je=int(sw/2), mask=mask)
    select case (sw)
    case (3, 7)
      call weno2d%add_point(x=-0.5d0, y= 0.001d0)
    case (5)
      call weno2d%add_point(x=-0.5d0, y= 0.0d0)
      call weno2d%add_point(x= 0.5d0, y= 0.0d0)
      call weno2d%add_point(x= 0.0d0, y=-0.5d0)
      call weno2d%add_point(x= 0.0d0, y= 0.5d0)
    end select
    call weno2d%calc_ideal_coefs(ierr)
    if (ierr /= 0) then
      write(*, *) 'Failed to calculate WENO ideal coefficients!'
      write(*, *)
      return
    end if

    call weno2d%release_unused_memory()

    write(*, *) 'Ideal sub-stencil coefficients:'
    do j = weno2d%sub_sw, 1, -1
      do i = 1, weno2d%sub_sw
        do k = 1, weno2d%ns
          if (weno2d%subs(k)%id == (j - 1) * weno2d%sub_sw + i) then
            write(*, '(F10.5)', advance='no') weno2d%gamma(k,check_point_idx)
            exit
          end if
        end do
        if (k == weno2d%ns + 1) write(*, '(A10)', advance='no') ' XXXXXXXXX'
      end do
      write(*, *)
    end do
    write(*, *)
    write(*, *)

    deallocate(weno2d)

  end subroutine test_weno_2d

  subroutine test_poly_1d(sw)

    integer, intent(in) :: sw

    type(poly_tensor_product_type), allocatable :: poly1d
    real(8) dx, fi(sw), fo(20)
    integer i, ierr

    allocate(poly1d)

    call poly1d%init(nd=1, sw=sw)
    dx = dble(sw) / (size(fo) - 1)
    do i = 1, size(fo)
      call poly1d%add_point(x=-sw/2.0d0+(i-1)*dx)
    end do
    call poly1d%calc_recon_matrix(ierr)
    if (ierr /= 0) then
      write(*, *) 'Failed to calculate polynomial coefficients!'
      return
    end if

    !dx = 1
    !do i = 1, sw
    !  fi(i) = sw / (2 * pi) * (cos(2 * pi / sw * poly1d%xc(i) - 0.5d0) - cos(2 * pi / sw * poly1d%xc(i) + 0.5d0))
    !end do
    !call poly1d%reconstruct(fi, fo, ierr)
    !do i = 1, size(fo)
    !  write(*, '(F10.5, A1)', advance='no') fo(i), ','
    !end do
    !write(*, *)

    deallocate(poly1d)

  end subroutine test_poly_1d

  subroutine test_poly_2d(sw)

    integer, intent(in) :: sw

    type(poly_tensor_product_type), allocatable :: poly2d
    integer ierr

    allocate(poly2d)

    call poly2d%init(nd=2, sw=sw)
    call poly2d%add_point(x=-0.5d0, y=0.0d0)
    call poly2d%add_point(x=-sw/2.0d0, y=-sw/2.0d0)
    call poly2d%calc_recon_matrix(ierr)
    if (ierr /= 0) then
      print *, ierr
      write(*, *) 'Failed to calculate polynomial coefficients!'
      write(*, *)
      return
    end if

    call poly2d%release_unused_memory()

    deallocate(poly2d)

  end subroutine test_poly_2d

  subroutine test_sine()

    type(weno_tensor_product_type), allocatable :: weno1d
    integer mask(5,1)
    real(8) dx, xi(5), fi(5), fo(2)
    integer sw, i, k, ierr

    allocate(weno1d)

    mask = 1
    call weno1d%init(nd=1, sw=5, mask=mask)
    call weno1d%add_point(x=-0.5d0)
    call weno1d%add_point(x= 0.5d0)
    call weno1d%calc_ideal_coefs(ierr)
    call weno1d%release_unused_memory()

    dx = 2 * pi / 160
    do i = 1, 5
      xi(i) = (30 + i - 1) * dx
    end do
    fi = sin(xi)
    call weno1d%reconstruct(fi, fo, ierr)

    write(*, *) cos(xi(3)), (fo(2) - fo(1)) / dx

    deallocate(weno1d)

  end subroutine test_sine

  subroutine test_acker2016()

    type(weno_tensor_product_type), allocatable :: weno1d
    integer mask(5,1)
    real(8) dx, fi(5), fo(2)
    integer sw, i, k, ierr

    allocate(weno1d)

    mask = 1
    call weno1d%init(nd=1, sw=5, mask=mask)
    call weno1d%add_point(x=-0.5d0)
    call weno1d%add_point(x= 0.5d0)
    call weno1d%calc_ideal_coefs(ierr)
    call weno1d%release_unused_memory()

    dx = 0.05
    fi = [0.0, 0.6015171347078453, 0.7544760960555312, 0.37714810606547045, -0.18077485348614492]
    call weno1d%reconstruct(fi, fo, ierr)

    write(*, *) 0.2231, weno1d%subs(1)%beta
    write(*, *) 0.3172, weno1d%subs(2)%beta
    write(*, *) 0.1177, weno1d%subs(3)%beta
    write(*, *) -2.7820226261239465, (fo(2) - fo(1)) / dx

    deallocate(weno1d)

  end subroutine test_acker2016

end program test
