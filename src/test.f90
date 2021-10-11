program test

  use weno

  implicit none

  integer mask1d(5,1)
  integer mask2d(5,5)

  mask2d = 1
  call test_2d_5x5(mask2d, 1)

  mask2d = 1
  mask2d(5,5) = 0
  call test_2d_5x5(mask2d, 1)

  mask2d = 1
  mask2d(1,1) = 0
  call test_2d_5x5(mask2d, 2)

  mask2d = 1
  mask2d(5,5) = 0; mask2d(4,5) = 0
  call test_2d_5x5(mask2d, 1)

  mask2d = 1
  mask2d(1,1) = 0; mask2d(2,1) = 0
  call test_2d_5x5(mask2d, 2)

  mask2d = 1
  mask2d(5,5) = 0; mask2d(5,4) = 0
  call test_2d_5x5(mask2d, 1)

  mask2d = 1
  mask2d(1,1) = 0; mask2d(1,2) = 0
  call test_2d_5x5(mask2d, 2)

  mask2d = 1
  mask2d(5,5) = 0; mask2d(4,5) = 0; mask2d(5,4) = 0; mask2d(4,4) = 0
  call test_2d_5x5(mask2d, 1)

  mask2d = 1
  mask2d(1,1) = 0; mask2d(1,2) = 0; mask2d(2,1) = 0; mask2d(2,2) = 0
  call test_2d_5x5(mask2d, 2)

  mask2d = 1
  mask2d(5,1) = 0; mask2d(5,2) = 0; mask2d(4,1) = 0; mask2d(4,2) = 0
  call test_2d_5x5(mask2d, 1)

  mask2d = 1
  mask2d(1,5) = 0; mask2d(2,5) = 0; mask2d(1,4) = 0; mask2d(2,4) = 0
  call test_2d_5x5(mask2d, 2)

  mask1d = 1
  call test_1d_5(mask1d, 1)

  mask1d = 1
  call test_1d_5(mask1d, 2)

  mask1d = 1
  mask1d(5,1) = 0
  call test_1d_5(mask1d, 2)

  mask1d = 1
  mask1d(1,1) = 0
  call test_1d_5(mask1d, 1)

  mask1d = 1
  mask1d(5,1) = 0
  call test_1d_5(mask1d, 1)

  mask1d = 1
  mask1d(1,1) = 0
  call test_1d_5(mask1d, 2)

  mask1d = 1
  mask1d(1,1) = 0; mask1d(2,1) = 0
  call test_1d_5(mask1d, 2)

  mask1d = 1
  mask1d(5,1) = 0; mask1d(4,1) = 0
  call test_1d_5(mask1d, 1)

contains

  subroutine test_1d_5(mask, check_point_idx)

    integer, intent(in) :: mask(5,1)
    integer, intent(in) :: check_point_idx

    type(weno_tensor_product_type), allocatable :: weno1d
    integer i, ierr

    write(*, *) 'Mask:'
    do i = 1, 5
      write(*, '(I5)', advance='no') mask(i,1)
    end do
    write(*, *)

    allocate(weno1d)

    !  _____________________________
    ! |     |     |     |     |     |
    ! |     |     x i0  x     |     |
    ! |_____|_____|_____|_____|_____|

    call weno1d%init(nd=1, sw=5, mask=mask)
    call weno1d%add_point(x=-0.5d0)
    call weno1d%add_point(x= 0.5d0)
    call weno1d%calc_ideal_coefs(ierr)
    if (ierr /= 0) then
      write(*, *) 'Failed to calculate WENO ideal coefficients!'
      stop 1
    end if

    write(*, *) 'Ideal sub-stencil coefficients:'
    do i = 1, weno1d%ns
      write(*, '(F10.5)', advance='no') weno1d%ic(i,check_point_idx)
    end do
    write(*, *)
    write(*, *)

    deallocate(weno1d)

  end subroutine test_1d_5

  subroutine test_2d_5x5(mask, check_point_idx)

    integer, intent(in) :: mask(5,5)
    integer, intent(in) :: check_point_idx

    type(weno_tensor_product_type), allocatable :: weno2d
    integer i, j, ierr

    write(*, *) 'Mask:'
    do j = 5, 1, -1
      do i = 1, 5
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

    call weno2d%init(nd=2, sw=5, mask=mask)
    call weno2d%add_point(x=-0.5d0, y= 0.0d0)
    call weno2d%add_point(x= 0.5d0, y= 0.0d0)
    call weno2d%add_point(x= 0.0d0, y=-0.5d0)
    call weno2d%add_point(x= 0.0d0, y= 0.5d0)
    call weno2d%calc_ideal_coefs(ierr)
    if (ierr /= 0) then
      write(*, *) 'Failed to calculate WENO ideal coefficients!'
      write(*, *)
      return
    end if

    write(*, *) 'Ideal sub-stencil coefficients:'
    do i = 1, weno2d%ns
      write(*, '(F10.5)', advance='no') weno2d%ic(i,check_point_idx)
    end do
    write(*, *)
    write(*, *)

    deallocate(weno2d)

  end subroutine test_2d_5x5

end program test
