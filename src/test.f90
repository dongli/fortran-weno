program test

  use weno

  implicit none

  integer mask_5(5,1)
  integer mask_3x3(3,3)
  integer mask_5x5(5,5)
  integer mask_7x7(7,7)

  mask_7x7 = 1
  call test_2d(mask_7x7, 1)

  mask_3x3 = 1
  call test_2d(mask_3x3, 1)

  mask_5x5 = 1
  call test_2d(mask_5x5, 1)

  mask_5x5 = 1
  mask_5x5(5,5) = 0
  call test_2d(mask_5x5, 1)

  mask_5x5 = 1
  mask_5x5(1,1) = 0
  call test_2d(mask_5x5, 2)

  mask_5x5 = 1
  mask_5x5(5,5) = 0; mask_5x5(4,5) = 0
  call test_2d(mask_5x5, 1)

  mask_5x5 = 1
  mask_5x5(1,1) = 0; mask_5x5(2,1) = 0
  call test_2d(mask_5x5, 2)

  mask_5x5 = 1
  mask_5x5(5,5) = 0; mask_5x5(5,4) = 0
  call test_2d(mask_5x5, 1)

  mask_5x5 = 1
  mask_5x5(1,1) = 0; mask_5x5(1,2) = 0
  call test_2d(mask_5x5, 2)

  mask_5x5 = 1
  mask_5x5(5,5) = 0; mask_5x5(4,5) = 0; mask_5x5(5,4) = 0; mask_5x5(4,4) = 0
  call test_2d(mask_5x5, 1)

  mask_5x5 = 1
  mask_5x5(1,1) = 0; mask_5x5(1,2) = 0; mask_5x5(2,1) = 0; mask_5x5(2,2) = 0
  call test_2d(mask_5x5, 2)

  mask_5x5 = 1
  mask_5x5(5,1) = 0; mask_5x5(5,2) = 0; mask_5x5(4,1) = 0; mask_5x5(4,2) = 0
  call test_2d(mask_5x5, 1)

  mask_5x5 = 1
  mask_5x5(1,5) = 0; mask_5x5(2,5) = 0; mask_5x5(1,4) = 0; mask_5x5(2,4) = 0
  call test_2d(mask_5x5, 2)

  mask_5 = 1
  call test_1d(mask_5, 1)

  mask_5 = 1
  call test_1d(mask_5, 2)

  mask_5 = 1
  mask_5(5,1) = 0
  call test_1d(mask_5, 2)

  mask_5 = 1
  mask_5(1,1) = 0
  call test_1d(mask_5, 1)

  mask_5 = 1
  mask_5(5,1) = 0
  call test_1d(mask_5, 1)

  mask_5 = 1
  mask_5(1,1) = 0
  call test_1d(mask_5, 2)

  mask_5 = 1
  mask_5(1,1) = 0; mask_5(2,1) = 0
  call test_1d(mask_5, 2)

  mask_5 = 1
  mask_5(5,1) = 0; mask_5(4,1) = 0
  call test_1d(mask_5, 1)

contains

  subroutine test_1d(mask, check_point_idx)

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
    ! |     |     x i0  x     |     |
    ! |_____|_____|_____|_____|_____|

    call weno1d%init(nd=1, sw=sw, is=-int(sw/2), ie=int(sw/2), mask=mask)
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
          write(*, '(F10.5)', advance='no') weno1d%ic(k,check_point_idx)
          exit
        end if
      end do
      if (k == weno1d%ns + 1) write(*, '(A10)', advance='no') ' XXXXXXXXX'
    end do
    write(*, *)
    write(*, *)

    deallocate(weno1d)

  end subroutine test_1d

  subroutine test_2d(mask, check_point_idx)

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
      print *, ierr
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
            write(*, '(F10.5)', advance='no') weno2d%ic(k,check_point_idx)
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

  end subroutine test_2d

end program test
