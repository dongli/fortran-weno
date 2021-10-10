program test

  use weno

  implicit none

  type(weno_tensor_product_type), allocatable :: weno2d
  integer mask(5,5)
  integer i, j, ierr

  mask = 1
  mask(5,5) = 0
  mask(4,5) = 0

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
    stop 1
  end if

  deallocate(weno2d)

end program test
