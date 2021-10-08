module poly_utils_mod

  ! This module provides some utility procedures for dealing with polynomials.

  implicit none

  private

  public calc_poly_tensor_product_integral_coefs
  public calc_poly_tensor_product_integral_coef_matrix

contains

  subroutine calc_poly_tensor_product_integral_coefs(nx, ny, x0, x1, y0, y1, coefs, mask)

    integer, intent(in ) :: nx
    integer, intent(in ) :: ny
    real(8), intent(in ) :: x0
    real(8), intent(in ) :: x1
    real(8), intent(in ) :: y0
    real(8), intent(in ) :: y1
    real(8), intent(out) :: coefs(nx*ny)
    integer, intent(in ), optional :: mask(nx,ny)

    integer i, j, k

    k = 1
    do j = 1, ny
      do i = 1, nx
        if (merge(mask(i,j) == 1, .true., present(mask))) then
          coefs(k) = (x1**i - x0**i) * (y1**j - y0**j) / (i * j)
        end if
        k = k + 1
      end do
    end do

  end subroutine calc_poly_tensor_product_integral_coefs

  subroutine calc_poly_tensor_product_integral_coef_matrix(nx, ny, A, mask)

    integer, intent(in ) :: nx
    integer, intent(in ) :: ny
    real(8), intent(out) :: A(nx*ny,nx*ny)
    integer, intent(in ), optional :: mask(nx,ny)

    real(8) x, y, y1, dx, dy, dx0p5, dy0p5
    integer i0, j0, i, j, k

    ! Assume dx and dy is unity.
    dx = 1
    dy = 1

    dx0p5 = 0.5d0 * dx
    dy0p5 = 0.5d0 * dy

    i0 = -int((nx - 1) * 0.5d0)
    j0 = -int((ny - 1) * 0.5d0)

    k = 1
    do j = 0, ny - 1
      y = (j0 + j) * dy
      do i = 0, nx - 1
        x = (i0 + i) * dx
        call calc_poly_tensor_product_integral_coefs( &
          nx, ny, x - dx0p5, x + dx0p5, y - dy0p5, y + dy0p5, A(:,k), mask)
        k = k + 1
      end do
    end do

  end subroutine calc_poly_tensor_product_integral_coef_matrix

end module poly_utils_mod
