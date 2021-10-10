module poly_utils_mod

  ! This module provides some utility procedures for dealing with polynomials.

  implicit none

  private

  public calc_poly_tensor_product_monomial
  public calc_poly_tensor_product_integral_coefs
  public calc_poly_tensor_product_integral_coef_matrix

contains

  subroutine calc_poly_tensor_product_monomial(x, y, i, j, val)

    real(8) , intent(in ) :: x
    real(8) , intent(in ) :: y
    integer , intent(in ) :: i
    integer , intent(in ) :: j
    real(16), intent(out) :: val

    val = x**i * y**j

  end subroutine calc_poly_tensor_product_monomial

  subroutine calc_poly_tensor_product_integral_coefs(nx, ny, x0, x1, y0, y1, coefs, mask)

    integer , intent(in ) :: nx
    integer , intent(in ) :: ny
    real(16), intent(in ) :: x0
    real(16), intent(in ) :: x1
    real(16), intent(in ) :: y0
    real(16), intent(in ) :: y1
    real(16), intent(out) :: coefs(:)
    integer , intent(in ), optional :: mask(nx,ny)

    integer i, j, k

    k = 1
    do j = 1, ny
      do i = 1, nx
        if (merge(mask(i,j) == 1, .true., present(mask))) then
          coefs(k) = (x1**i - x0**i) * (y1**j - y0**j) / (i * j)
          k = k + 1
        end if
      end do
    end do

  end subroutine calc_poly_tensor_product_integral_coefs

  subroutine calc_poly_tensor_product_integral_coef_matrix(nx, ny, A, mask)

    integer , intent(in ) :: nx
    integer , intent(in ) :: ny
    real(16), intent(out) :: A(:,:)
    integer , intent(in ), optional :: mask(nx,ny)

    integer i, j, k
    real(16) x, y

    k  = 1
    do j = 1, ny
      do i = 1, nx
        if (merge(mask(i,j) == 1, .true., present(mask))) then
          x = i - 1 - int((nx - 1) / 2)
          y = j - 1 - int((ny - 1) / 2)
          call calc_poly_tensor_product_integral_coefs( &
            nx, ny, x - 0.5_16, x + 0.5_16, y - 0.5_16, y + 0.5_16, A(:,k), mask)
          k = k + 1
        end if
      end do
    end do

  end subroutine calc_poly_tensor_product_integral_coef_matrix

end module poly_utils_mod
