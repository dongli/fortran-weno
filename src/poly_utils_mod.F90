module poly_utils_mod

  ! This module provides some utility procedures for dealing with polynomials.

  implicit none

  private

  public calc_monomial
  public calc_poly_tensor_product_integral_coefs
  public calc_poly_tensor_product_integral_matrix
  public calc_poly_integral_coefs
  public calc_poly_integral_matrix

  interface calc_monomial
    module procedure calc_monomial_1d_r8
    module procedure calc_monomial_1d_r16
    module procedure calc_monomial_2d_r8
    module procedure calc_monomial_2d_r16
  end interface calc_monomial

  interface calc_poly_tensor_product_integral_coefs
    module procedure calc_poly_tensor_product_integral_coefs_1d_r8
    module procedure calc_poly_tensor_product_integral_coefs_1d_r16
    module procedure calc_poly_tensor_product_integral_coefs_2d_r8
    module procedure calc_poly_tensor_product_integral_coefs_2d_r16
  end interface calc_poly_tensor_product_integral_coefs

  interface calc_poly_tensor_product_integral_matrix
    module procedure calc_poly_tensor_product_integral_matrix_1d_r8
    module procedure calc_poly_tensor_product_integral_matrix_1d_r16
    module procedure calc_poly_tensor_product_integral_matrix_2d_r8
    module procedure calc_poly_tensor_product_integral_matrix_2d_r16
  end interface calc_poly_tensor_product_integral_matrix

  interface calc_poly_integral_coefs
    module procedure calc_poly_integral_coefs_r8
    module procedure calc_poly_integral_coefs_r16
  end interface calc_poly_integral_coefs

  interface calc_poly_integral_matrix
    module procedure calc_poly_integral_matrix_r8
    module procedure calc_poly_integral_matrix_r16
  end interface calc_poly_integral_matrix

contains

  subroutine calc_monomial_1d_r8(x, i, val)

    real(8), intent(in ) :: x
    integer, intent(in ) :: i
    real(8), intent(out) :: val

    val = x**i

  end subroutine calc_monomial_1d_r8

  subroutine calc_monomial_1d_r16(x, i, val)

    real(16), intent(in ) :: x
    integer , intent(in ) :: i
    real(16), intent(out) :: val

    val = x**i

  end subroutine calc_monomial_1d_r16

  subroutine calc_monomial_2d_r8(x, y, i, j, val)

    real(8), intent(in ) :: x
    real(8), intent(in ) :: y
    integer, intent(in ) :: i
    integer, intent(in ) :: j
    real(8), intent(out) :: val

    val = x**i * y**j

  end subroutine calc_monomial_2d_r8

  subroutine calc_monomial_2d_r16(x, y, i, j, val)

    real(16), intent(in ) :: x
    real(16), intent(in ) :: y
    integer , intent(in ) :: i
    integer , intent(in ) :: j
    real(16), intent(out) :: val

    val = x**i * y**j

  end subroutine calc_monomial_2d_r16

  subroutine calc_poly_tensor_product_integral_coefs_1d_r8(nx, x0, x1, coefs, mask)

    integer , intent(in ) :: nx
    real( 8), intent(in ) :: x0
    real( 8), intent(in ) :: x1
    real(16), intent(out) :: coefs(nx)
    integer , intent(in ), optional :: mask(nx)

    integer i, k

    k = 1
    do i = 1, nx
      if (merge(mask(i) == 1, .true., present(mask))) then
        coefs(k) = (x1**i - x0**i) / i
        k = k + 1
      end if
    end do

  end subroutine calc_poly_tensor_product_integral_coefs_1d_r8

  subroutine calc_poly_tensor_product_integral_coefs_1d_r16(nx, x0, x1, coefs, mask)

    integer , intent(in ) :: nx
    real(16), intent(in ) :: x0
    real(16), intent(in ) :: x1
    real(16), intent(out) :: coefs(nx)
    integer , intent(in ), optional :: mask(nx)

    integer i, k

    k = 1
    do i = 1, nx
      if (merge(mask(i) == 1, .true., present(mask))) then
        coefs(k) = (x1**i - x0**i) / i
        k = k + 1
      end if
    end do

  end subroutine calc_poly_tensor_product_integral_coefs_1d_r16

  subroutine calc_poly_tensor_product_integral_coefs_2d_r8(nx, ny, x0, x1, y0, y1, coefs, mask)

    integer , intent(in ) :: nx
    integer , intent(in ) :: ny
    real( 8), intent(in ) :: x0
    real( 8), intent(in ) :: x1
    real( 8), intent(in ) :: y0
    real( 8), intent(in ) :: y1
    real(16), intent(out) :: coefs(nx*ny)
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

  end subroutine calc_poly_tensor_product_integral_coefs_2d_r8

  subroutine calc_poly_tensor_product_integral_coefs_2d_r16(nx, ny, x0, x1, y0, y1, coefs, mask)

    integer , intent(in ) :: nx
    integer , intent(in ) :: ny
    real(16), intent(in ) :: x0
    real(16), intent(in ) :: x1
    real(16), intent(in ) :: y0
    real(16), intent(in ) :: y1
    real(16), intent(out) :: coefs(nx*ny)
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

  end subroutine calc_poly_tensor_product_integral_coefs_2d_r16

  subroutine calc_poly_tensor_product_integral_matrix_1d_r8(nx, x, A, cell_mask, poly_mask)

    integer , intent(in ) :: nx
    real( 8), intent(in ) :: x(nx)
    real(16), intent(out) :: A(nx,nx)
    integer , intent(in ), optional :: cell_mask(nx)
    integer , intent(in ), optional :: poly_mask(nx)

    integer i, k

    k = 1
    do i = 1, nx
      if (merge(cell_mask(i) == 1, .true., present(cell_mask))) then
        call calc_poly_tensor_product_integral_coefs( &
          nx, x(i) - 0.5_8, x(i) + 0.5_8, A(:,k), poly_mask)
        k = k + 1
      end if
    end do

  end subroutine calc_poly_tensor_product_integral_matrix_1d_r8

  subroutine calc_poly_tensor_product_integral_matrix_1d_r16(nx, x, A, cell_mask, poly_mask)

    integer , intent(in ) :: nx
    real(16), intent(in ) :: x(nx)
    real(16), intent(out) :: A(nx,nx)
    integer , intent(in ), optional :: cell_mask(nx)
    integer , intent(in ), optional :: poly_mask(nx)

    integer i, k

    k = 1
    do i = 1, nx
      if (merge(cell_mask(i) == 1, .true., present(cell_mask))) then
        call calc_poly_tensor_product_integral_coefs( &
          nx, x(i) - 0.5_16, x(i) + 0.5_16, A(:,k), poly_mask)
        k = k + 1
      end if
    end do

  end subroutine calc_poly_tensor_product_integral_matrix_1d_r16

  subroutine calc_poly_tensor_product_integral_matrix_2d_r8(nx, ny, x, y, A, cell_mask, poly_mask)

    integer , intent(in ) :: nx
    integer , intent(in ) :: ny
    real( 8), intent(in ) :: x(nx)
    real( 8), intent(in ) :: y(ny)
    real(16), intent(out) :: A(nx*ny,nx*ny)
    integer , intent(in ), optional :: cell_mask(nx,ny)
    integer , intent(in ), optional :: poly_mask(nx,ny)

    integer i, j, k

    k  = 1
    do j = 1, ny
      do i = 1, nx
        if (merge(cell_mask(i,j) == 1, .true., present(cell_mask))) then
          call calc_poly_tensor_product_integral_coefs( &
            nx, ny, x(i) - 0.5_16, x(i) + 0.5_16, y(j) - 0.5_16, y(j) + 0.5_16, A(:,k), poly_mask)
          k = k + 1
        end if
      end do
    end do

  end subroutine calc_poly_tensor_product_integral_matrix_2d_r8

  subroutine calc_poly_tensor_product_integral_matrix_2d_r16(nx, ny, x, y, A, cell_mask, poly_mask)

    integer , intent(in ) :: nx
    integer , intent(in ) :: ny
    real(16), intent(in ) :: x(nx)
    real(16), intent(in ) :: y(ny)
    real(16), intent(out) :: A(nx*ny,nx*ny)
    integer , intent(in ), optional :: cell_mask(nx,ny)
    integer , intent(in ), optional :: poly_mask(nx,ny)

    integer i, j, k

    k  = 1
    do j = 1, ny
      do i = 1, nx
        if (merge(cell_mask(i,j) == 1, .true., present(cell_mask))) then
          call calc_poly_tensor_product_integral_coefs( &
            nx, ny, x(i) - 0.5_16, x(i) + 0.5_16, y(j) - 0.5_16, y(j) + 0.5_16, A(:,k), poly_mask)
          k = k + 1
        end if
      end do
    end do

  end subroutine calc_poly_tensor_product_integral_matrix_2d_r16

  subroutine calc_poly_integral_coefs_r8(np, poly, x0, x1, y0, y1, coefs)

    integer , intent(in ) :: np
    integer , intent(in ) :: poly(2,np)
    real( 8), intent(in ) :: x0
    real( 8), intent(in ) :: x1
    real( 8), intent(in ) :: y0
    real( 8), intent(in ) :: y1
    real(16), intent(out) :: coefs(np)

    integer i, j, k

    do k = 1, np
      i = poly(1,k) + 1; j = poly(2,k) + 1
      coefs(k) = (x1**i - x0**i) * (y1**j - y0**j) / (i * j)
    end do

  end subroutine calc_poly_integral_coefs_r8

  subroutine calc_poly_integral_coefs_r16(np, poly, x0, x1, y0, y1, coefs)

    integer , intent(in ) :: np
    integer , intent(in ) :: poly(2,np)
    real(16), intent(in ) :: x0
    real(16), intent(in ) :: x1
    real(16), intent(in ) :: y0
    real(16), intent(in ) :: y1
    real(16), intent(out) :: coefs(np)

    integer i, j, k

    do k = 1, np
      i = poly(1,k) + 1; j = poly(2,k) + 1
      coefs(k) = (x1**i - x0**i) * (y1**j - y0**j) / (i * j)
    end do

  end subroutine calc_poly_integral_coefs_r16

  subroutine calc_poly_integral_matrix_r8(nc, np, poly, x, y, A)

    integer , intent(in ) :: nc         ! Number of cells
    integer , intent(in ) :: np         ! Number of monomials
    integer , intent(in ) :: poly(2,np) ! Monomials degrees (x^?, y^?)
    real( 8), intent(in ) :: x(nc)
    real( 8), intent(in ) :: y(nc)
    real(16), intent(out) :: A(np,nc)

    integer i

    do i = 1, nc
      call calc_poly_integral_coefs(np, poly, x(i) - 0.5_16, x(i) + 0.5_16, y(i) - 0.5_16, y(i) + 0.5_16, A(:,i))
    end do

  end subroutine calc_poly_integral_matrix_r8

  subroutine calc_poly_integral_matrix_r16(nc, np, poly, x, y, A)

    integer , intent(in ) :: nc         ! Number of cells
    integer , intent(in ) :: np         ! Number of monomials
    integer , intent(in ) :: poly(2,np) ! Monomials degrees (x^?, y^?)
    real(16), intent(in ) :: x(nc)
    real(16), intent(in ) :: y(nc)
    real(16), intent(out) :: A(np,nc)

    integer i

    do i = 1, nc
      call calc_poly_integral_coefs(np, poly, x(i) - 0.5_16, x(i) + 0.5_16, y(i) - 0.5_16, y(i) + 0.5_16, A(:,i))
    end do

  end subroutine calc_poly_integral_matrix_r16

end module poly_utils_mod
