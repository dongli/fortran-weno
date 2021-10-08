program test

  use poly_utils_mod

  implicit none

  integer, parameter :: m = 3
  real(8) A(m*m,m*m)
  integer i, j

  call calc_poly_tensor_product_integral_coef_matrix(m, m, A)

  do j = 1, m * m
    do i = 1, m * m
      write(*, '(F10.5)', advance='no') A(i,j)
    end do
    write(*, *)
  end do

end program test
