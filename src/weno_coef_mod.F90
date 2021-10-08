module weno_coef_mod

  use poly_utils_mod
  use weno_types_mod

  implicit none

  private

  public calc_weno_tensor_product_ideal_coefs

contains

  subroutine calc_weno_tensor_product_ideal_coefs(weno, ierr)

    type(weno_tensor_product_type), intent(inout) :: weno
    integer, intent(in) :: ierr

    integer i, j

  end subroutine calc_weno_tensor_product_ideal_coefs

end module weno_coef_mod
