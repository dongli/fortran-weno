module weno_interp_mod

  implicit none

  private

	public weno_interp_5th_order

contains

  subroutine weno_interp_5th_order(x, xi)

		! Stencils:
		!
		!                     |-----------S3----------|
		!             |-----------S2----------|       |
		!     |-----------S1----------|       |       |
		!  ...|---o---|---o---|---o---x---o---|---o---|...
		!     |  i-2  |  i-1  |   i   |  i+1  |  i+2  |
		!                             |
		!                           i+1/2

		real(8), intent(in ) :: x(-1:3)
		real(8), intent(out) :: xi

    ! Polynomial coefficients
		real(8), parameter :: p1(3) = [ 1.0d0 / 3.0d0, -7.0d0 / 6.0d0, 11.0d0 / 6.0d0]
    real(8), parameter :: p2(3) = [-1.0d0 / 6.0d0,  5.0d0 / 6.0d0,  1.0d0 / 3.0d0]
    real(8), parameter :: p3(3) = [ 1.0d0 / 3.0d0,  5.0d0 / 6.0d0, -1.0d0 / 6.0d0]
    ! Stencils weights (combine 3rd order stencils into 5th order big stencil)
    real(8), parameter :: g (3) = [ 1.0d0 / 10.0d0, 3.0d0 / 5.0d0, 3.0d0 / 10.0d0]
    ! Smoothness indicator coefficients
    real(8), parameter :: s1 = 13.0d0 / 12.0d0
    real(8), parameter :: s2 = 0.25d0
    real(8), parameter :: s11(3) = [ 1.0d0, -2.0d0,  1.0d0]
    real(8), parameter :: s12(3) = [ 1.0d0, -4.0d0,  3.0d0]
    real(8), parameter :: s21(3) = [ 1.0d0, -2.0d0,  1.0d0]
    real(8), parameter :: s22(3) = [ 1.0d0,  0.0d0, -1.0d0]
    real(8), parameter :: s31(3) = [ 1.0d0, -2.0d0,  1.0d0]
    real(8), parameter :: s32(3) = [ 3.0d0, -4.0d0,  1.0d0]
    real(8), parameter :: eps = 1.0d-40

		real(8) b (3) ! Smooth indicators
    real(8) w (3) ! Combination weights
    real(8) xs(3) ! Interpolation values on each stencil

		! Calculate value at interfaces for each stencil.
		xs(1) = dot_product(p1, x(-1:1))
    xs(2) = dot_product(p2, x( 0:2))
    xs(3) = dot_product(p3, x( 1:3))
    ! Calculate smooth indicators for each stencil.
    b(1) = s1 * dot_product(s11, x(-1:1))**2 + s2 * dot_product(s12, x(-1:1))**2
    b(2) = s1 * dot_product(s21, x( 0:2))**2 + s2 * dot_product(s22, x( 0:2))**2
    b(3) = s1 * dot_product(s31, x( 1:3))**2 + s2 * dot_product(s32, x( 1:3))**2
		! Calculate stencil linear combination weights considering smooth indicators.
		w = g / (eps + b)**2; w = w / sum(w)
		xi = dot_product(w, xs)

  end subroutine weno_interp_5th_order

end module weno_interp_mod
