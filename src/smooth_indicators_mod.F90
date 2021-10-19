module smooth_indicators_mod

  implicit none

  interface
    pure real(8) function smooth_indicator_interface(a)
      real(8), intent(in) :: a(:)
    end function
  end interface

contains

  pure real(8) function smooth_indicator_3(a) result(res)

    real(8), intent(in) :: a(:) ! 3

    res = a(2) * a(2) + 13 * a(3) * a(3) / 3

  end function smooth_indicator_3

  pure real(8) function smooth_indicator_3x3(a) result(res)

    real(8), intent(in) :: a(:) ! 9

    res = (  720 * a(2) * a(2) + 3120 * a(3) * a(3) + 720  * a(4) * a(4) + 840   * a(5) * a(5) &
          +  120 * a(4) * a(6) + 3389 * a(6) * a(6) + 3120 * a(7) * a(7) + 120   * a(2) * a(8) &
          + 3389 * a(8) * a(8) + 520  * a(3) * a(9) + 520  * a(7) * a(9) + 13598 * a(9) * a(9) ) / 720

  end function smooth_indicator_3x3

end module smooth_indicators_mod
