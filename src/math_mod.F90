module math_mod

  implicit none

  private

  public inverse_matrix

  interface inverse_matrix
    module procedure inverse_matrix_r8
  end interface inverse_matrix

contains

  subroutine inverse_matrix_r8(A, B, ierr)

    real(8), intent(in ) :: A(:,:)
    real(8), intent(out) :: B(size(A,1),size(A,1))
    integer, intent(out) :: ierr

    real(8) tmp
    integer n, i, j, k
    integer imax(size(A,1)), jmax(size(A,1))

    if (size(A, 1) /= size(A, 2)) then
      ierr = 1
      return
    end if

    n = size(A, 1)

    B = A

    do k = 1, n
      tmp = 0
      do j = k, n
        do i = k, n
          if (abs(A(i,j)) > tmp) then
            tmp = abs(A(i,j)); imax(k) = i; jmax(k) = j
          end if
        end do
      end do
      if (tmp == 0) then
        ierr = 2
        return
      end if
      do j = 1, n
        tmp = B(k,j); B(k,j) = B(imax(k),j); B(imax(k),j) = tmp
      end do
      do i = 1, n
        tmp = B(i,k); B(i,k) = B(i,jmax(k)); B(i,jmax(k)) = tmp
      end do
      B(k,k) = 1.0d0 / B(k,k)
      do j = 1, n
        if (j /= k) then
          B(k,j) = B(k,j) * B(k,k)
        end if
      end do
      do i = 1, n
        if (i /= k) then
          do j = 1, n
            if (j /= k) then
              B(i,j) = B(i,j) - B(i,k) * B(k,j)
            end if
          end do
        end if
      end do
      do i = 1, n
        if (i /= k) then
          B(i,k) = -B(i,k) * B(k,k)
        end if
      end do
    end do

    do k = n, 1, -1
      do j = 1, n
        tmp = B(k,j); B(k,j) = B(jmax(k),j); B(jmax(k),j) = tmp
      end do
      do i = 1, n
        tmp = B(i,k); B(i,k) = B(i,imax(k)); B(i,imax(k)) = tmp
      end do
    end do

  end subroutine inverse_matrix_r8

end module math_mod
