subroutine tridag_array (a, b, c, r, u)
! sc: modified for complex!
use types, only: rprec
use param
implicit none

real(rprec), dimension(lh,ny,nz+1), intent(in) :: a, b, c
complex(rprec), dimension(lh,ny,nz+1), intent(in) :: r
complex(rprec), dimension(lh,ny,nz+1), intent(out):: u

integer, save :: n

logical, parameter :: DEBUG = .false.

character(64) :: fmt
integer :: jx, jy, j, j_min, j_max
real(rprec) :: bet(lh, ny)
real(rprec), dimension(lh,ny,nz+1) :: gam

! initialize n
n = nz + 1

!want to skip ny/2+1 and 1, 1

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then

  do jy = 1, ny
    do jx = 1, lh-1

      if (b(jx, jy, 1) == 0._rprec) then
        write (*, *) 'tridag_array: rewrite eqs, jx, jy= ', jx, jy
        stop
      end if

    end do
  end do

  bet = b(:, :, 1)
  u(:, :, 1) = r(:, :, 1) / bet

  j_min = 1  !--this is only for backward pass
else
  j_min = 2  !--this is only for backward pass
end if

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  j_max = n
else
  j_max = n-1
end if

do j = 2, j_max

  do jy = 1, ny

    if (jy == ny/2+1) cycle

    do jx = 1, lh-1

!     if (jx*jy == 1) cycle

      gam(jx, jy, j) = c(jx, jy, j-1) / bet(jx, jy)
      bet(jx, jy) = b(jx, jy, j) - a(jx, jy, j)*gam(jx, jy, j)

      if (bet(jx, jy) == 0._rprec) then
        write (*, *) 'tridag_array failed at jx,jy,j=', jx, jy, j
        write (*, *) 'a,b,c,gam,bet=', a(jx, jy, j), b(jx, jy, j),  &
                     c(jx, jy, j), gam(jx, jy, j), bet(jx, jy)
        stop
      end if

      u(jx, jy, j) = (r(jx, jy, j) - a(jx, jy, j) * u(jx, jy, j-1)) /  &
                     bet(jx, jy)

    end do

  end do

end do

do j = n-1, j_min, -1

  do jy = 1, ny

    if (jy == ny/2+1) cycle

    do jx = 1, lh-1

!     if (jx*jy == 1) cycle

      u(jx, jy, j) = u(jx, jy, j) - gam(jx, jy, j+1) * u(jx, jy, j+1)

    end do

  end do

end do


end subroutine tridag_array
