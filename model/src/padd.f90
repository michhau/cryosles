!!-----------------------------------------------------------------------------------
!!  OBJECT: subroutine padd (u_big, u)
!!-----------------------------------------------------------------------------------
!!
!!  LAST MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 30/09/2011
!!
!!  DESCRIPTION:
!!
!!  This subroutine puts arrays into larger, zero-padded arrays.
!!  It automatically zeroes the oddballs.
!!  The u_big is filled with the elements of u till nx/2 and ny/2 (oddballs),
!!  then they put all the rest of the elements at the end (in y sense) of the u_big file,
!!  many zeroes are therefore left in the center after the nyquist y position.
!!  
!!  This subroutine can be parallelized on the z externally
!!
!!-----------------------------------------------------------------------------------
SUBROUTINE PADD (U_BIG,U)

USE TYPES, ONLY: RPREC
USE PARAM, ONLY: LH,LH_BIG,NX,NY,NY2

IMPLICIT NONE

INTEGER :: JX,JY

!note we're dealing with 2D COMPLEX arrays
COMPLEX(RPREC),DIMENSION(LH,NY),INTENT(IN) :: U
COMPLEX(RPREC),DIMENSION(LH_BIG,NY2),INTENT(OUT) :: U_BIG

!make sure the big array is zeroed
U_BIG=0.0_RPREC

!skip the Nyquist frequency since it should be zero anyway
do jy=1,ny/2
  do jx=1,nx/2
    u_big(jx,jy)=u(jx,jy)
  end do
end do

!fills the remaining elements to obtain a [lh_big x ny2] matrix
do jy=1,ny/2-1
  do jx=1,nx/2
    !goes at the end y: ny2, substracts (ny/2-1) which are the remaining, and starts adding them from that point on
    u_big(jx,jy+ny2-ny/2+1)=u(jx,jy+ny/2+1)
  end do
end do

END SUBROUTINE PADD
