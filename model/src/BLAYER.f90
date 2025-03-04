!!--------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE BLAYER(U,V,W,L_X,LENGTH)
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 01/02/2013
!!
!!  DESCRIPTION
!!
!!  This subroutine applies the buffer layer at the end of the computational domain
!!
!!--------------------------------------------------------------------------------
SUBROUTINE BLAYER(U,V,W,LENGTH)

USE PARAM
USE TYPES, ONLY:RPREC

IMPLICIT NONE

REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(OUT) :: U,V,W
REAL(RPREC),INTENT(IN) :: LENGTH
INTEGER :: I,J,K
INTEGER :: BNX,JXI,JXE
REAL(RPREC) :: FACTOR

!LENGTH OF BL + GRID NODES IN BL
BNX=(L_X/LENGTH)/DX+1
JXI=NX+1-BNX
JXE=NX+1

!SKIP FIRST AND LAST SINCE WE KNOW ALREADY
DO I=JXI+1,JXE-1

   FACTOR=0.5_RPREC*(1.0_RPREC-COS(PI*(I-JXI)/REAL(JXE-JXI,RPREC)) )
   U(I,1:NY,0:NZ) = U(JXI,1:NY,0:NZ)+FACTOR*(U(1,1:NY,0:NZ)-U(JXI,1:NY,0:NZ))
   V(I,1:NY,0:NZ) = V(JXI,1:NY,0:NZ)+FACTOR*(V(1,1:NY,0:NZ)-V(JXI,1:NY,0:NZ))
   W(I,1:NY,0:NZ) = W(JXI,1:NY,0:NZ)+FACTOR*(W(1,1:NY,0:NZ)-W(JXI,1:NY,0:NZ))

ENDDO

END SUBROUTINE BLAYER
