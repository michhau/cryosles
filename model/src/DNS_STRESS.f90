!!--------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE DNS_MFLX
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 26/03/2013
!!
!!--------------------------------------------------------------------------------
SUBROUTINE DNS_STRESS (PARAM)

USE TYPES,ONLY:RPREC
USE SIM_PARAM,ONLY: DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ, &
                    TXX,TYY,TZZ,TXY,TXZ,TYZ

IMPLICIT NONE

REAL(RPREC),INTENT(IN) :: PARAM

!COMPUTE STRESSES ON NODES     (DO WE NEED THE NZ LEVEL?)
TXX=2.0_RPREC*DUDX*PARAM
TYY=2.0_RPREC*DVDY*PARAM
TZZ=2.0_RPREC*DWDZ*PARAM
TXY=(DUDY+DVDX)*PARAM             !UVP       -2.0_RPREC*NU*0.5_RPREC*(DUDY+DVDX)
TXZ=(DUDZ+DWDX)*PARAM             !W         -2.0_RPREC*NU*0.5_RPREC*(DUDZ+DWDX)
TYZ=(DVDZ+DWDY)*PARAM             !W         -2.0_RPREC*NU*0.5_RPREC*(DVDZ+DWDY)


END SUBROUTINE DNS_STRESS
