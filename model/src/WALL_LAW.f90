MODULE MWALL_LAW

USE SERVICE
USE TYPES,ONLY:RPREC
USE PARAM
USE SIM_PARAM
use stretch
USE MO_SIMILARITY, ONLY:U1,V1,U_AVG
!use ls_model
IMPLICIT NONE

!-------------------------------------------------------------------------------
CONTAINS

!!------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE WALL_LAW()
!!------------------------------------------------------------------------------
!!
!!  LAST MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 17/10/2012
!!
!!  DESCRIPTION:
!!
!!  calls the proper subroutine to compute stresses and dudz,dvdz.
!!
!!  Computes DUDZ(:,:,1),DVDZ(:,:,1),TXZ(:,:,1),TYZ(:,:,1). 
!!  Out of this subroutine the values are put at W nodes so they're
!!  shifted by DZ/2 compared to the position where U was read.
!!
!!------------------------------------------------------------------------------
SUBROUTINE WALL_LAW()

IMPLICIT NONE

INTEGER :: JX,JY
real(rprec) :: z_local
!real(rprec) :: taupx,taupy
IF(VERBOSE) PRINT*,'IN WALL_LAW'

!TRY VERSION
DO JY=1,NY
DO JX=1,NX
   !CHECK FOR NAN VALUES IN CASE OF 0 VELOCITY
   IF(U_AVG(JX,JY).EQ.0.0_RPREC)THEN
      TXZ(JX,JY,1)=0.0_RPREC
      TYZ(JX,JY,1)=0.0_RPREC
      DUDZ(JX,JY,1)=0.0_RPREC
      DVDZ(JX,JY,1)=0.0_RPREC
   ELSE

      z_local = get_z_local(1) !c_stretch*(exp(0.5_rprec*dz)-1.0_rprec)

      DUDZ(JX,JY,1)=PHI_M(JX,JY)*USTAR_AVG(JX,JY)/(z_local*VONK)   &
                                *U1(JX,JY)/U_AVG(JX,JY)
      DVDZ(JX,JY,1)=PHI_M(JX,JY)*USTAR_AVG(JX,JY)/(z_local*VONK)   &
                                *V1(JX,JY)/U_AVG(JX,JY)

      !!! ADDIDING A CORRECTION TO ACCOUNT FOR THE PARTICLE STRESS !!!
!      if(lsm) then
!         taupx = tpxz(jx,jy)
!         taupy = tpyz(jx,jy)
!      else
!         taupx = 0.0_rprec
!         taupy = 0.0_rprec
!      endif
!DBM: Particle stress close to the surface is already taken into account in the source terms Fx,Fy,Fz
      TXZ(JX,JY,1)=USTAR_AVG(JX,JY)**2/U_AVG(JX,JY)*U1(JX,JY) !- taupx
      TYZ(JX,JY,1)=USTAR_AVG(JX,JY)**2/U_AVG(JX,JY)*V1(JX,JY) !- taupy
   ENDIF
ENDDO
ENDDO

IF(VERBOSE) PRINT*,'OUT WALL_LAW'

END SUBROUTINE WALL_LAW

END MODULE MWALL_LAW
