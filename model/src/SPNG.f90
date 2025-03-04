MODULE SPNG

USE PARAM
USE TYPES, ONLY: RPREC
USE stretch
USE SIM_PARAM, ONLY: SPONGE

IMPLICIT NONE

!--------------------------------------------------------------------------------
CONTAINS


!!--------------------------------------------------------------------------------
!!  OBJECT: INIT_SPONGE()
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 25/07/2012
!!
!!  DESCRIPTION:
!!
!!  This subroutine computes the damp coefficient at each JZ of the domain
!!
!!--------------------------------------------------------------------------------
SUBROUTINE INIT_SPONGE()

USE PARAM
IMPLICIT NONE

REAL(RPREC) ::          Z_D,CFRDMP,Z_LOCAL,FACTOR,SPONGE_top
INTEGER ::              JZ,JZ_GLOBAL
integer :: kk

! Initial value
SPONGE = 0.0_rprec

if (SPONGE_TYPE==2) then

  Z_D = (1.0_rprec-SPONGE_H)*L_Z
  CFRDMP = 3.9_rprec
  SPONGE=0._rprec
  DO JZ=1,NZ
     !Z_LOCAL = (JZ-1._rprec)*DZ  !JF corrected to W node
!#IF(MPI)
     !Z_LOCAL = Z_LOCAL + ETA_DISPL(coord)*DZ
!#ENDIF
     kk = jz + eta_displ(coord)
     Z_LOCAL =  get_z_local(2*kk) !c_stretch*(exp(kk*dz)-1.0_rprec) !sum(dz(1:kk-1))
     IF(Z_LOCAL .ge. Z_D .and. Z_LOCAL .le. L_z)THEN
         SPONGE(JZ)=CFRDMP*0.5_rprec*(1.0_rprec-COS(PI*(Z_LOCAL-Z_D)/(L_Z-Z_D)))
     ELSE
         SPONGE(JZ)=0.0_rprec
     ENDIF
  ENDDO

elseif (SPONGE_TYPE==1) then

! sets relaxation term to vertical momentum equation in top quarter
! of domain, relaxation timescale on top is 50s with a FACTOR of 5 if we
! had Nz=40 for the layers 40...31, Nieuwstadt et al. 1991, turbulent shear 
! flows

#IF(MPI)
  !--the non-MPI recursive form is inconvenient, so instead replace with
  !  analytically evaluated version (should be this way anyway...)
  SPONGE=0._rprec
  FACTOR=9._rprec/(NZ_TOT - 3*(NZ_TOT-1)/4 - 1)  !JF made corrections
  SPONGE_top = L_Z / (50._rprec) !*U_STAR)
  DO JZ=1,NZ
    JZ_GLOBAL = JZ + eta_displ(coord)
    IF(JZ_GLOBAL .ge. 3*(NZ_TOT-1)/4 + 1)THEN  !JF made corrections
      SPONGE(JZ) = SPONGE_top * 5._rprec**((JZ_GLOBAL-NZ_TOT) * FACTOR)
    ENDIF
  ENDDO
#ELSE
  SPONGE=0._rprec
  FACTOR=9._rprec/(NZ-3*(NZ-1)/4-1)  !JF made corrections
  SPONGE(NZ)=L_Z/(50._rprec) !*U_STAR)
  do JZ=NZ-1,3*(NZ-1)/4+1,-1
     SPONGE(JZ)=SPONGE(JZ+1)/5._rprec**FACTOR
  end do
#ENDIF

end if

END SUBROUTINE INIT_SPONGE


!!--------------------------------------------------------------------------------
!!  OBJECT: SET_VSPONGE(U,V,W,RHSx,RHSy,RHSz,SPNGFLAG)
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 25/07/2012
!!
!!  DESCRIPTION:
!!
!!  This subroutine adds the sponge term to the RHS of the equations for the 
!!  velocity, it uses the SPONGE variable which must be initialized by calling 
!!  the SPONGE_INIT subroutine defined above.
!!
!!--------------------------------------------------------------------------------
SUBROUTINE SET_VSPONGE(U,V,W,RHSx,RHSy,RHSz)

IMPLICIT NONE

REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(IN) ::       U,V,W
REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(INOUT) ::    RHSx,RHSy,RHSz
INTEGER ::                                            JZ

   !ADD SPONGE LAYER (UVP for RHSx and RHSy need interpolation)
   IF(SPONGE_TYPE==2)THEN
      !$OMP PARALLEL DO
      DO JZ=1,NZ-1 
         RHSx(1:NX,1:NY,JZ) = RHSx(1:NX,1:NY,JZ)                        &
                            - 0.5_rprec*(SPONGE(JZ)+SPONGE(JZ+1))*      &
                              (U(1:NX,1:NY,JZ)-SUM(U(1:NX,1:NY,JZ))/(NX*NY))
         RHSy(1:NX,1:NY,JZ) = RHSy(1:NX,1:NY,JZ)                        &
                            - 0.5_rprec*(SPONGE(JZ)+SPONGE(JZ+1))*      &
                              (V(1:NX,1:NY,JZ)-SUM(V(1:NX,1:NY,JZ))/(NX*NY))
         RHSz(1:NX,1:NY,JZ) = RHSz(1:NX,1:NY,JZ)                        &
                            - 1.0_rprec*SPONGE(JZ)*                     &
                              (W(1:NX,1:NY,JZ)-SUM(W(1:NX,1:NY,JZ))/(NX*NY))
      ENDDO
      !$OMP END PARALLEL DO
   ELSEIF(SPONGE_TYPE==1)THEN
      !$OMP PARALLEL DO
      DO JZ=1,NZ-1 
         RHSz(1:NX,1:NY,JZ) = RHSz(1:NX,1:NY,JZ)-SPONGE(JZ)*W(1:NX,1:NY,JZ)
      ENDDO
      !$OMP END PARALLEL DO
   ENDIF
   
END SUBROUTINE SET_VSPONGE


!!--------------------------------------------------------------------------------
!!  OBJECT: SET_SCSPONGE(SC,RHSsc)
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 25/07/2012
!!
!!  DESCRIPTION:
!!
!!  This subroutine adds the sponge term to the RHS of the equations for the 
!!  scalar field, it uses the SPONGE variable which must be initialized by calling 
!!  the SPONGE_INIT subroutine defined above.
!!
!!--------------------------------------------------------------------------------
SUBROUTINE SET_SCSPONGE(SC,RHSsc)

USE PARAM
IMPLICIT NONE

REAL(RPREC),DIMENSION(LD,NY,NZ),INTENT(IN) ::       SC
REAL(RPREC),DIMENSION(LD,NY,NZ),INTENT(INOUT) ::    RHSsc
INTEGER ::                                          JZ

   !ADD SPONGE LAYER (UVP for RHSx and RHSy need interpolation)
   DO JZ=1,NZ-1 
      RHSsc(1:NX,1:NY,JZ) = RHSsc(1:NX,1:NY,JZ)                        &
                           - 0.5_rprec*(SPONGE(JZ)+SPONGE(JZ+1))*      &
                           (SC(1:NX,1:NY,JZ)-SUM(SC(1:NX,1:NY,JZ))/(NX*NY))
   ENDDO
   
END SUBROUTINE SET_SCSPONGE


END MODULE SPNG

