MODULE SGS_FLUXES

USE PARAM
USE TYPES,ONLY:RPREC
USE SERVICE
use stretch

!REAL(RPREC),DIMENSION(LD,NY,0:NZ) :: T_LAG

!--------------------------------------------------------------------------------
CONTAINS   


!!-------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE SGS_FLUX ()
!!-------------------------------------------------------------------------------
!!
!!  LAST MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 16/05/2013
!!
!!  DESCRIPTION:
!!
!!
!!-------------------------------------------------------------------------------
SUBROUTINE SGS_FLUX(THETA,DTDX,DTDY,DTDZ,CS,NU_T,NU_MLC,PR_SGS,DS,NU_S,TXTH,TYTH,TZTH,beta_S)

USE SGSMODULE,ONLY:U_LAG,V_LAG,W_LAG
USE SIM_PARAM,ONLY:FILTER_SIZE,DT
USE LASD 

IMPLICIT NONE

REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(IN) :: THETA,DTDX,DTDY,DTDZ,NU_T
REAL(RPREC),DIMENSION(LD,NY,NZ),  INTENT(IN) :: CS
REAL(RPREC),INTENT(IN) :: NU_MLC,PR_SGS
REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(OUT) :: NU_S,BETA_S
REAL(RPREC),DIMENSION(LD,NY,NZ),INTENT(OUT) :: DS
REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(OUT) :: TXTH,TYTH,TZTH
REAL(RPREC),DIMENSION(LD,NY,NZ) :: S11,S12,S13,S23,S22,S33,IDTDX,IDTDY,IDTDZ
REAL(RPREC),DIMENSION(NZ) :: DS_TMP,L
INTEGER :: JX,JY,JZ,JZ_MIN
REAL(RPREC) :: STMP
REAL(RPREC),SAVE :: LAGR_DT
real(rprec) :: ZZ,ZZ_LOCAL

IF(VERBOSE) PRINT*,"STARTED SGS_FLUX"

DO JZ=1,NZ
            IF((COORD==0).and.(JZ == 1))THEN
               ZZ=(JZ-0.5_RPREC)
               ZZ_LOCAL = get_z_local(int(2*zz)) !c_stretch*(exp(zz)-1.0_rprec)
             ELSE  ! W-NODES
               ZZ=(JZ-1)+ETA_DISPL(COORD)
               ZZ_LOCAL = get_z_local(int(2*zz)+1) - get_z_local(int(2*zz)-1) !c_stretch*(exp(zz)-1.0_rprec)
            ENDIF
            L(JZ)=FILTER_SIZE*(DX*DY*(ZZ_LOCAL))**(1._RPREC/3._RPREC)
ENDDO

S11=0.0_rprec
S12=0.0_rprec
S13=0.0_rprec
S23=0.0_rprec
S22=0.0_rprec
S33=0.0_rprec

!###########################################################################
!FILTER SIZE DEFINES MIXING LENGTH

!L=FILTER_SIZE*(DX*DY*DZ)**(1._RPREC/3._RPREC)       !NONDIMENSIONAL

!AVERAGE VEL_LAG/THETA_LAG DURING CS
IF((MODEL.EQ.7).AND.(JT_TOTAL.GE.DYN_INIT))THEN

   !INCREASE LAGR_DT
   LAGR_DT=LAGR_DT+DT

   !AT DYN_INIT DON'T HAVE YET CS_COUNT VALUES
   IF(JT_TOTAL==DYN_INIT)THEN
      LAGR_DT=LAGR_DT*CS_COUNT
   ENDIF

ENDIF

!#################################################################
!SGS COEFFICIENT (TO CHANGE IF MAKE DYNAMIC)
IF((JT.EQ. 1).OR.(MODEL==1).OR.(MODEL==2).OR.(MODEL==5))THEN 

   DS=CS/PR_SGS                  !fixed turbulent Prandtl number
   NU_S=NU_T/PR_SGS

ELSEIF((JT_TOTAL==DYN_INIT) .OR. ((JT_TOTAL.GT.DYN_INIT).AND.(MOD(JT_TOTAL,CS_COUNT).EQ.0))) THEN

   !SIJ ON UVP NODE AT JZ=1, ON W NODES FOR THE OTHERS (JZ>1)
   CALL CALC_SIJ(DTDX,DTDY,DTDZ,S11,S12,S13,S23,S22,S33,IDTDX,IDTDY,IDTDZ)
   IF(MODEL==6)THEN
      CALL STD_DYN_SC(THETA,IDTDX,IDTDY,IDTDZ,S11,S12,S13,S22,S23,S33,DS_TMP)
      DO JZ=1,NZ
         DS(:,:,JZ)=DS_TMP(JZ)
      ENDDO
   ELSEIF(MODEL==7)THEN
      call lagrange_Sdep_Scalar(S11,S12,S13,S22,S23,S33,theta,idTdx,idTdy,idTdz,Ds,beta_S)
      LAGR_DT=0.0_RPREC
      do jz = 0,nz
         u_lag(:,:,jz) = 0._rprec
         v_lag(:,:,jz) = 0._rprec
         w_lag(:,:,jz) = 0._rprec
      enddo
   ENDIF

   !EDDY VISCOSITY
   DO JZ=1,NZ
   DO JY=1,NY
   DO JX=1,NX
      STMP = SQRT(2._RPREC*(S11(JX,JY,JZ)**2 + S22(JX,JY,JZ)**2 +                &
                            S33(JX,JY,JZ)**2 + 2._RPREC*(S12(JX,JY,JZ)**2 +      &
                            S13(JX,JY,JZ)**2 + S23(JX,JY,JZ)**2)))
      NU_S(JX,JY,JZ)=STMP*DS(JX,JY,JZ)*L(JZ)**2
   ENDDO
   ENDDO
   ENDDO

ENDIF

!#################################################################
!SGS FLUXES

!JZ=1 SPECIAL TREATMENT (TZTH COMP. IN WALL_LAW ALTERNATIVELY)
IF(COORD==0)THEN
   TXTH(:,:,1) = (NU_S(:,:,1)+NU_MLC)*DTDX(:,:,1)        !UVP
   TYTH(:,:,1) = (NU_S(:,:,1)+NU_MLC)*DTDY(:,:,1)        !UVP
   IF(LLBC_SPECIAL.NE.'WALL_LAW')THEN
      TZTH(:,:,1)=NU_MLC*DTDZ(:,:,1) 
   ENDIF
   JZ_MIN=2
ELSE
   JZ_MIN=1
ENDIF

!$OMP PARALLEL DO
DO JZ=JZ_MIN,NZ-1 
   TXTH(:,:,JZ) = 0.5_RPREC*(NU_S(:,:,JZ)+NU_S(:,:,JZ+1)+2.0_RPREC*NU_MLC)*DTDX(:,:,JZ) !UVP
   TYTH(:,:,JZ) = 0.5_RPREC*(NU_S(:,:,JZ)+NU_S(:,:,JZ+1)+2.0_RPREC*NU_MLC)*DTDY(:,:,JZ) !UVP      
   TZTH(:,:,JZ) = (NU_S(:,:,JZ)+NU_MLC)*DTDZ(:,:,JZ)                            !W
ENDDO
!$OMP END PARALLEL DO

!MISSING TERM FOR W
IF((.NOT.USE_MPI).OR.(USE_MPI.AND.COORD==NPROC-1))THEN
   TZTH(:,:,NZ) = (NU_S(:,:,NZ)+NU_MLC)*DTDZ(:,:,NZ)  !W
ENDIF

!SEND RECEIVE
#IF(MPI)
   !--send scalar info down & recv scalar info from above
   call mpi_sendrecv (TZTH(1, 1, 1), ld*ny, MPI_RPREC, down, 26,  &
                      TZTH(1, 1, nz), ld*ny, MPI_RPREC, up, 26,   &
                      comm, status, ierr)
#ENDIF



IF(VERBOSE) PRINT*,"FINISHED SGS_FLUX"

END SUBROUTINE SGS_FLUX


!!------------------------------------------------------------------------------
!!  OBJECT: subroutine CALC_SIJ()
!!------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 27/03/2013
!!
!!  DESCRIPTION:
!!
!!  COMPUTES SIJ ON UVP NODE AT JZ=1, ON W NODES FOR THE OTHERS (JZ>1)
!!
!!------------------------------------------------------------------------------
SUBROUTINE CALC_SIJ(DTDX,DTDY,DTDZ,S11,S12,S13,S23,S22,S33,IDTDX,IDTDY,IDTDZ)

USE SIM_PARAM,ONLY:DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ 

IMPLICIT NONE

REAL(RPREC) :: UX,UY,UZ,VX,VY,VZ,WX,WY,WZ
REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(IN) :: DTDX,DTDY,DTDZ
REAL(RPREC),DIMENSION(LD,NY,NZ),INTENT(OUT) :: S11,S12,S13,S23,S22,S33,IDTDX,IDTDY,IDTDZ
INTEGER :: JX,JY,JZ,JZ_MIN
integer :: kk
real(rprec) :: fac_1,fac_2,fac_3

IF(COORD==0)THEN
   !FIRST LEVEL HAS EVERYTHING ON UVP NODES
   DO JY=1,NY
   DO JX=1,NX     

      UX=DUDX(JX,JY,1)       !UVP-NODE @ DZ/2
      UY=DUDY(JX,JY,1)       !UVP-NODE @ DZ/2
      IF(LLBC_SPECIAL.EQ.'WALL_LAW')THEN
         UZ=DUDZ(JX,JY,1)
         VZ=DVDZ(JX,JY,1)
      ELSE
         UZ=0.5_RPREC*(DUDZ(JX,JY,1)+DUDZ(JX,JY,2))
         VZ=0.5_RPREC*(DVDZ(JX,JY,1)+DVDZ(JX,JY,2))
      ENDIF
      VX=DVDX(JX,JY,1)                            !UVP-NODE @ DZ/2
      VY=DVDY(JX,JY,1)                            !UVP-NODE @ DZ/2

      fac_1 = get_z_local(2)
      fac_2 = get_z_local(0)
      fac_3 = get_z_local(1)

      ! interpolate from w to uvp nodes
      wx = dwdx(jx,jy,1) + ((dwdx(jx,jy,2)-dwdx(jx,jy,1))/(fac_1-fac_2))*(fac_3-fac_2)
      wy = dwdy(jx,jy,1) + ((dwdy(jx,jy,2)-dwdy(jx,jy,1))/(fac_1-fac_2))*(fac_3-fac_2)
      !WX=0.5_RPREC*(DWDX(JX,JY,1)+DWDX(JX,JY,2))  !UVP-NODE SO I NEED TO INTERP
      !WY=0.5_RPREC*(DWDY(JX,JY,1)+DWDY(JX,JY,2))  !UVP-NODE SO I NEED TO INTERP
      
      WZ=DWDZ(JX,JY,1)                            !UVP-NODE SO I NEED TO INTERP

      S11(JX,JY,1)=UX                             !UVP-NODE @ DZ/2
      S12(JX,JY,1)=0.5_RPREC*(UY+VX)              !UVP-NODE @ DZ/2
      S13(JX,JY,1)=0.5_RPREC*(UZ+WX)              !UVP-NODE @ DZ/2
      S22(JX,JY,1)=VY                             !UVP-NODE @ DZ/2
      S23(JX,JY,1)=0.5_RPREC*(VZ+WY)              !UVP-NODE @ DZ/2
      S33(JX,JY,1)=WZ                             !UVP-NODE @ DZ/2

      IDTDX(JX,JY,1) = DTDX(JX,JY,1)  ! UVP-NODE
      IDTDY(JX,JY,1) = DTDY(JX,JY,1)  ! UVP-NODE
      IDTDZ(JX,JY,1) = DTDZ(JX,JY,1)  ! UVP-NODE (CAREFUL IF NOT USING LAW OF THE  WALL!!)
   ENDDO
   ENDDO
  
   JZ_MIN=2

ELSE
   JZ_MIN=1
ENDIF


!necessary cause we pass just u,v,w and not dwdz,
#IF(MPI)
  call mpi_sendrecv (dwdz(1, 1, 1), ld*ny, MPI_RPREC, down, 25,  &
                     dwdz(1, 1, nz), ld*ny, MPI_RPREC, up, 25,   &
                     comm, status, ierr)
#ENDIF


!STANDARD NODES (PUT EVERYTHING ON W, WHERE CS AND NU_T ARE DEFINED)
DO JZ=JZ_MIN,NZ
DO JY=1,NY
DO JX=1,NX

   kk=jz + eta_displ(coord)

   fac_1=get_z_local(2*kk-1)
   fac_2=get_z_local(2*kk-3)
   fac_3=get_z_local(2*kk-2)

   ux=dudx(jx,jy,jz-1) + ((dudx(jx,jy,jz)-dudx(jx,jy,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)
   uy=dudy(jx,jy,jz-1) + ((dudy(jx,jy,jz)-dudy(jx,jy,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)
   vx=dvdx(jx,jy,jz-1) + ((dvdx(jx,jy,jz)-dvdx(jx,jy,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)
   vy=dvdy(jx,jy,jz-1) + ((dvdy(jx,jy,jz)-dvdy(jx,jy,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)

   wz=dwdz(jx,jy,jz-1) + ((dwdz(jx,jy,jz)-dwdz(jx,jy,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)

   !UX=0.5_RPREC*(DUDX(JX,JY,JZ) + DUDX(JX,JY,JZ-1))     !W-NODE
   !UY=0.5_RPREC*(DUDY(JX,JY,JZ) + DUDY(JX,JY,JZ-1))     !W-NODE
   UZ=DUDZ(JX,JY,JZ)                                    !W-NODE
   !VX=0.5_RPREC*(DVDX(JX,JY,JZ) + DVDX(JX,JY,JZ-1))     !W-NODE
   !VY=0.5_RPREC*(DVDY(JX,JY,JZ) + DVDY(JX,JY,JZ-1))     !W-NODE
   VZ=DVDZ(JX,JY,JZ)                                    !W-NODE
   WX=DWDX(JX,JY,JZ)                                    !W-NODE
   WY=DWDY(JX,JY,JZ)                                    !W-NODE
   !WZ=0.5_RPREC*(DWDZ(JX,JY,JZ) + DWDZ(JX,JY,JZ-1))     !W-NODE

   S11(JX,JY,JZ)=UX                                     !W-NODE
   S12(JX,JY,JZ)=0.5_RPREC*(UY+VX)                      !W-NODE
   S13(JX,JY,JZ)=0.5_RPREC*(UZ+WX)                      !W-NODE
   S22(JX,JY,JZ)=VY                                     !W-NODE
   S23(JX,JY,JZ)=0.5_RPREC*(VZ+WY)                      !W-NODE
   S33(JX,JY,JZ)=WZ                                     !W-NODE
    
   idtdx(jx,jy,jz) = dtdx(jx,jy,jz-1) + ((dtdx(jx,jy,jz)-dtdx(jx,jy,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)
   idtdy(jx,jy,jz) = dtdy(jx,jy,jz-1) + ((dtdy(jx,jy,jz)-dtdy(jx,jy,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)

   !IDTDX(JX,JY,JZ) = 0.5_RPREC*(DTDX(JX,JY,JZ) + DTDX(JX,JY,JZ-1)) !W-NODE
   !IDTDY(JX,JY,JZ) = 0.5_RPREC*(DTDY(JX,JY,JZ) + DTDY(JX,JY,JZ-1)) !W-NODE
   IDTDZ(JX,JY,JZ) = DTDZ(JX,JY,JZ) !W-NODE

ENDDO
ENDDO
ENDDO


END SUBROUTINE CALC_SIJ

END MODULE SGS_FLUXES
