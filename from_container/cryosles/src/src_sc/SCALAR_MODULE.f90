MODULE SCALARS_MODULE

USE TYPES,ONLY:RPREC
USE PARAM
USE SIM_PARAM
USE DERIV
USE SPNG
USE SERVICE,ONLY:Z_MPROFILE
USE SGS_FLUXES
USE LS_MODEL
use stretch
IMPLICIT NONE

INTEGER, PARAMETER:: TAG_COUNTER = 200

!------------------------------------------------------------------------------
! MCB. Added a variable value for dTdz_top, so the convective BL, do not grow,
! dTdz_top equivalent to the imposed surface heat flux.


!--------------------------------------------------------------------------------
CONTAINS

!!-------------------------------------------------------------------------------
!!  OBJECT: subroutine SCALAR_SOLVER
!!-------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 17/05/2013
!!
!!  DESCRIPTION:
!!
!!  Builds the Rhs for the scalar equation and solves it one step in time
!!  adopting the actual velocity field (GAUSS-SEIDEL approach).
!!  
!!  VARIABLES POSITIONING:
!!
!!  U,V				         -->	UVP
!!  W				            -->	W
!!  THETA		         	-->	UVP
!!  TXTH,TYTH,DTDX,DTDY		-->	UVP
!!  DZTH,DTDZ			      -->	W
!!  BETA_SCAL			      -->	UVP
!!
!!------------------------------------------------------------------------------
SUBROUTINE SCALAR_SOLVER ( T,DTDX,DTDY,DTDZ,TXTH,TYTH,TZTH,       &
                           RHS_T,RHS_TF,                          &
                           LLBC_T,DLBC_T,LBC_T,                   &
                           LUBC_T,DUBC_T,UBC_T,                   &
                           PR_T,PR_T_SGS,         &
                           CS,NU_T,DS,NU_S,T_STAR,PSI_T,PHI_T,opt_part_forc)

USE TEST_FILTERMODULE
USE BC 

IMPLICIT NONE

REAL(RPREC),DIMENSION(LD,NY),INTENT(IN) :: LBC_T,UBC_T
REAL(RPREC),INTENT(IN) ::  PR_T,PR_T_SGS
CHARACTER(*),INTENT(IN) :: LLBC_T,LUBC_T
REAL(RPREC),INTENT(IN) ::  DLBC_T,DUBC_T
REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(INOUT) :: DTDX,DTDY,DTDZ,TXTH,TYTH,TZTH,     &
                                                   T,RHS_T,RHS_TF

REAL(RPREC),DIMENSION(LD,NY,NZ),OPTIONAL,INTENT(OUT) :: DS
REAL(RPREC),DIMENSION(LD,NY,0:NZ),OPTIONAL,INTENT(OUT) :: NU_S
REAL(RPREC),DIMENSION(LD,NY,NZ),OPTIONAL,INTENT(IN) :: CS
REAL(RPREC),DIMENSION(LD,NY,0:NZ),OPTIONAL,INTENT(IN) :: NU_T
REAL(RPREC),DIMENSION(LD,NY,0:NZ) :: particle_forcing
REAL(RPREC),DIMENSION(LD,NY,0:NZ),optional :: opt_part_forc
REAL(RPREC),DIMENSION(LD,NY),OPTIONAL,INTENT(IN) :: T_STAR
REAL(RPREC),DIMENSION(LD,NY),OPTIONAL,INTENT(IN) :: PSI_T
REAL(RPREC),DIMENSION(LD,NY),OPTIONAL,INTENT(IN) :: PHI_T

INTEGER :: I,J,K
REAL(RPREC) :: A,B,C,P1,P2

INTEGER :: JZ,JZ_GLOBAL

IF(VERBOSE) PRINT*,'IN THETA_ALL_IN_ONE'

if(present(opt_part_forc)) then
    particle_forcing=opt_part_forc
else
    particle_forcing=0.0_rprec
endif
!###################################################################
!MAIN LOOP

!normalizatin factors
IF(SLOPE.and.sf_norm_1)THEN
    P1 = 1.0_RPREC/(RE*(DSIN(ALPHA*PI/180.0_RPREC)**0.5_RPREC))
    P2 = (DSIN(ALPHA*PI/180.0_RPREC)**0.5_RPREC)/RE
ELSE
    P1 = RI_SC1**2
    P2 = 1.0_RPREC/RE
ENDIF


!UPDATE RHS
RHS_TF=RHS_T

!SEND RECEIVE
#IF(MPI)
   !--send scalar info down & recv scalar info from above
   call mpi_sendrecv (T(1, 1, 1), ld*ny, MPI_RPREC, down, 20,  &
                      T(1, 1, nz), ld*ny, MPI_RPREC, up, 20,   &
                      comm, status, ierr)
   !--send scalar info up & recv scalar info from below
   call mpi_sendrecv (T(1, 1, nz-1), ld*ny, MPI_RPREC, up, 21,  &
                      T(1, 1, 0), ld*ny, MPI_RPREC, down, 21,   &
                      comm, status, ierr)
#ENDIF

!LBC
if(coord==0)       CALL SET_BC(0,T,LLBC_T,LBC_T,1)
!UBC
if(coord==nproc-1) CALL SET_BC(1,T,LUBC_T,UBC_T,1)

!DERIVATIVES (shift the order by JFANG)
CALL FILT_DA(T,DTDX,DTDY)
CALL DDZ_UV(DTDZ,T)

!stress calculations (need the right dtdz but attention not to overwrite ttz)
if(sgs)then
   call sgs_flux(t,dtdx,dtdy,dtdz,cs,nu_t,               &
                 p2/pr_t,pr_t_sgs,ds,nu_s,   &
                 txth,tyth,tzth,beta_sc1)
else
   txth=dtdx*p2/pr_t
   tyth=dtdy*p2/pr_t
   tzth=dtdz*p2/pr_t
endif

!overwrite deriv. and impose stresses at jz=1 if wall law
if((coord==0).and.(llbc_special.eq.'WALL_LAW'))then
   do j=1,ny
   do i=1,nx
      if(ustar_avg(i,j)==0.0_rprec)then
         tzth(i,j,1)=0.00000001_rprec
      else
         !remember that no - needed if i assume sc_star=-flux/ustar_avg
         dtdz(i,j,1) = PHI_T(i,j)*T_STAR(i,j)/(vonk*get_z_local(1))
         tzth(i,j,1) = T_STAR(i,j)*USTAR_AVG(i,j)
      endif
   enddo
   enddo
endif


!COMPUTE RHS FOR SCALAR EQ.
CALL CONVEC_SC(RHS_T,DTDX,DTDY,DTDZ)

!COMPUTE DERIVATIVES OF SGS FLUXES ONE PER TIME FOR CACHE+MEMORY
CALL DIVSTRESS(RHS_T,TXTH,TYTH,TZTH,.TRUE.)

IF(LSM) then
 DO JZ=0,NZ
    RHS_T(:,:,jz) = RHS_T(:,:,jz) + particle_forcing(:,:,jz)
 ENDDO
ENDIF

!###################################################################
!ADD MISSING TERM IN RHS IN CASE OF SLOPING SURFACES
IF(SLOPE)THEN
  
   A=DSIN(ALPHA*PI/180.0_RPREC)
   B=DCOS(ALPHA*PI/180.0_RPREC)
   !ADD NEW TERMS ON RHS FOR SC
   DO K=1,NZ-1
      RHS_T(:,:,K)=RHS_T(:,:,K)+P1*(U(:,:,K)*A-0.5_RPREC*(W(:,:,K)+W(:,:,K+1))*B)
   ENDDO
   
ENDIF

!ADD SPONGE LAYER (UVP for RHS need interpolation)
IF(T_SPONGE) CALL SET_SCSPONGE(T,RHS_T)

!CALL CALC_BETA(T) ! Done before advancement in time

!ADVANCE IN TIME, ABS2 (EXPLICIT)
CALL STEP_SCALAR(DT,T,RHS_T,RHS_TF)

!NEED TO RE-ENFORCE UBC/LBC (GOOD FOR CALCBETA AND OUTPUT)
if(coord==0)       CALL SET_BC(0,T,LLBC_T,LBC_T,1)
if(coord==nproc-1) CALL SET_BC(1,T,LUBC_T,UBC_T,1)

#IF(MPI)
   !--send scalar info down & recv scalar info from above
   call mpi_sendrecv (T(1, 1, 1), ld*ny, MPI_RPREC, down, 23,  &
                      T(1, 1, nz), ld*ny, MPI_RPREC, up, 23,   &
                      comm, status, ierr)
   !--send scalar info up & recv scalar info from below
   call mpi_sendrecv (T(1, 1, nz-1), ld*ny, MPI_RPREC, up, 24,  &
                      T(1, 1, 0), ld*ny, MPI_RPREC, down, 24,   &
                      comm, status, ierr)
#ENDIF


!###################################################################
!BOUYANCY TERM FOR VERT.MOM.EQUATION
!IF(.NOT.SLOPE)THEN
!   CALL CALCBETA(T,T_INIT,DTDZ_INIT,T_REF,              &
!                 ENT_T_INIT,ENT_DTDZ_INIT,Z_ENT_T,BETA_T)
!   !CALL CALCBETA(T,SC1_INIT/T_SCALE,DSC1DZ_INIT/T_SCALE*Z_I,SC1_REF/T_SCALE,        &
!   !              ENT_SC1_INIT/T_SCALE,ENT_DSC1DZ_INIT/T_SCALE*Z_I,Z_ENT_SC1*L_Z/Z_I,BETA_SCAL)
!ENDIF


IF(VERBOSE) PRINT*,'OUT THETA_ALL_IN_ONE'

END SUBROUTINE SCALAR_SOLVER



!!-----------------------------------------------------------------------------------
!!  OBJECT: subroutine INIT_SC
!!-----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 20/10/2012
!!
!!  DESCRIPTION:
!!
!!  this subroutine aims at providing the I.C. for the problem, it was set to be very
!!  general, should do the same also for u,v,w!
!!
!!-----------------------------------------------------------------------------------
SUBROUTINE INIT_SC(VAR,SC_INIT,DSCDZ_INIT,NOISE_FACTOR,Z_TURB,     &
                   Z_ENT,ENT_SC_INIT,ENT_DSCDZ_INIT,RHS_VAR,       &
                   SC_STAR_INIT, SC_Z0)

USE IC,ONLY:ADD_NOISE

IMPLICIT NONE

REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(OUT) :: VAR
REAL(RPREC),INTENT(IN) :: SC_INIT,DSCDZ_INIT,NOISE_FACTOR,Z_TURB,     &
                          Z_ENT,ENT_SC_INIT,ENT_DSCDZ_INIT,           &
                          SC_STAR_INIT, SC_Z0
REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(INOUT) :: RHS_VAR
REAL(RPREC) :: Z,W_STAR,T_STAR,T_AVG_LBC,WT_AVG_LBC
INTEGER :: JX,JY,JZ
REAL(RPREC) :: SIGMA_0,SIGMA_C,A
integer :: kk

PRINT*,'INSIDE INIT_SC'


IF(VI1)THEN !STANDARD LINEAR/CONST PROFILE

DO JZ=1,NZ
   !COMPUTE Z
   kk = jz+eta_displ(coord)
   Z = get_z_local(2*kk-1)
!   Z=(DBLE(JZ)-0.5_rprec)
!#IF(MPI)
!   !write(*,*) ' i am here: jz,z:',jz,z
!   Z= (DBLE(JZ)-0.5_rprec) +real(eta_displ(coord),rprec)
!#ENDIF
   !ASSIGN VALUE
   VAR(:,:,JZ)=SC_INIT+DSCDZ_INIT*Z
   !INVERSION LAYER
   IF(Z.GE.Z_ENT)THEN
      VAR(:,:,JZ)=VAR(:,:,JZ)+ENT_SC_INIT+ENT_DSCDZ_INIT*(Z-Z_ENT)
   ENDIF
ENDDO

ELSEIF(VI2)THEN !LOG PROFILE

  DO JZ=1,NZ
    kk = jz+eta_displ(coord)
    Z = get_z_local(2*kk-1)
    VAR(:,:,JZ) = SC_INIT + SC_STAR_INIT/VONK * DLOG(Z/SC_Z0)
  ENDDO

ELSEIF(VI3)THEN

   !KATABATIC PRANDTL NORMALIZED SOLUTION
   A=(DSIN(ALPHA*PI/180.0_RPREC)/2.0_rprec/(1.0_RPREC/RE))**0.5
   DO JZ=1,NZ
      Z=(JZ-0.5_RPREC)*DZ
#IF(MPI)
      Z=Z+eta_displ(coord)*DZ
#ENDIF
      SC1(:,:,JZ)=DLBC_SC1*DCOS(Z*A)*EXP(-Z*A) !NON-DIMENSIONAL
 ENDDO
   PRINT*,"SET IC VI3 FOR SCALAR FIELD"

ELSE
  PRINT*,'PLEASE SPECIFY A REASONABLE IC FOR SCALAR'
ENDIF

!ADD NOISE
CALL ADD_NOISE(VAR,NOISE_FACTOR,Z_TURB)

RHS_VAR=0.0_rprec

PRINT*, 'LEAVING INIT_SC'


END SUBROUTINE INIT_SC


!!--------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE CONVEC_SC(RHS,DTDX,DTDY,DTDZ,TXTH,TYTH,TZTH)
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 15/05/2013
!!
!!  DESCRIPTION: This subroutine computes the convective term for the scalar, it
!!               dealias with the 3/2 rule
!!
!!  TO BE IMPROVED: Parallelize accepting higher memory usage!
!!--------------------------------------------------------------------------------
SUBROUTINE CONVEC_SC(RHS,DTDX,DTDY,DTDZ)

USE DEALIAS
USE TEST_FILTERMODULE

IMPLICIT NONE

REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(IN) :: DTDX,DTDY,DTDZ
REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(OUT) :: RHS
INTEGER :: I,J,K
REAL(RPREC) :: A,B,C
!REAL(RPREC),DIMENSION(LD,NY,0:NZ) :: SCR_1,SCR_2
REAL(RPREC),DIMENSION(LD_BIG,NY2,0:NZ):: UI_PAD,DTDXI_PAD,RHS_M

!THIS I TRIED TO MAKE MEMORY EFFICIENT SO MANY SMALL LOOPS
!CHANGE OF VAR NOT TO DESTROY DTDXi AND Ui WITH DEALIAS1

!INITIALIZE RHS_M
RHS_M=0.0_RPREC

!X TERM
!PADDING
CALL DEALIAS1(U,UI_PAD)
CALL DEALIAS1(DTDX,DTDXI_PAD)
!ADD X TERM TO RHS_M
DO K=1,NZ-1
  RHS_M(:,:,K)=RHS_M(:,:,K)-UI_PAD(:,:,K)*DTDXI_PAD(:,:,K)
ENDDO

!Y TERM
!PADDING
CALL DEALIAS1(V,UI_PAD)
CALL DEALIAS1(DTDY,DTDXI_PAD)
DO K=1,NZ-1
  RHS_M(:,:,K)=RHS_M(:,:,K)-UI_PAD(:,:,K)*DTDXI_PAD(:,:,K)
ENDDO

!Z TERM
!PADDING
CALL DEALIAS1(W,UI_PAD)
CALL DEALIAS1(DTDZ,DTDXI_PAD)
DO K=1,NZ-1
  RHS_M(:,:,K)=RHS_M(:,:,K)-0.5_RPREC*(UI_PAD(:,:,K)*DTDXI_PAD(:,:,K)   &
              +UI_PAD(:,:,K+1)*DTDXI_PAD(:,:,K+1))
ENDDO

!TRUNCATION
CALL DEALIAS2(RHS,RHS_M)



END SUBROUTINE CONVEC_SC


!!--------------------------------------------------------------------------------
! !!  OBJECT: SUBROUTINE STEP_SCALAR(RHS,sgs_vert,dTdx,dTdy,dTdz,txTh,tyTh,tzTh)
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 30/07/2012
!!
!!  DESCRIPTION: This subroutine computes the buoyancy term given
!!
!!
!!--------------------------------------------------------------------------------
SUBROUTINE STEP_SCALAR(DT,SCALAR,RHS_1,RHS_2)

IMPLICIT NONE

INTEGER :: I,J,K
REAL(RPREC) :: DT
REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(IN) :: RHS_1,RHS_2
REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(OUT) :: SCALAR

!ABS2 TIME ADVANCEMENT
DO K=1,NZ-1
  SCALAR(:,:,K)=SCALAR(:,:,K)+DT*(1.5_rprec*RHS_1(:,:,K)-0.5_rprec*RHS_2(:,:,K))
ENDDO     

END SUBROUTINE STEP_SCALAR



! !!--------------------------------------------------------------------------------
! !!  OBJECT: SUBROUTINE CALCBETA(SCALAR,BETA_SCAL)
! !!--------------------------------------------------------------------------------
! !!
! !!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 19/03/2013
! !!
! !!  DESCRIPTION: This subroutine computes the buoyancy term
! !!
! !!  NOTES: -BETA is on UVP nodes as well as SC, will need to interpolate 
! !!          adeguately to W nodes when adding it to RHSz 
! !!          (bottom/top def by BC on scalars)
! !!
! !!--------------------------------------------------------------------------------
! SUBROUTINE CALCBETA(SC,SC_INIT,DSCDZ_INIT,SC_REF,ENT_SC_INIT,ENT_DSCDZ_INIT,Z_ENT,BETA_SCAL)
! 
! IMPLICIT NONE
! 
! REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(IN) :: SC
! REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(OUT) :: BETA_SCAL
! REAL(RPREC),INTENT(IN) :: SC_INIT,DSCDZ_INIT,SC_REF,ENT_SC_INIT,ENT_DSCDZ_INIT,Z_ENT
! REAL(RPREC),DIMENSION(LD,NY,0:NZ) :: SC_BAR,SC_TMP
! REAL(RPREC) :: G_HAT,ABOVE,BELOW
! INTEGER :: I,J,K,JZ_MIN
! REAL(RPREC) :: NN,NORM,Z
! 
! !INITIALIZE         @@@@ RECHECK CAREFULLY HOW THIS PART BECOMES WITH THETA
! G_HAT=10.0_RPREC    !@@@@ RANDOM NUMBER TO BE CHANGED!
! BETA_SCAL=0.0_rprec
! 
! !INITIALIZE
! SC_BAR(:,:,:)=0.0_RPREC
! NORM=1.0_RPREC/REAL(NX*NY,RPREC)
! 
! !PLANAR AVERAGE CASE
! IF(SC1_AVGTYPE==1)THEN
! 
! DO K=0,NZ
! DO J=1,NY
! DO I=1,NX
!     SC_BAR(1,1,K)=SC_BAR(1,1,K)+SC(I,J,K)*NORM
! ENDDO
! ENDDO
! ENDDO
! 
! !COMPUTE THE BUOYANCY TERM +g(theta'/<theta>) ON UVP NODES
! !$OMP PARALLEL DO
! DO K=1,NZ-1
! DO J=1,NY
! DO I=1,NX
!     BETA_SCAL(I,J,K)=G_HAT*((SC(I,J,K)-SC_BAR(1,1,K))/SC_BAR(1,1,K))
! ENDDO
! ENDDO
! ENDDO
! !$OMP END PARALLEL DO
! 
! !!VOLUME AVERAGE CASE
! !ELSEIF(SC1_AVGTYPE==2)THEN
! !
! !CALL TEST_FILTER_INIT (2._rprec*BETA_G1,G1,IFILTER_T)
! !
! !
! !!CREATE FILTERED THETA FIELD (UVP)
! !SC_TMP=SC
! !DO K=1,NZ
! !    CALL TEST_FILTER(SC_TMP(:,:,K),G1)
! !ENDDO
! !    
! !!COMPUTE FOR 1NODE UVP AND SET JZ_MIN 
! !!(HERE I SHOULD USE THE BOTTOM ACTUALLY TEMP SO NEED
! !!TO INCLUDE IT ALSO IN THE WALLSTRESS_SCALAR
! !SCALAR_BAR(1) IS NOT T_S OF COURSE, T_S IS A BC WE USE TO THEN COMPUTE 
! !THE DERIVATIVE AND T(1) ADJUST ACCORDINGLY
! !
! !DO J=1,NY
! !DO I=1,NX
! !  SC_BAR(I,J,1) = (SC_TMP(I,J,1)         &
! !                +  SC_TMP(I,J,2))*0.5_RPREC      !SHOULD ADD T_S AND SOMEHOW TAK 						 !INTO CONSIDERATION THAT IS CLOSER
! !ENDDO
! !ENDDO
! !JZ_MIN = 2
! !
! !COMPUTE THE MEAN SCALAR FROM 1 TO NZ-1 (UVP) (PROB NEED MPI SYNC)
! !NOTE THAT WE NEED TO HAVE SCALAR DEFINED AT NZ THROUGH EITHER BC OR SOMETHING ELSE
! !NN=1.0_RPREC/3.0_RPREC
! !!OMP PARALLEL DO
! !DO K=JZ_MIN,NZ-1
! !DO J=1,NY
! !DO I=1,NX
! !    SC_BAR(I,J,K) = (SC_TMP(I,J,K-1)          &
! !		  +  SC_TMP(I,J,K)            &
! !		  +  SC_TMP(I,J,K+1))*NN
! !ENDDO
! !ENDDO
! !ENDDO
! !!OMP END PARALLEL DO
! 
! !COMPUTE THE BUOYANCY TERM ON UVP NODES FROM 1 TO NZ-1
! !OMP PARALLEL DO
! !DO K=1,NZ-1
! !DO J=1,NY
! !DO I=1,NX
! !    BETA_SCAL(I,J,K)=G_HAT*(SC(I,J,K)-SC_BAR(I,J,K))/SC_BAR(I,J,K)
! !ENDDO
! !ENDDO
! !ENDDO
! !OMP END PARALLEL DO
! 
! !!TOPBC
! !IF((.NOT. USE_MPI) .OR. (USE_MPI .AND. COORD == NPROC-1))THEN    
! !    DO J=1,NY
! !    DO I=1,NX
! !      BETA_SCAL(I,J,NZ) = 0.0_RPREC
! !    ENDDO
! !    ENDDO
! !ELSE
! !    DO J=1,NY
! !    DO I=1,NX
! !      BETA_SCAL(I,J,NZ) = BOGUS
! !    ENDDO
! !    ENDDO
! !ENDIF   
! 
! !LBC
! !BETA_SCAL(:,:,0)=BOGUS
!     
! !REFERENCE TEMPERATURE CASE
! ELSEIF(SC1_AVGTYPE==3)THEN
! 
! !DEFINE CONSTANT REF. SC_BAR
! DO K=1,NZ-1
!     !COMPUTE Z
!     Z=(DBLE(K)-0.5_rprec)*DZ
!     SC_BAR(1,1,K)=SC_INIT+DSCDZ_INIT*Z
!     !INVERSION LAYER
!     IF(Z.GE.Z_ENT)THEN
!       SC_BAR(1,1,K)=SC_BAR(1,1,K)+ENT_SC_INIT+ENT_DSCDZ_INIT*(Z-Z_ENT)
!     ENDIF
! ENDDO
! 
! !COMPUTE THE BUOYANCY TERM +g(theta'/theta_ref) ON UVP NODES
! !$OMP PARALLEL DO
! DO K=1,NZ-1
!     BETA_SCAL(:,:,K)=G_HAT/SC_REF*(SC(:,:,K)-SC_BAR(1,1,K))
! ENDDO
! !$OMP END PARALLEL DO
! 
! ENDIF
! 
! 
! END SUBROUTINE CALCBETA

subroutine calc_beta()

use SP_SC1,ONLY:SC1
use SP_SC2,ONLY:SC2

implicit none
integer :: i,j,k,jz_min
real(rprec) :: g_hat,above,below
real(rprec),dimension(0:nz) :: scalar_bar
real(rprec),dimension(ld,ny,0:nz) :: scalar

g_hat = g*(z_i/u_star**2)
BUOYANCY = 0.0_rprec

   IF(TEMPERATURE) then
      IF(HUMIDITY) then
        !scalar=SC1*(1.0_rprec+(0.61_rprec*(SC2+SC3)))
         scalar=SC1*(1.0_rprec+(0.61_rprec*(SC2)))
      ELSE
        scalar=SC1
      ENDIF
   ENDIF

do k=0,nz
  scalar_bar(k) = 0.0_rprec
  do j=1,ny
     do i=1,nx
      scalar_bar(k) = scalar_bar(k) + scalar(i,j,k)
     enddo
  enddo
  scalar_bar(k) = scalar_bar(k) / (real(nx*ny,rprec))
enddo

if(coord==0) then
  jz_min = 2
else
  jz_min = 1
endif

do k=jz_min,nz
    !write(*,*) 'coord,k,scalar:',coord,k,scalar_bar(k)
     do j=1,ny
        do i=1,nx
           above=(scalar(i,j,k)-scalar_bar(k))/scalar_bar(k)
           below=(scalar(i,j,k-1)-scalar_bar(k-1))/scalar_bar(k-1)
           BUOYANCY(i,j,k) = g_hat*(above+below)/2.0_rprec
       enddo
    enddo
enddo
end subroutine calc_beta
END MODULE SCALARS_MODULE


