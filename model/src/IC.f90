MODULE IC

USE TYPES,ONLY:RPREC
USE SERVICE
USE PARAM
USE SIM_PARAM
use stretch

IMPLICIT NONE


! #####################################################################
CONTAINS


!!-----------------------------------------------------------------------------------
!!  OBJECT: subroutine INIT(u,v,w,theta,humid)
!!-----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 20/07/2012
!!
!!  DESCRIPTION:
!!
!!  this subroutine aims at providing the I.C. for the problem.
!!  It also enforces the boundary conditions so that they are satisfied by the I.C.
!!
!!  TO BE IMPROVED: Encapsulate dns_bc
!!
!!-----------------------------------------------------------------------------------
SUBROUTINE INIT(U,V,W,U_STAR,PFX,PFY,DSPL_H)


IMPLICIT NONE

REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(OUT) :: U,V,W
REAL(RPREC) :: Z,U_STAR,WT_AVG_LBC,T_AVG_LBC
REAL(RPREC),DIMENSION(LD,NY) :: ARG,ARG2
REAL(RPREC),INTENT(IN) :: PFX,PFY,DSPL_H
REAL(RPREC) :: DIR_X,DIR_Y,NORM,A
INTEGER :: JX,JZ
integer :: kk

 !FIX W_STAR, T_STAR AND Q_STAR (NOTE THAT THIS IS JUST FOR THE VERY 
!NOTE: WE TRY TO ASSIGN A NOISE PROP TO LOCAL REF VELOCITY, O.06 IS A REF FLUX
!IF(S_FLAG)THEN
!   IF(LLBC_SC1 .EQ. 'DRC')THEN
!      T_AVG_LBC=DLBC_SC1
!      U_STAR=(G/T_AVG_LBC*0.1_rprec)**(1._rprec/3._rprec) !It is dimensional: (m/s)
!   ELSEIF((LLBC_SC1.EQ.'FLX').OR.(LLBC_SC1.EQ.'NEU'))THEN
!      T_AVG_LBC=DLBC_SC1
!      WT_AVG_LBC = DLBC_SC1
!      U_STAR = SIGN((G/T_AVG_LBC*abs(WT_AVG_LBC))**(1._rprec/3._rprec), &
!               DLBC_SC1)   !It is dimensional: (m/s)
!   ENDIF
!   IF(U_STAR.EQ.0.0_RPREC)THEN
!      PRINT*,'ATTENTION W_STAR=0'
!      STOP
!   ENDIF
!ELSE
!    U_STAR=1.0_RPREC
!ENDIF


!NORMALIZE DIRECTION
NORM=(PFX**2+PFY**2)**0.5_RPREC
IF(NORM .GT. 0.0_RPREC)THEN
   DIR_X=PFX/NORM
   DIR_Y=PFY/NORM
ELSE
   DIR_X=1.0_RPREC
   DIR_Y=0.0_RPREC
ENDIF

IF(COORD==0)THEN
   print*,'DIR_X=',DIR_X
   print*,'DIR_Y=',DIR_Y
ENDIF

!MEAN VELOCITY PROFILE -------------------------------------------------
IF(VI1)THEN

   U=U0
   V=V0
   W=W0
   IF(COORD==0) PRINT*,'SET IC V_INIT1'
  
ELSEIF(VI2)THEN
   
   write(*,*) 'Displacement_height:',dspl_h 
   DO JZ=1,NZ
      !DEFINE HEIGHT
      kk = jz + eta_displ(coord)
      Z = get_z_local(2*kk-1) !c_stretch*(exp(kk*dz)-1.0_rprec)
      !Z=sum(dz(1:kk))-(dz(kk)/2.0_rprec)
      IF(Z.GE.DSPL_H)THEN
         !IMPOSE LOG LAW (VECTORS AS LBC_ZO MIGHT VARY IN SPACE)
         ARG2(1:NX,1:NY)=(Z-DSPL_H)/LBC_Z0(1:NX,1:NY)
         ARG=0.0_RPREC
         ARG(1:NX,1:NY)=(1._RPREC/VONK)*DLOG(ARG2(1:NX,1:NY))
         U(:,:,JZ)= DIR_X*SQRT(PFX*get_z_local(2*nz_tot-2))*ARG(:,:)
         V(:,:,JZ)= DIR_Y*SQRT(PFX*get_z_local(2*nz_tot-2))*ARG(:,:)
         !U(:,:,JZ) = DIR_X*PFX*ARG(:,:)
         !V(:,:,JZ) = DIR_Y*PFX*ARG(:,:)        
         W(:,:,JZ)= 0.00_RPREC
      ELSE
         U(:,:,JZ)=0.0_RPREC
         V(:,:,JZ)=0.0_RPREC
         W(:,:,JZ)=0.0_RPREC
      ENDIF
   ENDDO
   IF(COORD==0) PRINT*,'SET IC V_INIT2'
   
ELSEIF(VI3)THEN   

   !KAT_START
   !U=0.0_RPREC
   !V=0.0_RPREC
   !W=0.0_RPREC
   !IF(COORD==0.or.coord==1.or.coord==2)THEN
   !DO JX=1,NX
   !   U(JX,:,1:NZ) = 0.1*DSIN(10.0_RPREC*2.0_RPREC*PI*(JX-1)/NX)
   !ENDDO
   !ENDIF
   !POISEUILLE
   !DO JZ=1,NZ
   !   Z=(JZ-0.5_RPREC)*DZ+ETA_DISPL(COORD)*DZ
   !   U(:,:,JZ)=RE*Z*(1.0_RPREC-0.5_RPREC*Z)  !NON-DIMENSIONAL
   !ENDDO
   !V=0.0_RPREC
   !W=0.0_RPREC
   
   !KATABATIC PRANDTL NORMALIZED SOLUTION
   A=(DSIN(ALPHA*PI/180.0_RPREC)/2.0_rprec/(1.0_RPREC/RE))**0.5
   DO JZ=1,NZ

      kk = jz + eta_displ(coord)
!      Z=sum(dz(1:kk))-(dz(kk)/2.0_rprec)
      Z = get_z_local(2*kk) !c_stretch*(exp(kk*dz)-1.0_rprec)
      U(:,:,JZ)=-0.5*DLBC_SC1*DSIN(Z*A)*EXP(-Z*A) !NON-DIMENSIONAL
 ENDDO
   V=0.0_RPREC
   W=0.0_RPREC
   IF(COORD==0) PRINT*,"SET IC V_INIT3"

ELSEIF(VI4)THEN 
   
   DO JZ=1,NZ

      kk = jz + eta_displ(coord)
      !Z=sum(dz(1:kk))-(dz(kk)/2.0_rprec)
      Z = get_z_local(2*kk)
      !Z = c_stretch*(exp(kk*dz)-1.0_rprec)
      U(:,:,JZ)=(U0+DUDZ0*Z)  !CONST. GRADIENT
   ENDDO
   V=V0
   W=W0
   IF(COORD==0) PRINT*,'SET IC V_INIT4'

ELSEIF(VI5)THEN
   
   write(*,*) 'Using composite profile with Z_BLEND=',Z_BLEND
   DO JZ=1,NZ
      kk = jz + eta_displ(coord)
      Z = get_z_local(2*kk-1)
      
      IF(Z.GE.DSPL_H)THEN
         IF(Z.LE.Z_BLEND) THEN
            !IMPOSE LOG LAW (VECTORS AS LBC_ZO MIGHT VARY IN SPACE)
            ARG2(1:NX,1:NY)=(Z-DSPL_H)/LBC_Z0(1:NX,1:NY)
            ARG=0.0_RPREC
            ARG(1:NX,1:NY)=(1._RPREC/VONK)*DLOG(ARG2(1:NX,1:NY))
            U(:,:,JZ)= DIR_X*SQRT(PFX*get_z_local(2*nz_tot-2))*ARG(:,:)
            V(:,:,JZ)= DIR_Y*SQRT(PFX*get_z_local(2*nz_tot-2))*ARG(:,:)
            W(:,:,JZ)= 0.00_RPREC
         ELSE
            ! Constant profile above Z_BLEND
            ! Get the velocity at Z_BLEND
            !ARG2(1:NX,1:NY)=(Z_BLEND-DSPL_H)/LBC_Z0(1:NX,1:NY)
            !ARG=(1._RPREC/VONK)*DLOG(ARG2(1,1))
            !U(:,:,JZ)= DIR_X*SQRT(PFX*get_z_local(2*nz_tot-2))*ARG(:,:)
            !V(:,:,JZ)= DIR_Y*SQRT(PFX*get_z_local(2*nz_tot-2))*ARG(:,:)
            !W(:,:,JZ)= 0.00_RPREC
            ARG2(1:NX,1:NY)=(Z-DSPL_H)/LBC_Z0(1:NX,1:NY)
            ARG=0.0_RPREC
            ARG(1:NX,1:NY)=(1._RPREC/VONK)*DLOG(ARG2(1:NX,1:NY))
            U(:,:,JZ)= DIR_X*SQRT(PFX*get_z_local(2*nz_tot-2))*ARG(:,:)
            V(:,:,JZ)= DIR_Y*SQRT(PFX*get_z_local(2*nz_tot-2))*ARG(:,:)
            W(:,:,JZ)= 0.00_RPREC
         ENDIF
      ELSE
         U(:,:,JZ)=0.0_RPREC
         V(:,:,JZ)=0.0_RPREC
         W(:,:,JZ)=0.0_RPREC
      ENDIF
   ENDDO
   IF(COORD==0) PRINT*,'SET IC V_INIT5 (COMPOSITE PROFILE)'

ENDIF

!NOISE -----------------------------------------------------
IF(COORD==0) PRINT*,'ADDING NOISE TO MOMENTUM'
IF(COORD==0) PRINT*,'U_STAR=',U_STAR
CALL MPI_BARRIER(COMM,IERR)
CALL ADD_NOISE(U,U_STAR,Z_TURB*L_Z)
CALL MPI_BARRIER(COMM,IERR)
CALL ADD_NOISE(V,U_STAR,Z_TURB*L_Z)
CALL MPI_BARRIER(COMM,IERR)
CALL ADD_NOISE(W,U_STAR,Z_TURB*L_Z)

!!! SETTING RHS=0
RHSX=0.0_rprec
RHSy=0.0_rprec
RHSz=0.0_rprec

END SUBROUTINE INIT


!MARCO: SUPPORT SUBROUTINE, ASSUME UNIFORM DISTRIBUTION
SUBROUTINE ADD_NOISE(X,U_STAR,ZTURB)

   IMPLICIT NONE
   
   REAL(RPREC),DIMENSION(:,:,:),INTENT(INOUT) :: X
   REAL(RPREC),INTENT(IN) :: U_STAR,ZTURB
   REAL(RPREC) :: NOISE,Z,RAN3,VARIANCE_FACT,AVG_NOISE
   INTEGER :: JX,JY,JZ,SEED,ID
   INTEGER,SAVE :: RND=0
integer :: kk

   IF(U_STAR.GT.0.0_RPREC)THEN
      RND=RND+4
      ID=0
      AVG_NOISE=0
      DO JZ=1,NZ

         SEED=-80-JZ-RND
#IF(MPI)
         SEED=SEED-eta_displ(coord)
#ENDIF
      kk = jz + eta_displ(coord)
!      Z=sum(dz(1:kk))-(dz(kk)/2.0_rprec)
Z = get_z_local(2*kk) !c_stretch*(exp(kk*dz)-1.0_rprec)
         DO JY=1,NY
         DO JX=1,NX
            IF (Z .LE. ZTURB) THEN
            NOISE=RAN3(SEED)*2.0_RPREC-1.0_RPREC
            AVG_NOISE=AVG_NOISE+NOISE
            ID=ID+1
            VARIANCE_FACT=3.0_RPREC*VNF*U_STAR**2
            X(JX,JY,JZ)=X(JX,JY,JZ)+VARIANCE_FACT*(NOISE*(get_z_local(2*nz_tot-2)-Z*0.9_RPREC)/get_z_local(2*nz_tot-2))  
            ENDIF
         ENDDO
         ENDDO
         ENDDO
   !IF(COORD==0) PRINT*,'ID=',ID
   !IF(ID.GT.0)THEN 
   !   PRINT*,'AVG_NOISE = ',AVG_NOISE/ID,' FOR COORD=',COORD
   !ELSE
   !   PRINT*,'AVG_NOISE = 0 FOR COORD ',COORD
   !ENDIF
   ENDIF
END SUBROUTINE ADD_NOISE   


END MODULE IC
