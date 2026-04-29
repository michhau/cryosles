MODULE BC

USE TYPES
USE PARAM
USE SERVICE
USE SIM_PARAM
use stretch

IMPLICIT NONE

!---------------------------------------------------------------------------------
CONTAINS


!!--------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE SET_ZO()
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 16/10/2012
!!
!!  ITS DONE LIKE THIS SO THAT CAN BE NICELY CONNECTED WITH REMOTE OR DIFF VALUES
!!
!!--------------------------------------------------------------------------------
SUBROUTINE SET_LBC_ZO(LBC_ZO)

IMPLICIT NONE

REAL(RPREC),DIMENSION(:,:),INTENT(INOUT) :: LBC_ZO 

LBC_ZO=LBC_DZ0

!!!xxxnew implementation for roughness length for heterogeneous surfaces

!LBC_ZO(16:48,:) = LBC_DZ0 * 10.0_rprec
 
END SUBROUTINE SET_LBC_ZO



!!--------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE SET_LBC_ZO_HETERO()
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Michi Haugeneder 29.04.2026
!!
!!  Set lower BC for aerodynamic roughness length for heterogeneous surfaces.
!!
!!--------------------------------------------------------------------------------
SUBROUTINE SET_LBC_ZO_HETERO(LBC_ZO, PATCHES, NX, NY, Z0_0, Z0_1)
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: NX, NY
   LOGICAL, INTENT(IN) :: PATCHES(NX,NY)
   REAL(RPREC), INTENT(IN) :: Z0_0
   REAL(RPREC), INTENT(IN) :: Z0_1

   REAL(RPREC),DIMENSION(:,:),INTENT(OUT) :: LBC_ZO

   INTEGER :: I, J, PATCH_I

   ! Active physical x-cells are 1:NX.  Extra LD padding is filled
   ! periodically only to keep full-array operations initialized.
   DO J = 1, NY
     DO I = 1, SIZE(LBC_ZO,1)
       PATCH_I = MOD(I-1, NX) + 1
       IF (PATCHES(PATCH_I,J)) THEN
         LBC_ZO(I,J) = Z0_1
       ELSE
         LBC_ZO(I,J) = Z0_0
       END IF
     END DO
   END DO

END SUBROUTINE SET_LBC_ZO_HETERO



!!--------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE SET_LBC_SC_ZO()
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 16/10/2012
!!
!!  ITS DONE LIKE THIS SO THAT CAN BE NICELY CONNECTED WITH REMOTE OR DIFF VALUES
!!
!!--------------------------------------------------------------------------------
SUBROUTINE SET_LBC_SC_ZO(LBC_SC_ZO,SCALAR)

IMPLICIT NONE

REAL(RPREC),DIMENSION(:,:),INTENT(INOUT) :: LBC_SC_ZO
INTEGER,INTENT(IN) :: SCALAR

IF(SCALAR .EQ. 1) then
 ! z0h = 0.05*z0m according to kB^(-1)=3 (see Obsidian Michi Haugeneder 29.04.2026)
 LBC_SC_ZO=LBC_DZ0*0.05_rprec
ELSEIF(SCALAR .EQ. 2) then
 LBC_SC_ZO=LBC_DZ0*0.1_rprec
ENDIF   

END SUBROUTINE SET_LBC_SC_ZO



!!--------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE SET_LBC_SC_ZO_HETERO()
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Michi Haugeneder 02.04.2025
!!
!!  Set lower BC for Z0 for heterogeneous surfaces
!!
!!--------------------------------------------------------------------------------
SUBROUTINE SET_LBC_SC_ZO_HETERO(LBC_SC_ZO, PATCHES, NX, NY, Z0_0, Z0_1)
   IMPLICIT NONE
   
   ! Input variables
   INTEGER, INTENT(IN) :: NX, NY
   LOGICAL, INTENT(IN) :: PATCHES(NX,NY)
   REAL(RPREC), INTENT(IN) :: Z0_0      ! Value for false elements
   REAL(RPREC), INTENT(IN) :: Z0_1      ! Value for true elements
   
   ! Output variable
   REAL(RPREC),DIMENSION(:,:),INTENT(OUT) :: LBC_SC_ZO
   
   ! Local variables
   INTEGER :: I, J, PATCH_I
   
   ! Active physical x-cells are 1:NX.  Extra LD padding is filled
   ! periodically only to keep full-array operations initialized.
   ! z0h = 0.05*z0m according to kB^(-1)=3 (see Obsidian Michi Haugeneder 29.04.2026)
   DO J = 1, NY
     DO I = 1, SIZE(LBC_SC_ZO,1)
       PATCH_I = MOD(I-1, NX) + 1
       IF (PATCHES(PATCH_I,J)) THEN
         LBC_SC_ZO(I,J) = Z0_1*0.05_rprec
       ELSE
         LBC_SC_ZO(I,J) = Z0_0*0.05_rprec
       END IF
     END DO
   END DO
   
 END SUBROUTINE SET_LBC_SC_ZO_HETERO



!!--------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE SET_INLET(U,V,W,BC_INLET,U_REF,Z_REF,ALPHA)
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 16/10/2012
!!
!!  U/V/W(LD,NY,0:NZ)
!!
!!--------------------------------------------------------------------------------
SUBROUTINE SET_INLET()

use SP_SC1,only:SC1
use SP_SC2,only:SC2
!use SP_SC3,only:SC3
IMPLICIT NONE

character(256) :: opath,filename

IF(VERBOSE) WRITE (*,*) 'IN SET_INLET'

write(opath,'(a)') trim(out_path)//'inlet_slices/'

write(filename,'(a,i0)') trim(opath)//'u_YZ.c',coord
open(111,file=trim(filename),status='unknown',access='direct',recl=8*ny*(nz-1))
read(111,rec=jt) u((nx/2),1:ny,1:nz-1)
close(111)

write(filename,'(a,i0)') trim(opath)//'v_YZ.c',coord
open(111,file=trim(filename),status='unknown',access='direct',recl=8*ny*(nz-1))
read(111,rec=jt) v((nx/2),1:ny,1:nz-1)
close(111)

write(filename,'(a,i0)') trim(opath)//'w_YZ.c',coord
open(111,file=trim(filename),status='unknown',access='direct',recl=8*ny*(nz-1))
read(111,rec=jt) w((nx/2),1:ny,1:nz-1)
close(111)

if(temperature) then
  write(filename,'(a,i0)') trim(opath)//'t_YZ.c',coord
  open(111,file=trim(filename),status='unknown',access='direct',recl=8*ny*(nz-1))
  read(111,rec=jt) sc1((nx/2),1:ny,1:nz-1)
  close(111)
endif

if(humidity) then
  write(filename,'(a,i0)') trim(opath)//'q_YZ.c',coord
  open(111,file=trim(filename),status='unknown',access='direct',recl=8*ny*(nz-1))
  read(111,rec=jt) sc2((nx/2),1:ny,1:nz-1)
  close(111)
endif
!if(humidity) then
!  write(filename,'(a,i0)') trim(opath)//'q_part_YZ.c',coord
!  open(111,file=trim(filename),status='unknown',access='direct',recl=8*ny*(nz-1))
!  read(111,rec=jt) sc3((nx/2),1:ny,1:nz-1)
!  close(111)
!endif

CALL MPI_SENDRECV (U(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 1,  &
                U(1, 1, NZ), LD*NY, MPI_RPREC, UP, 1,   &
                COMM, STATUS, IERR)
CALL MPI_SENDRECV (V(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 2,  &
                V(1, 1, NZ), LD*NY, MPI_RPREC, UP, 2,   &
                COMM, STATUS, IERR)
CALL MPI_SENDRECV (W(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 3,  &
                W(1, 1, NZ), LD*NY, MPI_RPREC, UP, 3,   &
                COMM, STATUS, IERR)
CALL MPI_SENDRECV (U(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 4,  &
                U(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 4,   &
                COMM, STATUS, IERR)
CALL MPI_SENDRECV (V(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 5,  &
                V(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 5,   &
                COMM, STATUS, IERR)
CALL MPI_SENDRECV (W(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 6,  &
                W(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 6,   &
                COMM, STATUS, IERR)

IF(TEMPERATURE) then
CALL MPI_SENDRECV (SC1(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 300,  &
                SC1(1, 1, NZ), LD*NY, MPI_RPREC, UP, 300,   &
                COMM, STATUS, IERR)
CALL MPI_SENDRECV (SC1(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 400,  &
                SC1(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 400,   &
                COMM, STATUS, IERR)
ENDIF

IF(HUMIDITY) then
CALL MPI_SENDRECV (SC2(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 300,  &
                SC2(1, 1, NZ), LD*NY, MPI_RPREC, UP, 300,   &
                COMM, STATUS, IERR)
CALL MPI_SENDRECV (SC2(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 400,  &
                SC2(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 400,   &
                COMM, STATUS, IERR)
ENDIF

!IF(HUMIDITY) then
!CALL MPI_SENDRECV (SC3(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 300,  &
!                SC3(1, 1, NZ), LD*NY, MPI_RPREC, UP, 300,   &
!                COMM, STATUS, IERR)
!CALL MPI_SENDRECV (SC3(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 400,  &
!                SC3(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 400,   &
!                COMM, STATUS, IERR)
!ENDIF

IF(VERBOSE) WRITE (*,*) 'OUT SET_INLET'

END SUBROUTINE SET_INLET



!!--------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE COMPUTE_LBC(LBC_MOM)
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 06/02/2013
!!
!!  ITS DONE LIKE THIS SO THAT CAN BE NICELY CONNECTED WITH REMOTE OR DIFF VALUES
!!
!!--------------------------------------------------------------------------------
SUBROUTINE COMPUTE_LBC(LBC_U,LBC_V,LBC_W)

IMPLICIT NONE

REAL(RPREC),DIMENSION(:,:),INTENT(OUT) :: LBC_U,LBC_V,LBC_W	!(LD,NY)

IF(VERBOSE) WRITE (*,*) 'IN COMPUTE_LBC'

!STANDARD CASES
IF(LLBC_U.EQ.'DRC') LBC_U=DLBC_U
IF(LLBC_U.EQ.'NEU') LBC_U=DLBC_U
IF(LLBC_V.EQ.'DRC') LBC_V=DLBC_U
IF(LLBC_V.EQ.'NEU') LBC_V=DLBC_U
IF(LLBC_W.EQ.'DRC') LBC_W=DLBC_U
IF(LLBC_W.EQ.'NEU') LBC_W=DLBC_U

!SPECIAL CASES

IF(VERBOSE) WRITE (*,*) 'OUT COMPUTE_LBC'

END SUBROUTINE COMPUTE_LBC



!!--------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE COMPUTE_LBC(LBC_MOM)
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 06/02/2013
!!
!!  ITS DONE LIKE THIS SO THAT CAN BE NICELY CONNECTED WITH REMOTE OR DIFF VALUES
!!
!!--------------------------------------------------------------------------------
SUBROUTINE COMPUTE_UBC(UBC_U,UBC_V,UBC_W)

IMPLICIT NONE

REAL(RPREC),DIMENSION(:,:),INTENT(OUT) :: UBC_U,UBC_V,UBC_W

IF(VERBOSE) WRITE (*,*) 'IN COMPUTE_UBC'

!STANDARD CASES
IF(LUBC_U.EQ.'DRC') UBC_U=DUBC_U
IF(LUBC_U.EQ.'NEU') UBC_U=DUBC_U
IF(LUBC_V.EQ.'DRC') UBC_V=DUBC_V
IF(LUBC_V.EQ.'NEU') UBC_V=DUBC_V
IF(LUBC_W.EQ.'DRC') UBC_W=DUBC_W
IF(LUBC_W.EQ.'NEU') UBC_W=DUBC_W

!SPECIAL CASES

IF(VERBOSE) WRITE (*,*) 'OUT COMPUTE_UBC'

END SUBROUTINE COMPUTE_UBC



!!--------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE COMPUTE_LBC_SC(LBC_SC,SCALAR)
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 16/10/2012
!!
!!  ITS DONE LIKE THIS SO THAT CAN BE NICELY CONNECTED WITH REMOTE OR DIFF VALUES
!!
!!--------------------------------------------------------------------------------
SUBROUTINE COMPUTE_LBC_SC(LBC_SC,SCALAR)

IMPLICIT NONE

INTEGER,INTENT(IN) :: SCALAR
REAL(RPREC),DIMENSION(:,:),INTENT(OUT) :: LBC_SC

IF(VERBOSE) WRITE (*,*) 'IN COMPUTE_LBC_SC'

!VERY STRAIGHT FOR NOW
IF(SCALAR.EQ.1)THEN
   IF(LLBC_SC1.EQ.'DRC') LBC_SC=DLBC_SC1
   IF(LLBC_SC1.EQ.'NEU') LBC_SC=DLBC_SC1
   IF(LLBC_SC1.EQ.'FLX') LBC_SC=DLBC_SC1
ELSEIF(SCALAR.EQ.2)THEN
   IF(LLBC_SC2.EQ.'DRC') LBC_SC=DLBC_SC2
   IF(LLBC_SC2.EQ.'NEU') LBC_SC=DLBC_SC2
   IF(LLBC_SC2.EQ.'FLX') LBC_SC=DLBC_SC2
ENDIF

!xxximplement LBC for temperature here for heterogeneous surfaces

IF(VERBOSE) WRITE (*,*) 'OUT COMPUTE_LBC_SC'

END SUBROUTINE COMPUTE_LBC_SC



!!--------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE COMPUTE_LBC_SC_HETERO()
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Michi Haugeneder 02.04.2025
!!
!!  Set lower BC for temperature (SC) for heterogeneous surfaces
!!
!!--------------------------------------------------------------------------------
SUBROUTINE COMPUTE_LBC_SC_HETERO(LBC_SC, PATCHES, NX, NY, T_0, T_1)
   IMPLICIT NONE
   
   ! Input variables
   INTEGER, INTENT(IN) :: NX, NY
   LOGICAL, INTENT(IN) :: PATCHES(NX,NY)
   REAL(RPREC), INTENT(IN) :: T_0      ! Value for false elements
   REAL(RPREC), INTENT(IN) :: T_1      ! Value for true elements
   
   ! Output variable
   REAL(RPREC),DIMENSION(:,:),INTENT(OUT) :: LBC_SC
   
   ! Local variables
   INTEGER :: I, J, PATCH_I
   REAL(RPREC) :: T_0_NORM, T_1_NORM
   
   ! Normalize the temperature values by dividing by T_SCALE
   T_0_NORM = T_0 / T_SCALE
   T_1_NORM = T_1 / T_SCALE
   
   ! Active physical x-cells are 1:NX.  Extra LD padding is filled
   ! periodically only to keep full-array operations initialized.
   DO J = 1, NY
     DO I = 1, SIZE(LBC_SC,1)
       PATCH_I = MOD(I-1, NX) + 1
       IF (PATCHES(PATCH_I,J)) THEN
         LBC_SC(I,J) = T_1_NORM
       ELSE
         LBC_SC(I,J) = T_0_NORM
       END IF
     END DO
   END DO
   
 END SUBROUTINE COMPUTE_LBC_SC_HETERO



!!--------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE COMPUTE_UBC_SC(UBC_SC,SCALAR)
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 16/10/2012
!!
!!  ITS DONE LIKE THIS SO THAT CAN BE NICELY CONNECTED WITH REMOTE OR DIFF VALUES
!!  
!!  TO IMPROVE:
!!
!!  Perform the adimensionalization directly in PARAM.F90
!!
!!--------------------------------------------------------------------------------
SUBROUTINE COMPUTE_UBC_SC(UBC_SC,SCALAR)

IMPLICIT NONE

INTEGER,INTENT(IN) :: SCALAR
REAL(RPREC),DIMENSION(:,:),INTENT(OUT) :: UBC_SC

IF(VERBOSE) WRITE (*,*) 'IN COMPUTE_UBC_SC'

!STANDARD CASES FOR NOW
IF(SCALAR.EQ.1)THEN
   IF(LUBC_SC1.EQ.'DRC') UBC_SC=DUBC_SC1
   IF(LUBC_SC1.EQ.'NEU') UBC_SC=DUBC_SC1
   IF(LUBC_SC1.EQ.'FLX') UBC_SC=DUBC_SC1
ELSEIF(SCALAR.EQ.2)THEN
   IF(LUBC_SC2.EQ.'DRC') UBC_SC=DUBC_SC2
   IF(LUBC_SC2.EQ.'NEU') UBC_SC=DUBC_SC2
   IF(LUBC_SC2.EQ.'FLX') UBC_SC=DUBC_SC2
ENDIF

IF(VERBOSE) WRITE (*,*) 'OUT COMPUTE_UBC_SC'

END SUBROUTINE COMPUTE_UBC_SC



!!--------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE SET_SC_LBC(VAR,LLBC_SC,LBC_SC)
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 16/10/2012
!!
!!  NOTE: VARIABLES SHOULD BE ALREADY NORMALIZED
!!
!!  BC(LD,NY)
!!  VAR(LD,NY,0:NZ)
!!
!!  BC --> 0 if LBC, 1 if UBC
!!--------------------------------------------------------------------------------
SUBROUTINE SET_BC(IND,VAR,BC_TYPE,BC,STAG)

IMPLICIT NONE

INTEGER,INTENT(IN) :: IND,STAG
REAL(RPREC),DIMENSION(:,:),INTENT(IN) :: BC
REAL(RPREC),DIMENSION(:,:,0:),INTENT(INOUT) :: VAR
CHARACTER(*),INTENT(IN) :: BC_TYPE
INTEGER :: JX,JY

IF(VERBOSE) WRITE (*,*) 'IN SET_BC'

IF(IND.EQ.0)THEN

   IF(STAG.EQ.1)THEN
      IF(BC_TYPE .EQ. "DRC") VAR(:,:,0)=2.0_RPREC*BC(:,:)-VAR(:,:,1)
      IF(BC_TYPE .EQ. "NEU") VAR(:,:,0)=VAR(:,:,1)-BC(:,:)*get_dz_local(1)*2.0_rprec
   ELSEIF(STAG.EQ.0)THEN
      IF(BC_TYPE .EQ. "DRC") VAR(:,:,1)=BC(:,:)
      IF(BC_TYPE .EQ. "NEU") PRINT*,'DO NOT NEED (SET_BC1)'
   ENDIF

ELSEIF(IND.EQ.1)THEN

   IF(STAG.EQ.1)THEN
      IF(BC_TYPE .EQ. "DRC") VAR(:,:,NZ)=2.0_RPREC*BC(:,:)-VAR(:,:,NZ-1)
      IF(BC_TYPE .EQ. "NEU") VAR(:,:,NZ)=VAR(:,:,NZ-1)+BC(:,:)*get_dz_local(2*nz_tot)
   ELSEIF(STAG.EQ.0)THEN
      IF(BC_TYPE .EQ. "DRC") VAR(:,:,NZ)=BC(:,:)
      IF(BC_TYPE .EQ. "NEU") PRINT*,'DO NOT NEED (SET_BC)'
   ENDIF

ENDIF

IF(VERBOSE) WRITE (*,*) 'OUT SET_BC'

END SUBROUTINE SET_BC



!!--------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE SET_VALUE(VAR,LAYER)
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 31/07/2012
!!
!!--------------------------------------------------------------------------------
SUBROUTINE SET_VALUE(VAR,LAYER,VALUE)

IMPLICIT NONE

INTEGER,INTENT(IN) :: LAYER
REAL(RPREC),INTENT(IN) :: VALUE
REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(OUT) :: VAR

IF(VERBOSE) WRITE (*,*) 'IN SET_VALUE'

IF(LAYER==NZ .OR. LAYER==0)THEN
   IF((.NOT. USE_MPI) .OR. (USE_MPI .AND. COORD==NPROC-1))THEN
      VAR(:,:,LAYER) = VALUE
   ENDIF
ELSE
   VAR(:,:,LAYER) = VALUE
ENDIF

IF(VERBOSE) WRITE (*,*) 'OUT SET_VALUE'

END SUBROUTINE SET_VALUE


END MODULE BC
