MODULE CNVC

CONTAINS

!!--------------------------------------------------------------------------------
!!  OBJECT: subroutine convec (cx, cy, cz)
!!-----------------------------------------------------------------------------------
!!
!!  LAST CEHCKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 02/10/2011
!!
!!  DESCRIPTION:
!!
!!  This subroutine provides cx, cy, cz at 1:nz-1;
!!  cx,y,z = - (u X vorticity);
!!  Compute rotational convective term in physical space  (X-mom eq.);
!!  w, dwdx, dwdy, dudz and dvdz are on W nodes, rest on UVP nodes;
!!  (Fixed k=Nz plane on 1/24)   -->  ??
!!  corrected near wall on 4/14/96 {w(DZ/2)=0.5*w(DZ)};
!!  Note from someone: better to have 0.25*w(dz), since really parabolic;
!!
!!-----------------------------------------------------------------------------------
SUBROUTINE CONVEC(CX,CY,CZ)
USE TYPES,ONLY:RPREC
USE SERVICE,ONLY:Z_MPROFILE
USE PARAM
USE SIM_PARAM, ONLY : U,V,W,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY,T3A,T3B
USE FFT

IMPLICIT NONE

!NOTE THAT CX,CY,CZ HAVE 2(NX/2+1) AS DIMENTION AND Y (SYMM OF FOURIER TRANSF OF REAL SIGNAL)
REAL(RPREC),DIMENSION(LD,NY,0:NZ) :: CX,CY,CZ
LOGICAL, PARAMETER :: DEBUG = .FALSE.
INTEGER :: JZ
INTEGER :: JZ_MIN, JZ_MAX
REAL(RPREC),DIMENSION(LD_BIG,NY2,NZ) :: CC_BIG
REAL(RPREC),DIMENSION(LD_BIG,NY2,0:NZ) :: U_BIG, V_BIG, W_BIG
REAL(RPREC),DIMENSION(LD_BIG,NY2,NZ) :: VORT1_BIG, VORT2_BIG, VORT3_BIG
REAL(RPREC) :: IGNORE_ME, CONST, FACTOR

IF(VERBOSE) WRITE (*,*) 'STARTED CONVEC'

!RECALL DUDZ, AND DVDZ ARE ON UVP NODE FOR K=1 ONLY
   !USE CX,CY,CZ FOR TEMP STORAGE HERE

!normalization constant
CONST=1.0_RPREC/(NX*NY)

! ####  COULD POTENTIALLY MERGE ALSO THESE NEXT TWO LOOPS
! ####  CHECK IF THE IF CLAUSE IS SLOWING DOWN COMPUATATIONS OF NOT! IT MIGHT
!FOUNDAMENTAL TO HAVE U_BIG(JZ=0) E V_BIG(JZ=0) AND W_BIG(NZ) FOR INTERPOLATION
!$OMP PARALLEL DO
DO JZ=0,NZ

   IF(JZ.LT.NZ)THEN
      !USE CX,CY,CZ FOR TEMP STORAGE HERE
      CX(:,:,JZ)=-CONST*U(:,:,JZ)
      CY(:,:,JZ)=-CONST*V(:,:,JZ)

      !DO FORWARD FFT ON NORMAL-SIZE ARRAYS
      CALL DFFTW_EXECUTE_DFT_R2C(FORW,CX(:,:,JZ),CX(:,:,JZ))
      CALL DFFTW_EXECUTE_DFT_R2C(FORW,CY(:,:,JZ),CY(:,:,JZ))
  
      !ZERO PAD: PADD TAKES ALSO CARE OF THE ODDBALLS
      CALL PADD(U_BIG(:,:,JZ),CX(:,:,JZ))
      CALL PADD(V_BIG(:,:,JZ),CY(:,:,JZ))
    
      !BACK TO PHYSICAL SPACE (THE NORMALIZATION SHOULD BE OK)
      CALL DFFTW_EXECUTE_DFT_C2R(BACK_BIG,U_BIG(:,:,JZ),U_BIG(:,:,JZ))
      CALL DFFTW_EXECUTE_DFT_C2R(BACK_BIG,V_BIG(:,:,JZ),V_BIG(:,:,JZ))
   ENDIF

   !THESE OPERATIONS DO NOT NEED LEVEL 0
   IF(JZ.GT.0)THEN
      CZ(:,:,JZ)=-CONST*W(:,:,JZ)
      CALL DFFTW_EXECUTE_DFT_R2C(FORW,CZ(:,:,JZ),CZ(:,:,JZ))
      CALL PADD(W_BIG(:,:,JZ),CZ(:,:,JZ))
      CALL DFFTW_EXECUTE_DFT_C2R(BACK_BIG,W_BIG(:,:,JZ),W_BIG(:,:,JZ))

!SPECIAL TREATMENT JZ=1 CAUSE VORT1,VORT2 IS ON UVP (VORT3 ALWAYS UVP)
!IF((SGS.EQV..TRUE.).AND.(LLBC_SPECIAL.EQ.'WALL_LAW').and.(COORD==0))THEN
!      CX(:,:,1) = CONST*(0.5_RPREC*(DWDY(:,:,1)+DWDY(:,:,2))-DVDZ(:,:,1)) !SPECIAL UVP
!      CY(:,:,1) = CONST*(DUDZ(:,:,1)-0.5_RPREC*(DWDX(:,:,1)+DWDX(:,:,2))) !SPECIAL UVP
!      CZ(:,:,1) = CONST*(DVDX(:,:,1)-DUDY(:,:,1))                         !NORMAL UVP
!      JZ_MIN=2
!ELSE
!      JZ_MIN=1
!ENDIF

      !COMPUTE VORTICITY AND PADD
      !SINGLE COMPONENTS OF VORTICITY FIELD
      CX(:,:,JZ) = CONST*(DWDY(:,:,JZ)-DVDZ(:,:,JZ))     !W NODES
      CY(:,:,JZ) = CONST*(DUDZ(:,:,JZ)-DWDX(:,:,JZ))     !W NODES
      CZ(:,:,JZ) = CONST*(DVDX(:,:,JZ)-DUDY(:,:,JZ))     !UVP NODES
      CALL DFFTW_EXECUTE_DFT_R2C(FORW,CX(:,:,JZ),CX(:,:,JZ))
      CALL DFFTW_EXECUTE_DFT_R2C(FORW,CY(:,:,JZ),CY(:,:,JZ))
      !ZERO PAD: PADD TAKES CARE OF THE ODDBALLS
      CALL PADD(VORT1_BIG(:,:,JZ),CX(:,:,JZ))
      CALL PADD(VORT2_BIG(:,:,JZ),CY(:,:,JZ))
      !BACK TO PHYSICAL SPACE
      CALL DFFTW_EXECUTE_DFT_C2R(BACK_BIG,VORT1_BIG(:,:,JZ),VORT1_BIG(:,:,JZ))
      CALL DFFTW_EXECUTE_DFT_C2R(BACK_BIG,VORT2_BIG(:,:,JZ),VORT2_BIG(:,:,JZ))
      IF(JZ.LT.NZ)THEN
         CALL DFFTW_EXECUTE_DFT_R2C(FORW,CZ(:,:,JZ),CZ(:,:,JZ))
         CALL PADD(VORT3_BIG(:,:,JZ),CZ(:,:,JZ))
         CALL DFFTW_EXECUTE_DFT_C2R(BACK_BIG,VORT3_BIG(:,:,JZ),VORT3_BIG(:,:,JZ))
      ENDIF
   
   ENDIF
ENDDO
!$OMP END PARALLEL DO


! ##################################################################################
! MG:  UP TO HERE ALL SEEMS PRETTY MUCH FINE!

!REDEFINITION OF CONST
CONST=1.0_RPREC/(NX2*NY2)

!CX CONVECTIVE ACCELERATION TERM CALCULATION (UVP NODES)
!OPERATION: V * (DU/DY-DV/DX) + INTERP_W*(DV/DW-DW/DY)  !CHECKED TWICE ITS OK SIGNS
!SPECIAL TREATMENT JZ=1 CAUSE VORT1,VORT2 IS ON UVP (VORT3 ALWAYS UVP)
! JZ_MIN=1
! IF((.NOT.USE_MPI).OR.(USE_MPI.AND.COORD==0))THEN
!    IF((SGS.EQV..TRUE.).AND.(LBC_MOM.EQ.'WALL_LAW'))THEN
!       CC_BIG(:,:,1) = CONST*(V_BIG(:,:,1)*(-VORT3_BIG(:,:,1))      &       !UVP NODES
!                   + 0.5_RPREC*W_BIG(:,:,2)*VORT2_BIG(:,:,2))
!       JZ_MIN=2
!    ENDIF
! ENDIF
!$OMP PARALLEL DO
DO JZ=1,NZ-1       !UVP NODES SO NZ-1
   CC_BIG(:,:,JZ) = CONST*(V_BIG(:,:,JZ)*(-VORT3_BIG(:,:,JZ))           &
                  + 0.5_RPREC*(W_BIG(:,:,JZ+1)*VORT2_BIG(:,:,JZ+1)      &
                  + W_BIG(:,:,JZ)*VORT2_BIG(:,:,JZ))) !here should take away this JZ cause vort2(1) is on uvp not on w nodes!!
   CALL DFFTW_EXECUTE_DFT_R2C(FORW_BIG,CC_BIG(:,:,JZ),CC_BIG(:,:,JZ)) 
   CALL UNPADD(CX(:,:,JZ),CC_BIG(:,:,JZ))                   !UN-ZERO PAD
   CALL DFFTW_EXECUTE_DFT_C2R(BACK,CX(:,:,JZ),CX(:,:,JZ))   !BACK TO PHYS SPACE 
!END DO
!!$OMP END PARALLEL DO

!CY CONVECTIVE ACCELERATION TERM CALCULATION
!OPERATION: U*(DU/DY-DV/DX) + INTERP_W*(DV/DZ-DW/DV)
!SPECIAL TREATMENT JZ=1 CAUSE VORT1,VORT2 IS ON UVP IN THIS CASE (SO WE WOULD DO WRONG)
! JZ_MIN=1
! IF((.NOT. USE_MPI).OR.(USE_MPI.AND.COORD == 0))THEN
!    IF((SGS.EQV..TRUE.).AND.(LBC_MOM.EQ.'WALL_LAW'))THEN      
!       CC_BIG(:,:,1) = CONST*(U_BIG(:,:,1)*(VORT3_BIG(:,:,1))         & !UVP NODES
!                      + 0.5_RPREC*W_BIG(:,:,2)*(-VORT1_BIG(:,:,2)))
!       JZ_MIN=2
!    ENDIF
! ENDIF
!!$OMP PARALLEL DO
!DO JZ = 1,NZ-1       !UVP NODES SO NZ-1
   CC_BIG(:,:,JZ)=CONST*(U_BIG(:,:,JZ)*VORT3_BIG(:,:,JZ)&
        +0.5_RPREC*(W_BIG(:,:,JZ)*(-VORT1_BIG(:,:,JZ))&
        +W_BIG(:,:,JZ+1)*(-VORT1_BIG(:,:,JZ+1))))
   CALL DFFTW_EXECUTE_DFT_R2C(FORW_BIG,CC_BIG(:,:,JZ),CC_BIG(:,:,JZ))
   CALL UNPADD(CY(:,:,JZ),CC_BIG(:,:,JZ))
   CALL DFFTW_EXECUTE_DFT_C2R(BACK,CY(:,:,JZ),CY(:,:,JZ))
!END DO
!!$OMP END PARALLEL DO


!CZ CALCULATION
!if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
!  cc_big(:,:,1)=0._rprec
!  jz_min = 2
!else
!  jz_min = 1
!end if
!if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
!   cc_big(:,:,nz)=0._rprec ! according to JDA paper p.242
!end if
!OPERATION: INT_U * (DW/DX-DU/DZ) + INTERP_V*(DW/DY-DV/DW)
!HERE WE USE U_BIG(JZ=0) AND V_BIG(JZ=0)
!!$OMP PARALLEL DO
!DO JZ=1,NZ
  CC_BIG(:,:,JZ) = CONST*0.5_RPREC*(                            &
    (U_BIG(:,:,JZ)+U_BIG(:,:,JZ-1))*(-VORT2_BIG(:,:,JZ))      &
    + (V_BIG(:,:,JZ)+V_BIG(:,:,JZ-1))*VORT1_BIG(:,:,JZ))
   CALL DFFTW_EXECUTE_DFT_R2C(FORW_BIG,CC_BIG(:,:,JZ),CC_BIG(:,:,JZ))
   CALL UNPADD(CZ(:,:,JZ),CC_BIG(:,:,JZ))
   CALL DFFTW_EXECUTE_DFT_C2R(BACK,CZ(:,:,JZ),CZ(:,:,JZ))
END DO
!$OMP END PARALLEL DO


IF (VERBOSE) WRITE (*,*) 'FINISHED CONVEC'

END SUBROUTINE CONVEC

!-----------------------------------------------------------------------
! CONVEC DESCRIPTION:

! First it saves into cx,y,z the values of const*u1,2,3;
! Then it calculates the FFT of each one;
! Padds the vectors for the 3/2 dealiasing saving them into u1_big, u2_big, u3_big;
! Goes back to physical space for u1_big, u2_big, u3_big;
! 
! Saves into cx,y,z the vorticity terms es for cx: dw/dy-dv/dz;
! Then it calculates the FFT of each one;
! Padds the vectors for the 3/2 dealiasing saving them into vort1_big, vort2_big, vort3_big;
! Goes back to physical space for vort1_big, vort2_big, vort3_big;
! 
! Calculates the convective terms u_big*vort;
! Then it calculates the FFT of each one;
! Unpadd the vectors in the wavenumber space;
! Goes back in physical space for cx,cy,cz;       FINAL RESULT



END MODULE CNVC
