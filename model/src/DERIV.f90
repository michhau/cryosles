MODULE DERIV

USE TYPES
USE PARAM
USE SERVICE
use stretch
IMPLICIT NONE


!##############################################################################
CONTAINS


!!------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE DDZ_UV (DFDZ, F)
!!------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 5/10/2011
!!
!!  DESCRIPTION:
!!
!!  This sub computes the z derivative dfdz[ld,ny,nz] of a function f[lh,ny,nz].
!!  It reads the value of f at jz and jz-1 and calc dfdz as 1/jz*[f(jz)-f(jz-1)]
!!
!!
!!------------------------------------------------------------------------------
SUBROUTINE DDZ_UV(DFDZ,F)


IMPLICIT NONE

!F AND DFDX ARE REAL ARRAYS AS THEY DON'T NEED TO GO THROUGH THE FFT
REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(IN) :: F
REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(INOUT) :: DFDZ
INTEGER :: K,KK
REAL(RPREC) :: CONST


IF(VERBOSE) PRINT*,'IN DDZ_UV'

!VERTICAL DISTANCE NEEDED TO CALCULATE DERIVATIVES
!OPERATE FROM 0-NZ UVP --> 0-NZ W
!$OMP PARALLEL DO
DO K=1,NZ
     KK = K + eta_displ(coord)
     CONST= 1.0_rprec/(get_dz_local(2*kk-2))! (1._RPREC/dz)*(1.0_rprec/dexp(KK*DZ)) * (1.0_rprec/c_stretch)
   DFDZ(:,:,K)=CONST*(F(:,:,K)-F(:,:,K-1))   !UVP I TAKE ADV OF LEVEL 0
ENDDO
!$OMP END PARALLEL DO

IF(VERBOSE) PRINT*,'OUT DDZ_UV'

END SUBROUTINE DDZ_UV



!!------------------------------------------------------------------------------
!!  OBJECT: subroutine ddz_w (dfdz, f)
!!------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 5/10/2011
!!
!!  DESCRIPTION:
!!
!!  This subroutine computes the first derivative for the w variable on z direction
!!
!!
!!------------------------------------------------------------------------------
SUBROUTINE DDZ_W (DFDZ,F)

IMPLICIT NONE

!F AND DFDX ARE REAL ARRAYS AS THEY DON'T NEED TO GO THROUGH THE FFT
REAL(RPREC), DIMENSION(LD,NY,0:NZ), INTENT(IN) :: F
REAL(RPREC), DIMENSION(LD,NY,0:NZ), INTENT(OUT) :: DFDZ
REAL(RPREC) :: CONST
INTEGER :: K,KK

IF(VERBOSE) PRINT*,'IN DDZ_W'

!$OMP PARALLEL DO
DO K=0,NZ-1    !1 NOT 0 CAUSE ANYWAY THE 0 FOR UVP IS BELOW COMP. DOMAIN
   KK = K + eta_displ(coord)
   CONST=1.0_rprec/get_dz_local(2*kk-1)   !(1.0_RPREC/DZ)*(1.0_rprec/exp((real(KK,rprec)+0.5_rprec)*DZ)) * (1.0_rprec/c_stretch)
   DFDZ(:,:,K)=CONST*(F(:,:,K+1)-F(:,:,K))
ENDDO
!$OMP END PARALLEL DO

IF(VERBOSE) PRINT*,'IN DDZ_W'

END SUBROUTINE DDZ_W



!!------------------------------------------------------------------------------
!!  OBJECT: subroutine ddx (dfdx, f)
!!------------------------------------------------------------------------------
!!
!!  LAST MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 30/09/2011
!!
!!  DESCRIPTION:
!!
!!  This subroutine computes the x der dfdx[lh,ny,nz] of a function f[lh,ny,nz].
!!  oriz slices calc the FFT mult by ik and IFFT back physical space.
!!
!!------------------------------------------------------------------------------
! subroutine ddx (dfdx, f)        !f remains untouched
! 
! use fft
! 
! implicit none
! 
! integer :: jz
! 
! 
! !f and dfdx are COMPLEX arrays as they need to go through the FFT
! REAL(RPREC),DIMENSION(:,:,0:),INTENT(IN) :: F
! REAL(RPREC),DIMENSION(:,:,0:),INTENT(OUT) :: DFDX
! 
! !complex(rprec), dimension(lh,ny,0:nz), intent(in) :: f
! !complex(rprec), dimension(lh,ny,0:nz), intent(out) :: dfdx
! 
! real(rprec) :: const, ignore_me
! 
! !normalization constant
! CONST=1._RPREC/(NX*NY)
! 
! !$OMP PARALLEL DO
! DO JZ=1,NZ
! 
!    !TEMPORARY STORAGE IN DFDX
!    DFDX(:,:,JZ)=CONST*F(:,:,JZ)
! 
!    !FFT TO GO IN WAVENUMBER SPACE (AND COMPUTE THERE THE DERIVATIVES TO HAVE HIGH ACCURACY)
!    CALL DFFTW_EXECUTE_DFT_R2C( FORW, DFDX(:,:,JZ), DFDX(:,:,JZ) )
! 
! 
!    !KILLS ODDBALLS
!    DFDX(LH,:,JZ)=0._RPREC
!    DFDX(:,NY/2+1,JZ)=0._RPREC
! 
!    !COMPUTE COEFFICIENTS FOR PSEUDOSPECTRAL DERIVATIVE CALCULATION
!    DFDX(:,:,JZ)=EYE*KX(:,:)*DFDX(:,:,JZ)
! 
!    !INVERSE TRANSFORM TO GET PSEUDOSPECTRAL DERIVATIVE (VALUES BACK IN PHYSICAL SPACE)
!    CALL DFFTW_EXECUTE_DFT_C2R( BACK, DFDX(:,:,JZ), DFDX(:,:,JZ) )      
! 
! ENDDO
! !$OMP END PARALLEL DO
! 
! end subroutine ddx
! 
! 
! !!------------------------------------------------------------------------------
! !!  OBJECT: subroutine ddy (dfdy, f)
! !!------------------------------------------------------------------------------
! !!
! !!  LAST MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 30/09/2011
! !!
! !!  DESCRIPTION:
! !!
! !!  This subroutine computes the x der dfdx[lh,ny,nz] of a function f[lh,ny,nz].
! !!  oriz slices calc the FFT mult by ik and IFFT back physical space.
! !!
! !!------------------------------------------------------------------------------
! subroutine ddy (dfdy, f)
! ! f remains untouched
! use fft
! 
! implicit none
! integer :: jz
! !complex(rprec), dimension(lh,ny,0:nz), intent(in) :: f
! !complex(rprec), dimension(lh,ny,0:nz), intent(inout) :: dfdy
! REAL(RPREC),DIMENSION(:,:,0:),INTENT(IN) :: F
! REAL(RPREC),DIMENSION(:,:,0:),INTENT(OUT) :: DFDY
! real(rprec) :: const, ignore_me
! 
! CONST=1._RPREC/(NX*NY)
! 
! ! Loop through horizontal slices
! !$OMP PARALLEL DO
! DO JZ=1,NZ  !--CAN BE NZ-1 EXCEPT FOR DTYZDY OF THE TOP PROCESS
! 
!    ! TEMPORARY STORAGE IN DFDX
!    DFDY(:,:,JZ)=CONST*F(:,:,JZ)
! 
!    CALL DFFTW_EXECUTE_DFT_R2C( FORW, DFDY(:,:,JZ), DFDY(:,:,JZ) )
! 
!    ! ODDBALLS
!    DFDY(LH,:,JZ)=0._RPREC
!    DFDY(:,NY/2+1,JZ)=0._RPREC
!    
!    ! COMPUTE COEFFICIENTS FOR PSEUDOSPECTRAL DERIVATIVE CALCULATION.
!    DFDY(:,:,JZ)=EYE*KY(:,:)*DFDY(:,:,JZ)
!    ! INVERSE TRANSFORM TO GET PSEUDOSPECTRAL DERIVATIVE.
! 
!    CALL DFFTW_EXECUTE_DFT_C2R( BACK, DFDY(:,:,JZ), DFDY(:,:,JZ) )    
! 
! END DO
! !$END OMP PARALLEL DO
! 
! 
! end subroutine ddy



!!--------------------------------------------------------------------------------
!!  OBJECT: subroutine divstress (rhs,tx,ty,tz)
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 15/05/2013
!!
!!  DESCRIPTION:
!!
!!  Provides divt, jz=1:nz-1
!!
!!  NOTE: THIS WAY COULD PARALLELIZE THIS 3 DERIVATIVE CALCULATION AT TH
!!        PRICE OF A LITTLE MORE MEMORY, ALTERNATIVE USE A SCR FOR DTXDX
!!        ETC AND UPDATE DIRECTLY THE RHS
!!
!!  NOTE2: WE COMPUTE THE STRESS DIVERGENCE FROM 1 TO NZ-1 SINCE FOR UVP
!!         THE NZ LEVEL IS OUTSIDE AND SINCE WE PRESCRIBE A BC FOR W(NZ)
!!
!!--------------------------------------------------------------------------------
SUBROUTINE DIVSTRESS (RHS,TX,TY,TZ,STAG)

USE TYPES, ONLY: RPREC
USE PARAM, ONLY: LD,NY,NZ,BOGUS,VERBOSE


IMPLICIT NONE

REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(INOUT) :: RHS
REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(IN) :: TX,TY,TZ
REAL(RPREC),DIMENSION(LD,NY,0:NZ) :: DTXDX,DTYDY,DTZDZ!,TEST
LOGICAL,INTENT(IN) :: STAG
INTEGER I

IF(VERBOSE) WRITE (*,*) 'STARTED DIVSTRESS_UV'

!COMPUTE STRESS GRADIENTS:
CALL DDX(DTXDX,TX)
CALL DDY(DTYDY,TY)
IF(STAG) CALL DDZ_W(DTZDZ,TZ)
IF(.NOT.STAG) CALL DDZ_UV(DTZDZ,TZ)

!$OMP PARALLEL DO
DO I=0,NZ
   RHS(:,:,I) = RHS(:,:,I)+DTXDX(:,:,I)+DTYDY(:,:,I)+DTZDZ(:,:,I)
ENDDO
!$OMP END PARALLEL DO

IF(VERBOSE) WRITE (*,*) 'FINISHED DIVSTRESS_UV'

END SUBROUTINE DIVSTRESS


END MODULE DERIV

