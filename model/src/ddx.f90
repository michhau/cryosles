!!------------------------------------------------------------------------------
!!  OBJECT: subroutine ddx (dfdx, f)
!!------------------------------------------------------------------------------
!!
!!  LAST MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 30/09/2011
!!
!!  DESCRIPTION:
!!
!!  This subroutine computes the x derivative dfdx[lh,ny,nz] of a function f[lh,ny,nz].
!!  It loops on horizontal slices calculating the FFT mult by ik and IFFT back physical space.
!!
!!
!!
!!
!!
!!------------------------------------------------------------------------------
subroutine ddx (dfdx, f)        !f remains untouched

use types, only : rprec
use param, only : lh, nx, ny, nz
use fft

implicit none

integer :: jz


!f and dfdx are COMPLEX arrays as they need to go through the FFT
complex(rprec), dimension(lh,ny,0:nz), intent(in) :: f
complex(rprec), dimension(lh,ny,0:nz), intent(inout) :: dfdx

real(rprec) :: const, ignore_me

!normalization constant
CONST=1._RPREC/(NX*NY)

!$OMP PARALLEL DO
DO JZ=1,NZ

   !TEMPORARY STORAGE IN DFDX
   DFDX(:,:,JZ)=CONST*F(:,:,JZ)

   !FFT TO GO IN WAVENUMBER SPACE (AND COMPUTE THERE THE DERIVATIVES TO HAVE HIGH ACCURACY)
   CALL DFFTW_EXECUTE_DFT_R2C( FORW, DFDX(:,:,JZ), DFDX(:,:,JZ) )


   !KILLS ODDBALLS
   DFDX(LH,:,JZ)=0._RPREC
   DFDX(:,NY/2+1,JZ)=0._RPREC

   !COMPUTE COEFFICIENTS FOR PSEUDOSPECTRAL DERIVATIVE CALCULATION
   DFDX(:,:,JZ)=EYE*KX(:,:)*DFDX(:,:,JZ)

   !INVERSE TRANSFORM TO GET PSEUDOSPECTRAL DERIVATIVE (VALUES BACK IN PHYSICAL SPACE)
   CALL DFFTW_EXECUTE_DFT_C2R( BACK, DFDX(:,:,JZ), DFDX(:,:,JZ) )      

ENDDO
!$OMP END PARALLEL DO

end subroutine ddx
