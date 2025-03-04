!!-----------------------------------------------------------------------------------
!!  OBJECT: subroutine filt_da ()
!!-----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 01/10/2011
!!
!!  DESCRIPTION:
!!
!!      - reads the f_c, dfdx_c, dfdy_c that are passed with dimention (nx/2+1, ny, nz-1)
!!      - calculates the Fourier coefficients along x and y on slices at different jz;
!!      - kills the odd-balls;
!!      - moltiplies the Fourier coefficients by i through the "eye" function;
!!      - multiplies the Fourier coefficients by k;
!!      - calculates the inverse FFT to go back to physical space.
!!
!!  NOTES:
!!
!!  The dimention nx/2+1 comes from how the dfftw3 is conceived.
!!
!!-----------------------------------------------------------------------------------

! kills the oddball components in f and calculates in plane derivatives
SUBROUTINE FILT_DA(F_C,DFDX_C,DFDY_C)
USE TYPES,ONLY:RPREC
USE PARAM,ONLY:LH,NX,NY,NZ,PROFILING
USE SIM_PARAM,ONLY:T1A,T1B
USE FFT

IMPLICIT NONE

INTEGER :: JZ
COMPLEX(RPREC),DIMENSION(lh,ny,0:nz),INTENT(INOUT) :: F_C,DFDX_C,DFDY_C
REAL(RPREC) :: CONST,IGNORE_ME


IF(PROFILING) CALL DATE_AND_TIME(VALUES=T1A)


const=1._rprec/(nx*ny)          !note that is not a comment!!!

!loop through horizontal slices (slice means all the x-y domain)
!$OMP PARALLEL DO 
do jz=0,nz

   !normalize all nodes at defined jz (otherwise forward and back fft increments by nx*ny the result)
   f_c(:,:,jz)=const*f_c(:,:,jz)

   !calculate the Fourier coefficients, gives forw or back and the input-output (input=output)
   call dfftw_execute_dft_r2c( forw, f_c(:,:,jz), f_c(:,:,jz) )

   !kill the oddballs inside f_c 
   f_c(lh,:,jz)=0._rprec
   f_c(:,ny/2+1,jz)=0._rprec

   !in-plane wavenumber-space derivative = molt by i every element and then by k
   dfdy_c(:,:,jz)=eye*ky(:,:)*f_c(:,:,jz)
   dfdx_c(:,:,jz)=eye*kx(:,:)*f_c(:,:,jz)

   !inverse transform 
   call dfftw_execute_dft_c2r( back, f_c(:,:,jz), f_c(:,:,jz) )    
   call dfftw_execute_dft_c2r( back, dfdx_c(:,:,jz), dfdx_c(:,:,jz) )
   call dfftw_execute_dft_c2r( back, dfdy_c(:,:,jz), dfdy_c(:,:,jz) )      

end do
!$OMP END PARALLEL DO


IF(PROFILING) CALL DATE_AND_TIME(VALUES=T1B)

end subroutine filt_da
