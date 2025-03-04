module test_filtermodule
use types, only: rprec
use param, only: lh,ny

private lh,ny
!TS Truely grid refinement test needs to keep the filter_size
!TS the same as that in coarse grid (Double grid size: filter_size=2. etc)

end module test_filtermodule            !@@@@@@@@@THE MODULE HERE ENDS

!THIS IS NOT IN THE MODULE!

SUBROUTINE TEST_FILTER(F,G_TEST)
!NOTE: THIS FILTERS IN-PLACE, SO INPUT IS RUINED
USE TYPES, ONLY: RPREC
USE PARAM, ONLY: LH,NY
USE FFT

IMPLICIT NONE

!NOTE WE'RE TREATING AS COMPLEX HERE
COMPLEX(RPREC), DIMENSION(LH,NY), INTENT(INOUT) ::F
REAL(RPREC), DIMENSION(LH,NY), INTENT(IN) :: G_TEST
REAL(RPREC) :: IGNORE_ME


CALL DFFTW_EXECUTE_DFT_R2C( FORW, F, F )     

!THE NORMALIZATION IS "IN" G_TEST
F = G_TEST*F  !NYQUIST TAKEN CARE OF WITH G_TEST


CALL DFFTW_EXECUTE_DFT_C2R( BACK, F, F )     

END SUBROUTINE TEST_FILTER



!!-----------------------------------------------------------------------------------
!!  OBJECT: test_filter_init (alpha,G_test)
!!-----------------------------------------------------------------------------------
!!
!!  LAST MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 30/09/2011
!!
!!  DESCRIPTION:
!!
!!  This subroutine has different behaviour depending on the filter tipe (1,2,3).
!!  It creates a matrix "G_test" of rank 2 with (lh x ny) extensions and fills in it 
!!  the values of k2, previously calculated in the subroutine "init_wavenumber (fft.f90)", 
!!  opportunely modified (reduced) depending on which filter is used:
!!
!!      - with the spectral cut-off filter it assumes G_test = 0 where k2 > than the kc2 defined by the filter size;
!!      - with gaussian filter it assumes G_test to be a Gaussian function of k2 (transfer function);
!!      - with top hat filter  it assumes G_test to be a Sinc function of k2 (transfer function);
!!
!!-----------------------------------------------------------------------------------
subroutine test_filter_init (ALPHA,G_TEST,IFILTER)
! spectral cutoff filter at width alpha*delta
! note the normalization for FFT's is already in G! (see 1/(nx*ny))

use types, only: rprec
use param, only: lh,nx,ny,dx,dy,pi,model
use fft
implicit none
real(rprec) :: alpha, delta, kc2
real(rprec), dimension(lh,ny) :: G_test
INTEGER,INTENT(IN) :: IFILTER
G_test=1._rprec/(nx*ny)         !normalization for the forward FFT
delta = alpha*sqrt(dx*dy)       !2D-delta, not full 3D one, and assumes a mean value between dx and dy

!spectral cutoff filter case
if(IFILTER==1) then
  kc2 = (pi/(delta))**2                         !(max kc)^2 with alpha delta filter size
  where (k2 >= kc2) G_test = 0._rprec           !k2 = kx*kx + ky*ky from fft.f90

  !little check
  !if (model==6.OR.model==7) then
  !  print *, 'Use Gaussian or Top-hat filter for mixed models'
  !  stop
  !indif

!Gaussian filter case
else if(IFILTER==2) then
  G_test=exp(-(2._rprec*delta)**2*k2/(4._rprec*6._rprec))*G_test

!top-hat filter case
else if(IFILTER==3) then
  G_test= (sin(kx*delta/2._rprec)*sin(ky*delta/2._rprec)+1E-8)/&
           (kx*delta/2._rprec*ky*delta/2._rprec+1E-8)*G_test
endif

!since our k2 has zero at nyquist, we have to do this by hand
G_test(lh,:) = 0._rprec
G_test(:,ny/2+1) = 0._rprec

end subroutine test_filter_init
