!!-----------------------------------------------------------------------------------
!!  OBJECT: module elementar
!!-----------------------------------------------------------------------------------
!!
!!  LAST MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 23/10/2011
!!
!!  DESCRIPTION:
!!
!!  This module provides elementar useful functions and subroutines
!!  
!!
!!
!!
!!
!!
!!-----------------------------------------------------------------------------------

module elementar

use types,only:rprec

use param

use sim_param, only: u,v,w,p,DT


implicit none



!------------------------------------------------------------------------------------
CONTAINS  !--------------------------------------------------------------------------
!------------------------------------------------------------------------------------





!!-----------------------------------------------------------------------------------
!!  OBJECT: subroutine fftx(in, out)
!!-----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 10/11/2011
!!
!!  DESCRIPTION:
!!
!!  it calculates the fft on the x direction for the given variable
!!
!!-----------------------------------------------------------------------------------
subroutine fftx(in, out)
use param
use fft

implicit none

real(kind=rprec),dimension(nx) :: in
double complex,dimension(nx/2+1) :: out
integer :: i
integer*8, save :: plan
logical, save :: init = .false.

!make the plan for the FFT
if (.not. init) then
    call dfftw_plan_dft_r2c_1d(plan,nx,in,out,FFTW_ESTIMATE)
    init = .true.
end if

do i=1,nx
    in(i) = in(i)/nx
enddo

call dfftw_execute_dft_r2c(plan, in, out)

end subroutine fftx


!OTHER VERSION FOR AUTOCORRELATION
subroutine fftx2(in, out)
use param
use fft

implicit none

real(kind=rprec),dimension(2*nx) :: in
double complex,dimension(nx+1) :: out
integer :: i
integer*8, save :: plan
logical, save :: init = .false.

!make the plan for the FFT
if (.not. init) then
    call dfftw_plan_dft_r2c_1d(plan,2*nx,in,out,FFTW_ESTIMATE)
    init = .true.
end if

do i=1,2*nx
    in(i) = in(i)/2*nx
enddo

!make the plan for the FFT
call dfftw_execute_dft_r2c(plan, in, out)

end subroutine fftx2



!!-----------------------------------------------------------------------------------
!!  OBJECT: subroutine ifftx(out, in)
!!-----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 10/11/2011
!!
!!  DESCRIPTION:
!!
!!  it calculates the ifft on the x direction for the given variable
!!
!!-----------------------------------------------------------------------------------
subroutine ifftx(out, in)
use param
use fft

implicit none

real(kind=rprec),dimension(nx) :: in
double complex,dimension(nx/2+1) :: out

integer*8, save :: plan
logical, save :: init_inv = .false.

!make the plan for the FFT
if (.not. init_inv) then
    call dfftw_plan_dft_c2r_1d(plan,nx,out,in,FFTW_ESTIMATE)
    init_inv = .true.
end if

call dfftw_execute_dft_c2r(plan, out, in)

end subroutine ifftx


!OTHER VERSION FOR AUTOCORRELATION
subroutine ifftx2(out, in)
use param
use fft

implicit none

real(kind=rprec),dimension(2*nx) :: in
double complex,dimension(nx+1) :: out

integer*8, save :: plan
logical, save :: init_inv = .false.

!make the plan for the FFT
if (.not. init_inv) then
    call dfftw_plan_dft_c2r_1d(plan,2*nx,out,in,FFTW_ESTIMATE)
    init_inv = .true.
end if

call dfftw_execute_dft_c2r(plan, out, in)

end subroutine ifftx2



!!-----------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE PSD(IN, OUT)
!!-----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: MARCO GIOMETTO ( MGIOMETTO@GMAIL.COM ) ON 10/02/2013
!!
!!  DESCRIPTION:
!!
!!  IT CALCULATES THE NON NORMALIZED PSD FOR THE GIVEN SIGNAL (THE LENGTH OF THE SIGNAL
!!  MUST BE AN EVEN NUMBER.
!!
!!-----------------------------------------------------------------------------------
SUBROUTINE PSD(IN,OUT,N)

USE FFT

IMPLICIT NONE

REAL(RPREC),DIMENSION(N),INTENT(INOUT) :: IN
DOUBLE COMPLEX,DIMENSION(N/2+1) :: IN_TMP
REAL(RPREC),DIMENSION(N/2+1),INTENT(OUT) :: OUT
INTEGER,INTENT(IN) :: N

INTEGER :: I
INTEGER*8,SAVE :: PLAN
LOGICAL,SAVE :: INIT=.FALSE.

IF(MOD(SIZE(IN),2).NE.0)THEN
   PRINT*,"ATTENTION THE GIVEN SIGNAL IS NOT EVEN"
   STOP
ENDIF

!MAKE THE PLAN FOR THE FFT
IF(.NOT.INIT)THEN
    CALL DFFTW_PLAN_DFT_R2C_1D(PLAN,N,IN,IN_TMP,FFTW_ESTIMATE)
    INIT=.TRUE.
ENDIF

!NORMALIZATION FOR FFT
DO I=1,N
    IN(I)=IN(I)/N
ENDDO

!COMPUTE FFT
CALL DFFTW_EXECUTE_DFT_R2C(PLAN,IN,IN_TMP)

!PSD
DO I=1,N/2+1
    OUT(I) = REAL(REAL(IN_TMP(I)*CONJG(IN_TMP(I))),RPREC)
ENDDO

END SUBROUTINE PSD



!!-----------------------------------------------------------------------------------
!!  OBJECT: subroutine spectr_x(in, out)
!!-----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 10/11/2011
!!
!!  DESCRIPTION:
!!
!!  it calculates the energy spectrum on the x direction for the given variable
!!
!!-----------------------------------------------------------------------------------
subroutine spectr_x(in, out)
use param
use fft

implicit none

real(kind=rprec),dimension(nx) :: in
double complex,dimension(nx/2+1) :: out
integer :: i

call fftx(in, out)

!out(1) = 0.5*(out(1)*CONJG(out(1)))
! out(nx/2+1) = 0.0
do i=1,nx/2
    out(i) = (out(i)*CONJG(out(i)))
end do

end subroutine spectr_x


!OTHER VERSION FOR AUTOCORRELATION
subroutine spectr_2x(in, out)
use param
use fft

implicit none

real(kind=rprec),dimension(2*nx) :: in
double complex,dimension(nx+1) :: out
integer :: i

call fftx2(in, out)

!out(1) = 0.5*(out(1)*CONJG(out(1)))
do i=1,nx+1
    out(i) = (out(i)*CONJG(out(i)))
end do

end subroutine spectr_2x



!!-----------------------------------------------------------------------------------
!!  OBJECT: subroutine autocorr_x(in, out)
!!-----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 10/11/2011
!!
!!  DESCRIPTION:
!!
!!  it calculates the autocorreation function throuth FFT on the x direction for the given variable
!!
!!-----------------------------------------------------------------------------------
subroutine autocorr_x(in, out)
use param
use fft

implicit none

real(kind=rprec),dimension(nx) :: in
real(kind=rprec),dimension(2*nx) :: in_padd
real(kind=rprec),dimension(2*nx) :: out
double complex,dimension(nx+1) :: out_temp

call acorrpadd_x(in, in_padd)
call spectr_2x(in_padd, out_temp)
call ifftx2(out_temp, out)

end subroutine autocorr_x



!!-----------------------------------------------------------------------------------
!!  OBJECT: subroutine acorrpadd(in, out_padd)
!!-----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 10/11/2011
!!
!!  DESCRIPTION:
!!
!!  padds the in variables with zeros so the autocorrelation will be on n+1 and not just n/2+1
!!
!!-----------------------------------------------------------------------------------
subroutine acorrpadd_x(in, out_padd)
use param
use fft

real(kind=rprec),dimension(nx), intent(in) :: in
real(rprec), dimension(2*nx), intent(out) :: out_padd
integer :: jx

!make sure the big array is zeroed
out_padd = 0.0

do jx=1,nx
    out_padd(jx)=in(jx)
end do

end subroutine acorrpadd_x



end module elementar
