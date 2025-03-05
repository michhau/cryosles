!---------------------------------------------------------------------------
!fftw (both 2.1.X and 3.2.X versions)
!---------------------------------------------------------------------------

module fft
use types, only: rprec
use param, only: lh, ny
implicit none

!save means that what happens here to the variables is permanent through calls
save

public :: kx, ky, k2, eye, forw, back, forw_big, back_big, init_fft

!plans for FFTW3, integer*8 as recommended in the manual
integer*8 :: forw, back, forw_big, back_big

!other variables
real(rprec), dimension(lh,ny) :: kx, ky, k2

!flag variables for FFTW3
integer :: iflag, iflag2

!note that is a COMPLEX parameter
complex(rprec), parameter :: eye = (0._rprec,1._rprec)

include "fftw3.f"


!---------------------------------------------------------------------------
contains  !-----------------------------------------------------------------
!---------------------------------------------------------------------------

!!-----------------------------------------------------------------------------------
!!  OBJECT: subroutine init_fft ()
!!-----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 30/09/2011
!!
!!  DESCRIPTION:
!!
!!      -calls the fft subroutine to compute the Fourier coefficients for:
!!      -forward transform f(x_i) -> a(k_i)
!!      -inverse transform a(k_i) -> f(x_i)
!!      -3/2 forward
!!      -3/2 backward
!!      -calls the init_wavenumber to fill and renormalize kx,ky,k2!!
!!
!!-----------------------------------------------------------------------------------
subroutine init_fft ()
use param, only: nx, ny, nx2, ny2, lh, lh_big, ld, ld_big
implicit none

integer :: ierr, iret

!creates the output vectors of size 2(nx/2+1)
real(rprec) :: rarr1(ld,ny), rarr2(ld_big,ny2)

!comment out the following lines to use the multi-threaded fftw and specify nthreads
!call dfftw_init_threads(iret)
!call dfftw_plan_with_nthreads(8)
!print *,"iret =      ",iret

!define flags
iflag  = FFTW_MEASURE + FFTW_UNALIGNED + FFTW_DESTROY_INPUT
iflag2 = FFTW_MEASURE + FFTW_UNALIGNED + FFTW_DESTROY_INPUT

!formulate the fft plans (forw, n0, n1, double in, complex out, flags)
call dfftw_plan_dft_r2c_2d( forw, nx, ny, rarr1, rarr1, iflag )
call dfftw_plan_dft_c2r_2d( back, nx, ny, rarr1, rarr1, iflag )
call dfftw_plan_dft_r2c_2d( forw_big, nx2, ny2, rarr2, rarr2, iflag2 )
call dfftw_plan_dft_c2r_2d( back_big, nx2, ny2, rarr2, rarr2, iflag2 )
call init_wavenumber ()

end subroutine init_fft


!!-----------------------------------------------------------------------------------
!!  OBJECT: subroutine init_wavenumber ()
!!-----------------------------------------------------------------------------------
!!
!!  LAST MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 30/09/2011
!!
!!  DESCRIPTION:
!!
!!  This subroutine considers 3 matrices for kx, ky and k2 of dimension (nx/2 x ny)
!!  each filled with vales of k from 0 to nx/2, putting the nyquist value nx/2+1 = 0 for
!!  the first derivative in both x and y direction. The operations are:
!!
!!      -fills the kx,ky,k2 matrices with their values;
!!      -adds at nx,y/2+1 the nyquist values in kx, ky = 0 for the derivative reason;
!!      -renormalizes the wavenumbers to take into account the real domain dimention (and not 2pix2pi)
!!
!!  EXAMPLE OF KX FOR NX=NY=6 ( which means kc = +-3 ):
!!
!!  0 1 2         0 1 2 0          0 1 2 0
!!  0 1 2         0 1 2 0          0 1 2 0
!!  0 1 2   -->   0 1 2 0   -->    0 1 2 0
!!  0 1 2         0 1 2 0          0 0 0 0
!!  0 1 2         0 1 2 0          0 1 2 0
!!  0 1 2         0 1 2 0          0 1 2 0
!!
!!-----------------------------------------------------------------------------------
subroutine init_wavenumber ()
use param,only: lh,nx,ny,L_x,L_y,pi       !lh=nx/2+1
implicit none
integer :: jx,jy

!fills kx starting from jx=0 to jx=nx/2 (the freq nx/2+1 is filled with 0 later)
do jx=1,lh-1
   kx(jx,:) = real(jx-1,kind=rprec)
end do

do jy=1,ny
   ky(:,jy) = real(modulo(jy - 1 + ny/2,ny) - ny/2, kind=rprec)
end do

!eliminates the nyquist value to calculate a correct first derivative later
  kx(lh,:)=0._rprec
  ky(lh,:)=0._rprec
  kx(:,ny/2+1)=0._rprec
  ky(:,ny/2+1)=0._rprec

!renormalization for the aspect ratio change
  kx=2._rprec*pi/L_x*kx
  ky=2._rprec*pi/L_y*ky

!magnitude squared:
  k2 = kx*kx + ky*ky

end subroutine init_wavenumber


end module fft
