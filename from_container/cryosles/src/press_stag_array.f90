!!----------------------------------------------------------------------------------
!  OBJECT: subroutine press_stag_array (p_hat, dfdx, dfdy)
!!-----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 06/10/2011
!!
!!  DESCRIPTION:
!!
!!
!!
!!
!!
!!-----------------------------------------------------------------------------------

subroutine press_stag_array (p_hat,dfdx,dfdy)

!p_hat contains the physical space pressure on exit
!provides p_hat 0:nz-1, dfdx, dfdy 1:nz-1
!----------------------------------------------------------------------------
! Boundary Layer version with 4th order derivs in vertical.
! 04 December 1995
! Modifications:
! 12/6:  added upper and lower boundary conditions.     
! 12/8:  Corrected sign on imag. x and y deriv of RHS.
! 12/17: Forcing pressure equal zero at wall 
!        Removed forcing of <P>=0 for all z.
!        Will need to change BC when hetero. surface stress.
! 12/17: Added ficticious node below wall (Note solver on Nz+1 x Nz+1)
! 12/18: Revised 2st deriv stencil at wall (avg of deriv at -1 and at 1)
! 12/21: Redid FDD to 2nd order accurate.
! 12/22  Now value of P(wall) diagnosed from prev P, to match gradient BC
! 1/13:  major changes! Broke out mean pressure for separate solution
! 1/20:  back to tridag for natrix solution (same sol'n as LUDCMP...A Keeper!)
! 1/23:  Staggered solution
!        Using Nz+1 levels for computing P, but tossing out level below ground
! 4/1:   Changed sign on Div T_iz at wall and lid (five places)
! 16 March 2011
! Mods by JF:
! 3/16:  added fx,fy,fz to rH for iterative sol. of IB force and pressure
! 3/16:  changed bottom and top bc. using w^(*)
!----------------------------------------------------------------------------
use types, only: rprec
use param
use sim_param, only: u, v, w, DT
use fft
use service,only:Z_MPROFILE
use stretch

implicit none

complex(rprec), dimension(lh,ny,0:nz) :: p_hat
real(rprec), dimension(ld,ny,0:nz_max) :: rH_x, rH_y, rH_z
complex(rprec), dimension(lh,ny,0:nz_max) :: H_x, H_y, H_z
equivalence (rH_x,H_x),(rH_y,H_y),(rH_z,H_z) ! equivalence statement only works with dimensions declared as parameters
real(rprec), dimension(ld,ny) :: rtopw, rbottomw
complex(rprec), dimension(lh,ny) :: topw, bottomw
equivalence (rtopw,topw),(rbottomw,bottomw)
complex(rprec),dimension(lh,ny,nz),intent(out) :: dfdx, dfdy

real(rprec) :: const, ignore_me
integer :: jx, jy, jz, k, jz_min
integer :: kk
character(64) :: fname
logical, parameter :: DEBUG = .false.

complex(rprec), dimension(lh,ny,nz+1) :: RHS_col
real(rprec), dimension(lh,ny,nz+1) :: a, b, c
REAL(RPREC) :: NRM,SCR1,SCR2
real(rprec) :: f1,f2,f3
!---------------------------------------------------------------------

if (VERBOSE) write (*, *) 'started press_stag_array'

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  p_hat(:, :, 0) = (0._rprec, 0._rprec)
else
  p_hat(:, :, 0) = BOGUS
end if

!=====================================================================
!get the right hand side ready
!loop over levels

!ADIMENTIONALIZATION CONSTANT
NRM=1.0_RPREC/(NX*NY)
SCR1=NRM/TADV1
SCR2=SCR1/DT
!$OMP PARALLEL DO
DO JZ=1,NZ-1

   !TEMP STORAGE FOR SUM OF RHS TERMS NORMALIZED FOR FFT
   !THE OLD TIMESTEP ALREADY CONTAINS THE PRESSURE TERM
   RH_X(:,:,JZ) = SCR2*U(:,:,JZ)    !UVP-NODE
   RH_Y(:,:,JZ) = SCR2*V(:,:,JZ)    !UVP-NODE
   RH_Z(:,:,JZ) = SCR2*W(:,:,JZ)    !W-NODE

   !FFT
   CALL DFFTW_EXECUTE_DFT_R2C( FORW, RH_X(:,:,JZ), RH_X(:,:,JZ) )
   CALL DFFTW_EXECUTE_DFT_R2C( FORW, RH_Y(:,:,JZ), RH_Y(:,:,JZ) )
   CALL DFFTW_EXECUTE_DFT_R2C( FORW, RH_Z(:,:,JZ), RH_Z(:,:,JZ) )

ENDDO
!$OMP END PARALLEL DO

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
rbottomw(:, :) = const / tadv1 * (w(:, :, 1)/dt)
call dfftw_execute_dft_r2c( forw, rbottomw(:,:), rbottomw(:,:) )
H_z(:, :, 1) = bottomw(:, :)  ! JF added              ! MG: isn't it dangerous as it was a pointer?
end if

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
rtopw(:, :) = const / tadv1 * (w(:, :, nz)/dt)     !I have the nz as I passed the velocities before
call dfftw_execute_dft_r2c( forw, rtopw(:,:), rtopw(:,:) )
H_z(:, :, nz) = topw(:, :)  !JF added
end if

!kill all the oddballs
!$OMP PARALLEL DO
do jz = 1, nz
  H_x(lh, :, jz)=0._rprec
  H_y(lh, :, jz)=0._rprec
  H_z(lh, :, jz)=0._rprec
  H_x(:, ny/2+1, jz)=0._rprec
  H_y(:, ny/2+1, jz)=0._rprec
  H_z(:, ny/2+1, jz)=0._rprec
end do
!$OMP END PARALLEL DO

!with MPI topw and bottomw are only on top & bottom processes
!kill top-bottom oddballs
topw(lh, :)=0._rprec
topw(:, ny/2+1)=0._rprec
bottomw(lh, :)=0._rprec
bottomw(:, ny/2+1)=0._rprec

!==========================================================================
! Loop over (Kx,Ky) to solve for Pressure amplitudes

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
  !bottom nodes (is at the bottom to obtain the pressure p_0)
  a(:, :, 1) = BOGUS                            !was 0._rprec
  b(:, :, 1) = -1._rprec
  c(:, :, 1) = 1._rprec
  RHS_col(:, :, 1) = bottomw(:, :) * get_z_local(1) * 2.0_rprec         !JF: no - sign as above
end if

jz_min = 2

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) then
  !top nodes (nz+1 is actually nz as jz_min before was set to 2)
  a(:, :, nz+1) = -1._rprec          !(MG: points at nz not nz+1 because of numbering purposes)
  b(:, :, nz+1) = 1._rprec
  c(:, :, nz+1) = BOGUS  !was 0._rprec
  RHS_col(:, :, nz+1) = topw(:, :) * (get_z_local(2*nz_tot-1) - get_z_local(2*nz_tot-3)) !fills the last element of the array with topw*dz
end if

#IF(MPI)
  !--send 1->nz for H_z only
  call mpi_sendrecv (H_z(1, 1, 1), lh*ny, MPI_CPREC, down, 9,  &
                     H_z(1, 1, nz), lh*ny, MPI_CPREC, up, 9,   &
                     comm, status, ierr)

  call mpi_sendrecv (H_z(1, 1, nz-1), lh*ny, MPI_CPREC, up, 9,  &
                     H_z(1, 1, 0), lh*ny, MPI_CPREC, down, 9,   &
                     comm, status, ierr)


#ENDIF

!this is the calculation of the second derivative on x,y (see coeff k**2 in b) and z (coefficients a,b,c)
!!$OMP PARALLEL DO
do jz = jz_min, nz         !MG: its from 2 to nz but see in the RHS_col i put values jz-1 so indexes at the end are from 1 to nz-1)
  do jy = 1, ny

    if (jy == ny/2 + 1) cycle           !this skips the oddballs which are killed later in the p_hat
    do jx=1, lh-1                       !skipping the oddball at lh
      !if (jx*jy == 1) cycle

      kk = jz + eta_displ(coord)
      f3 = get_dz_local(2*kk-2) !exp(kk*dz)*C_STRETCH
      f1 = get_dz_local(2*kk-4) !exp((kk-1)*dz)*C_STRETCH
      f2 = get_dz_local(2*kk-3) ! exp((kk-0.5_rprec)*dz)*C_STRETCH

      a(jx, jy, jz) = (1.0_rprec/(dz*dz)) * (1.0_rprec/f2) * (1.0_rprec/f1)
      b(jx, jy, jz) = -(kx(jx, jy)**2 + ky(jx, jy)**2) - (1.0_rprec/(dz*dz))*(1.0_rprec/f2)*((1.0_rprec/f1) + (1.0_rprec/f3))
      c(jx, jy, jz) = (1.0_rprec / (dz*dz)) * (1.0_rprec/f2) * (1.0_rprec/f3)

      RHS_col(jx, jy, jz) = eye * (kx(jx, jy) * H_x(jx, jy, jz-1) +   &    !values at jz-1
                                   ky(jx, jy) * H_y(jx, jy, jz-1)) +  &
                            ((H_z(jx, jy, jz) - H_z(jx, jy, jz-1))/dz)*(1.0_rprec/f2)     !W interp with level above
    end do
  end do
end do
!!$OMP END PARALLEL DO

!JF: set p_hat(1,1,0) = 0
!first component of Fourier its the 0 frequency
if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
a(1, 1, 1) = 0._rprec
b(1, 1, 1) = 1._rprec
c(1, 1, 1) = 0._rprec
RHS_col(1, 1, 1) = (0._rprec, 0._rprec)
end if


!JF: this does not skip zero wavenumber solution, nyquist freqs
#IF(MPI)
  call tridag_array_pipelined (20, a, b, c, RHS_col, p_hat)
#ELSE
  call tridag_array (a, b, c, RHS_col, p_hat)
#ENDIF



#IF(MPI)
  !--make sure 0 <-> nz-1 are syncronized
  !-- 1 <-> nz should be in sync already
  call mpi_sendrecv (p_hat(1, 1, nz-1), lh*ny, MPI_CPREC, up, 10,  &
                     p_hat(1, 1, 0), lh*ny, MPI_CPREC, down, 10,   &
                     comm, status, ierr)

  call mpi_sendrecv (p_hat(1, 1, 1), lh*ny, MPI_CPREC, down, 10,  &
                     p_hat(1, 1, nz), lh*ny, MPI_CPREC, up, 10,   &
                     comm, status, ierr)

#ENDIF

!zero the nyquist freqs
p_hat(lh, :, :) = 0._rprec
p_hat(:, ny/2+1, :) = 0._rprec


!now need to get p_hat(wave,level) to physical p(jx,jy,jz), loop over height levels
call dfftw_execute_dft_c2r( back, p_hat(:,:,0), p_hat(:,:,0) )        
!IFFT for all remaining levels
!$omp parallel do
do jz=1,nz-1            !used to be nz
   do jy=1,ny
   do jx=1,lh
      !complex (both dfdx and p_hat were defined as complex)
      dfdx(jx,jy,jz)=eye*kx(jx,jy)*p_hat(jx,jy,jz)
      dfdy(jx,jy,jz)=eye*ky(jx,jy)*p_hat(jx,jy,jz)
      !note the oddballs of p_hat are already 0, so we should be OK here
   end do
   end do
   call dfftw_execute_dft_c2r( back, dfdx(:,:,jz), dfdx(:,:,jz) )    
   call dfftw_execute_dft_c2r( back, dfdy(:,:,jz), dfdy(:,:,jz) )    
   call dfftw_execute_dft_c2r( back, p_hat(:,:,jz), p_hat(:,:,jz) )  
enddo
!$omp end parallel do

#IF(MPI)
  !--make sure 0 <-> nz-1 are syncronized
  !-- 1 <-> nz should be in sync already
  !call mpi_sendrecv (p_hat(1, 1, nz-1), lh*ny, MPI_CPREC, up, 10,  &
  !                   p_hat(1, 1, 0), lh*ny, MPI_CPREC, down, 10,   &
  !                   comm, status, ierr)

  !call mpi_sendrecv (p_hat(1, 1, 1), lh*ny, MPI_CPREC, down, 10,  &
  !                   p_hat(1, 1, nz), lh*ny, MPI_CPREC, up, 10,   &
  !                   comm, status, ierr)

#ENDIF



!nz level is not needed elsewhere (although its valid)
dfdx(:, :, nz) = BOGUS
dfdy(:, :, nz) = BOGUS
!p_hat(:, :, nz) = p_hat(:,:,nz-1)

if (VERBOSE) write (*, *) 'finished press_stag_array'

end subroutine press_stag_array
