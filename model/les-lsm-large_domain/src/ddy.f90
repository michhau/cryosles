subroutine ddy (dfdy, f)
! f remains untouched
use types, only: rprec
use param, only: lh, nx, ny, nz
use fft

implicit none
integer :: jz
complex(rprec), dimension(lh,ny,0:nz), intent(in) :: f
complex(rprec), dimension(lh,ny,0:nz), intent(inout) :: dfdy
real(rprec) :: const, ignore_me

CONST=1._RPREC/(NX*NY)

! Loop through horizontal slices
!$OMP PARALLEL DO
DO JZ=1,NZ  !--CAN BE NZ-1 EXCEPT FOR DTYZDY OF THE TOP PROCESS

   ! TEMPORARY STORAGE IN DFDX
   DFDY(:,:,JZ)=CONST*F(:,:,JZ)

   CALL DFFTW_EXECUTE_DFT_R2C( FORW, DFDY(:,:,JZ), DFDY(:,:,JZ) )

   ! ODDBALLS
   DFDY(LH,:,JZ)=0._RPREC
   DFDY(:,NY/2+1,JZ)=0._RPREC
   
   ! COMPUTE COEFFICIENTS FOR PSEUDOSPECTRAL DERIVATIVE CALCULATION.
   DFDY(:,:,JZ)=EYE*KY(:,:)*DFDY(:,:,JZ)
   ! INVERSE TRANSFORM TO GET PSEUDOSPECTRAL DERIVATIVE.

   CALL DFFTW_EXECUTE_DFT_C2R( BACK, DFDY(:,:,JZ), DFDY(:,:,JZ) )    

END DO
!$END OMP PARALLEL DO


end subroutine ddy
