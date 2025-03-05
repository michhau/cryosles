!!-----------------------------------------------------------------------------------
!!  OBJECT: module sgsmodule
!!-----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 14/10/2011
!!
!!  DESCRIPTION:
!!
!!
!!
!!
!!
!!-----------------------------------------------------------------------------------
module sgsmodule
use types, only: rprec
use param, only: LD,NY,NZ,MODEL,verbose 

IMPLICIT NONE

REAL(RPREC),DIMENSION(:,:,:),ALLOCATABLE :: NU_T,     &
                                            Beta,    &
                                            U_LAG,V_LAG,W_LAG
REAL(RPREC),DIMENSION(:,:,:),ALLOCATABLE :: F_LM,F_MM,F_QN,F_NN,Betaclip, CS_OPT2
!REAL(RPREC),DIMENSION(LD,NY,0:NZ) :: NU_T

real(rprec), parameter :: opftime = 1.5_rprec
!real(rprec), dimension(ld,ny,nz) :: F_LM,F_MM,F_QN,F_NN,Beta,Betaclip
real(rprec), dimension(:), allocatable :: Beta_avg,Betaclip_avg,Beta_avg_s,Betaclip_avg_s

!xxxx----- Added by Vij - 04/14/04--xxxx----------------
! Nu_t is needed for scalar sgs
! dissip is dissipation calculated in sgs_stag and outputted when needed
! For more details look into scalars_module.f90
!real (rprec), dimension (ld,ny,$lbz:nz) :: Nu_t,dissip,Nu_s
!real(kind=rprec),dimension(ld,ny,nz)::Cs_opt2,Cs_opt2_avg,Ds_opt2
!real(kind=rprec),dimension(ld,ny,nz)::F_KX_SC1,F_XX_SC1,F_PY_SC1,F_YY_SC1,Beta_s,Betaclip_s
!real(kind=rprec),dimension(ld,ny,nz)::F_KX_SC2,F_XX_SC2,F_PY_SC2,F_YY_SC2
!xxxx---- Vij change ends here --xxxx-------------
!real(kind=rprec),dimension(ld,ny,0:nz) :: u_lag=0.0_RPREC,v_lag=0.0_RPREC,w_lag=0.0_RPREC
!real(kind=rprec),dimension(ld,ny,0:nz) :: T_lag=0.0_RPREC
!real(kind=rprec),dimension(ld,ny,$lbz:nz)::u_lag1,v_lag1,w_lag1
integer ::jt_count

!real(rprec), dimension(ld,ny,nz) :: Cs_Ssim

!-----------------------------------------------------------------------------------
CONTAINS


 SUBROUTINE ALLOCATE_SGS()
 
 IMPLICIT NONE
 
 IF(VERBOSE) PRINT*,'IN ALLOCATE_SGS'
 
 ALLOCATE(CS_OPT2(LD,NY,NZ))
 ALLOCATE(NU_T(LD,NY,0:NZ))
 
 !IF((MODEL==5).OR.(MODEL==7))THEN
    
    ALLOCATE(F_LM(LD,NY,NZ))
    ALLOCATE(F_MM(LD,NY,NZ))
    ALLOCATE(F_QN(LD,NY,NZ))
    ALLOCATE(F_NN(LD,NY,NZ))
    ALLOCATE(BETA(LD,NY,0:NZ))
    ALLOCATE(BETACLIP(LD,NY,NZ))
    ALLOCATE(U_LAG(LD,NY,0:NZ))
    ALLOCATE(V_LAG(LD,NY,0:NZ))
    ALLOCATE(W_LAG(LD,NY,0:NZ))
 
 !ENDIF
 
 ALLOCATE(Beta_avg(nz))
 ALLOCATE(Betaclip_avg(nz))
 ALLOCATE(Beta_avg_s(nz))
 ALLOCATE(Betaclip_avg_s(nz))
 
 IF(VERBOSE) PRINT*,'OUT ALLOCATE_SGS'
 
 END SUBROUTINE ALLOCATE_SGS
 
 
 SUBROUTINE DEALLOCATE_SGS()
 
 IMPLICIT NONE
 
 DEALLOCATE(CS_OPT2)
 DEALLOCATE(NU_T)
 
 !IF((MODEL==5).OR.(MODEL==7))THEN
 
    DEALLOCATE(F_LM)
    DEALLOCATE(F_MM)
    DEALLOCATE(F_QN)
    DEALLOCATE(F_NN)
    DEALLOCATE(BETA)
    DEALLOCATE(BETACLIP)
    DEALLOCATE(U_LAG)
    DEALLOCATE(V_LAG)
    DEALLOCATE(W_LAG)
 
 !ENDIF
 
 DEALLOCATE(Beta_avg)
 DEALLOCATE(Betaclip_avg)
 DEALLOCATE(Beta_avg_s)
 DEALLOCATE(Betaclip_avg_s)
 
 END SUBROUTINE DEALLOCATE_SGS



real(rprec) function rtnewt(A, jz)

use types, only: rprec

integer, parameter :: jmax=100
real(rprec) :: x1,x2,xacc
integer :: j, jz
real(rprec) :: df,dx,f
real(rprec), dimension(0:5) :: A

!setting parameters
x1 = 0._rprec
x2 = 15._rprec          !try to find the largest root first
xacc = 0.001_rprec      !doesn't need to be that accurate
rtnewt = 0.5_rprec*(x1+x2)

!Newton-Rapson loop to find root
do j=1,jmax

   !function
   f = A(0)+rtnewt*(A(1)+rtnewt*(A(2)+rtnewt*(A(3)+rtnewt*(A(4)+rtnewt*A(5)))))

   !first derivative of the function
   df = A(1) + rtnewt*(2._rprec*A(2) + rtnewt*(3._rprec*A(3) +&
        rtnewt*(4._rprec*A(4) + rtnewt*(5._rprec*A(5)))))

   !Newton-Rapson advancement step
   dx=f/df
   rtnewt = rtnewt - dx

   !check convergence
   if (abs(dx) < xacc) return
end do

!if didn't converge then put this guess value (scale sim)
rtnewt = 1._rprec
write(6,*) 'using beta=1 at jz= ', jz

end function rtnewt


END MODULE SGSMODULE
