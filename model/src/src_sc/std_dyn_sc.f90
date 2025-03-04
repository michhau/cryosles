!!-----------------------------------------------------------------------------------
!!  OBJECT: subroutine std_dynamic(CS_Z,S11,S12,S13,S22,S23,S33)
!!-----------------------------------------------------------------------------------
!!
!!  LAST MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 14/10/2011
!!
!!  DESCRIPTION:
!!
!!  this is the w-node version
!!  provides CS_Z 1:nz
!!  this subroutine computes the value of Cs for every node in the simulation
!!  assuming scale invariance and averaging over planes the Cs coefficient (x-y)
!!
!!
!!-----------------------------------------------------------------------------------
!MPI: requires u,v 0:nz, except bottom process only 1:nz
subroutine std_dyn_sc(s,idsdx,idsdy,idsdz,s11,s12,s13,s22,s23,s33,cs_z)

!can save more mem if necessary. mem requirement ~ n^2, not n^3
use types,only:rprec
use param,only:ld,ny,nz,dx,dy,dz, USE_MPI, coord
USE SIM_PARAM
use test_filtermodule

implicit none

integer :: jz
real(kind=rprec), dimension(nz),intent(out):: cs_z
real(kind=rprec), dimension(ld,ny,0:nz),intent(in) :: s,idsdx,idsdy,idsdz
real(kind=rprec),dimension(ld,ny,nz),intent(in)::s11,s12,s13,s22,s23,s33
real(kind=rprec), dimension(ld,ny) :: idsdx_bar,idsdy_bar,idsdz_bar,                &
                                      ss_bar,s11_bar,s12_bar,s13_bar,s22_bar,       &
                                      s23_bar,s33_bar,            &
                                      s_idsdx_bar, s_idsdy_bar, s_idsdz_bar
real(kind=rprec), dimension(ld,ny) :: ss,u_bar,v_bar,w_bar,s_bar,k1,k2,k3,x1,x2,x3
real(kind=rprec) :: delta,const,den,num


!filter size as function of the 3 dx
delta = filter_size*(dx*dy*dz)**(1._rprec/3._rprec)


!THINK IT NICELY FOR PARALLELIZATION CAUSE HERE IS OF CRUCIAL IMPORTANCE
do jz=1,nz

! using L_ij as temp storage here
! watch the 0.25's:  recall w = c*z^2 close to wall, so get 0.25
   if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and.  &
        (jz == 1) ) then

     !calculations on uvp-nodes
     k1(:,:) = s(:,:,1)*u(:,:,1)                !uv-node
     k2(:,:) = s(:,:,1)*v(:,:,1)                !uv-node
     k3(:,:) = s(:,:,1)*0.25_rprec*w(:,:,2)     !parabolic interp.

     s_bar(:,:) = s(:,:,1)                       !uv-node
     u_bar(:,:) = u(:,:,1)                       !uv-node
     v_bar(:,:) = v(:,:,1)                       !uv-node
     w_bar(:,:) = 0.25_rprec*w(:,:,2)            !parabolic interp.

   else  

     !put everything (L_ij and u_bar) on w-nodes
     k1(:,:) = 0.5_rprec*(s(:,:,jz)+s(:,:,jz-1))*0.5_rprec*(u(:,:,jz)+u(:,:,jz-1))
     k2(:,:) = 0.5_rprec*(s(:,:,jz)+s(:,:,jz-1))*0.5_rprec*(v(:,:,jz)+v(:,:,jz-1))
     k3(:,:) = 0.5_rprec*(s(:,:,jz)+s(:,:,jz-1))*w(:,:,jz)

     s_bar(:,:) = 0.5_rprec*(s(:,:,jz) + s(:,:,jz-1))
     u_bar(:,:) = u(:,:,jz)
     v_bar(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
     w_bar(:,:) = w(:,:,jz)

   end if

   !calculating the second part for the Germano identity
   call test_filter(s_bar,G_test)
   call test_filter(u_bar,G_test)                        !in-place filtering
   call test_filter(v_bar,G_test)
   call test_filter(w_bar,G_test)

   !calculating the L_ij term for the Germano identity
   call test_filter(k1,G_test)                          !in-place filtering
   k1 = k1 - s_bar*u_bar
   call test_filter(k2,G_test)
   k2 = k2 - s_bar*v_bar
   call test_filter(k3,G_test)
   k3 = k3 - s_bar*w_bar

   !calculate |S| (w-nodes)
   ss(:,:) = sqrt(2._rprec*(s11(:,:,jz)**2 + s22(:,:,jz)**2 + &
        s33(:,:,jz)**2 + 2._rprec*(s12(:,:,jz)**2 + &
        s13(:,:,jz)**2 + s23(:,:,jz)**2)))

   !initializing the s_ij doubly filtered (w-nodes)
   s11_bar(:,:) = s11(:,:,jz)
   s12_bar(:,:) = s12(:,:,jz)
   s13_bar(:,:) = s13(:,:,jz)
   s22_bar(:,:) = s22(:,:,jz)
   s23_bar(:,:) = s23(:,:,jz)
   s33_bar(:,:) = s33(:,:,jz)

   !filtering the gradients if S_ij for M_ij calculation
   call test_filter(S11_bar,G_test)
   call test_filter(S12_bar,G_test)
   call test_filter(S13_bar,G_test)
   call test_filter(S22_bar,G_test)
   call test_filter(S23_bar,G_test)
   call test_filter(S33_bar,G_test)

   !calculating the |S| doubly filtered (w-nodes)
   SS_bar = sqrt(2._rprec*(S11_bar**2 + S22_bar**2 + S33_bar**2 + &
          2._rprec*(S12_bar**2 + S13_bar**2 + S23_bar**2)))

   !calculating the product of |S| and S_ij not test-filtered
   S_idsdx_bar(:,:) = SS(:,:)*idsdx(:,:,jz)
   S_idsdy_bar(:,:) = SS(:,:)*idsdy(:,:,jz)
   S_idsdz_bar(:,:) = SS(:,:)*idsdz(:,:,jz)

   !test-filtering the product of |S| and S_ij (first part of M_ij)
   call test_filter(S_idsdx_bar,G_test)
   call test_filter(S_idsdy_bar,G_test)
   call test_filter(S_idsdz_bar,G_test)

   !create the M_ij with Smagorinsky coefficient taken away (scale similarity)
   const = delta**2
   x1 = const*(S_idsdx_bar - 4._rprec*SS_bar*idsdx_bar)
   x2 = const*(S_idsdy_bar - 4._rprec*SS_bar*idsdy_bar)
   x3 = const*(S_idsdz_bar - 4._rprec*SS_bar*idsdz_bar)

   NUM=sum(k1*x1+k2*x2+k3*x3) 
   DEN=sum(x1**2 + x2**2 + x3**2) 
   
   !calculating the Smagorinsky coefficient by plane averaging
   IF(DEN.NE.0.0_RPREC)THEN
      CS_Z(jz)=NUM/DEN
   ELSE
      CS_Z(jz)=0.0_RPREC   !CLIP TO ZERO
   ENDIF

   !clipping the negative values
   CS_Z(jz) = max(0._rprec, CS_Z(jz))

   !HANDLE 0 GRADIENT
   !IF(DEN.EQ.0.0_RPREC)THEN   
   !   CS_Z(jz)=0._rprec
   !   print*,'CLIPPED TO ZERO!!'
   !ENDIF

enddo


end subroutine std_dyn_sc
