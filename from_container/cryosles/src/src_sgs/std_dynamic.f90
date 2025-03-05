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
subroutine std_dynamic(CS_Z,S11,S12,S13,S22,S23,S33)

!can save more mem if necessary. mem requirement ~ n^2, not n^3
use types,only:rprec
use param,only:ld,ny,nz,dx,dy,dz, USE_MPI, coord
USE SIM_PARAM
use test_filtermodule

implicit none

integer :: jz
real(kind=rprec), dimension(nz),intent(out):: CS_Z
real(kind=rprec),dimension(ld,ny,nz),intent(in)::S11,S12,S13,S22,S23,S33
real(kind=rprec), dimension(ld,ny)::L11,L12,L13,L22,L23,L33,&
     M11,M12,M13,M22,M23,M33
real(kind=rprec), dimension(ld,ny) :: S_bar,S11_bar,S12_bar,&
     S13_bar,S22_bar,S23_bar,S33_bar,S_S11_bar, S_S12_bar,&
     S_S13_bar, S_S22_bar, S_S23_bar, S_S33_bar
real(kind=rprec), dimension(ld,ny) :: u_bar,v_bar,w_bar
real(kind=rprec), dimension(ld,ny) :: S
real(kind=rprec) :: delta, const,DEN,NUM

!filter size as function of the 3 dx
delta = filter_size*(dx*dy*dz)**(1._rprec/3._rprec)

do jz=1,nz
! using L_ij as temp storage here
! watch the 0.25's:  recall w = c*z^2 close to wall, so get 0.25
   if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and.  &
        (jz == 1) ) then

     !calculations on uvp-nodes
     L11(:,:) = u(:,:,1)*u(:,:,1)                !uv-node
     L12(:,:) = u(:,:,1)*v(:,:,1)                !uv-node
     L13(:,:) = u(:,:,1)*0.25_rprec*w(:,:,2)     !parabolic interp.
     L22(:,:) = v(:,:,1)*v(:,:,1)                !uv-node
     L23(:,:) = v(:,:,jz)*0.25_rprec*w(:,:,2)    !parabolic interp.
     L33(:,:) = (0.25_rprec*w(:,:,2))**2         !parabolic interp.

     u_bar(:,:) = u(:,:,1)                       !uv-node
     v_bar(:,:) = v(:,:,1)                       !uv-node
     w_bar(:,:) = 0.25_rprec*w(:,:,2)            !parabolic interp.

   else  

     !put everything (L_ij and u_bar) on w-nodes
     L11(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))*&
                0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))
     L12(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))*&
                0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
     L13(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))*w(:,:,jz)
     L22(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))*&
                0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
     L23(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))*w(:,:,jz)
     L33(:,:) = w(:,:,jz)*w(:,:,jz)

     u_bar(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))
     v_bar(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
     w_bar(:,:) = w(:,:,jz)

   end if

   !calculating the second part for the Germano identity
   call test_filter(u_bar,G_test)                        !in-place filtering
   call test_filter(v_bar,G_test)
   call test_filter(w_bar,G_test)

   !calculating the L_ij term for the Germano identity
   call test_filter(L11,G_test)                          !in-place filtering
   L11 = L11 - u_bar*u_bar
   call test_filter(L12,G_test)
   L12 = L12 - u_bar*v_bar
   call test_filter(L13,G_test)
   L13 = L13 - u_bar*w_bar
   call test_filter(L22,G_test)
   L22 = L22 - v_bar*v_bar
   call test_filter(L23,G_test)
   L23 = L23 - v_bar*w_bar
   call test_filter(L33,G_test)
   L33 = L33 - w_bar*w_bar

   !calculate |S| (w-nodes)
   S(:,:) = sqrt(2._rprec*(S11(:,:,jz)**2 + S22(:,:,jz)**2 + &
        S33(:,:,jz)**2 + 2._rprec*(S12(:,:,jz)**2 + &
        S13(:,:,jz)**2 + S23(:,:,jz)**2)))

   !initializing the S_ij doubly filtered (w-nodes)
   S11_bar(:,:) = S11(:,:,jz)
   S12_bar(:,:) = S12(:,:,jz)
   S13_bar(:,:) = S13(:,:,jz)
   S22_bar(:,:) = S22(:,:,jz)
   S23_bar(:,:) = S23(:,:,jz)
   S33_bar(:,:) = S33(:,:,jz)

   !filtering the gradients if S_ij for M_ij calculation
   call test_filter(S11_bar,G_test)
   call test_filter(S12_bar,G_test)
   call test_filter(S13_bar,G_test)
   call test_filter(S22_bar,G_test)
   call test_filter(S23_bar,G_test)
   call test_filter(S33_bar,G_test)

   !calculating the |S| doubly filtered (w-nodes)
   S_bar = sqrt(2._rprec*(S11_bar**2 + S22_bar**2 + S33_bar**2 + &
          2._rprec*(S12_bar**2 + S13_bar**2 + S23_bar**2)))

   !calculating the product of |S| and S_ij not test-filtered
   S_S11_bar(:,:) = S(:,:)*S11(:,:,jz)
   S_S12_bar(:,:) = S(:,:)*S12(:,:,jz)
   S_S13_bar(:,:) = S(:,:)*S13(:,:,jz)
   S_S22_bar(:,:) = S(:,:)*S22(:,:,jz)
   S_S23_bar(:,:) = S(:,:)*S23(:,:,jz)
   S_S33_bar(:,:) = S(:,:)*S33(:,:,jz)

   !test-filtering the product of |S| and S_ij (first part of M_ij)
   call test_filter(S_S11_bar,G_test)
   call test_filter(S_S12_bar,G_test)
   call test_filter(S_S13_bar,G_test)
   call test_filter(S_S22_bar,G_test)
   call test_filter(S_S23_bar,G_test)
   call test_filter(S_S33_bar,G_test)

   !create the M_ij with Smagorinsky coefficient taken away (scale similarity)
   const = 2._rprec*delta**2
   M11 = const*(S_S11_bar - 4._rprec*S_bar*S11_bar)
   M12 = const*(S_S12_bar - 4._rprec*S_bar*S12_bar)
   M13 = const*(S_S13_bar - 4._rprec*S_bar*S13_bar)
   M22 = const*(S_S22_bar - 4._rprec*S_bar*S22_bar)
   M23 = const*(S_S23_bar - 4._rprec*S_bar*S23_bar)
   M33 = const*(S_S33_bar - 4._rprec*S_bar*S33_bar)

   NUM=sum(L11*M11+L22*M22+L33*M33+2._rprec*(L12*M12+L13*M13+L23*M23)) 
   DEN=sum(M11**2 + M22**2 + M33**2 + 2._rprec*(M12**2 + M13**2 + M23**2)) 
   
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

end do


end subroutine std_dynamic
