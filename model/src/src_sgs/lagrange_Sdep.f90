!!-----------------------------------------------------------------------------------
!!  OBJECT: subroutine lagrange_Sdep (S11,S12,S13,S22,S23,S33)
!!-----------------------------------------------------------------------------------
!!
!!  LAST MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 14/10/2011
!!
!!  DESCRIPTION:
!!
!!  this is the w-node version
!!  provides Cs_opt2 1:nz
!!  MPI: required u,v on 0:nz, except bottom node 1:nz
!!
!!
!!
!!
!!
!!
!!
!!-----------------------------------------------------------------------------------
subroutine lagrange_Sdep (S11,S12,S13,S22,S23,S33,LAGRAN_DT)
use types, only: rprec
use param
use sim_param
use sgsmodule, only: F_LM,F_MM,F_QN,F_NN,NU_T,CS_OPT2,beta,opftime,Beta_avg,Betaclip_avg
use test_filtermodule
use stretch

implicit none

integer :: jx, jy, jz
integer :: i, px, py, lx, ly, lz
integer :: counter1,counter2,counter3,counter4,counter5
integer :: istart, iend

real(rprec) :: tf1,tf2,tf1_2,tf2_2     !size of the second test filter
real(rprec) :: fractus
real(rprec), dimension(ld,ny,nz),intent(in) :: S11,S12,S13,S22,S23,S33
real(rprec), dimension(ld,ny) :: Cs_opt2_2d,Cs_opt2_4d      !2d to save mem.
real(rprec) :: Betaclip                !made a scalar to save mem.

real(rprec), dimension(ld,ny) :: S,tempos
real(rprec), dimension(ld,ny) :: L11,L12,L13,L22,L23,L33
real(rprec), dimension(ld,ny) :: Q11,Q12,Q13,Q22,Q23,Q33
real(rprec), dimension(ld,ny) :: M11,M12,M13,M22,M23,M33
real(rprec), dimension(ld,ny) :: N11,N12,N13,N22,N23,N33

real(rprec), dimension(nz) :: LMvert,MMvert,QNvert,NNvert
real(rprec), dimension(ld,ny) :: LM,MM,QN,NN,Tn,epsi,dumfac

real(rprec), dimension(ld,ny) :: S_bar,S11_bar,S12_bar,&
     S13_bar,S22_bar,S23_bar,S33_bar,S_S11_bar, S_S12_bar,&
     S_S13_bar, S_S22_bar, S_S23_bar, S_S33_bar
real(rprec), dimension(ld,ny) :: S_hat,S11_hat,S12_hat,&
     S13_hat,S22_hat,S23_hat,S33_hat,S_S11_hat, S_S12_hat,&
     S_S13_hat, S_S22_hat, S_S23_hat, S_S33_hat

real(rprec), dimension(ld,ny) :: u_bar,v_bar,w_bar
real(rprec), dimension(ld,ny) :: u_hat,v_hat,w_hat

real(rprec) :: delta,const
real(rprec) :: opftdelta,powcoeff
REAL(RPREC),INTENT(IN) :: LAGRAN_DT
!logical, save :: F_LM_MM_init = .false.
!logical, save :: F_QN_NN_init = .false.
integer :: kk
real(rprec) :: fac_1,fac_2,fac_3

!---------------------------------------------------------------------
if (VERBOSE) write (*, *) 'started lagrange_Sdep'

!fix the delta for Smagorinsky
!test filter sizes
tf1=2._rprec
tf2=4._rprec   !usually equal to tf1**2, can be different

!square of test filter sizes
tf1_2=tf1**2.0_rprec
tf2_2=tf2**2.0_rprec

!MS: don't need probably
!TS Add the opftdelta
!TS FOR BUILDING (opftime=1*1.5)
powcoeff = -1._rprec/8._rprec

call interpolag_Sdep ()

!time step used to advance the subgrid transport equation
fractus= 1.0_rprec/REAL(ny*nx,kind=rprec)

do jz = 1,nz
  kk = jz + eta_displ(coord)
  delta = filter_size*(dx*dy*(get_z_local(2*kk-1)-get_z_local(2*kk-3)))**(1._rprec/3._rprec)
  opftdelta = opftime*delta
  const = 2.0_rprec*(delta**2)

  if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and. (jz == 1) ) then
    !!! watch the 0.25's:  recall w = c*z^2 close to wall, so get 0.25
    ! put on uvp-nodes
    u_bar(:,:) = u(:,:,1)  ! first uv-node
    v_bar(:,:) = v(:,:,1)  ! first uv-node
    w_bar(:,:) = 0.25_rprec*w(:,:,2)  ! first uv-node
  else
    ! w-nodes
 
   fac_1=get_z_local(2*kk-1)
   fac_2=get_z_local(2*kk-3)
   fac_3=get_z_local(2*kk-2)

   u_bar(:,:) = u(:,:,jz-1) + ((u(:,:,jz)-u(:,:,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)
   v_bar(:,:) = v(:,:,jz-1) + ((v(:,:,jz)-v(:,:,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)
  
    !u_bar(:,:) = 0.5_rprec*(u(:,:,jz) + u(:,:,jz-1))
    !v_bar(:,:) = 0.5_rprec*(v(:,:,jz) + v(:,:,jz-1))
    w_bar(:,:) = w(:,:,jz)
  end if

  u_hat = u_bar
  v_hat = v_bar
  w_hat = w_bar
  L11 = u_bar*u_bar
  L12 = u_bar*v_bar
  L13 = u_bar*w_bar
  L22 = v_bar*v_bar
  L23 = v_bar*w_bar
  L33 = w_bar*w_bar
  Q11 = u_bar*u_bar
  Q12 = u_bar*v_bar
  Q13 = u_bar*w_bar
  Q22 = v_bar*v_bar
  Q23 = v_bar*w_bar
  Q33 = w_bar*w_bar

  call test_filter(u_bar,G_test)
  call test_filter(v_bar,G_test)
  call test_filter(w_bar,G_test)
  call test_filter(u_hat,G_test_test)
  call test_filter(v_hat,G_test_test)
  call test_filter(w_hat,G_test_test)

  call test_filter(L11,G_test)  ! in-place filtering
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
  call test_filter(Q11,G_test_test)
  Q11 = Q11 - u_hat*u_hat
  call test_filter(Q12,G_test_test)
  Q12 = Q12 - u_hat*v_hat
  call test_filter(Q13,G_test_test)
  Q13 = Q13 - u_hat*w_hat
  call test_filter(Q22,G_test_test)
  Q22 = Q22 - v_hat*v_hat
  call test_filter(Q23,G_test_test)
  Q23 = Q23 - v_hat*w_hat
  call test_filter(Q33,G_test_test)
  Q33 = Q33 - w_hat*w_hat


  !do jx=1,nx
  !   do jy=1,ny
  !      write(*,*) coord,jx,jy,jz,S11(jx,jy,jz),S22(jx,jy,jz),S33(jx,jy,jz)
  !   enddo
  !enddo

  ! calculate |S|
  S(:,:) = sqrt(2._rprec*(S11(:,:,jz)**2.0_rprec + S22(:,:,jz)**2.0_rprec + S33(:,:,jz)**2.0_rprec +&
                2._rprec*(S12(:,:,jz)**2.0_rprec + S13(:,:,jz)**2.0_rprec + S23(:,:,jz)**2.0_rprec)))

  S11_bar(:,:) = S11(:,:,jz)  ! S_ij already on w-nodes
  S12_bar(:,:) = S12(:,:,jz)
  S13_bar(:,:) = S13(:,:,jz)
  S22_bar(:,:) = S22(:,:,jz)
  S23_bar(:,:) = S23(:,:,jz)
  S33_bar(:,:) = S33(:,:,jz)

  S11_hat = S11_bar
  S12_hat = S12_bar
  S13_hat = S13_bar
  S22_hat = S22_bar
  S23_hat = S23_bar
  S33_hat = S33_bar

  call test_filter(S11_bar,G_test)
  call test_filter(S12_bar,G_test)
  call test_filter(S13_bar,G_test)
  call test_filter(S22_bar,G_test)
  call test_filter(S23_bar,G_test)
  call test_filter(S33_bar,G_test)

  call test_filter(S11_hat,G_test_test)
  call test_filter(S12_hat,G_test_test)
  call test_filter(S13_hat,G_test_test)
  call test_filter(S22_hat,G_test_test)
  call test_filter(S23_hat,G_test_test)
  call test_filter(S33_hat,G_test_test)

  S_bar = sqrt(2._rprec*(S11_bar**2 + S22_bar**2 + S33_bar**2 +&
               2._rprec*(S12_bar**2 + S13_bar**2 + S23_bar**2)))
  S_hat = sqrt(2._rprec*(S11_hat**2 + S22_hat**2 + S33_hat**2 +&
               2._rprec*(S12_hat**2 + S13_hat**2 + S23_hat**2)))
  S_S11_bar(:,:) = S(:,:)*S11(:,:,jz)
  S_S12_bar(:,:) = S(:,:)*S12(:,:,jz)
  S_S13_bar(:,:) = S(:,:)*S13(:,:,jz)
  S_S22_bar(:,:) = S(:,:)*S22(:,:,jz)
  S_S23_bar(:,:) = S(:,:)*S23(:,:,jz)
  S_S33_bar(:,:) = S(:,:)*S33(:,:,jz)

  S_S11_hat(:,:) = S_S11_bar(:,:)
  S_S12_hat(:,:) = S_S12_bar(:,:)
  S_S13_hat(:,:) = S_S13_bar(:,:)
  S_S22_hat(:,:) = S_S22_bar(:,:)
  S_S23_hat(:,:) = S_S23_bar(:,:)
  S_S33_hat(:,:) = S_S33_bar(:,:)
  
  call test_filter(S_S11_bar,G_test)
  call test_filter(S_S12_bar,G_test)
  call test_filter(S_S13_bar,G_test)
  call test_filter(S_S22_bar,G_test)
  call test_filter(S_S23_bar,G_test)
  call test_filter(S_S33_bar,G_test)     

  call test_filter(S_S11_hat,G_test_test)
  call test_filter(S_S12_hat,G_test_test)
  call test_filter(S_S13_hat,G_test_test)
  call test_filter(S_S22_hat,G_test_test)
  call test_filter(S_S23_hat,G_test_test)
  call test_filter(S_S33_hat,G_test_test)  
  
  M11 = const*(S_S11_bar - tf1_2*S_bar*S11_bar)  ! assume beta=1
  M12 = const*(S_S12_bar - tf1_2*S_bar*S12_bar)
  M13 = const*(S_S13_bar - tf1_2*S_bar*S13_bar)
  M22 = const*(S_S22_bar - tf1_2*S_bar*S22_bar)
  M23 = const*(S_S23_bar - tf1_2*S_bar*S23_bar)
  M33 = const*(S_S33_bar - tf1_2*S_bar*S33_bar)

  N11 = const*(S_S11_hat - tf2_2*S_hat*S11_hat)  ! assume beta=1, so beta^2=1
  N12 = const*(S_S12_hat - tf2_2*S_hat*S12_hat)
  N13 = const*(S_S13_hat - tf2_2*S_hat*S13_hat)
  N22 = const*(S_S22_hat - tf2_2*S_hat*S22_hat)
  N23 = const*(S_S23_hat - tf2_2*S_hat*S23_hat)
  N33 = const*(S_S33_hat - tf2_2*S_hat*S33_hat)
  
  LM = L11*M11+L22*M22+L33*M33+2._rprec*(L12*M12+L13*M13+L23*M23)
  MM = M11**2+M22**2+M33**2+2._rprec*(M12**2+M13**2+M23**2)
  QN = Q11*N11+Q22*N22+Q33*N33+2._rprec*(Q12*N12+Q13*N13+Q23*N23)
  NN = N11**2+N22**2+N33**2+2._rprec*(N12**2+N13**2+N23**2)

!  if (inilag) then
!   if ((.not. F_LM_MM_init) .and. (jt == cs_count .or. jt == DYN_init)) then
    if (jt_total .eq. dyn_init) then
      print *,'F_MM and F_LM initialized' 
      F_MM (:,:,jz) = MM
      F_LM (:,:,jz) = 0.03_rprec*MM
      F_MM(ld-1:ld,:,jz)=1._rprec
      F_LM(ld-1:ld,:,jz)=1._rprec

      !if (jz == nz) F_LM_MM_init = .true.
    end if
!  end if

!   if (inflow) then
!     iend = floor (buff_end * nx + 1._rprec)
!     iend = modulo (iend - 1, nx) + 1
!     istart = floor ((buff_end - buff_len) * nx + 1._rprec)
!     istart = modulo (istart - 1, nx) + 1
!       
!     Tn=merge(.1_rprec*const*S**2,MM,MM.le..1_rprec*const*S**2)
!     MM=Tn
!     LM(istart + 1:iend, 1:ny) = 0._rprec
!     F_LM(istart + 1:iend, 1:ny, jz) = 0._rprec
!     Tn=merge(.1_rprec*const*S**2,NN,NN.le..1_rprec*const*S**2)
!     NN=Tn
!     QN(istart + 1:iend, 1:ny) = 0._rprec
!     F_QN(istart + 1:iend, 1:ny, jz) = 0._rprec
!   endif

  !TS FIX THE BUG IN COMPAQ COMPILER
  !TS ADD the following
  Tn = max(real(F_LM(:,:,jz)*F_MM(:,:,jz)),real(1E-24))
  Tn = opftdelta*(Tn**powcoeff)
  Tn(:,:) = max(real(1E-24),real(Tn(:,:)))
  dumfac = lagran_dt/Tn 
  epsi = dumfac / (1._rprec+dumfac)
  
  F_LM(:,:,jz) = (epsi*LM + (1._rprec-epsi)*F_LM(:,:,jz))
  F_MM(:,:,jz) = (epsi*MM + (1._rprec-epsi)*F_MM(:,:,jz))
  F_LM(:,:,jz) = max(real(1E-24),real(F_LM(:,:,jz)))
  F_MM(:,:,jz) = max(real(1E-24),real(F_MM(:,:,jz)))
  
  Cs_opt2_2d(:,:) = F_LM(:,:,jz)/F_MM(:,:,jz)
  Cs_opt2_2d(ld,:) = 1E-24
  Cs_opt2_2d(ld-1,:) = 1E-24
  Cs_opt2_2d(:,:) = max(real(1E-24),real(Cs_opt2_2d(:,:)))

!  if (inilag) then
!   if ((.not. F_QN_NN_init) .and. (jt == cs_count .or. jt == DYN_init)) then
    if (jt_total .eq. dyn_init) then
      print *,'F_NN and F_QN initialized'
      F_NN (:,:,jz) = NN
      F_QN (:,:,jz) = 0.03_rprec*NN
      F_NN(ld-1:ld,:,jz)=1._rprec
      F_QN(ld-1:ld,:,jz)=1._rprec
  
      !if (jz == nz) F_QN_NN_init = .true.
    end if
!  end if

  !TS FIX THE BUG IN COMPAQ COMPILER
  !TS ADD the following
  Tn = max(real(F_QN(:,:,jz)*F_NN(:,:,jz)),real(1E-24))
  Tn = opftdelta*(Tn**powcoeff)
  Tn(:,:) = max(real(1E-24),real(Tn(:,:)))
  dumfac = lagran_dt/Tn
  epsi = dumfac / (1._rprec+dumfac)
  
  F_QN(:,:,jz) = (epsi*QN + (1._rprec-epsi)*F_QN(:,:,jz))
  F_NN(:,:,jz) = (epsi*NN + (1._rprec-epsi)*F_NN(:,:,jz))
  F_QN(:,:,jz) = max(real(1E-24),real(F_QN(:,:,jz)))
  F_NN(:,:,jz) = max(real(1E-24),real(F_NN(:,:,jz)))

  !write(*,*) F_NN(:,:,jz)
 
  Cs_opt2_4d(:,:) = F_QN(:,:,jz)/F_NN(:,:,jz)
  Cs_opt2_4d(ld,:) = 1E-24
  Cs_opt2_4d(ld-1,:) = 1E-24
  Cs_opt2_4d(:,:) = max(real(1E-24),real(Cs_opt2_4d(:,:)))
  
  !--changed to save mem.
  Beta(:,:,jz)=(Cs_opt2_4d(:,:)/Cs_opt2_2d(:,:))&
               **(log(tf1)/(log(tf2)-log(tf1)))
  counter1=0
  counter2=0

  do jy=1,Ny
    do jx=1,Nx
      if (Beta(jx,jy,jz).LE.1/(tf1*tf2)) counter1=counter1+1
    end do
  end do

!  !--MPI: this is not valid
!  if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) .and.  &
!       (jz == nz) ) then
!    Beta(:,:,jz)=1._rprec  !JF: why to do this here???
!  end if

  do jy = 1, ny
    do jx = 1, nx  !--nx was ld, changed by JF
      Betaclip=max(real(Beta(jx,jy,jz)),real(1._rprec/(tf1*tf2)))
      Cs_opt2(jx,jy,jz)=Cs_opt2_2d(jx,jy)/Betaclip  !--2d array
    end do
  end do
  
  Cs_opt2(ld,:,jz) = 1E-24
  Cs_opt2(ld-1,:,jz) = 1E-24
  Cs_opt2(:,:,jz) = max(real(1E-24),real(Cs_opt2(:,:,jz)))
  
  !JF: which is the clipped one???
  Beta_avg(jz)=sum(Cs_opt2_2d(1:nx,1:ny))/sum(Cs_opt2(1:nx,1:ny,jz))
  Betaclip_avg(jz)=sum(Cs_opt2_4d(1:nx,1:ny))/sum(Cs_opt2_2d(1:nx,1:ny))

   !FIXING TO 0 THE CS COEFFICIENT
   IF(CLIP_SGS)THEN
      DO JX=1,LD
      DO JY=1,NY
         IF((F_NN(JX,JY,JZ).lt.1E-24_rprec).OR.(F_MM(JX,JY,JZ).lt.1E-24_rprec) )THEN 
            Cs_opt2(JX,JY,JZ)=Co**2
         ENDIF
      ENDDO
      ENDDO
   ENDIF  
  
!   if (calc_lag_Ssim) then !calc_lag_Ssim is specified in param.f90
!     ! This basically saves the Cs calculated using the scale-invariance
!     ! assumption within the scale-dependent model to Cs_Ssim
! 
!     Cs_Ssim(:,:,jz)=Cs_opt2_2d
!     if (mod(jt,100) == 0) then
!       if (jz==1.or.jz==nz/4.or.jz==nz/2) then
!         !--2d array to save mem.
!         write(91)real(jt*dt),real(jz),real(Cs_opt2_2d(1:NX,1:NY))
!         write(92)real(jt*dt),real(jz),real(Cs_opt2_4d(1:NX,1:NY))
!       end if
!     end if
! 
!     if (mod(jt,200) == 0) then
!       LMvert(jz) = sum(sum(LM,DIM=1),DIM=1)/ny/nx
!       MMvert(jz) = sum(sum(MM,DIM=1),DIM=1)/ny/nx
!       QNvert(jz) = sum(sum(QN,DIM=1),DIM=1)/ny/nx
!       NNvert(jz) = sum(sum(NN,DIM=1),DIM=1)/ny/nx
!     end if
! 
!   end if

end do  ! this ends the main jz=1,nz loop          

!TS
!if(use_bldg)then
!  do i=1,n_bldg
!    px=bldg_pts(1,i)
!    py=bldg_pts(2,i)
!    lx=bldg_pts(3,i)
!    ly=bldg_pts(4,i)
!    lz=bldg_pts(5,i)
!    Cs_opt2(px:px+lx,py:py+ly,1:lz)=1.E-24
!  enddo
!endif
!

if (VERBOSE) write (*, *) 'finished lagrange_Sdep'

end subroutine lagrange_Sdep
