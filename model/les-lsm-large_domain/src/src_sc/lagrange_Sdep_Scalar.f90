MODULE LASD

USE TYPES
USE PARAM
USE SIM_PARAM, ONLY: DS_OPT2_2D,DS_OPT2_4D, F_KX,F_XX,F_PY,F_YY
use stretch

IMPLICIT NONE



!----------------------------------------------------------------------------------
CONTAINS

! this is the w-node version
!--provides Cs_opt2 1:nz
!--MPI: required u,v on 0:nz, except bottom node 1:nz
subroutine lagrange_Sdep_Scalar(S11,S12,S13,S22,S23,S33,theta,dThdx,dThdy,dThdz,Ds_opt2,beta_s)
! standard dynamic model to calculate the Smagorinsky coefficient
! this is done layer-by-layer to save memory
! everything is done to be on uv-nodes
! -note: we need to calculate |S| here, too.
! stuff is done on uv-nodes
! can save more mem if necessary. mem requirement ~ n^2, not n^3
use service
use sim_param
use sgsmodule,only:F_LM,F_MM,F_QN,F_NN,beta,opftime,Beta_avg,Betaclip_avg, &
    Beta_avg_s,Betaclip_avg_s
use test_filtermodule

implicit none


REAL(KIND=RPREC),DIMENSION(LD,NY,0:NZ),INTENT(IN) :: THETA
REAL(KIND=RPREC),DIMENSION(LD,NY,NZ),INTENT(OUT) ::   DS_OPT2
REAL(KIND=RPREC),DIMENSION(LD,NY,NZ),INTENT(OUT) ::   BETA_S
REAL(KIND=RPREC), DIMENSION(LD,NY,NZ),INTENT(IN) ::  DTHDX,DTHDY,DTHDZ

integer :: jx,jy,jz
!TS
integer :: i,j,k, px, pyy, lx, ly, lz
integer :: counter1,counter2,counter3,counter4,counter5
real(kind=rprec):: tf1,tf2,tf1_2,tf2_2 ! Size if the second test filter
real(kind=rprec) , dimension(ld,ny,nz) :: visc,dissip
real(kind=rprec) :: fractus
real(kind=rprec), dimension(ld,ny,nz),intent(in) :: S11,S12,S13,S22,S23,S33
real(kind=rprec), dimension(ld,ny,nz) :: Betaclip,Betaclip_s
real(kind=rprec), dimension(ld,ny,nz) :: T_avg, T2_avg, SIGMA_T
integer, dimension(nz) :: added

real(kind=rprec), dimension(ld,ny) :: S,tempos
real(kind=rprec), dimension(ld,ny) :: L11,L12,L13,L22,L23,L33
real(kind=rprec), dimension(ld,ny) :: Q11,Q12,Q13,Q22,Q23,Q33
real(kind=rprec), dimension(ld,ny) :: M11,M12,M13,M22,M23,M33
real(kind=rprec), dimension(ld,ny) :: N11,N12,N13,N22,N23,N33
real(kind=rprec), dimension(ld,ny) :: K1,K2,K3
real(kind=rprec), dimension(ld,ny) :: X1,X2,X3
real(kind=rprec), dimension(ld,ny) :: P1,P2,P3
real(kind=rprec), dimension(ld,ny) :: Y1,Y2,Y3

real(kind=rprec), dimension(nz) :: LMvert,MMvert,QNvert,NNvert,KXvert,XXvert,PYvert,YYvert
real(kind=rprec), dimension(ld,ny) :: LM,MM,QN,NN,Tn,Tn_s,epsi,epsi_s,dumfac,dumfac_s,KX,XX,PY,YY
real(kind=rprec), dimension(ld,ny,nz) :: Tn1_s,Tn2_s,Tn1,Tn2

real(kind=rprec), dimension(ld,ny) :: S_bar,S11_bar,S12_bar,&
     S13_bar,S22_bar,S23_bar,S33_bar,S_S11_bar, S_S12_bar,&
     S_S13_bar, S_S22_bar, S_S23_bar, S_S33_bar
real(kind=rprec), dimension(ld,ny) :: S_hat,S11_hat,S12_hat,&
     S13_hat,S22_hat,S23_hat,S33_hat,S_S11_hat, S_S12_hat,&
     S_S13_hat, S_S22_hat, S_S23_hat, S_S33_hat
real(kind=rprec), dimension(ld,ny) :: dThdx_bar,dThdy_bar,dThdz_bar
real(kind=rprec), dimension(ld,ny) :: dThdx_hat,dThdy_hat,dThdz_hat
real(kind=rprec), dimension(ld,ny) :: S_dThdx_bar,S_dThdy_bar,S_dThdz_bar
real(kind=rprec), dimension(ld,ny) :: S_dThdx_hat,S_dThdy_hat,S_dThdz_hat

real(kind=rprec), dimension(ld,ny) :: u_bar,v_bar,w_bar,theta_bar
real(kind=rprec), dimension(ld,ny) :: u_hat,v_hat,w_hat,theta_hat

real(kind=rprec) :: delta,const,const2,lagran_dt
real(kind=rprec) :: opftdelta,powcoeff,powcoeff_s

!real (rprec), dimension (ld, ny, 0:nz) :: u_temp, v_temp, theta_temp

logical, save :: F_LM_MM_init = .false.
logical, save :: F_QN_NN_init = .false.

logical, save :: F_KX_XX_init = .false.
logical, save :: F_PY_YY_init = .false.

integer :: kk
real(rprec) :: fac_1,fac_2,fac_3

!---------------------------------------------------------------------

if (VERBOSE) write (*, *) 'started lagrange_Sdep'



lagran_dt=dt*real(cs_count,kind=rprec)


tf1=2._rprec
tf2=4._rprec
tf1_2=tf1**2
tf2_2=tf2**2

!TS Add the opftdelta
!TS FOR BUILDING (opftime=1*1.5)
powcoeff = -1._rprec/8._rprec
powcoeff_s = -1._rprec/8._rprec

call interpolag_Sdep_Scalar(F_KX,F_XX,F_PY,F_YY)

!! (this is the time step used in the lagrangian computations)
fractus= 1._rprec/real(ny*nx,kind=rprec)

do jz = 1,nz
   kk = jz+eta_displ(coord)
   delta = filter_size*(dx*dy*(get_z_local(2*kk-1)-get_z_local(2*kk-3)))**(1._rprec/3._rprec)
   opftdelta = opftime*delta
   const = 2._rprec*(delta**2)
   const2 = (delta**2)

   !For computing SIGMA_T = sqrt( avg_(T**2) - (avg_T*avgT) )
   !----------------------------------------------------------
   !T_avg(:,:,jz) = T_avg(:,:,jz) + theta(:,:,jz)
   !T2_avg(:,:,jz) = T2_avg(:,:,jz) + (theta(:,:,jz)**2)
   !added(jz) = added(jz) + 1
   !SIGMA_T(:,:,jz) = sqrt( (T2_avg(:,:,jz)/added(jz)) - ((T_avg(:,:,jz)/added(jz))**2.d0)  )

   SIGMA_T(:,:,jz)=1.0_rprec !d0
  
   !----------------------------------------------------------

! using L_ij as temp storage here
   if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) .and. (jz == 1) ) then
      !!! watch the 0.25's:  recall w = c*z^2 close to wall, so get 0.25
      ! put on uvp-nodes
      u_bar(:,:) = u(:,:,1)! first uv-node
      v_bar(:,:) = v(:,:,1)! first uv-node
      w_bar(:,:) = 0.25_rprec*w(:,:,2)! first uv-node  

      theta_bar(:,:) = theta(:,:,1)! first uv-node

   else ! w-nodes

  fac_1=get_z_local(2*kk-1)
  fac_2=get_z_local(2*kk-3)
  fac_3=get_z_local(2*kk-2)

   u_bar(:,:) = u(:,:,jz-1) + ((u(:,:,jz)-u(:,:,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)
   v_bar(:,:) = v(:,:,jz-1) + ((v(:,:,jz)-v(:,:,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)
   theta_bar(:,:) = theta(:,:,jz-1) + ((theta(:,:,jz)-theta(:,:,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)

      w_bar(:,:) = w(:,:,jz)


   end if

   u_hat = u_bar
   v_hat = v_bar
   w_hat = w_bar
     

   theta_hat = theta_bar

   K1 = theta_bar*u_bar
   K2 = theta_bar*v_bar
   K3 = theta_bar*w_bar

   P1 = theta_bar*u_bar
   P2 = theta_bar*v_bar
   P3 = theta_bar*w_bar

   call test_filter(u_bar,G_test)
   call test_filter(v_bar,G_test)
   call test_filter(w_bar,G_test)
   call test_filter(theta_bar,G_test)

   call test_filter(K1,G_test)
   K1 = K1 - theta_bar*u_bar
   call test_filter(K2,G_test)
   K2 = K2 - theta_bar*v_bar
   call test_filter(K3,G_test)
   K3 = K3 - theta_bar*w_bar

   call test_filter(u_hat,G_test_test)
   call test_filter(v_hat,G_test_test)
   call test_filter(w_hat,G_test_test)
   call test_filter(theta_hat,G_test_test)

   call test_filter(P1,G_test_test)
   P1 = P1 - theta_hat*u_hat
   call test_filter(P2,G_test_test)
   P2 = P2 - theta_hat*v_hat
   call test_filter(P3,G_test_test)
   P3 = P3 - theta_hat*w_hat

   ! calculate |S|
   S(:,:) = sqrt(2._rprec*(S11(:,:,jz)**2 + S22(:,:,jz)**2 + S33(:,:,jz)**2 +&
        2._rprec*(S12(:,:,jz)**2 + S13(:,:,jz)**2 + S23(:,:,jz)**2)))

   ! S_ij already on w-nodes; this is the way it is imported from "sgs_stag()"
   ! Same for dThdx, dThdy, dThdz
   S11_bar(:,:) = S11(:,:,jz)  
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

   dThdx_bar = dThdx(:,:,jz)
   dThdy_bar = dThdy(:,:,jz)
   dThdz_bar = dThdz(:,:,jz)

   dThdx_hat = dThdx_bar
   dThdy_hat = dThdy_bar
   dThdz_hat = dThdz_bar

   call test_filter(S11_bar,G_test)
   call test_filter(S12_bar,G_test)
   call test_filter(S13_bar,G_test)
   call test_filter(S22_bar,G_test)
   call test_filter(S23_bar,G_test)
   call test_filter(S33_bar,G_test)

   call test_filter(dThdx_bar,G_test)
   call test_filter(dThdy_bar,G_test)
   call test_filter(dThdz_bar,G_test)

   call test_filter(S11_hat,G_test_test)
   call test_filter(S12_hat,G_test_test)
   call test_filter(S13_hat,G_test_test)
   call test_filter(S22_hat,G_test_test)
   call test_filter(S23_hat,G_test_test)
   call test_filter(S33_hat,G_test_test)

   call test_filter(dThdx_hat,G_test_test)
   call test_filter(dThdy_hat,G_test_test)
   call test_filter(dThdz_hat,G_test_test)

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

   S_dThdx_bar(:,:) = S(:,:)*dThdx(:,:,jz)
   S_dThdy_bar(:,:) = S(:,:)*dThdy(:,:,jz)
   S_dThdz_bar(:,:) = S(:,:)*dThdz(:,:,jz)

   call test_filter(S_S11_bar,G_test)
   call test_filter(S_S12_bar,G_test)
   call test_filter(S_S13_bar,G_test)
   call test_filter(S_S22_bar,G_test)
   call test_filter(S_S23_bar,G_test)
   call test_filter(S_S33_bar,G_test)

   call test_filter(S_dThdx_bar,G_test)
   call test_filter(S_dThdy_bar,G_test)
   call test_filter(S_dThdz_bar,G_test)
 

   X1 = const2*(S_dThdx_bar - tf1_2*S_bar*dThdx_bar)
   X2 = const2*(S_dThdy_bar - tf1_2*S_bar*dThdy_bar)
   X3 = const2*(S_dThdz_bar - tf1_2*S_bar*dThdz_bar)  
   
   S_S11_hat(:,:) = S_S11_bar(:,:)
   S_S12_hat(:,:) = S_S12_bar(:,:)
   S_S13_hat(:,:) = S_S13_bar(:,:)
   S_S22_hat(:,:) = S_S22_bar(:,:)
   S_S23_hat(:,:) = S_S23_bar(:,:)
   S_S33_hat(:,:) = S_S33_bar(:,:)

   S_dThdx_hat(:,:) = S_dThdx_bar(:,:)
   S_dThdy_hat(:,:) = S_dThdy_bar(:,:)
   S_dThdz_hat(:,:) = S_dThdz_bar(:,:)

   call test_filter(S_S11_hat,G_test_test)
   call test_filter(S_S12_hat,G_test_test)
   call test_filter(S_S13_hat,G_test_test)
   call test_filter(S_S22_hat,G_test_test)
   call test_filter(S_S23_hat,G_test_test)
   call test_filter(S_S33_hat,G_test_test)  

   call test_filter(S_dThdx_hat,G_test_test)
   call test_filter(S_dThdy_hat,G_test_test)
   call test_filter(S_dThdz_hat,G_test_test)

   Y1 = const2*(S_dThdx_hat - tf2_2*S_hat*dThdx_hat)
   Y2 = const2*(S_dThdy_hat - tf2_2*S_hat*dThdy_hat)
   Y3 = const2*(S_dThdz_hat - tf2_2*S_hat*dThdz_hat)

   KX = K1*X1 + K2*X2 + K3*X3 
   XX = X1**2 + X2**2 + X3**2
   PY = P1*Y1 + P2*Y2 + P3*Y3 
   YY = Y1**2 + Y2**2 + Y3**2
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Lagrangian average Part !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   !A) For the first filtering Scale:
   
   if (inilag) then
        if ((.not. F_KX_XX_init) .and. (jt == cs_count .or. jt == DYN_init)) then 
   	  
        !For the scalar case, we just copy the exact same initialization
        !as for the momentum case, above. MCB 24/7/10
        print *,'F_XX and F_KX initialized' 
!          print*, 'MM=',MM
!          print*, 'XX=',XX
        F_XX (:,:,jz) = XX
        F_KX (:,:,jz) = 0.03_rprec*XX
        F_XX(ld-1:ld,:,jz)=1._rprec
        F_KX(ld-1:ld,:,jz)=1._rprec        

        if (jz == nz) F_KX_XX_init = .true.
         
        endif
   endif
  
      
   !This section computes the "EPSILON", for the Lagrangian time integration. MCB 24/7/10
   !-----------------------------------------------------------------------------------
   
   Tn =max(real(F_LM(:,:,jz)*F_MM(:,:,jz),rprec),real(1E-24,rprec))
   Tn=opftdelta*(Tn**powcoeff)
   Tn(:,:) = max(real(1E-24,rprec),real(Tn(:,:),rprec))
   dumfac = lagran_dt/Tn 
   epsi = dumfac / (1._rprec+dumfac)

   Tn1(:,:,jz) = epsi(:,:)


   !------------------------------------------------------------------------------------


   !For the Scalar case:
   F_KX(:,:,jz)=(epsi(:,:)*KX + (1._rprec-epsi(:,:))*F_KX(:,:,jz))
   F_XX(:,:,jz)=(epsi(:,:)*XX + (1._rprec-epsi(:,:))*F_XX(:,:,jz))

   F_KX(:,:,jz)= max(real(1E-24,rprec),real(F_KX(:,:,jz),rprec)) !We insure no-zero values with "clipping"
   F_XX(:,:,jz)= max(real(1E-24,rprec),real(F_XX(:,:,jz),rprec)) !We insure no-zero values with "clipping"  
   
   !---------------------------------------------------
   !Computation of the "Cs" constant (at the first Filtering scale):
   !----------------------------------------------------------------
   
   !SCALAR PART
   Ds_opt2_2d(:,:,jz) = F_KX(:,:,jz)/F_XX(:,:,jz)
   Ds_opt2_2d(ld,:,jz) = 1E-24
   Ds_opt2_2d(ld-1,:,jz) = 1E-24
   Ds_opt2_2d(:,:,jz)=max(real(1E-24,rprec),real(Ds_opt2_2d(:,:,jz),rprec))
   IF(VERBOSE) print*,'Ds_opt2_2d(1:NX/2,1:NY/2,jz)=',Ds_opt2_2d(1:NX/2,1:NY/2,jz)

   
   !B) For the second filtering Scale:
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   
   if (inilag) then
   	  if ((.not. F_PY_YY_init) .and. (jt == cs_count .or. jt == DYN_init)) then 
   	  
         !For the scalar case, we just copy the exact same initialization
         !as for the momentum case, above. MCB 24/7/10
         F_YY (:,:,jz) = YY
         F_PY (:,:,jz) = 0.03_rprec*YY
         F_YY(ld-1:ld,:,jz)=1._rprec
         F_PY(ld-1:ld,:,jz)=1._rprec

         if (jz == nz) F_PY_YY_init = .true.
         
      end if
   end if



   !This section computes the "EPSILON", for the Lagrangian time integration. MCB 24/7/10
   !-----------------------------------------------------------------------------------

   !TS FIX THE BUG IN COMPAQ COMPILER
   !TS ADD the following
   Tn =max(real(F_QN(:,:,jz)*F_NN(:,:,jz),rprec),real(1E-24,rprec))
   Tn=opftdelta*(Tn**powcoeff)
   Tn(:,:) = max(real(1E-24,rprec),real(Tn(:,:),rprec))
   dumfac = lagran_dt/Tn
   epsi = dumfac / (1._rprec+dumfac)

   !Tn2(:,:,jz) = epsi(:,:)

   !For the scalar:
   Tn_s =max(real(F_PY(:,:,jz)*F_YY(:,:,jz),rprec),real(1E-24,rprec))
   Tn_s=opftdelta*SIGMA_T(:,:,jz)*(Tn_s**powcoeff_s)
   Tn_s(:,:) = max(real(1E-24,rprec),real(Tn_s(:,:),rprec))
   dumfac_s = lagran_dt/Tn_s
   epsi_s = dumfac_s / (1._rprec+dumfac_s) 

   Tn2_s(:,:,jz) = epsi_s(:,:)

   !------------------------------------------------------------------------------------

   !Here we apply the actual "Lagrangian time average":
   !For the Scalar case:
   F_PY(:,:,jz)=(epsi(:,:)*PY + (1._rprec-epsi(:,:))*F_PY(:,:,jz))
   F_YY(:,:,jz)=(epsi(:,:)*YY + (1._rprec-epsi(:,:))*F_YY(:,:,jz))

   F_PY(:,:,jz)= max(real(1E-24,rprec),real(F_PY(:,:,jz),rprec)) !We insure no-zero values with "clipping"
   F_YY(:,:,jz)= max(real(1E-24,rprec),real(F_YY(:,:,jz),rprec)) !We insure no-zero values with "clipping"

   !---------------------------------------------------
   !Computation of the "Cs" constant (at the first Filtering scale):
   !----------------------------------------------------------------
   
   !For the scalar: (We basically copy the same as we had for momentum,...)
   Ds_opt2_4d(:,:,jz) = F_PY(:,:,jz)/F_YY(:,:,jz)
   Ds_opt2_4d(ld,:,jz) = 1E-24
   Ds_opt2_4d(ld-1,:,jz) = 1E-24
   Ds_opt2_4d(:,:,jz)=max(real(1E-24,rprec),real(Ds_opt2_4d(:,:,jz),rprec))

   Beta_s(:,:,jz)=(Ds_opt2_4d(:,:,jz)/Ds_opt2_2d(:,:,jz))**(log(tf1)/(log(tf2)-log(tf1)))
   
   counter1=0      
   counter2=0
   
   !--MPI: this is not valid
   if ( ((.not. USE_MPI) .or. (USE_MPI .and. coord == nproc-1)) .and. (jz == nz) ) then
      Beta_s(:,:,jz)=1._rprec 
   end if

   !COMPUTING Ds_opt2
   Betaclip_s(:,:,jz)=max(real(Beta_s(:,:,jz),rprec),real(1._rprec/(tf1*tf2),rprec))
   Ds_opt2(:,:,jz)=Ds_opt2_2d(:,:,jz)/Betaclip_s(:,:,jz)
   Ds_opt2(ld,:,jz) = 1E-24
   Ds_opt2(ld-1,:,jz) = 1E-24
   Ds_opt2(:,:,jz)=max(real(1E-24,rprec),real(Ds_opt2(:,:,jz),rprec))

   Beta_avg_s(jz)=sum(Ds_opt2_2d(1:nx,1:ny,jz))/sum(Ds_opt2(1:nx,1:ny,jz))
   Betaclip_avg_s(jz)=sum(Ds_opt2_4d(1:nx,1:ny,jz))/sum(Ds_opt2_2d(1:nx,1:ny,jz))

   !FIXING TO 0 THE DS COEFFICIENT
   IF(CLIP_SGS)THEN
      DO JX=1,LD
      DO JY=1,NY
         IF((F_NN(JX,JY,JZ).lt.1E-24_rprec).OR.(F_YY(JX,JY,JZ).lt.1E-24_rprec) &
                                           .OR.(F_XX(JX,JY,JZ).lt.1E-24_rprec) &
                                           .OR.(F_MM(JX,JY,JZ).lt.1E-24_rprec) )THEN 
            Ds_opt2(JX,JY,JZ)=Co**2/1.0_rprec !(was /pr)
         ENDIF
      ENDDO
      ENDDO
   ENDIF
   
   
   
end do ! this ends the main jz=1-nz loop      

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

if (VERBOSE) write (*, *) 'finished lagrange_Sdep'

end subroutine lagrange_Sdep_Scalar


END MODULE LASD
