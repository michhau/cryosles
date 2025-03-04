! this is the w-node version
subroutine interpolag_Sdep_Scalar(F_KX,F_XX,F_PY,F_YY)
! This subroutine computes the values of F_LM and F_MM, F_KX and F_XX 
! at positions (x-u*dt) for use in the subroutine lagrangian
use types,only:rprec
use param
use sgsmodule
use sim_param,only:dt
use stretch

implicit none
integer :: jx, jy, jz 
integer :: jjx, jjy, jjz
integer :: jz_min


!--saves here: force heap storage
real(kind=rprec), save, dimension(ld,ny,0:nz_max) :: xp, yp, zp ! save attribute only works with dimensions declared as parameters
real(kind=rprec), save, dimension(ld,ny,0:nz_max) :: u_temp, v_temp !,theta_temp
real(kind=rprec) :: frac_x,frac_y,frac_z
real(kind=rprec) :: comp_x,comp_y,comp_z
real(kind=rprec), save, dimension(nx+2,ny+2,nz_max+2) :: FF_KX,FF_XX
real(kind=rprec), save, dimension(nx+2,ny+2,nz_max+2) :: FF_PY,FF_YY
integer :: addx, addy, addz
integer :: jxaddx, jyaddy, jzaddz
!real(kind=rprec), dimension(nx+2,ny+2,nz+2) :: Beta_t
real(kind=rprec),dimension(ld,ny,nz)::F_KX,F_XX,F_PY,F_YY
integer :: kk
real(rprec) :: fac_1,fac_2,fac_3
!---------------------------------------------------------------------

!     creates dummy arrays FF_LM and FF_MM to use in the subroutine
 
do jz=1,nz
   do jy=1,ny
      do jx=1,nx
         FF_KX(jx+1,jy+1,jz+1) = F_KX(jx,jy,jz)
         FF_XX(jx+1,jy+1,jz+1) = F_XX(jx,jy,jz)
         FF_PY(jx+1,jy+1,jz+1) = F_PY(jx,jy,jz)
         FF_YY(jx+1,jy+1,jz+1) = F_YY(jx,jy,jz)
      end do
   end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This is a bit like witch craft but the following lines do take care
!of all the edges including the corners

!a) For FF_KX:
!-----------------------

FF_KX(1,:,:) = FF_KX(nx+1,:,:)
FF_KX(nx+2,:,:) = FF_KX(2,:,:)

FF_KX(:,1,:) = FF_KX(:,ny+1,:) 
FF_KX(:,ny+2,:) = FF_KX(:,2,:)
 
  
#if(MPI)
   call mpi_sendrecv (FF_KX(1, 1, nz), (nx+2)*(ny+2), MPI_RPREC, up, 1,   &
                     FF_KX(1, 1, 1), (nx+2)*(ny+2), MPI_RPREC, down, 1,  &
                     comm, status, ierr)
#endif
if(coord==0) FF_KX(:,:,1)=FF_KX(:,:,2)
  
#IF(MPI)
call mpi_sendrecv (FF_KX(1, 1, 3), (nx+2)*(ny+2), MPI_RPREC, down, 2,   &
                     FF_KX(1, 1, nz+2), (nx+2)*(ny+2), MPI_RPREC, up, 2,  &
                     comm, status, ierr)
#endif
if((.not.USE_MPI).or.(USE_MPI.and.coord==nproc-1)) FF_KX(:,:,nz+2)=FF_KX(:,:,nz+1)

!b) For FF_XX:
!-----------------------

FF_XX(1,:,:) = FF_XX(nx+1,:,:)
FF_XX(nx+2,:,:) = FF_XX(2,:,:)

FF_XX(:,1,:) = FF_XX(:,ny+1,:) 
FF_XX(:,ny+2,:) = FF_XX(:,2,:) 

#if(MPI)
   call mpi_sendrecv (FF_XX(1, 1, nz), (nx+2)*(ny+2), MPI_RPREC, up, 3,   &
                     FF_XX(1, 1, 1), (nx+2)*(ny+2), MPI_RPREC, down, 3,  &
                     comm, status, ierr)
#endif
if(coord==0) FF_XX(:,:,1)=FF_XX(:,:,2)
#IF(MPI)
call mpi_sendrecv (FF_XX(1, 1, 3), (nx+2)*(ny+2), MPI_RPREC, down, 4,   &
                     FF_XX(1, 1, nz+2), (nx+2)*(ny+2), MPI_RPREC, up, 4,  &
                     comm, status, ierr)
#endif
if((.not.USE_MPI).or.(USE_MPI.and.coord==nproc-1)) FF_XX(:,:,nz+2)=FF_XX(:,:,nz+1)


!C) For FF_PY:
!-----------------------
FF_PY(1,:,:) = FF_PY(nx+1,:,:)
FF_PY(nx+2,:,:) = FF_PY(2,:,:)

FF_PY(:,1,:) = FF_PY(:,ny+1,:) 
FF_PY(:,ny+2,:) = FF_PY(:,2,:)

#if(MPI)
   call mpi_sendrecv (FF_PY(1, 1, nz), (nx+2)*(ny+2), MPI_RPREC, up, 5,   &
                     FF_PY(1, 1, 1), (nx+2)*(ny+2), MPI_RPREC, down, 5,  &
                     comm, status, ierr)
#endif
if(coord==0) FF_PY(:,:,1)=FF_PY(:,:,2)
#IF(MPI)
call mpi_sendrecv (FF_PY(1, 1, 3), (nx+2)*(ny+2), MPI_RPREC, down, 6,   &
                     FF_PY(1, 1, nz+2), (nx+2)*(ny+2), MPI_RPREC, up, 6,  &
                     comm, status, ierr)
#endif
if((.not.USE_MPI).or.(USE_MPI.and.coord==nproc-1)) FF_PY(:,:,nz+2)=FF_PY(:,:,nz+1)


!d) For FF_YY:
!-----------------------
FF_YY(1,:,:) = FF_YY(nx+1,:,:)
FF_YY(nx+2,:,:) = FF_YY(2,:,:)

FF_YY(:,1,:) = FF_YY(:,ny+1,:) 
FF_YY(:,ny+2,:) = FF_YY(:,2,:) 

#if(MPI)
   call mpi_sendrecv (FF_YY(1, 1, nz), (nx+2)*(ny+2), MPI_RPREC, up, 7,   &
                     FF_YY(1, 1, 1), (nx+2)*(ny+2), MPI_RPREC, down, 7,  &
                     comm, status, ierr)
#endif
if(coord==0) FF_PY(:,:,1)=FF_YY(:,:,2)
#IF(MPI)
call mpi_sendrecv (FF_YY(1, 1, 3), (nx+2)*(ny+2), MPI_RPREC, down, 8,   &
                     FF_YY(1, 1, nz+2), (nx+2)*(ny+2), MPI_RPREC, up, 8,  &
                     comm, status, ierr)
#endif
if((.not.USE_MPI).or.(USE_MPI.and.coord==nproc-1)) FF_YY(:,:,nz+2)=FF_YY(:,:,nz+1)

! end of witch craft
! END of Boundary Conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!puts u_lag and and v_lag on cs nodes
u_temp(:,:,0:nz) = u_lag/real(cs_count,kind=rprec)
v_temp(:,:,0:nz) = v_lag/real(cs_count,kind=rprec)

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
   u_lag (:,:,1) = u_temp(:,:,1)
   v_lag (:,:,1) = v_temp(:,:,1)
   jz_min = 2
else
   jz_min=1
endif

!$OMP PARALLEL DO
do jz = jz_min, nz  !--uv_lag at 0 level not needed (W-node)
   kk = jz + eta_displ(coord)

  fac_1=get_z_local(2*kk-1)
  fac_2=get_z_local(2*kk-3)
  fac_3=get_z_local(2*kk-2)

   u_lag(1:nx,:,jz)=u_temp(1:nx,:,jz-1) + ((u_temp(1:nx,:,jz)-u_temp(1:nx,:,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)
   v_lag(1:nx,:,jz)=v_temp(1:nx,:,jz-1) + ((v_temp(1:nx,:,jz)-v_temp(1:nx,:,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)

   !u_lag(1:nx,:,jz) = 0.5_rprec*(u_temp(1:nx,:,jz) + u_temp(1:nx,:,jz-1))
   !v_lag(1:nx,:,jz) = 0.5_rprec*(v_temp(1:nx,:,jz) + v_temp(1:nx,:,jz-1))
end do
!$OMP END PARALLEL DO

w_lag = w_lag/real(cs_count,kind=rprec)
if(coord==0) w_lag (:,:,1) = 0.25_rprec*w_lag (:,:,2)

!     computes the 3-D inverse displacement arrays that describe
!     the location where the point was at the previous step

do jz = 1, nz
  kk = jz+eta_displ(coord)
  xp(1:nx,:,jz) = -u_lag(1:nx,:,jz)*dt*real(cs_count,kind=rprec)/dx  !! is u already normalized
  yp(1:nx,:,jz) = -v_lag(1:nx,:,jz)*dt*real(cs_count,kind=rprec)/dy
  zp(1:nx,:,jz) = -w_lag(1:nx,:,jz)*dt*real(cs_count,kind=rprec)/(get_z_local(2*kk-1)-get_z_local(2*kk-3))
end do

!xp = -u_lag*DT*real(cs_count,kind=rprec)/dx    !! is u already normalized
!yp = -v_lag*DT*real(cs_count,kind=rprec)/dy   
!zp = -w_lag*DT*real(cs_count,kind=rprec)/dz


  !	because first plane is on u,v,p nodes
  !	this corrects for the fact that the first cell
  !	in the z direction has height dz/2
  !	it doubles the zp fraction if this fraction relates to the cell
  !     beneath it
if(coord==0)then
do jx= 1,nx
do jy=1,ny
   zp(jx,jy,2) = zp(jx,jy,2) + min (zp(jx,jy,2), 0._rprec)
   zp(jx,jy,1) = zp(jx,jy,1) + max (zp(jx,jy,1), 0._rprec)
end do
end do

zp(:,:,2) = min (1._rprec, zp(:,:,2))
zp(:,:,1) = min (1._rprec, zp(:,:,2))
zp(:,:,2) = max (-1._rprec, zp(:,:,2))
zp(:,:,1) = max (-1._rprec, zp(:,:,2))

endif

!if(mod(jt,wbase).eq.0) print*,'Lagrangian CFL condition for coord= ',coord,'=',  &
!                         maxval(abs(xp(1:nx, :, 1:nz)))


do jz=1,nz
   jjz = jz+1
   do jy=1,ny
      jjy = jy+1
      do jx=1,nx
         jjx = jx+1
         !     the are the values to add to the indices jx, jy, and jz
         !     the are +1 or -1 depending on what cube should be used for interpolation
         addx = int(sign(1._rprec,xp(jx,jy,jz)))
         addy = int(sign(1._rprec,yp(jx,jy,jz)))
         addz = int(sign(1._rprec,zp(jx,jy,jz)))
         jxaddx = jjx + addx 
         jyaddy = jjy + addy
         jzaddz = jjz + addz
         !     computes the relative weights given to F_** in the cube depending on point location
         comp_x = abs(xp(jx,jy,jz))
         comp_y = abs(yp(jx,jy,jz))
         comp_z = abs(zp(jx,jy,jz)) 
         frac_x = 1._rprec - comp_x
         frac_y = 1._rprec - comp_y
         frac_z = 1._rprec - comp_z

         ! The interpolation works in such a way that takes the values of the F_** function
         ! around the actual point (i,j,k), and multiplies it with a weight function. This weight function
         ! is the result of integrating over a streamline considering Taylor's Hypothesis. (Depending on how
         ! far is your previous function value, the weight is bigger or smaller. The weight is based in spatial
         ! distance).


         !For the Scalars:
         !----------------

         !     computes interpolated F_KX

         F_KX(jx,jy,jz)=frac_x*frac_y*&
              (FF_KX(jjx,jjy,jjz)*frac_z+FF_KX(jjx,jjy,jzaddz)*comp_z)&
              + frac_x*comp_y*&
              (FF_KX(jjx,jyaddy,jjz)*frac_z+FF_KX(jjx,jyaddy,jzaddz)*comp_z)&
              + comp_x*frac_y*&
              (FF_KX(jxaddx,jjy,jjz)*frac_z+FF_KX(jxaddx,jjy,jzaddz)*comp_z)&
              + comp_x*comp_y*&
              (FF_KX(jxaddx,jyaddy,jjz)*frac_z+FF_KX(jxaddx,jyaddy,jzaddz)*comp_z)
 
         !     computes interpolated F_XX
         F_XX(jx,jy,jz)=frac_x*frac_y*&
              (FF_XX(jjx,jjy,jjz)*frac_z+FF_XX(jjx,jjy,jzaddz)*comp_z)&
              + frac_x*comp_y*&
              (FF_XX(jjx,jyaddy,jjz)*frac_z+FF_XX(jjx,jyaddy,jzaddz)*comp_z)&
              + comp_x*frac_y*&
              (FF_XX(jxaddx,jjy,jjz)*frac_z+FF_XX(jxaddx,jjy,jzaddz)*comp_z)&
              + comp_x*comp_y*&
              (FF_XX(jxaddx,jyaddy,jjz)*frac_z+FF_XX(jxaddx,jyaddy,jzaddz)*comp_z)

         !     computes interpolated F_PY
         F_PY(jx,jy,jz)=frac_x*frac_y*&
              (FF_PY(jjx,jjy,jjz)*frac_z+FF_PY(jjx,jjy,jzaddz)*comp_z)&
              + frac_x*comp_y*&
              (FF_PY(jjx,jyaddy,jjz)*frac_z+FF_PY(jjx,jyaddy,jzaddz)*comp_z)&
              + comp_x*frac_y*&
              (FF_PY(jxaddx,jjy,jjz)*frac_z+FF_PY(jxaddx,jjy,jzaddz)*comp_z)&
              + comp_x*comp_y*&
              (FF_PY(jxaddx,jyaddy,jjz)*frac_z+FF_PY(jxaddx,jyaddy,jzaddz)*comp_z)
 
         !     computes interpolated F_YY
         F_YY(jx,jy,jz)=frac_x*frac_y*&
              (FF_YY(jjx,jjy,jjz)*frac_z+FF_YY(jjx,jjy,jzaddz)*comp_z)&
              + frac_x*comp_y*&
              (FF_YY(jjx,jyaddy,jjz)*frac_z+FF_YY(jjx,jyaddy,jzaddz)*comp_z)&
              + comp_x*frac_y*&
              (FF_YY(jxaddx,jjy,jjz)*frac_z+FF_YY(jxaddx,jjy,jzaddz)*comp_z)&
              + comp_x*comp_y*&
              (FF_YY(jxaddx,jyaddy,jjz)*frac_z+FF_YY(jxaddx,jyaddy,jzaddz)*comp_z)

      end do
   end do
end do

u_lag = 0._rprec 
v_lag = 0._rprec
w_lag = 0._rprec


end subroutine interpolag_Sdep_Scalar
