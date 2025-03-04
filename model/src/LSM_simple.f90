!------------------------------------------------------------------------------!
! MODULE: ls_model
!
! MAIN DEVELOPERS: Francesco Comola (francesco.comola.gmail.com)
!                  Varun Sharma (varun.sharma@epfl.ch)
! CONTRIBUTORS: Armin Sigmund (armin.sigmund@epfl.ch)
!               Daniela Brito Melo (daniela.britomelo@epfl.ch)
!
! LAST UPDATE: May 19, 2021
!
! DESCRIPTION: This module contains the Lagrangian Stochastic Model (LSM)
!
!------------------------------------------------------------------------------! 

module ls_model
!------------------------------------------------------------------------------!
use types
use param
use sim_param
use sgsmodule, only: cs_opt2
use sgs_stresses, only: calc_Sij
use service,only: mag  !, interp_var
use sp_sc1
use sp_sc2
!use sp_sc3
use iso_fortran_env, only: int64
use ifport
use stretch
!------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------!
integer                                    :: np,lsm_step
real(rprec)                                :: lsm_dt, disp_length, delta, entrain_dt
real(rprec)                                :: ve_flux, vd_flux, vs_flux, av_conc
real(rprec), dimension(:), allocatable     :: u_sgs, v_sgs, w_sgs,    &
                                              x_p, y_p, z_p,          &
                                              u_p, v_p, w_p,          &
                                              ax_p, ay_p, az_p,       &
                                              mass_p, sigma2, d_p,T_p, Q_p
real(rprec), dimension(:),allocatable :: re_p,sh_p,nu_p,timescale_stokes
real(rprec), dimension(:),allocatable :: delta_r,delta_m,delta_T
real(rprec), dimension(:,:), allocatable   :: ais, lsm_surface, tau, entrainment
real(rprec), dimension(:,:,:), allocatable :: sum_strain,sum_cs,sum_diss, &
                                              Fx_lsm, Fy_lsm, Fz_lsm,Ft_lsm_1,Ft_lsm_2,Ft_lsm,Fq_lsm,      &
                                              conc, xis, nis, sis, h_flux, v_flux, ppn, c_h_flux, c_v_flux
!real(rprec), dimension(:,:,:), allocatable :: u_regroup,v_regroup,w_regroup,T_regroup,Q_regroup
!real(rprec), dimension(:,:,:), allocatable :: diss,sgstke,restke
logical, parameter                         :: screen = .false.
logical, dimension(:), allocatable         :: p_exist, p_exchange
real(rprec), dimension(:,:,:),allocatable :: sum_u2,sum_u,sum_v2,sum_v,tke_u,tke_v
real(rprec), dimension(:,:,:),allocatable :: sum_w2,sum_w,tke_w
integer(kind=int64),dimension(:),allocatable :: p_ip
integer(kind=int64) :: np_global
real(rprec),allocatable,dimension(:),public :: exchange_x_p,exchange_y_p,exchange_z_p
real(rprec),allocatable,dimension(:),public :: exchange_u_p,exchange_v_p,exchange_w_p
real(rprec),allocatable,dimension(:),public :: exchange_u_sgs_p,exchange_v_sgs_p,exchange_w_sgs_p,exchange_sigma_p
real(rprec),allocatable,dimension(:),public :: exchange_ax_p,exchange_ay_p,exchange_az_p
real(rprec),allocatable,dimension(:),public :: exchange_d_p,exchange_mass_p
real(rprec),allocatable,dimension(:),public :: exchange_T_p
integer(kind=int64),allocatable,dimension(:),public :: exchange_np_global
integer :: ene_cons_ok,ene_cons_no,mom_cons_ok,mom_cons_no
real(rprec),dimension(:,:),allocatable,public :: tpxz,tpyz
real(rprec) :: averaging_weight
real(rprec),parameter :: averaging_time = 1000.0
real(rprec) :: sum_delta_m,sum_delta_q
!------------------------------------------------------------------------------!
contains

!------------------------------------------------------------------------------!
!                            ALLOCATE_LSM                                      !
!------------------------------------------------------------------------------!
subroutine allocate_lsm()

implicit none

!Forcings onto the LES flow field ( NOTE that the variables are 3D global arrays )
allocate(Fx_lsm(ld, ny,0:nz)) ! force in the x direction
allocate(Fy_lsm(ld, ny,0:nz)) ! force in the y direction
allocate(Fz_lsm(ld, ny,0:nz)) ! force in the z direction
IF(TEMPERATURE) then
  allocate(Ft_lsm(ld,ny,0:nz)) ! thermal forcing
  allocate(Ft_lsm_1(ld,ny,0:nz))
  allocate(Ft_lsm_2(ld,ny,0:nz))
ENDIF
IF(HUMIDITY) then
  allocate(Fq_lsm(ld,ny,0:nz)) ! vapor forcing
ENDIF

! Building three dimensional arrays for all principal flow variables
!allocate(u_regroup(ld,ny,0:nz)) ! u vel
!allocate(v_regroup(ld,ny,0:nz)) ! v vel
!allocate(w_regroup(ld,ny,0:nz)) ! w vel
allocate(diss(ld,ny,0:nz)) ! dissipation - from average_les()
allocate(restke(ld,ny,0:nz)) ! resolved tke - from average_les()
allocate(sgstke(ld,ny,0:nz)) ! subgrid-scale tke - from average_les()
!allocate(phi(ld,ny,0:nz)) ! phi function for immersed boundary
!IF(TEMPERATURE) THEN
! allocate(T_regroup(ld,ny,0:nz)) ! temperature
!endif
!IF(HUMIDITY) THEN
! allocate(Q_regroup(ld,ny,0:nz)) ! specific humidity
!endif
allocate(h_flux(ld,ny,0:nz))
allocate(v_flux(ld,ny,0:nz))
allocate(ppn(ld,ny,0:nz))
allocate(conc(ld, ny, 0:nz))
allocate(c_h_flux(ld,ny,0:nz))
allocate(c_v_flux(ld,ny,0:nz))

! Variables needed for average_les()
allocate(sum_u2(ld,ny,0:nz))
allocate(sum_v2(ld,ny,0:nz))
allocate(sum_u(ld,ny,0:nz))
allocate(sum_v(ld,ny,0:nz))
allocate(tke_u(ld,ny,0:nz))
allocate(tke_v(ld,ny,0:nz))
allocate(sum_w2(ld,ny,0:nz))
allocate(sum_w(ld,ny,0:nz))
allocate(tke_w(ld,ny,0:nz))
!allocate(diss(ld, ny, 0:nz))
!allocate(restke(ld, ny, 0:nz))
!allocate(sgstke(ld, ny, 0:nz))
allocate(sum_strain(ld, ny, 0:nz))
allocate(sum_cs(ld,ny,0:nz))
allocate(sum_diss(ld,ny,0:nz))

! variables needed for surface processes
allocate(lsm_surface(ld, ny))
allocate(entrainment(ld, ny))
allocate(tau(ld, ny))
allocate(xis(3,nx,ny))
allocate(nis(3,nx,ny))
allocate(sis(2,nx,ny))
allocate(ais(nx,ny))

allocate(tpxz(ld,ny))
allocate(tpyz(ld,ny))

! parcticle quantities
allocate(u_sgs(np_max)) ! sgs velocity in x
allocate(v_sgs(np_max)) ! sgs velocity in y
allocate(w_sgs(np_max)) ! sgs velocity in z
allocate(x_p(np_max)) ! x loc (global)
allocate(y_p(np_max)) ! y loc (global)
allocate(z_p(np_max)) ! z loc (global)
allocate(u_p(np_max)) ! vel in x direction
allocate(v_p(np_max)) ! vel in y direction
allocate(w_p(np_max)) ! vel in z direction
allocate(ax_p(np_max)) ! acc. in x direction
allocate(ay_p(np_max)) ! acc. in y direction
allocate(az_p(np_max)) ! acc. in z direction
allocate(mass_p(np_max)) ! mass of particle
allocate(d_p(np_max)) ! diameter of particle
allocate(sigma2(np_max)) ! noise for sgs vel of the particle
IF(TEMPERATURE) then
 allocate(T_p(np_max)) ! temperature of the particle
ENDIF
IF(HUMIDITY) then
 allocate(Q_p(np_max))
ENdIF

allocate(p_ip(np_max)) ! particle id
allocate(p_exist(np_max)) ! a check for the particle
allocate(p_exchange(np_max))

! non_dimensional numbers
allocate(re_p(np_max))
allocate(sh_p(np_max))
allocate(nu_p(np_max))
! stokes timescale
allocate(timescale_stokes(np_max))
!increments in mass and temperature of drop
allocate(delta_r(np_max))
allocate(delta_m(np_max))
allocate(delta_T(np_max))

end subroutine allocate_lsm

!------------------------------------------------------------------------------!
!                             DEALLOCATE_LSM                                   !
!------------------------------------------------------------------------------!
subroutine deallocate_lsm()

implicit none

deallocate(Fx_lsm) ! force in the x direction
deallocate(Fy_lsm) ! force in the y direction
deallocate(Fz_lsm) ! force in the z direction
IF(TEMPERATURE) then
  deallocate(Ft_lsm)! thermal forcing
  deallocate(Ft_lsm_1)
  deallocate(Ft_lsm_2)
ENDIF
IF(HUMIDITY) then
  deallocate(Fq_lsm) ! vapor forcing
ENDIF

! Building three dimensional arrays for all principal flow variables
!deallocate(u_regroup) ! u vel
!deallocate(v_regroup) ! v vel
!deallocate(w_regroup) ! w vel
deallocate(diss) ! dissipation - from average_les()
deallocate(restke) ! resolved tke - from average_les()
deallocate(sgstke) ! subgrid-scale tke - from average_les()
!!deallocate(phi) ! phi function for immersed boundary
!IF(TEMPERATURE) THEN
! deallocate(T_regroup) ! temperature
!endif
!IF(HUMIDITY) THEN
! deallocate(Q_regroup) ! specific humidity
!endif

! Variables needed for average_les()
deallocate(sum_u2)
deallocate(sum_v2)
deallocate(sum_u)
deallocate(sum_v)
deallocate(tke_u)
deallocate(tke_v)
deallocate(sum_w2)
deallocate(sum_w)
deallocate(tke_w)
!deallocate(diss)
!deallocate(restke)
!deallocate(sgstke)
deallocate(sum_strain)
deallocate(sum_cs)
deallocate(sum_diss)
! variables needed for surface processes
deallocate(lsm_surface)
deallocate(entrainment)
deallocate(tau)
deallocate(h_flux)
deallocate(v_flux)
deallocate(c_h_flux)
deallocate(c_v_flux)

deallocate(xis)
deallocate(nis)
deallocate(sis)
deallocate(ais)

deallocate(conc)

! parcticle quantities
deallocate(u_sgs) ! sgs velocity in x
deallocate(v_sgs) ! sgs velocity in y
deallocate(w_sgs) ! sgs velocity in z
deallocate(x_p) ! x loc (global)
deallocate(y_p) ! y loc (global)
deallocate(z_p) ! z loc (global)
deallocate(u_p) ! vel in x direction
deallocate(v_p) ! vel in y direction
deallocate(w_p) ! vel in z direction
deallocate(ax_p) ! acc. in x direction
deallocate(ay_p) ! acc. in y direction
deallocate(az_p) ! acc. in z direction
deallocate(mass_p) ! mass of particle
deallocate(d_p) ! diameter of particle
deallocate(sigma2) ! noise for sgs vel of the particle
IF(TEMPERATURE) then
 deallocate(T_p) ! temperature of the particle
ENDIF
IF(HUMIDITY) then
 deallocate(Q_p)
ENDIF

deallocate(p_ip) ! particle id
deallocate(p_exist) ! a check for the particle
deallocate(p_exchange)

! non_dimensional numbers
deallocate(re_p)
deallocate(sh_p)
deallocate(nu_p)
! stokes timescale
deallocate(timescale_stokes)
!increments in mass and temperature of drop
deallocate(delta_r)
deallocate(delta_m)
deallocate(delta_T)

deallocate(tpxz)
deallocate(tpyz)

end subroutine deallocate_lsm

!------------------------------------------------------------------------------!
!                            INITIALIZE_LSM                                    !
!------------------------------------------------------------------------------!
subroutine initialize_lsm()

implicit none

integer                                :: i_seed, ix, iy, iz
integer, dimension(:), allocatable     :: a_seed
integer, dimension(1:8)                :: dt_seed
integer :: folder_exists, file_exists
character(64) :: filename

lsm_dt      = ref_dt*lsm_freq

if(jt_total .lt. lsm_init) then
lsm_step    = 0
elseif( jt_total .gt. lsm_init) then
lsm_step = (jt_total-lsm_init)/lsm_freq + 1
endif

lsm_path=trim(out_path)//'/lsm/'

write(*,*) 'PPP_MIN,PPP_MAX:',PPP_MIN,PPP_MAX

np          = 0
!np_global   = 0
!delta       = (dx*dy*dz)**(1._rprec/3._rprec)

!delta = filter_size*(dx*dy*dz(jz+eta_displ(coord)))**(1._rprec/3._rprec)

!if (llbc_special .eq. 'IBM_DFA') then
!   disp_length = dz/2.0_rprec
!else
   disp_length = 4.0_rprec * d_m * (1.0_rprec/z_i)
!end if

!disp_length = dz/2._rprec
!disp_length = 4._rprec*d_m

write(*,*) disp_length,lbc_dz0

ve_flux = 0._rprec
vd_flux = 0._rprec
vs_flux = 0._rprec
av_conc = 0._rprec


if(coord==0)then
      !create output folder
      print*,'Creating folder for output...'
      folder_exists=system('mkdir '//trim(lsm_path))
      if(folder_exists.eq.0)then
         print*,'Created output folder!'
      else
         print*,'Output already exists!'
         print*,'Result=',folder_exists
      endif

      inquire(file=trim(lsm_path)//'mass_balance.txt', EXIST=file_exists)
      if (.not. file_exists) then
         open(300,file=trim(lsm_path)//'mass_balance.txt',status='new')
         write(300,'(a115)') "Time(s)   u_star2(m2/s2)   entrainment(kg/m2/s)   deposition(kg/m2/s)   ejection(kg/m2/s)   suspension(kg/m3)"
         close(300)
      endif
endif
call mpi_barrier(comm,ierr)

write(filename,'(a,i0,a)') trim(lsm_path)//'mom_balance_c',coord,'.txt'
inquire(file=trim(filename), EXIST=file_exists)
if (.not. file_exists) then
   open(300,file=trim(filename),status='new')
   write(300,'(a70)') "Time(s)   px_forcing(m2/s2)   pz_forcing(m2/s2)   u_star2(m2/s2)"
   close(300)
endif

print*, "MODE=", mode
print*, "PARTICLE TYPE=", p_type
print*, "INITIALIZING LSM. DT_LSM=", lsm_dt,"s"
print*, "LSM EVERY", lsm_freq, "STEPS"

! conc, Fx, Fy, Fz are set to 0 at every step, so no need to initialize here
u_sgs(:)          = 0._rprec
v_sgs(:)          = 0._rprec
w_sgs(:)          = 0._rprec
x_p(:)            = 0._rprec
y_p(:)            = 0._rprec
z_p(:)            = 0._rprec
u_p(:)            = 0._rprec
v_p(:)            = 0._rprec
w_p(:)            = 0._rprec
ax_p(:)           = 0._rprec
ay_p(:)           = 0._rprec
az_p(:)           = 0._rprec
mass_p(:)         = 0._rprec
d_p(:)            = 0._rprec
sigma2(:)         = 0._rprec
IF(TEMPERATURE) then
 T_p(:) = 0.0_rprec
ENDIF
IF(HUMIDITY) then
 Q_p(:) = 0.0_rprec
ENDIF

p_exist(:)        = .false.
p_ip(:) = 0

sum_strain(:,:,:) = 0._rprec
sum_cs(:,:,:) = 0.0_rprec
sum_diss(:,:,:)= 0.0_rprec
entrainment(:,:)  = 0._rprec
h_flux(:,:,:)       = 0._rprec
v_flux(:,:,:)       = 0._rprec
ppn(:,:,:) = 1.E-30
c_h_flux(:,:,:) = 0.0_rprec
c_v_flux(:,:,:) = 0.0_rprec

!u_regroup(:,:,:) = 0.0_rprec
!v_regroup(:,:,:) = 0.0_rprec
!w_regroup(:,:,:) = 0.0_rprec
diss(:,:,:) = 0.0_rprec
restke(:,:,:) = 0.0_rprec
sgstke(:,:,:) = 0.0_rprec
phi(:,:,:) = 0.0_rprec

Fx_lsm(:,:,:) = 0.0_rprec
Fy_lsm(:,:,:) = 0.0_rprec
Fz_lsm(:,:,:) = 0.0_rprec

IF(TEMPERATURE) then
! T_regroup(:,:,:) = 0.0_rprec
 Ft_lsm(:,:,:) = 0.0_rprec
 Ft_lsm_1(:,:,:) = 0.0_rprec
 Ft_lsm_2(:,:,:) = 0.0_rprec
ENDIF
IF(HUMIDITY) then
! Q_regroup(:,:,:) = 0.0_rprec
 Fq_lsm(:,:,:) = 0.0_rprec
ENDIF

tpxz(:,:) = 0.0_rprec
tpyz(:,:) = 0.0_rprec

sum_u2 = 0.0_rprec
sum_v2 = 0.0_rprec
sum_u = 0.0_rprec
sum_v = 0.0_rprec
sum_w2 = 0.0_rprec
sum_w = 0.0_rprec
tke_u = 0.0_rprec
tke_v = 0.0_rprec
tke_w = 0.0_rprec

diss = 0.0_rprec
restke = 0.0_rprec
sgstke = 0.0_rprec

re_p = 0.0_rprec
sh_p = 0.0_rprec
nu_p = 0.0_rprec
timescale_stokes = 0.0_rprec
delta_r = 0.0_rprec
delta_m = 0.0_rprec
delta_T = 0.0_rprec

averaging_weight = 1.0_rprec / (1.0_rprec + averaging_time)

call random_seed(size=i_seed)
allocate(a_seed(1:i_seed))
call random_seed(get=a_seed)
call date_and_time(values=dt_seed)
a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
call random_seed(put=a_seed)
deallocate(a_seed)

end subroutine initialize_lsm

!------------------------------------------------------------------------------!
!                           INITIALIZE_LSM_SURFACE                             !
!------------------------------------------------------------------------------!
subroutine initialize_lsm_surface

implicit none

integer        :: ix,iy,iz,ixm1,ixp1,iym1,iyp1
real(rprec)    :: temp, min_height
character(100) :: file_path
logical :: file_exists

do iz=0,nz
     phi(:,:,iz) = get_z_local(2*((iz-1) + eta_displ(coord))) !(iz-1)*dz + real(coord*(nz-1)*dz,rprec)
     norm(1,:,:,iz) = 0._rprec
     norm(2,:,:,iz) = 0._rprec
     norm(3,:,:,iz) = 1._rprec
end do

do iy=1,ny
   do ix=1,nx
      xis(1,ix,iy) = (ix-1)*dx
      xis(2,ix,iy) = (iy-1)*dy
      xis(3,ix,iy) = 0.0_rprec

      nis(1,ix,iy) = 0.0_rprec
      nis(2,ix,iy) = 0.0_rprec
      nis(3,ix,iy) = 1.0_rprec
   enddo
enddo

ais = dx*dy*z_i*z_i

do iy=1,ny
    do ix=1,nx
        ixm1 = modulo(ix-2,nx)+1
        iym1 = modulo(iy-2,ny)+1
        ixp1 = modulo(ix,nx)+1
        iyp1 = modulo(iy,ny)+1
        sis(1,ix,iy) = atan((xis(3,ixp1,iy)-xis(3,ixm1,iy))/(2._rprec*dx))
        sis(2,ix,iy) = atan((xis(3,ix,iyp1)-xis(3,ix,iym1))/(2._rprec*dy))
    end do
end do

if (r_dep) then
    file_path = rst_path//'/'//file_dep
    open(100,file=trim(file_path),status='unknown')
    read(100,*)
    do iy=1,ny
        do ix=1,nx
            read(100,*) temp, temp, lsm_surface(ix,iy)
        end do
    end do
    close(100)
else
        lsm_surface(:,:)  = init_dep
end if

if(coord .eq. 0) then
   write (file_path, '(a,a,i0.8,a)') trim(lsm_path)//'lsm_surface','_',jt_total,'.out'
   inquire(file=file_path,EXIST=file_exists)
   if(file_exists) then
     open(111,file=trim(file_path),status='unknown')
     write(*,*) 'Reading from:',trim(file_path)
     do ix=1,nx
       do iy=1,ny
        read(111,*) lsm_surface(ix,iy)
       enddo
     enddo
       close(111)
   endif
endif

if(mode .eq. 'saltation') then
  if(coord .eq. 0) then
    write (file_path, '(a,a,i0.8,a)') trim(lsm_path)//'entrainment','_',jt_total,'.out'
    inquire(file=file_path,EXIST=file_exists)
    if(file_exists) then
      open(111,file=trim(file_path),status='unknown')
      write(*,*) 'Reading from:',trim(file_path)
      do ix=1,nx
       do iy=1,ny
         read(111,*) entrainment(ix,iy)
       enddo
      enddo
        close(111)
    endif
  endif
endif

end subroutine initialize_lsm_surface

!------------------------------------------------------------------------------!
!                           INITIALIZE_PARCELS                                 !
!------------------------------------------------------------------------------!
subroutine initialize_parcels

implicit none

character(256) :: fname
integer :: dummy
integer :: i
logical :: file_exists
character(16) :: coord_str

  write(*,*) 'Inside initialize_parcels!!!'
  write(coord_str,'(i0)') coord

  write (fname, '(a,i0,a,i0,a,i0.8,a)') '/particles_',coord,'_',colour,'_',jt_total,'.out'
  fname = trim(lsm_path)//trim(fname)
  inquire(file=fname,EXIST=file_exists)
  if(file_exists) then
    write(*,*) 'READING PARTICLES FROM: ',fname
    call execute_command_line('cat '//fname//' | wc -l > tmp.txt'//coord_str)
    open(111,file='tmp.txt'//coord_str)
    read(111,*) np
    close(111)
    open(200,file=fname,form='formatted',action='read')
    do i=1,np
      IF(sublimation) then
        read(200,*) dummy,p_ip(i),x_p(i),y_p(i),z_p(i),u_p(i),v_p(i),w_p(i),T_p(i),d_p(i), &
                    re_p(i),nu_p(i),sh_p(i),mass_p(i),delta_m(i),delta_r(i),delta_T(i),timescale_stokes(i)
      else
        read(200,*) dummy,p_ip(i),x_p(i),y_p(i),z_p(i),u_p(i),v_p(i),w_p(i),ax_p(i),ay_p(i),az_p(i),d_p(i), &
                 re_p(i),mass_p(i),delta_m(i),delta_r(i),delta_T(i),timescale_stokes(i)
      endif
      p_exist(i) = .true.
      p_exchange(i) = .false.
    enddo
    close(200)
    call execute_command_line('rm tmp.txt'//coord_str)
  else
   np=0
  endif
  ! Exchange parcels between ranks after np is known for all ranks because particles files from a previous simulation were written before exchanging parcels
  call exchange_parcels()

end subroutine initialize_parcels

!------------------------------------------------------------------------------!
!                            ENTRAIN_PARCELS                                   !
!------------------------------------------------------------------------------!
subroutine entrain_parcels()

implicit none

integer                   :: ix, iy
real(rprec)               :: tau_thresh, tau_excess,        &
                             n_entrain, d_g
real(rprec)               :: mean, std,                     &
                             v_ang, h_ang,                  &
                             val_u, val_v, val_w, vel_fric, &
                             e_vel
real(rprec)               :: slope, temp_mass, dep_mass
real(rprec), dimension(3) :: coord_p
!real(rprec) :: rand_theta
!------------------------------------------------------------------------------!

if (screen) print*, "in entrain_parcels"
! Aerodynamic entrainment is only computed by the lower rank
if (coord .eq. 0) then
   ! Loop at the surface nodes
   do iy=1,ny
      do ix=1,nx
         ! Fluid threshold shear stress from Bagnold's expression
         tau_thresh              = 0.01_rprec*g*d_m*(p_rho-f_rho) ![Pa]
         ! Compute number of entrained particles from Anderson and Haff 1991
         tau_excess              = max(0._rprec, tau(ix,iy)-tau_thresh) ![Pa]
         n_entrain               = 1.5_rprec*tau_excess/(8._rprec*pi*(d_m**2.0_rprec)) ![Number of entrained particles/m2/s]: Entrainment coefficient from Doorschot and Lehning, 2002
         n_entrain               = (n_entrain*lsm_dt*ais(ix,iy)*z_i) !Number of entrained particles
         entrainment(ix,iy)      = entrainment(ix,iy) + n_entrain

         if (entrainment(ix,iy).gt.ppp_min) then
            ! Compute particle diameter in the parcel from distribution
            call                sample_distribution(d_m,d_s,'lognormal',d_g)
            d_g                 = max(d_min, min(d_max, d_g))            
            ! Eroded mass
            temp_mass           = entrainment(ix,iy)*p_rho*pi*(d_g**3.0_rprec)/6._rprec ![kg] 
            ! Deposited mass (available to be eroded)
            dep_mass            = lsm_surface(ix,iy)*ais(ix,iy) ![kg]
            ! Constrain temp_mass = min(temp_mass,dep_mass)
            if ((dep_mass.gt.0._rprec).and.(dep_mass.lt.temp_mass)) then
               temp_mass       = dep_mass
            end if
            if (temp_mass.le.dep_mass) then
               ! ENTRAIN PARCEL
               np = np + 1               ! Define new parcel index
               np_global = np_global + 1 ! Update total number of parcels
               if (np>np_max) then
                  print*,     "INCREASE NUMBER OF PARCELS (A)!"
                  stop
               end if
               ! Save new parcel data
               d_p(np)         = d_g
               mass_p(np)      = temp_mass
               p_exist(np)     = .true.
               u_sgs(np)       = 0._rprec
               v_sgs(np)       = 0._rprec
               w_sgs(np)       = 0._rprec
               x_p(np)         = xis(1,ix,iy) + disp_length*nis(1,ix,iy)
               y_p(np)         = xis(2,ix,iy) + disp_length*nis(2,ix,iy)
               z_p(np)         = xis(3,ix,iy) + disp_length*nis(3,ix,iy)
               p_ip(np)        = np_global
               ! Save mass of particles aloft
               av_conc         = av_conc+mass_p(np)
               ! Check if the domain borders were crossed
               call check_boundaries(np)
               if (p_exist(np)) then
                  ! Parcel location ([m], [m], global index)
                  coord_p     = (/ x_p(np), y_p(np), get_eta_local(z_p(np)) /)
                  ! Get flow velocity at particle location
                  call        interp_vel(coord_p,val_u,val_v,val_w,np)
                  ! Update deposited mass at the surface
                  call        update_lsm_surface(np,-1._rprec)
                  ! Save mass of entrained particles
                  ve_flux     = ve_flux+mass_p(np)               
               ! Parcel ejection angle
               h_ang       = atan2(val_v,val_u) ! Horizontal angle of ejection = horizontal flow direction
               mean        = (75._rprec-55._rprec*(1._rprec-exp(-d_p(np)/(175.E-06))))/180._rprec*pi
               std         = 15._rprec/180._rprec*pi
               call        sample_distribution(mean, std, 'lognormal', v_ang) ! Get vertical ejection angle (from Clifton and Lehning, 2008)
               slope       = sis(1,ix,iy)*cos(h_ang)+sis(2,ix,iy)*sin(h_ang)
               v_ang       = min(0.5_rprec*pi, max(-0.5_rprec*pi, v_ang+slope))
               ! Parcel ejection velocity
               vel_fric    = sqrt(tau(ix,iy)/f_rho)
               mean        = 3.5_rprec*vel_fric
               std         = 2.5_rprec*vel_fric
               call        sample_distribution(mean, std, 'lognormal', e_vel) ! Get ejection velocity (from Clifton and Lehning, 2008)
               u_p(np)     = e_vel*cos(v_ang)*cos(h_ang)
               v_p(np)     = e_vel*cos(v_ang)*sin(h_ang)
               w_p(np)     = e_vel*sin(v_ang)
               ! Parcel initial aceleration
               ax_p(np)    = 0._rprec
               ay_p(np)    = 0._rprec
               az_p(np)    = 0._rprec
               sigma2(np)  = 0._rprec ! (used in compute_sgs)
               IF (TEMPERATURE) then
                  !call random_number(rand_theta)
                  if (LLBC_SC1 .eq. 'DRC') then
                     T_p(np)     = (DLBC_SC1*T_SCALE)/T_SCALE !+ 0.002_rprec*DLBC_SC1*(2.0*rand_theta-1.0_rprec) !0.9_rprec
                  else
                     T_p(np) = 1.0_rprec
                  endif
                  !T_p(np) = DLBC_SC1 !0.9_rprec
               ENDIF
               endif !(p_exist(np))
            end if !(temp_mass.le.dep_mass)
            entrainment(ix,iy)  = 0._rprec
         end if !(entrainment(ix,iy).gt.ppp_min)
      end do !(ix=1,nx)
   end do !(iy=1,ny)
endif !(coord .eq. 0)

end subroutine entrain_parcels

!------------------------------------------------------------------------------!
!                               LINE_RELEASE                                   !
!------------------------------------------------------------------------------!
!subroutine line_release()
!
!implicit none
!
!real(rprec),parameter :: height_rel=0.40000
!real(rprec) :: g_coord_uvp
!integer :: proc_z
!real(rprec), dimension(3) :: point_rel, coord_p
!real(rprec)               :: d_g, mass_rel, mass_count
!integer :: i,j
!
!if(mod(jt,2000) .eq. 0) then
!
!   g_coord_uvp = (height_rel/dz) + 0.5_rprec
!   proc_z = ceiling(g_coord_uvp/real((nz-1),rprec)) - 1
!
!if(coord .eq. proc_z) then
!
! do j=1,ny
!   do i=1,nx
!    np=np+1
!    np_global = np_global+1
!    if (np.gt.np_max) then
!        print*, "INCREASE NUMBER OF PARCELS"
!        stop
!    end if
!
!    p_exist(np) = .true.
!    call sample_distribution(d_m, d_s, 'lognormal', d_g)
!    d_p(np) = max(min(d_g,d_max),d_min)
!    mass_p(np) = ppp_min*(p_rho * pi*(d_p(np)**3)/6._rprec)
!
!    x_p(np) = real(i,rprec)*dx ! 2.0_rprec*dx
!    y_p(np) = real(j,rprec)*dy  !point_rel(2)
!    z_p(np) = height_rel
!
!    call check_boundaries(np)
!    if (p_exist(np)) then
!
!        coord_p(1) = x_p(np)
!        coord_p(2) = y_p(np)
!        coord_p(3) = z_p(np) - (real(eta_displ(coord)*dz,rprec))
!
!        call interp_vel(coord_p,u_p(np),v_p(np),w_p(np),np)
!    end if
!    u_sgs(np) = 0._rprec
!    v_sgs(np) = 0._rprec
!    w_sgs(np) = 0._rprec
!
!    ax_p(np) = 0._rprec
!    ay_p(np) = 0._rprec
!    az_p(np) = 0._rprec
!
!    p_ip(np) = np_global
!    sigma2(np) = 0.0_rprec
!    if(temperature) then
!     T_p(np) = 0.0_rprec
!    endif
!
!  end do
! enddo
!endif ! MPI branch
!
!endif ! jt
!
!end subroutine line_release

!------------------------------------------------------------------------------!
!                               POINT_RELEASE                                  !
!------------------------------------------------------------------------------!
!subroutine point_release()
!
!implicit none
!
!real(rprec), parameter    :: mass_rate=6.1875*1.0E-01 !kg/s
!real(rprec), dimension(3) :: point_rel, coord_p
!real(rprec)               :: d_g, mass_rel, mass_count
!integer :: i,j
!
!do i=-3,3
!  do j=-3,3
!point_rel=(/ 0.5_rprec*l_x+real(i,rprec), 0.5_rprec*l_y+real(j,rprec), 1._rprec /)
!if (point_rel(3).gt.l_z) then
!    print*, "SOURCE POINT TOO HIGH"
!    stop
!end if
!mass_rel = mass_rate*lsm_dt ! kg released every lsm_dt
!
!mass_count=0._rprec
!do while (mass_count.lt.mass_rel)
!    np=np+1
!    if (np.gt.np_max) then
!        print*, "INCREASE NUMBER OF PARCELS"
!        stop
!    end if
!
!    p_exist(np) = .true.
!
!    call sample_distribution(d_m, d_s, 'lognormal', d_g)
!    d_p(np) = max(min(d_g,d_max),d_min)
!    mass_p(np) = ppp_min*(p_rho * pi*(d_p(np)**3)/6._rprec)
!    mass_count = mass_count+mass_p(np)
!
!    x_p(np) = point_rel(1)
!    y_p(np) = point_rel(2)
!    z_p(np) = point_rel(3)
!
!    call check_boundaries(np)
!    if (p_exist(np)) then
!        coord_p = (/x_p(np), y_p(np), z_p(np) /)
!        call interp_vel(coord_p,u_p(np),v_p(np),w_p(np),np)
!        u_p(np) = 0.0_rprec
!        v_p(np) = 0.0_rprec
!        w_p(np) = 4.0_rprec
!
!    end if
!    u_sgs(np) = 0._rprec
!    v_sgs(np) = 0._rprec
!    w_sgs(np) = 0._rprec
!
!    ax_p(np) = 0._rprec
!    ay_p(np) = 0._rprec
!    az_p(np) = 0._rprec
!end do
!
!enddo ! end of j loop
!enddo ! end of i loop
!
!end subroutine point_release

!------------------------------------------------------------------------------!
!                                    SNOWFALL                                  !
!------------------------------------------------------------------------------!
subroutine snowfall()

implicit none

integer                   :: ix, iy
real(rprec), parameter    :: precip=1.000_rprec ! mm/h
real(rprec), parameter    :: height_rel=800.00_rprec !dimension-less
real(rprec)               :: d_g, mass_rel, mass_g
real(rprec), dimension(3) :: coord_p
real(rprec) :: g_coord_uvp
integer :: proc_z

if (height_rel.gt.l_z) then
    print*, "SNOWFALL TOO HIGH"
    stop
end if

g_coord_uvp = (get_eta_local(height_rel)/dz) + 0.5_rprec
proc_z = count( REAL(NZ_CUMUL, RPREC) < g_coord_uvp )

if(proc_z .eq. coord) then

mass_rel = (precip/(3600._rprec*1000._rprec))*(dx*z_i)*(dy*z_i)*(lsm_dt)*p_rho ! kg released every lsm step from each node

do iy=1,ny
    do ix=1,nx
        mass_g = p_rho*pi*(d_m**3)/6._rprec

        entrainment(ix,iy) = entrainment(ix,iy)+mass_rel/mass_g
        if (entrainment(ix,iy).ge.ppp_min) then

            np=np+1
            np_global = np_global+1
            if (np.gt.np_max) then
                print*, "INCREASE NUMBER OF PARCELS (B)"
                stop
            end if

            p_exist(np) = .true.

            call sample_distribution(d_m,d_s,'lognormal',d_g)
            d_g = min(d_max,max(d_min,d_g))
            d_p(np) = d_g

            mass_p(np) = ppp_min*(p_rho*pi*(d_p(np)**3)/6._rprec)
            av_conc = av_conc+mass_p(np)
            ve_flux = ve_flux+mass_p(np)

            x_p(np) = real(ix-1)*dx
            y_p(np) = real(iy-1)*dy
            z_p(np) = height_rel

            p_ip(np) = np_global ! Define p_ip here because it is used in interp_vel for debugging reasons

            if (p_exist(np)) then
                  coord_p(1) = x_p(np)
                  coord_p(2) = y_p(np)
                  coord_p(3) = get_eta_local(z_p(np)) - (real(eta_displ(coord)*dz,rprec))
                  call interp_vel(coord_p,u_p(np),v_p(np),w_p(np),np)
            end if

            u_sgs(np) = 0._rprec
            v_sgs(np) = 0._rprec
            w_sgs(np) = 0._rprec

            sigma2(np) = 0.0_rprec

            ax_p(np) = 0._rprec
            ay_p(np) = 0._rprec
            az_p(np) = 0._rprec

            IF(TEMPERATURE) then
              T_p(np) = DLBC_SC1 !0.9_rprec
            ENDIF
            entrainment(ix,iy) = 0._rprec
        end if
    end do
end do

endif ! to place snowfall in the right coordinate
end subroutine snowfall

!------------------------------------------------------------------------------!
!                               COMPUTE_SGS                                    !
!------------------------------------------------------------------------------!
subroutine compute_sgs()

implicit none

integer                     :: ip
real(rprec), parameter      :: c_0 = 6._rprec ! source: Weil 2004
real(rprec)                 :: mean, std
real(rprec)                 :: dcsi, f_s
real(rprec)                 :: ds2dx, ds2dy, ds2dz, ds2dt
real(rprec)                 :: sgstke_p, restke_p, diss_p
real(rprec), dimension(3)   :: coord_p, x_diff
real(rprec), dimension(2,2,2) :: var_section
real(rprec) :: debug_a
real(rprec) :: gamma_p,gamma_l
integer :: i,j,k
integer, dimension(2) :: ix_close,iy_close,iz_close ! indexes for interpolation
logical :: at_surf

if (screen) print*, "In compute sgs"

mean = 0.0_rprec
std = sqrt(lsm_dt)

!$omp parallel do private(coord_p,f_s,diss_p,restke_p,sgstke_p,ds2dt,ds2dx,ds2dy,ds2dz,gamma_l,gamma_p,dcsi), &
!$omp& private(var_section,ix_close,iy_close,iz_close,x_diff,at_surf) &
!$omp& shared(sigma2,x_p,y_p,z_p,u_p,v_p,w_p,u_sgs,v_sgs,w_sgs,diss,restke,sgstke,np)
do ip=1,np
    !coord_p = (/x_p(ip),y_p(ip),z_p(ip)/)

     coord_p(1) = x_p(ip)
     coord_p(2) = y_p(ip)
     coord_p(3) = get_eta_local(z_p(ip)) - (real(eta_displ(coord)*dz,rprec))

    ! get indexes of closest W nodes
    call get_indexes_3d(coord_p,1._rprec,ix_close,iy_close,iz_close,x_diff,at_surf,1,ip)

    ! Interpolate between the 8 closest grid nodes
    ! diss
    var_section = diss(ix_close,iy_close,iz_close)
    call interp_var_lsm(var_section,diss_p,x_diff(1),x_diff(2),x_diff(3))
    ! restke
    var_section = restke(ix_close,iy_close,iz_close)
    call interp_var_lsm(var_section,restke_p,x_diff(1),x_diff(2),x_diff(3))
    ! sgstke
    var_section = sgstke(ix_close,iy_close,iz_close)
    call interp_var_lsm(var_section,sgstke_p,x_diff(1),x_diff(2),x_diff(3))

    f_s = sgstke_p/(restke_p+sgstke_p)
    ! compute stress and derivative of the stress
    ds2dt = sigma2(ip)
    sigma2(ip) = 2._rprec/3._rprec*sgstke_p
    ds2dt = (sigma2(ip)-ds2dt)/lsm_dt

    call interp_grad(coord_p,1.0_rprec,sgstke,ds2dx,ds2dy,ds2dz)
    ds2dx = 2._rprec/3._rprec * ds2dx
    ds2dy = 2._rprec/3._rprec * ds2dy
    ds2dz = 2._rprec/3._rprec * ds2dz

    gamma_l = 2.0_rprec * sigma2(ip)/(c_0*diss_p)
    if (trim(p_type).eq.'inertial') then
        gamma_p = gamma_l/sqrt(1.0_rprec + ((2.0_rprec*w_p(ip))**2.0_rprec)/sigma2(ip))
    else
        gamma_p = gamma_l
    end if
    ! update the sgs velocities of the flow
    call sample_distribution(mean, std, 'normal', dcsi)
    u_sgs(ip) = u_sgs(ip) - f_s*u_sgs(ip)*lsm_dt/gamma_p &
        + 0.5_rprec*(ds2dt*u_sgs(ip)/sigma2(ip)+ds2dx)*lsm_dt + sqrt(f_s*2.0_rprec*sigma2(ip)/gamma_p)*dcsi

    call sample_distribution(mean, std, 'normal', dcsi)
    v_sgs(ip) = v_sgs(ip) - f_s*v_sgs(ip)*lsm_dt/gamma_p &
        + 0.5_rprec*(ds2dt*v_sgs(ip)/sigma2(ip)+ds2dy)*lsm_dt + sqrt(f_s*2.0_rprec*sigma2(ip)/gamma_p)*dcsi

    call sample_distribution(mean, std, 'normal', dcsi)
    w_sgs(ip) = w_sgs(ip) - f_s*w_sgs(ip)*lsm_dt/gamma_p &
        + 0.5_rprec*(ds2dt*w_sgs(ip)/sigma2(ip)+ds2dz)*lsm_dt + sqrt(f_s*2.0_rprec*sigma2(ip)/gamma_p)*dcsi

    !write(*,*) 'compute_sgs,u_sgs,v_sgs,w_sgs: ',jt,u_sgs(ip),v_sgs(ip),w_sgs(ip)

    if (isnan(u_sgs(ip))) then
         print*, 'nan value u sgs'
         write(*,*) 'ip,fs,gamma_p,ds2dt,ds2dx,sigma2(ip),dcsi: ',p_ip(ip),f_s,gamma_p,ds2dt,ds2dx,sigma2(ip),dcsi
         write(*,*) 'u,v,w: ',u_p(ip),v_p(ip),w_p(ip)
         write(*,*) 'u_sgs,v_sgs,w_sgs: ', u_sgs(ip),v_sgs(ip),w_sgs(ip)
         write(*,*) 'coord,x,y,z: ',coord,x_p(ip),y_p(ip),z_p(ip)
         write(*,*) 'diss,restke,sgstke: ',diss_p,restke_p,sgstke_p
         write(*,*) 'cs_opt2:'
         do i=1,nx
            do j=1,ny
               do k=1,nz-1
                  write(*,*) coord,i,j,k,cs_opt2(i,j,k)
               enddo
            enddo
         enddo
         write(*,*) diss
         write(*,*) '============================================'
         write(*,*) restke
         write(*,*) '============================================'
         write(*,*) sgstke
         write(*,*) '============================================'
         stop
    endif

    if (isnan(v_sgs(ip))) then
         print*, 'nan value v sgs'
         write(*,*) 'fs,gamma_p,ds2dt,ds2dx,sigma2(ip),dcsi',f_s,gamma_p,ds2dt,ds2dx,sigma2(ip),dcsi
         stop
    endif

    if (isnan(w_sgs(ip))) then
         print*, 'nan value w sgs'
         write(*,*) 'fs,gamma_p,ds2dt,ds2dx,sigma2(ip),dcsi',f_s,gamma_p,ds2dt,ds2dx,sigma2(ip),dcsi
         stop
    endif

end do
!$omp end parallel do

end subroutine compute_sgs

!------------------------------------------------------------------------------!
!                             COMPUTE_ADVECTION                                !
!------------------------------------------------------------------------------!
subroutine compute_advection()

implicit none

integer                   :: ip
real(rprec)               :: val_u,val_v,val_w
real(rprec), dimension(3) :: coord_p
real(rprec) :: ppp_local
!------------------------------------------------------------------------------!

if (screen) print*, "In compute advection"

if (trim(p_type).eq.'passive') then
    do ip=1,np
        ! Get flow velocity at particle location
        coord_p = (/x_p(ip), y_p(ip), z_p(ip)/) ! Shouldn't z component be in eta units?
        call    interp_vel(coord_p,val_u,val_v,val_w,ip)

        ! Update velocity
        u_p(ip) = val_u + u_sgs(ip)
        v_p(ip) = val_v + v_sgs(ip)
        w_p(ip) = val_w + w_sgs(ip)

        ! Update position
        x_p(ip) = x_p(ip) + u_p(ip)*lsm_dt
        y_p(ip) = y_p(ip) + v_p(ip)*lsm_dt
        z_p(ip) = z_p(ip) + w_p(ip)*lsm_dt

        ! Check if the domain borders were crossed
        call check_boundaries(ip)
    end do
else if (trim(p_type).eq.'settling') then
    do ip = 1,np
        ! Get flow velocity at particle location
        coord_p = (/x_p(ip), y_p(ip), z_p(ip)/) ! Shouldn't z component be in eta units?
        call    interp_vel(coord_p,val_u,val_v,val_w,ip)

        ! Update velocity
        u_p(ip) = val_u + u_sgs(ip)
        v_p(ip) = val_v + v_sgs(ip)
        w_p(ip) = val_w + w_sgs(ip) - settling_vel

        ! Update position
        x_p(ip) = x_p(ip) + u_p(ip)*lsm_dt
        y_p(ip) = y_p(ip) + v_p(ip)*lsm_dt
        z_p(ip) = z_p(ip) + w_p(ip)*lsm_dt

    end do
else if (trim(p_type).eq.'inertial') then
    !$omp parallel do shared(x_p,y_p,z_p,ax_p,ay_p,az_p,T_p,d_p,mass_p,np,delta_T,delta_m,delta_r,ppp_local)
    do ip=1,np
        ! Update position
        x_p(ip) = x_p(ip) + u_p(ip)*lsm_dt + 0.5_rprec*ax_p(ip)*lsm_dt*lsm_dt
        y_p(ip) = y_p(ip) + v_p(ip)*lsm_dt + 0.5_rprec*ay_p(ip)*lsm_dt*lsm_dt
        z_p(ip) = z_p(ip) + w_p(ip)*lsm_dt + 0.5_rprec*az_p(ip)*lsm_dt*lsm_dt

        ! Update velocity
        u_p(ip) = u_p(ip) + ax_p(ip)*lsm_dt
        v_p(ip) = v_p(ip) + ay_p(ip)*lsm_dt
        w_p(ip) = w_p(ip) + az_p(ip)*lsm_dt

        ! Update temperature
        if(temperature) then
          T_p(ip) = T_p(ip) + lsm_dt*delta_T(ip)
        endif

        ! Update diameter
        if(sublimation) then
          ppp_local = mass_p(ip) / (p_rho * pi*((d_p(ip))**3.0_rprec)/6.0_rprec) !Number of particles in the parcel
          d_p(ip) = d_p(ip) + lsm_dt*delta_r(ip)*2.0_rprec
          mass_p(ip) = mass_p(ip) + lsm_dt*delta_m(ip)*ppp_local
        endif

        ! Check if the domain borders were crossed
        call check_boundaries(ip)
    enddo
    !$omp end parallel do
end if

end subroutine compute_advection

!------------------------------------------------------------------------------!
!                             CHECK_BOUNDARIES                                 !
!------------------------------------------------------------------------------!
subroutine check_boundaries(ip)

implicit none

integer, intent(in)         :: ip
!------------------------------------------------------------------------------!

if (outlet) then ! No periodic BC for particles
    if ((x_p(ip).gt.l_x).or.(x_p(ip).lt.0._rprec)) then
        p_exist(ip) = .false.
    end if
    if ((y_p(ip).gt.l_y).or.(y_p(ip).lt.0._rprec)) then
        p_exist(ip) = .false.
    end if
else ! Periodic BC for particles
!   if(colour .eq. 0 ) then
!    x_p(ip) = modulo(x_p(ip),0.900_rprec*l_x)
!   else
    x_p(ip) = modulo(x_p(ip),l_x)
!   endif
    y_p(ip) = modulo(y_p(ip),l_y)
end if

if (z_p(ip) .gt. L_Z) then
    p_exist(ip) = .false.
    av_conc = av_conc - mass_p(ip)
end if

if(trim(mode) .ne. 'saltation') then
  if (z_p(ip) .lt. get_z_local(2)) then
    if (llbc_special.eq.'IBM_DFA') then ! this shouldn't happen, get rid of the parcel in case
        p_exist(ip) = .false.
        av_conc = av_conc-mass_p(ip)
    else ! move it above the lsm_surface, but close enough to compute lsm_surface processes
        !z_p(ip) = 0.5_rprec*disp_length
        p_exist(ip) = .false.
        !av_conc = av_conc-mass_p(ip)
        call update_lsm_surface(ip,1._rprec)
        vd_flux = vd_flux-mass_p(ip)
        av_conc = av_conc-mass_p(ip)
    end if
  end if
endif

if (trim(mode) .eq. 'saltation') then
   if (z_p(ip) .lt. 0.0_rprec) then
      if (llbc_special.eq.'IBM_DFA') then ! this shouldn't happen, get rid of the parcel in case
         p_exist(ip) = .false.
         av_conc = av_conc - mass_p(ip)
      else ! move the parcel above the lsm_surface, but close enough to compute lsm_surface processes
         z_p(ip) = 0.5_rprec*disp_length
      end if
   end if
endif

!if (d_p(ip) .lt. real(1e-5,rprec)) then
   !p_exist(ip) = .false.
   !write(*,*) 'now passive particle id:',p_ip(ip),x_p(ip),y_p(ip),z_p(ip),d_p(ip),colour,coord
!endif

end subroutine check_boundaries

!------------------------------------------------------------------------------!
!                             CHECK_LSM_SURFACE                                !
!------------------------------------------------------------------------------!
subroutine check_lsm_surface()

implicit none

integer                     :: ip,np1
real(rprec)                 :: phi_p, cosang
real(rprec), dimension(3)   :: coord_p, norm_s, dir_p, x_diff
real(rprec), dimension(2,2,2) :: var_section
integer, dimension(2)       :: ix_close,iy_close,iz_close ! indexes for interpolation
logical                     :: at_surf
!------------------------------------------------------------------------------!

! Check mass and momentum conservation
ene_cons_ok=0
mom_cons_ok=0
ene_cons_no=0
mom_cons_no=0

! Save current number of parcels
np1 = np

!!!$omp parallel do private(coord_p,phi_p,dir_p,norm_s,cosang,nz)
do ip=1,np1 ! Loop on all parcels
   ! Compute parcel normal height above the surface
   coord_p = (/ x_p(ip), y_p(ip), get_eta_local(z_p(ip))-(real(eta_displ(coord)*dz,rprec)) /) ! Parcel location ([m], [m], local index)

   ! Parcel height above surface [m]
   if (llbc_special.ne.'IBM_DFA') then
     phi_p = z_p(ip)
   else
     ! get indexes of closest W nodes
     call get_indexes_3d(coord_p,1._rprec,ix_close,iy_close,iz_close,x_diff,at_surf,2,ip)
     ! Interpolate phi between the 8 closest grid nodes
     var_section = phi(ix_close,iy_close,iz_close)
     call interp_var_lsm(var_section, phi_p, x_diff(1), x_diff(2), x_diff(3))
   endif

   ! Check if parcel is close to the surface
   if (phi_p.lt.disp_length) then
       ! Components of parcel velocity versor
       dir_p(1) = u_p(ip)/sqrt(u_p(ip)**2+v_p(ip)**2+w_p(ip)**2)
       dir_p(2) = v_p(ip)/sqrt(u_p(ip)**2+v_p(ip)**2+w_p(ip)**2)
       dir_p(3) = w_p(ip)/sqrt(u_p(ip)**2+v_p(ip)**2+w_p(ip)**2)

       ! Get surface normal vector at point of impact
       call interp_2dh(coord_p, nis(1,:,:), norm_s(1))
       call interp_2dh(coord_p, nis(2,:,:), norm_s(2))
       call interp_2dh(coord_p, nis(3,:,:), norm_s(3))

       norm_s = norm_s/mag(norm_s)
       cosang = dot_product(dir_p,norm_s)

       ! Check if parcel is aproaching the surface
       if ((cosang.lt.0._rprec).and.(trim(mode).eq.'saltation')) then
          ! Move parcel closer to the surface
          x_p(ip) = x_p(ip) - (disp_length-phi_p)*norm_s(1)
          y_p(ip) = y_p(ip) - (disp_length-phi_p)*norm_s(2)
          z_p(ip) = z_p(ip) - (disp_length-phi_p)*norm_s(3)

          ! Compute rebound and splash
          call lsm_surface_processes(ip)

       else if ((cosang.lt.0._rprec).and.(trim(mode).ne.'saltation')) then
          ! Parcel gets deposited
          p_exist(ip) = .false.
          call update_lsm_surface(ip,1._rprec)

          ! Update deposited mass and mass aloft
          vd_flux = vd_flux-mass_p(ip)
          av_conc = av_conc-mass_p(ip)
       end if
   end if
end do
!!!$omp end parallel do

!if(np .gt. 0) then
! write(*,*) ene_cons_ok,mom_cons_ok,ene_cons_no,mom_cons_no
!endif

end subroutine check_lsm_surface

!------------------------------------------------------------------------------!
!                             LSM_SURFACE_PROCESSES                            !
!------------------------------------------------------------------------------!
subroutine lsm_surface_processes(ip)

implicit none

integer, intent(in)         :: ip
integer                     :: ip1,np1
real(rprec) :: rand_theta
real(rprec)                 :: prob_reb, rand, d_g,             &
                               v_ang, h_ang, slope,             &
                               xsis_p, ysis_p
real(rprec), dimension(2)   :: i_ang
real(rprec), dimension(3)   :: coord_p
real(rprec), parameter      :: corr_e = 0._rprec, corr_m = 0._rprec
real(rprec)                 :: cos_a, cos_b, cos_i,             &
                               i_vel, r_vel, e_vel,             &
                               i_ene, i_mom,                    &
                               epsilon_r, epsilon_f, mu_r, mu_f,&
                               temp_mass, dep_mass, mass_g,     &
                               mean, std, n_splash, ais_p,      &
                               n_splash1, n_splash2, n_impact,  &
                               av_mass, sd_mass, av_vel,        &
                               av_vel2, sd_vel2, av_massvel2,   &
                               sd_vel, av_massvel, av_d3, sd_d3, &
                               out_ene,out_mom,tmp_out_vel
!------------------------------------------------------------------------------!

! Total energy and momentum after impact
out_ene = 0.0_rprec
out_mom = 0.0_rprec

! IMPACT PROPERTIES
coord_p         = (/x_p(ip), y_p(ip), get_eta_local(z_p(ip))- (real(eta_displ(coord)*dz,rprec))/) ! Parcel location ([m], [m], local index)
i_vel           = sqrt(u_p(ip)**2.0_rprec+v_p(ip)**2.0_rprec+w_p(ip)**2.0_rprec) ! Parcel impact velocity [m/s]
n_impact        = mass_p(ip)/(p_rho*pi*(d_p(ip)**3)/6._rprec)                    ! Number of impacting particles
i_ene           = 0.5_rprec*mass_p(ip)*(i_vel**2.0_rprec)                        ! Impact parcel kinetic energy [J]
i_mom           = mass_p(ip)*i_vel                                               ! Impact parcel momentum [kg.m/s]
i_ang(1)        = -asin(w_p(ip)/i_vel)                                           ! Vertical angle of impact (w < 0)
i_ang(2)        = atan2(v_p(ip),u_p(ip))                                         ! Horizontal angle of impact
! Get slope at impact location
call            interp_2dh(coord_p,sis(1,:,:),xsis_p)
call            interp_2dh(coord_p,sis(2,:,:),ysis_p)
slope           = xsis_p*cos(i_ang(2))+ysis_p*sin(i_ang(2))

! REBOUND
prob_reb        = 0.9_rprec*(1._rprec-exp(-2._rprec*i_vel)) ! Probability of rebound from Anderson and Haff, 1991 (constant from Groot Zwaaftink et al 2011)
call            random_number(rand)

if (rand.lt.prob_reb) then ! Rebound occurs
    r_vel       = 0.5_rprec*i_vel ! Velocity after rebound (from Doorschot and Lehning, 2002)
    ! Compute vertical angle of rebound
    mean        = 45._rprec/180._rprec*pi
    call        sample_distribution(mean, 0._rprec, 'exponential', v_ang) ! From Kok and Renno, 2009
    v_ang       = min(0.5_rprec*pi, max(-0.5_rprec*pi, v_ang+slope))
    ! Horizontal angle of rebound (aligned with impact)
    h_ang       = i_ang(2)
    ! Update parcel velocity
    u_p(ip)     = r_vel * cos(v_ang)*cos(h_ang)
    v_p(ip)     = r_vel * cos(v_ang)*sin(h_ang)
    w_p(ip)     = r_vel * sin(v_ang)
    ! Update parcel position
    x_p(ip)     = coord_p(1)+u_p(ip)*lsm_dt
    y_p(ip)     = coord_p(2)+v_p(ip)*lsm_dt
    z_p(ip)     = z_p(ip) + w_p(ip)*lsm_dt
    ! Check if the domain borders were crossed
    call        check_boundaries(ip)
    ! Set sgs velocities after rebound
    u_sgs(ip)   = 0._rprec
    v_sgs(ip)   = 0._rprec
    w_sgs(ip)   = 0._rprec
    ! Set aceleration after rebound
    ax_p(ip)    = 0._rprec
    ay_p(ip)    = 0._rprec
    az_p(ip)    = 0._rprec

else ! No rebound
    ! Parcel gets deposited
    p_exist(ip) = .false.
    call        update_lsm_surface(ip,1._rprec)

    ! Update deposited mass and mass aloft
    vd_flux     = vd_flux-mass_p(ip)
    av_conc     = av_conc-mass_p(ip)
end if

! Compute outgoing energy and momentum (from rebound)
tmp_out_vel = sqrt(u_p(ip)**2.0_rprec+v_p(ip)**2.0_rprec+w_p(ip)**2.0_rprec)
out_ene     = out_ene + 0.5_rprec*mass_p(ip)*(tmp_out_vel**2.0_rprec)
out_mom     = out_mom + mass_p(ip)*tmp_out_vel

! NUMBER OF EJECTIONS
! Set splash parameters (from Comola and Lehning 2017)
mu_r            = 0.5_rprec                                 ! Fraction of impact momentum retained by the rebounding particle
epsilon_r       = 0.25_rprec                                ! Fraction of impact energy retained by the rebounding particle
mu_f            = 0.4_rprec                                 ! Fraction of impact momentum lost to the bed
epsilon_f       = 0.96_rprec*(1._rprec-prob_reb*epsilon_r)  ! Fraction of impact energy lost to the bed

av_d3           = (d_m+(d_s**2)/d_m)**3                           ! mean(d**3), d - diameter of ejected particle
sd_d3           = av_d3*sqrt((1._rprec+(d_s/d_m)**2)**9-1._rprec) ! std(d**3)
av_vel          = 0.25_rprec*(i_vel**0.3_rprec)                   ! mean(v), v - ejection velocity
av_vel2         = 2._rprec*(av_vel**2)                            ! mean(v**2)
av_mass         = p_rho*pi/6._rprec*av_d3                         ! mean(mass_e), mass_e - mass of ejected particle
sd_vel          = av_vel                                          ! std(v)
sd_vel2         = 2._rprec*sqrt(5._rprec)*(av_vel**2)             ! std(v**2)
sd_mass         = p_rho*pi/6._rprec*sd_d3                         ! std(mass_e), mass_e - mass of ejected particle

if (llbc_special.ne.'IBM_DFA') then
    cos_a       = 0.75_rprec    ! mean(cos(alpha)), alpha - vertical angle of ejection
    cos_b       = 0.96_rprec    ! mean(cos(beta)), beta - horizontal angle of ejection
    cos_i       = cos(i_ang(1))
else
    cos_a       = 1._rprec
    cos_b       = 1._rprec
    cos_i       = 1._rprec
end if

av_massvel      = av_mass*av_vel*cos_a*cos_b  + corr_m*sd_mass*sd_vel                             ! mean(mass_e * v)
av_massvel2     = av_mass*av_vel2 + corr_e*sd_mass*sd_vel2                                        ! mean(mass_e * v**2)
n_splash1       = i_ene*(1._rprec - prob_reb*epsilon_r - epsilon_f)/(0.5_rprec*av_massvel2+b_ene) ! Number of ejected grains (from energy conservation)
n_splash2       = i_mom*cos_i*(1._rprec - prob_reb*mu_r - mu_f)/av_massvel                        ! Number of ejected grains (from momentum conservation)
! Number of ejected grains after impact
n_splash        = min(n_splash1,n_splash2)

! SAMPLING NEW EJECTED PARCELS
if (n_splash.ge.n_impact) then
    np1                 = floor(n_splash/ppp_max)+1 ! Minimum number of parcels needed to eject n_splash
    ! Loop on ejected parcels
    do ip1=1,np1
        ! Compute particle diameter in the parcel from distribution
        call            sample_distribution(d_m,d_s,'lognormal',d_g)
        d_g             = min(d_max,max(d_g,d_min))
        ! Mass of one particle
        mass_g          = p_rho*pi*(d_g**3)/6._rprec
        ! Compute splashed mass
        if              (ip1.ne.np1) then
            temp_mass   = ppp_max*mass_g
        else
            temp_mass   = (n_splash-(np1-1)*ppp_max)*mass_g
        end if
        ! Compute deposited mass at impact location
        call            interp_2dh(coord_p,lsm_surface,dep_mass) ! Get deposited mass at impact location
        call            interp_2dh(coord_p,ais,ais_p)            ! Get grid surface area at impact location
        dep_mass        = dep_mass*ais_p
        ! Constrain temp_mass = min(temp_mass,dep_mass)
        if              ((dep_mass.gt.0._rprec).and.(dep_mass.lt.temp_mass)) then
            temp_mass   = dep_mass
        end if
        if              (temp_mass.le.dep_mass) then
            ! EJECT PARCEL
            np          = np + 1        ! Define new parcel index
            np_global   = np_global + 1 ! Update total number of parcels
            if          (np.gt.np_max) then
                print*, "INCREASE NUMBER OF PARCELS (C)"
                write(*,*) x_p(ip),y_p(ip),z_p(ip),u_p(ip),v_p(ip),w_p(ip),n_impact,n_splash,temp_mass,dep_mass,d_g,d_p(ip),i_vel,ais_p
                stop
            end if
            ! Compute parcel ejection velocity from distribution
            call        sample_distribution(av_vel,0._rprec,'exponential',e_vel)
            ! Save new parcel data
            d_p(np)     = d_g
            mass_p(np)  = temp_mass
            p_exist(np) = .true.
            p_ip(np)    = np_global
            ! Compute ejection vertical angle from distribution
            mean        = 50._rprec/180._rprec*pi
            call        sample_distribution(mean, 0._rprec, 'exponential', v_ang)
            v_ang       = min(0.5_rprec*pi, max(-0.5_rprec*pi, v_ang+slope))
            ! Compute ejection horizontal angle
            if          (llbc_special.ne.'IBM_DFA') then
                mean    = i_ang(2)
                std     = 15._rprec/180._rprec*pi
                call    sample_distribution(mean, std, 'normal', h_ang)
                h_ang   = min(pi, max(-pi, h_ang))
            else
                h_ang   = i_ang(2)
            end if
            ! Parcel initial velocity
            u_p(np)     = e_vel * cos(v_ang)*cos(h_ang)
            v_p(np)     = e_vel * cos(v_ang)*sin(h_ang)
            w_p(np)     = e_vel * sin(v_ang)
            ! Parcel initial position
            x_p(np)     = coord_p(1)+u_p(np)*lsm_dt
            y_p(np)     = coord_p(2)+v_p(np)*lsm_dt
            z_p(np)     = z_p(ip)+w_p(np)*lsm_dt
            ! Check if the domain borders were crossed
            call        check_boundaries(np)
            ! Save cummulative data
            vs_flux     = vs_flux+mass_p(np)
            av_conc     = av_conc+mass_p(np)
            ! Update deposited mass at the surface
            call        update_lsm_surface(np,-1._rprec)
            ! Parcel initial sgs velocities
            u_sgs(np)   = 0._rprec
            v_sgs(np)   = 0._rprec
            w_sgs(np)   = 0._rprec
            sigma2(np)  = 0._rprec ! (used in compute_sgs)
            ! Parcel initial aceleration
            ax_p(np)    = 0._rprec
            ay_p(np)    = 0._rprec
            az_p(np)    = 0._rprec
            IF (TEMPERATURE) then
               ! call random_number(rand_theta)
               if (LLBC_SC1 .eq. 'DRC') then
                  T_p(np)     = (DLBC_SC1*T_SCALE)/T_SCALE !+ 0.002_rprec*DLBC_SC1*(2.0*rand_theta-1.0_rprec) !0.9_rprec
               else
                  T_p(np) = 1.0_rprec !*T_SCALE
               endif
            ENDIF
            ! Update outgoing energy and momentum (add splashed particles)
            tmp_out_vel = sqrt(u_p(np)**2.0_rprec+v_p(np)**2.0_rprec+w_p(np)**2.0_rprec)
            out_ene     = out_ene + 0.5_rprec*mass_p(np)*(tmp_out_vel**2.0_rprec)
            out_mom     = out_mom + mass_p(np)*tmp_out_vel
        end if !(temp_mass.le.dep_mass)
    end do !ip1=1,np1
end if !(n_splash.ge.n_impact)

  !if(out_ene .gt. i_ene) then
  !  write(*,*) 'OH NO - more outgoing energy than incoming energy !'
  !  write(*,*) jt,p_ip(ip),out_ene,i_ene
  !  ene_cons_no=ene_cons_no+1
  !else
  !  ene_cons_ok=ene_cons_ok+1
  !endif

   !if(out_mom .gt. i_mom) then
   !  write(*,*) 'OH NO - more outgoing mom than incoming mom !'
   !  write(*,*) jt,p_ip(ip),out_mom,i_mom
   !  mom_cons_no=mom_cons_no+1
   !else
   !  mom_cons_ok=mom_cons_ok+1
   ! endif

end subroutine lsm_surface_processes

!------------------------------------------------------------------------------!
!                               REMOVE_PARCELS                                 !
!------------------------------------------------------------------------------!
subroutine remove_parcels()

implicit none

integer :: np_old
!------------------------------------------------------------------------------!

if (screen) print*, "In remove parcels"

! Save current number of parcels
np_old = np

! Calculate the new number of parcels
np = count(p_exist(1:np_old))

if(np < np_old) then

  ! pack the existing parcels in the first np positions
  mass_p(1:np) = pack(mass_p(1:np_old), p_exist(1:np_old))
  d_p(1:np)    = pack(d_p(1:np_old), p_exist(1:np_old))
  x_p(1:np)    = pack(x_p(1:np_old), p_exist(1:np_old))
  y_p(1:np)    = pack(y_p(1:np_old), p_exist(1:np_old))
  z_p(1:np)    = pack(z_p(1:np_old), p_exist(1:np_old))
  u_p(1:np)    = pack(u_p(1:np_old), p_exist(1:np_old))
  v_p(1:np)    = pack(v_p(1:np_old), p_exist(1:np_old))
  w_p(1:np)    = pack(w_p(1:np_old), p_exist(1:np_old))
  ax_p(1:np)   = pack(ax_p(1:np_old), p_exist(1:np_old))
  ay_p(1:np)   = pack(ay_p(1:np_old), p_exist(1:np_old))
  az_p(1:np)   = pack(az_p(1:np_old), p_exist(1:np_old))
  u_sgs(1:np)  = pack(u_sgs(1:np_old), p_exist(1:np_old))
  v_sgs(1:np)  = pack(v_sgs(1:np_old), p_exist(1:np_old))
  w_sgs(1:np)  = pack(w_sgs(1:np_old), p_exist(1:np_old))
  sigma2(1:np) = pack(sigma2(1:np_old), p_exist(1:np_old))
  IF(TEMPERATURE) then
    T_p(1:np) = pack(T_p(1:np_old),p_exist(1:np_old))
  endif
  IF(HUMIDITY) then
    Q_p(1:np) = pack(Q_p(1:np_old),p_exist(1:np_old))
  endif
  p_ip(1:np) = pack(p_ip(1:np_old),p_exist(1:np_old))

  re_p(1:np) = pack(re_p(1:np_old),p_exist(1:np_old))
  sh_p(1:np) = pack(sh_p(1:np_old),p_exist(1:np_old))
  nu_p(1:np) = pack(nu_p(1:np_old),p_exist(1:np_old))
  timescale_stokes(1:np) = pack(timescale_stokes(1:np_old),p_exist(1:np_old))
  delta_r(1:np) = pack(delta_r(1:np_old),p_exist(1:np_old))
  delta_m(1:np) = pack(delta_m(1:np_old),p_exist(1:np_old))
  delta_T(1:np) = pack(delta_T(1:np_old),p_exist(1:np_old))

  ! reset the remaining positions
  x_p(np+1:np_old)    = 0._rprec
  y_p(np+1:np_old)    = 0._rprec
  z_p(np+1:np_old)    = 0._rprec
  u_p(np+1:np_old)    = 0._rprec
  v_p(np+1:np_old)    = 0._rprec
  w_p(np+1:np_old)    = 0._rprec
  ax_p(np+1:np_old)   = 0._rprec
  ay_p(np+1:np_old)   = 0._rprec
  az_p(np+1:np_old)   = 0._rprec
  u_sgs(np+1:np_old)  = 0._rprec
  v_sgs(np+1:np_old)  = 0._rprec
  w_sgs(np+1:np_old)  = 0._rprec
  mass_p(np+1:np_old) = 0._rprec
  d_p(np+1:np_old)    = 0._rprec
  sigma2(np+1:np_old) = 0._rprec
  IF(TEMPERATURE) then
   T_p(np+1:np_old) = 0.0_rprec
  ENDIF
  IF(HUMIDITY) then
   Q_p(np+1:np_old) = 0.0_rprec
  ENDIF
  p_ip(np+1:np_old) = 0

  re_p(np+1:np_old)    = 0._rprec
  sh_p(np+1:np_old)    = 0._rprec
  nu_p(np+1:np_old)    = 0._rprec
  timescale_stokes(np+1:np_old) = 0._rprec
  delta_r(np+1:np_old) = 0._rprec
  delta_m(np+1:np_old) = 0._rprec
  delta_T(np+1:np_old) = 0._rprec

  ! now that parcels have been removed, the p_exist vector has only true entries
  p_exist(1:np)        = .true.
  p_exist(np+1:np_old) = .false.

endif

end subroutine remove_parcels

!------------------------------------------------------------------------------!
!                               COMPUTE_FORCES                                 !
!------------------------------------------------------------------------------!
subroutine compute_forces()

implicit none

integer                   :: ip !,iz
real(rprec)               :: Fx_p,Fy_p,Fz_p,Ft_p_1,Ft_p_2,FQ_p,n_g,mass_g,val_u,val_v,val_w,mag_vel_rel
real(rprec), dimension(3) :: coord_p,vel_rel,x_diff_uvp,x_diff_w
real(rprec), dimension(2,2,2) :: var_section
real(rprec) :: loc_T,loc_Q
real(rprec) :: Lv,sat_vap_p_1,sat_vap_p_2,q_conv,q_evap
real(rprec) :: temp_T
real(rprec) :: g_hat
real(rprec) :: drag_c
real(rprec),dimension(3) :: drag_f
real(rprec) :: diffusion_coeff, tm_coeff, conductivity
real(rprec) :: sat_spec_humidity_1,sat_spec_humidity_2
real(rprec) :: thickness, eta_loc, mass_air_norm
integer :: min_eta,max_eta
integer, dimension(2) :: ix_uvp,iy_uvp,iz_uvp, ix_w,iy_w,iz_w ! indexes for extrapolation
logical               :: at_surf_uvp, at_surf_w
!------------------------------------------------------------------------------!

g_hat = g*(z_i/u_star**2.0_rprec)

!$omp parallel do private(coord_p,val_u,val_v,val_w,loc_T,loc_Q,vel_rel,mag_vel_rel,Lv,sat_vap_p_1,sat_vap_p_2), &
!$omp& private(drag_c,drag_f,mass_g,n_g,Fx_p,Fy_p,Fz_p,FQ_p,q_conv,q_evap,Ft_p_1,Ft_p_2,eta_loc,thickness,min_eta,max_eta), &
!$omp& private(var_section,ix_uvp,iy_uvp,iz_uvp,ix_w,iy_w,iz_w,x_diff_uvp,x_diff_w,at_surf_uvp,at_surf_w,mass_air_norm) &
!$omp& reduction(+:sum_delta_m,sum_delta_q)
do ip=1,np ! Loop on parcels
    ! Parcel location ([m], [m], local index)
    coord_p = (/ x_p(ip),y_p(ip),get_eta_local(z_p(ip))-(real(eta_displ(coord)*dz,rprec)) /)

    ! Get LES variables to particle location (flow velocity, temperature and humidity)
    call interp_vel(coord_p,val_u,val_v,val_w,ip)
    IF(TEMPERATURE) then
      !call interp_var_lsm(coord_p,T_regroup,0.5_rprec,loc_T)
       call interp_T(coord_p,loc_T,ip)
    ENDIF
    IF(HUMIDITY) then
      !call interp_var_lsm(coord_p,Q_regroup,0.5_rprec,loc_Q)
       call interp_Q(coord_p,loc_Q,ip) ! specific humidity (kg kg^{-1})
    ENDIF

    ! Compute relative velocity between fluid and parcels
    vel_rel(1) = val_u + u_sgs(ip) - u_p(ip)
    vel_rel(2) = val_v + v_sgs(ip) - v_p(ip)
    vel_rel(3) = val_w + w_sgs(ip) - w_p(ip)
    mag_vel_rel = mag(vel_rel)

    !!!!!
    !!! A SIMPLE SANITY CHECK !!!!!
    if(mag_vel_rel .gt. 100.0_rprec) then
       write(*,*) ' A particle has gone insane: ',x_p(ip),y_p(ip),z_p(ip),get_eta_local(z_p(ip)),real(eta_displ(coord)*dz,rprec),u_p(ip),v_p(ip),w_p(ip),u_sgs(ip),v_sgs(ip),w_sgs(ip),val_u,val_v,val_w,mag_vel_rel
       STOP
       !p_exist(ip) = .false.
       !vel_rel(1) = 0.0_rprec
       !vel_rel(2) = 0.0_rprec
       !vel_rel(3) = 0.0_rprec
       !ax_p(ip)=0.0_rprec
       !ay_p(ip)=0.0_rprec
       !az_p(ip)=0.0_rprec
       !delta_m(ip)=0.0_rprec
       !delta_r(ip)=0.0_rprec
       !delta_T(ip)=0.0_rprec
       !cycle
       !continue
    endif
    !!!!!!!

    ! computing the non-dimensional numbers for momentum, mass and heat transfer
    re_p(ip) = (u_star)*mag_vel_rel*d_p(ip) / nu_molec ! Reynolds number
    sh_p(ip) = 1.79_rprec + 0.606_rprec*(re_p(ip)**(0.5_rprec))*((PR_SC2)**(1.0_rprec/3.0_rprec)) ! Sherwood number
    nu_p(ip) = 1.79_rprec + 0.606_rprec*(re_p(ip)**(0.5_rprec))*((PR_SC1)**(1.0_rprec/3.0_rprec)) ! Nusselt number
    ! timescale for particle
    timescale_stokes(ip) = p_rho*(d_p(ip)**2.0_rprec)/(18.0_rprec*f_rho*nu_molec) ! the stokes timescale [s]

    if(sublimation) then
    ! latent heat of vaporization / condensation
     !Lv = (25.0_rprec - 0.02274_rprec*(loc_T*288.0_rprec-273.15_rprec))*((10.0_rprec)**(5.0_rprec)) ! units J Kg^{-1} ! for water drop
     Lv = (2834.1_rprec - 0.29_rprec*(loc_T*T_SCALE-273.15_rprec) - 0.004_rprec*(loc_T*T_SCALE-273.15_rprec)**(2.0_rprec))*((10.0_rprec)**(3.0_rprec))
     ! Lv = (2500.8_rprec-2.36_rprec*(loc_T*T_SCALE-273.15_rprec) + 0.0016_rprec*((loc_T*T_SCALE-273.15_rprec)**(2.0_rprec)) + &
     !              0.00006_rprec*((loc_T*T_SCALE-273.15_rprec)**(3.0_rprec)))*((10.0_rprec)**(3.0_rprec))
    !!! COMPUTE ALL THE INCREMENTS FOR THE PARTICLES

    ! (A): MASS TRANSFER
    ! (A.1): Saturation vapor pressure
    !sat_vap_p = 610.94_rprec * exp( (17.625_rprec*(loc_T*288.0_rprec-273.15_rprec)) / (243.04_rprec+loc_T*288.0_rprec-273.15_rprec) )
    !sat_vap_p = 611.21_rprec * exp( (18.678_rprec - (loc_T*T_SCALE - 273.15_rprec)/234.5_rprec) * ((loc_T*T_SCALE-273.15_rprec)/(257.14_rprec + loc_T*T_SCALE-273.15)))
    sat_vap_p_1 = 611.15_rprec * exp( (23.036_rprec - (T_p(ip)*T_SCALE - 273.15_rprec)/333.7_rprec) * ((T_p(ip)*T_SCALE-273.15_rprec)/(279.82_rprec + T_p(ip)*T_SCALE-273.15)))
    sat_vap_p_2 = 611.15_rprec * exp( (23.036_rprec - (loc_T*T_SCALE - 273.15_rprec)/333.7_rprec) * ((loc_T*T_SCALE-273.15_rprec)/(279.82_rprec + loc_T*T_SCALE-273.15)))
    ! saturation vapor density (kg m^{-3}) at drop (although named spec_humidity)
    sat_spec_humidity_1 = (molec_weight_water * sat_vap_p_1)/(R_universal * T_p(ip) * T_SCALE)

    sat_spec_humidity_2 = (molec_weight_water * sat_vap_p_2)/(R_universal * loc_T * T_SCALE)

    ! (A.2): saturation vapor density (kg m^{-3}) at drop (although named spec_humidity) ! USING ARDEN BUCK EQUATION / BUCK (1996)
    Q_p(ip) = sat_spec_humidity_1 !(molec_weight_water/(R_universal*(T_p(ip)*T_SCALE)*f_rho))*sat_vap_p* &
                                  !exp((Lv*molec_weight_water/R_universal)*(1.0_rprec/(loc_T*T_SCALE) - 1.0_rprec/(T_p(ip)*T_SCALE))) ! + &
                                  ! 2.0_rprec*molec_weight_water*surface_tension/(R_universal*water_rho*(d_p(ip)*0.5_rprec)*(T_p(ip)*288.0_rprec)) )
    ! (A.3): compute change in radius and mass of the drop (Helgans and Richter, 2016)
    delta_r(ip) = (1.0_rprec/9.0_rprec)*(sh_p(ip)/PR_SC2)*((d_p(ip)*0.5_rprec)/timescale_stokes(ip))*(loc_Q-(Q_p(ip)/f_rho))*(z_i/u_star)! [m] if considering normalization, otherwise [m s-1]
    delta_m(ip) = pi*(d_p(ip)**(2.0_rprec))*p_rho*delta_r(ip) ! [kg] if considering normalization, otherwise [kg s-1]

    ! (B): HEAT TRANSFER
     delta_T(ip) =((-1.0_rprec/3.0_rprec)*(nu_p(ip)/PR_SC1)*(c_p_a/c_L)*(1.0_rprec/timescale_stokes(ip))*(T_p(ip)-loc_T)*T_SCALE* &
                  (z_i/u_star) + (3.0_rprec*Lv*(2.0_rprec/(d_p(ip)*c_L))*delta_r(ip)))*(1.0_rprec/T_SCALE)

    if(thorpe_mason) then
     diffusion_coeff = nu_molec/PR_SC2 ! [m^2 s^(-1)]
     conductivity = nu_molec * c_p_a * f_rho / PR_SC1
     tm_coeff = (Lv*molec_weight_water/(R_universal*loc_T*T_SCALE)) - 1.0_rprec
     delta_m(ip) = (pi*d_p(ip)*(loc_Q*f_rho/sat_spec_humidity_2 - 1.0_rprec))/( (Lv/(conductivity*loc_T*T_SCALE)) * (1.0_rprec/nu_p(ip)) * tm_coeff + &
                       (1.0_rprec/(diffusion_coeff*sat_spec_humidity_2*sh_p(ip))) ) * (z_i/u_star) ! [kg]
     delta_T(ip) = 0.0_rprec

     delta_r(ip) = delta_m(ip) / (p_rho * pi * (d_p(ip)**2.0_rprec))

    endif
   endif

   ! (C): MOMENTUM TRANSFER
   IF (Schiller) then ! Expression for Cd from Schiller and Nauman (in Clift et al, 1978 - page 111, Table 5.1)
      ! Parcel acceleration
      ax_p(ip) = ((1.0_rprec + 0.15_rprec*(re_p(ip)**(0.687_rprec)))*(1.0_rprec/timescale_stokes(ip))*(vel_rel(1))) * (z_i/u_star) ! dimension-less (D-L)
      ay_p(ip) = ((1.0_rprec + 0.15_rprec*(re_p(ip)**(0.687_rprec)))*(1.0_rprec/timescale_stokes(ip))*(vel_rel(2))) * (z_i/u_star) ! D-L
      az_p(ip) = ((1.0_rprec + 0.15_rprec*(re_p(ip)**(0.687_rprec)))*(1.0_rprec/timescale_stokes(ip))*(vel_rel(3))) * (z_i/u_star) - g_hat ! D-L
   endif

   IF (Mason) then
      if (re_p(ip) .ne. 0._rprec) then
         drag_c    = 24._rprec/re_p(ip) + 6._rprec/(1._rprec+sqrt(re_p(ip))) + 0.4_rprec ! Sphere drag coefficient, Cd
         drag_f(1) = (0.5_rprec*f_rho*pi*(d_p(ip)*d_p(ip)/4._rprec)*drag_c*vel_rel(1)*mag_vel_rel)*(u_star**2.0_rprec) ![N]
         drag_f(2) = (0.5_rprec*f_rho*pi*(d_p(ip)*d_p(ip)/4._rprec)*drag_c*vel_rel(2)*mag_vel_rel)*(u_star**2.0_rprec) ![N]
         drag_f(3) = (0.5_rprec*f_rho*pi*(d_p(ip)*d_p(ip)/4._rprec)*drag_c*vel_rel(3)*mag_vel_rel)*(u_star**2.0_rprec) ![N]
      else
         drag_f(1:3) = 0._rprec
      end if
      ! Parcel acceleration
      mass_g   = p_rho*pi*(d_p(ip)**3)/6._rprec ![kg]
      ax_p(ip) = (drag_f(1)/mass_g)*(z_i/(u_star**2.0_rprec))         ! D-L
      ay_p(ip) = (drag_f(2)/mass_g)*(z_i/(u_star**2.0_rprec))         ! D-L
      az_p(ip) = (drag_f(3)/mass_g - g_hat)*(z_i/(u_star**2.0_rprec)) ! D-L
   endif

    ! COMPUTE ALL THE COUPLING TERMS
    mass_g   = p_rho*pi*(d_p(ip)**3.0_rprec)/6._rprec ! Mass of one paerticle
    n_g      = mass_p(ip)/mass_g                      ! Number of particles on parcel

    eta_loc = get_eta_local(z_p(ip))                            ! Global index of particle position
    min_eta = floor(eta_loc)                                    ! Closest W-node below
    max_eta = min_eta+1                                         ! Closest W-node above
    thickness = get_z_local(2*max_eta) - get_z_local(2*min_eta) ! dz between the two neighboring W-nodes
    mass_air_norm = f_rho*(dx*dy*thickness)*(z_i**3.0_rprec)    ! mass of air (normalized)

    ! (A) MOMENTUM (Newton's second law)
    Fx_p = -1.0_rprec * n_g * (1.0_rprec/mass_air_norm) * (mass_g*ax_p(ip) + delta_m(ip)*u_p(ip)) ! m s-2
    Fy_p = -1.0_rprec * n_g * (1.0_rprec/mass_air_norm) * (mass_g*ay_p(ip) + delta_m(ip)*v_p(ip)) ! m s-2
    Fz_p = -1.0_rprec * n_g * (1.0_rprec/mass_air_norm) * (mass_g*(az_p(ip)+g_hat) + delta_m(ip)*w_p(ip)) ! m s-2

    ! (B) VAPOR MASS
    IF (sublimation) then
       FQ_p = -1.0_rprec * n_g * (1.0_rprec/mass_air_norm) * (delta_m(ip)) ! kg kg-1 s-1

       ! (C) HEAT
       q_conv = (-1.0_rprec/3.0_rprec)*(nu_p(ip)/PR_SC1)*(c_p_a)*(p_rho)*((pi*d_p(ip)**3.0_rprec)/(6.0_rprec*timescale_stokes(ip)))*(T_p(ip)-loc_T)*(z_i/u_star) ! J s-1
       q_evap = delta_m(ip)*c_p_v*(loc_T-T_p(ip)) ! J s-1
       if (thorpe_mason) then
          q_conv = -1.0_rprec * Lv * delta_m(ip) * (1.0_rprec/T_SCALE)
          q_evap = 0.0_rprec
       endif

       Ft_p_1 = -1.0_rprec * n_g * (1.0_rprec/(mass_air_norm*c_p_a)) * (q_conv) ! K s-1
       Ft_p_2 = -1.0_rprec * n_g * (1.0_rprec/(mass_air_norm*c_p_a)) * (q_evap) ! K s-1
    endif

    if (d_p(ip) .le. real(1e-5,rprec))  then
       ! Assume parcel to sublimate completely and immediately 
       !write(*,*) 'Particle: ',p_ip(ip),'d_p: ',d_p(ip),' sublimates completely and disappears!'
       p_exist(ip) = .false.
       Fx_p=0.0_rprec
       Fy_p=0.0_rprec
       Fz_p=0.0_rprec
       if(sublimation) then
          ! complete sublimation: consider the whole mass of the particle
          delta_m(ip) = mass_g
          FQ_p        = -1.0_rprec * n_g * (1.0_rprec/mass_air_norm) * mass_g
          ! energy for complete sublimation is taken from the atmosphere (same for thorpe_mason)
          q_conv = -1.0_rprec * Lv * mass_g * (1.0_rprec/T_SCALE)
          Ft_p_1 = -1.0_rprec * n_g * (1.0_rprec/(mass_air_norm*c_p_a)) * (q_conv) ! K s-1
          Ft_p_2 = 0.0_rprec
       endif
    endif

    if(sublimation) then
      sum_delta_m = sum_delta_m+n_g*delta_m(ip)
      sum_delta_q = sum_delta_q+n_g*q_conv
    endif

    ! Get indexes of closest W nodes
    call get_indexes_3d(coord_p,1._rprec,ix_w,iy_w,iz_w,x_diff_w,at_surf_w,17,ip)
    ! Get indexes of the closest UVP nodes
    call get_indexes_3d(coord_p,0.5_rprec,ix_uvp,iy_uvp,iz_uvp,x_diff_uvp,at_surf_uvp,18,ip)
    
    ! EXTRAPOLATE THE COUPLING TERMS TO THE LES MATRICES (to the nearst 8 nodes)
    ! Fx_p
    var_section                  = Fx_lsm(ix_uvp,iy_uvp,iz_uvp)
    call extrap_3d(Fx_p, var_section, x_diff_uvp(1),x_diff_uvp(2),x_diff_uvp(3), at_surf_uvp)
    Fx_lsm(ix_uvp,iy_uvp,iz_uvp) = var_section
    ! Fy_p
    var_section                  = Fy_lsm(ix_uvp,iy_uvp,iz_uvp)
    call extrap_3d(Fy_p, var_section, x_diff_uvp(1),x_diff_uvp(2),x_diff_uvp(3), at_surf_uvp)
    Fy_lsm(ix_uvp,iy_uvp,iz_uvp) = var_section
    ! Fz_p
    var_section            = Fz_lsm(ix_w,iy_w,iz_w)
    call extrap_3d(Fz_p, var_section, x_diff_w(1),x_diff_w(2),x_diff_w(3), at_surf_w)
    Fz_lsm(ix_w,iy_w,iz_w) = var_section
    if (sublimation) then
      IF (TEMPERATURE) then
          ! Ft_p_1
          var_section                    = Ft_lsm_1(ix_uvp,iy_uvp,iz_uvp)
          call extrap_3d(Ft_p_1, var_section, x_diff_uvp(1),x_diff_uvp(2),x_diff_uvp(3), at_surf_uvp)
          Ft_lsm_1(ix_uvp,iy_uvp,iz_uvp) = var_section
          ! Ft_p_2
          var_section                    = Ft_lsm_2(ix_uvp,iy_uvp,iz_uvp)
          call extrap_3d(Ft_p_2, var_section, x_diff_uvp(1),x_diff_uvp(2),x_diff_uvp(3), at_surf_uvp)
          Ft_lsm_2(ix_uvp,iy_uvp,iz_uvp) = var_section
      ENDIF
      IF (HUMIDITY) then
         ! FQ_p
         var_section                  = Fq_lsm(ix_uvp,iy_uvp,iz_uvp)
         call extrap_3d(FQ_p, var_section, x_diff_uvp(1),x_diff_uvp(2),x_diff_uvp(3), at_surf_uvp)
         Fq_lsm(ix_uvp,iy_uvp,iz_uvp) = var_section
      ENDIF
    endif

end do !ip=1,np
!$omp end parallel do

If (sublimation) then
   Ft_lsm = Ft_lsm_1 + Ft_lsm_2
endif

! compute the particle stress
!do iz=1,nz
!    if ((coord .eq. 0) .and. (iz.eq.1)) then
!        tpxz(:,:) = Fx_lsm(:,:,iz)*1.5_rprec*dz
!        tpyz(:,:) = Fy_lsm(:,:,iz)*1.5_rprec*dz
!    else
!DMB: Particle stresses close to the surface is already taken into account in source terms Fx,Fy,Fz
!        tpxz(:,:) = Fx_lsm(:,:,1)*get_z_local(1)
!        tpyz(:,:) = Fy_lsm(:,:,1)*get_z_local(1)
!    end if
!end do

!call mpi_allreduce(mpi_in_place,tpxz,ld*ny,mpi_rprec,mpi_sum,comm,ierr)
!call mpi_allreduce(mpi_in_place,tpyz,ld*ny,mpi_rprec,mpi_sum,comm,ierr)
!call mpi_barrier(comm,ierr)

end subroutine compute_forces

!------------------------------------------------------------------------------!
!                            COMPUTE_CONCENTRATION                             !
!------------------------------------------------------------------------------!
subroutine compute_concentration()

implicit none

integer                   :: ip
real(rprec)               :: conc_p
real(rprec), dimension(3) :: coord_p, x_diff
real(rprec), dimension(2,2,2) :: var_section
real(rprec) :: eta_loc,thickness
integer :: min_eta,max_eta
integer, dimension(2) :: ix_close,iy_close,iz_close ! indexes for extrapolation
logical               :: at_surf
!------------------------------------------------------------------------------!

if (screen) print*, "In compute concentration"

! Loop on parcels
do ip=1,np
    ! Parcel location ([m], [m], local index)
    coord_p(1) = x_p(ip)
    coord_p(2) = y_p(ip)
    coord_p(3) = get_eta_local(z_p(ip)) - (real(eta_displ(coord)*dz,rprec))

    eta_loc   = get_eta_local(z_p(ip))                          ! Global index of parcel location
    min_eta   = floor(eta_loc)                                  ! Closest W-node below
    max_eta   = min_eta+1                                       ! Closest W-node above
    thickness = get_z_local(2*max_eta) - get_z_local(2*min_eta) ! dz between the two neighboring W-nodes

    conc_p  = mass_p(ip)/(dx*dy*thickness*(z_i**3.0_rprec))     ! Parcel contribution to concentration [kg/m^3]

    ! Get indexes of closest W nodes
    call get_indexes_3d(coord_p,1._rprec,ix_close,iy_close,iz_close,x_diff,at_surf,19,ip)
    
    ! Extrapolate concentration to the nearst 8 W-nodes
    var_section                      = conc(ix_close,iy_close,iz_close)
    call extrap_3d_conc(conc_p, var_section, x_diff(1),x_diff(2),x_diff(3))
    conc(ix_close,iy_close,iz_close) = var_section
end do

end subroutine compute_concentration

!------------------------------------------------------------------------------!
!                              UPDATE_LSM_SURFACE                              !
!------------------------------------------------------------------------------!
subroutine update_lsm_surface(ip, norm)

implicit none

integer, intent(in)       :: ip
real(rprec), intent(in)   :: norm
integer                   :: ix,ix1,iy,iy1,nn,jx,jy
integer,dimension(4,2)    :: cn
real(rprec)               :: mass,x1,x2
real(rprec), dimension(3) :: coord_p
real(rprec), dimension(4) :: wg
!------------------------------------------------------------------------------!

! Update 'birth'/'death' files with ejected/deposited parcels
if (mod(p_ip(ip),1).eq.0) then
   if (norm .lt. 0.0_rprec) then ! norm = -1
      open(200,file=trim(lsm_path)//'/birth.log',form='formatted',action='write',position='append')
      write(200,*) p_ip(ip),d_p(ip),jt_total
      close(200)
   else ! norm = 1
      open(200,file=trim(lsm_path)//'/death.log',form='formatted',action='write',position='append')
      write(200,*) p_ip(ip),d_p(ip),jt_total
      close(200)
   endif
endif

! Parcel location coordinates [m]
coord_p = (/x_p(ip),y_p(ip),z_p(ip)/)

! Parcel mass (positive/negative if it is being added/ejected from the surface)
mass = norm*mass_p(ip)

! Indexes of the closest surface nodes
ix  = modulo(floor(coord_p(1)/dx),nx)+1
iy  = modulo(floor(coord_p(2)/dy),ny)+1
ix1 = modulo(ix,nx)+1
iy1 = modulo(iy,ny)+1

! Coordinates of the 4 closest surface nodes
cn(1,:) = (/ix,iy/)
cn(2,:) = (/ix1,iy/)
cn(3,:) = (/ix,iy1/)
cn(4,:) = (/ix1,iy1/)

! Fraction of the distance from parcel location to the bottom-left neighboring node
x1 = modulo(coord_p(1),dx)/dx
x2 = modulo(coord_p(2),dy)/dy

! Compute extrapolation weights
wg(1)   = (1.-x1)*(1.-x2)
wg(2)   = x1     *(1.-x2)
wg(3)   = (1.-x1)*x2
wg(4)   = x1     *x2

! Update lsm_surface at the closest surface nodes
do nn=1,4
    jx = cn(nn,1)
    jy = cn(nn,2)
    lsm_surface(jx,jy) = lsm_surface(jx,jy) + wg(nn)*mass/ais(jx,jy)
    if (lsm_surface(jx,jy).lt.0._rprec) then
        lsm_surface(jx,jy) = 0._rprec
    end if
end do

end subroutine update_lsm_surface

!------------------------------------------------------------------------------!
!                                 COMPUTE_FLUX                                 !
!------------------------------------------------------------------------------!
subroutine compute_flux()

implicit none

integer                             :: ip
real(rprec)                         :: n_g, m_g
real(rprec), dimension(3)           :: coord_p, x_diff
real(rprec), dimension(2,2,2)       :: var_section
integer, dimension(2)               :: ix_close,iy_close,iz_close ! indexes for extrapolation
logical                             :: at_surf
!------------------------------------------------------------------------------!

if (screen) print*, "In compute flux"

h_flux(:,:,:)   = 0._rprec
v_flux(:,:,:)   = 0._rprec
ppn(:,:,:)      = 1.E-24
c_h_flux(:,:,:) = 0.0_rprec
c_v_flux(:,:,:) = 0.0_rprec

! Loop on parcels
do ip=1,np
   coord_p = (/x_p(ip),y_p(ip),get_eta_local(z_p(ip))-(real(eta_displ(coord)*dz,rprec)) /) ! Parcel location ([m], [m], local index)
   m_g     = (pi/6._rprec)*p_rho*(d_p(ip)**3)                                          ! Mass of one particle
   n_g     = mass_p(ip)/m_g                                                            ! Number of particles in the parcels

   ! Get indexes of closest W nodes
   call get_indexes_3d(coord_p,1._rprec,ix_close,iy_close,iz_close,x_diff,at_surf,20,ip)

   ! Extrapolate parcel velocities and number of particles to the nearest 8 W-nodes
   ! number of particles times u_p
   var_section                        = h_flux(ix_close,iy_close,iz_close)
   call extrap_3d(n_g*u_p(ip), var_section, x_diff(1),x_diff(2),x_diff(3), at_surf)
   h_flux(ix_close,iy_close,iz_close) = var_section
   ! number of particles times w_p
   var_section                        = v_flux(ix_close,iy_close,iz_close)
   call extrap_3d(n_g*w_p(ip), var_section, x_diff(1),x_diff(2),x_diff(3), at_surf)
   v_flux(ix_close,iy_close,iz_close) = var_section
   ! n_g
   var_section                        = ppn(ix_close,iy_close,iz_close)
   call extrap_3d(n_g, var_section, x_diff(1),x_diff(2),x_diff(3), at_surf)
   ppn(ix_close,iy_close,iz_close)    = var_section
end do

! x and z parcel velocities weighted by number of particles at each W-node
h_flux = h_flux/ppn
v_flux = v_flux/ppn

! Horizontal/vertical flux at each W-node
c_h_flux = conc*h_flux
c_v_flux = conc*v_flux

end subroutine compute_flux

!------------------------------------------------------------------------------!
!                                 INTERP_GRAD                                  !
!------------------------------------------------------------------------------!
subroutine interp_grad(coord_p,stag,var,ddx,ddy,ddz)

implicit none

real(rprec), intent(in)                        :: stag
real(rprec), dimension(3), intent(in)          :: coord_p
real(rprec), dimension(ld,ny,0:nz), intent(in) :: var
real(rprec), intent(out)                       :: ddx, ddy, ddz
integer                                        :: ix,ix1,iy,iy1,iz,iz1,nn,jx,jy,jz
integer, dimension(8,3)                        :: cn
real(rprec)                                    :: x1,x2,x3,x_i,y_i,z_i,phi_i, &
                                                  var_ix,var_ix1,var_iy,var_iy1,var_iz,var_iz1
real(rprec), dimension(8)                      :: wg

ix  = modulo(floor(coord_p(1)/dx),nx)+1
iy  = modulo(floor(coord_p(2)/dy),ny)+1
ix1 = modulo(ix,nx)+1
iy1 = modulo(iy,ny)+1
iz  = floor(coord_p(3)/dz+stag)
iz1 = iz+1

cn(1,:) = (/ix,iy,iz/)
cn(2,:) = (/ix1,iy,iz/)
cn(3,:) = (/ix,iy1,iz/)
cn(4,:) = (/ix1,iy1,iz/)
cn(5,:) = (/ix,iy,iz1/)
cn(6,:) = (/ix1,iy,iz1/)
cn(7,:) = (/ix,iy1,iz1/)
cn(8,:) = (/ix1,iy1,iz1/)

! calculate weights
x1 = modulo(coord_p(1),dx)/dx
x2 = modulo(coord_p(2),dy)/dy
x3 = modulo(coord_p(3)-stag*dz,dz)/dz
wg(1)   = (1.-x1)*(1.-x2)*(1.-x3)
wg(2)   = x1     *(1.-x2)*(1.-x3)
wg(3)   = (1.-x1)*x2     *(1.-x3)
wg(4)   = x1     *x2     *(1.-x3)
wg(5)   = (1.-x1)*(1.-x2)*x3
wg(6)   = x1     *(1.-x2)*x3
wg(7)   = (1.-x1)*x2     *x3
wg(8)   = x1     *x2     *x3

var_ix  = var(ix,iy,iz)*wg(1)  + var(ix,iy1,iz)*wg(3)  + var(ix,iy,iz1)*wg(5)  + var(ix,iy1,iz1)*wg(7)
var_ix1 = var(ix1,iy,iz)*wg(2) + var(ix1,iy1,iz)*wg(4) + var(ix1,iy,iz1)*wg(6) + var(ix1,iy1,iz1)*wg(8)
ddx = (var_ix1-var_ix)/dx

var_iy  = var(ix,iy,iz)*wg(1)  + var(ix1,iy,iz)*wg(2)  + var(ix1,iy,iz1)*wg(6)  + var(ix,iy,iz1)*wg(5)
var_iy1 = var(ix,iy1,iz)*wg(3) + var(ix1,iy1,iz)*wg(4) + var(ix1,iy1,iz1)*wg(8) + var(ix,iy1,iz1)*wg(7)
ddy = (var_iy1-var_iy)/dy

var_iz  = var(ix,iy,iz)*wg(1)  + var(ix1,iy,iz)*wg(2)  + var(ix1,iy1,iz)*wg(4)  + var(ix,iy1,iz)*wg(3)
var_iz1 = var(ix,iy,iz1)*wg(5) + var(ix1,iy,iz1)*wg(6) + var(ix1,iy1,iz1)*wg(8) + var(ix,iy1,iz1)*wg(7)
ddz = (var_iz1-var_iz)/dz

end subroutine interp_grad

!------------------------------------------------------------------------------!
!                                EXTRAP_3D_CONC                                !
!------------------------------------------------------------------------------!
subroutine extrap_3d_conc(var_p, var, x1, x2, x3)

implicit none

real(rprec), intent(in)                     :: var_p
!real(rprec) :: tmp_p
real(rprec), dimension(2,2,2), intent(inout)  :: var
real(rprec), intent(in)                       :: x1,x2,x3
real(rprec), dimension(2,2,2)               :: wg
!------------------------------------------------------------------------------!

! Compute extrapolation weights
wg(1,1,1)   = (1.-x1)*(1.-x2)*(1.-x3)
wg(2,1,1)   = x1     *(1.-x2)*(1.-x3)
wg(1,2,1)   = (1.-x1)*x2     *(1.-x3)
wg(2,2,1)   = x1     *x2     *(1.-x3)
wg(1,1,2)   = (1.-x1)*(1.-x2)*x3
wg(2,1,2)   = x1     *(1.-x2)*x3
wg(1,2,2)   = (1.-x1)*x2     *x3
wg(2,2,2)   = x1     *x2     *x3

! Update VAR with the contribution of VAR_P
var = var + wg*var_p

end subroutine extrap_3d_conc

!------------------------------------------------------------------------------!
!                                  EXTRAP_3D                                   !
!------------------------------------------------------------------------------!
subroutine extrap_3d(var_p, var, x1, x2, x3, at_surface)

implicit none

real(rprec), intent(in)                     :: var_p
!real(rprec) :: tmp_p
real(rprec), dimension(2,2,2), intent(inout)  :: var
logical, intent(in)                         :: at_surface
real(rprec), intent(in)                     :: x1,x2,x3
real(rprec), dimension(2,2,2)               :: wg
!------------------------------------------------------------------------------!

! Compute extrapolation weights
wg(1,1,1)   = (1.-x1)*(1.-x2)*(1.-x3)
wg(2,1,1)   = x1     *(1.-x2)*(1.-x3)
wg(1,2,1)   = (1.-x1)*x2     *(1.-x3)
wg(2,2,1)   = x1     *x2     *(1.-x3)
wg(1,1,2)   = (1.-x1)*(1.-x2)*x3
wg(2,1,2)   = x1     *(1.-x2)*x3
wg(1,2,2)   = (1.-x1)*x2     *x3
wg(2,2,2)   = x1     *x2     *x3

! Save variable to be extrapolated
!tmp_p = var_p

! Take care of nodes in the boundary
! If the bottom nodes correspond to the 1st W nodes (surface) or
! if the bottom nodes correspond to the ghost UVP nodes (below the surface)
if ( at_surface ) then
!  tmp_p = var_p !/ 1.5_rprec
  wg(:,:,1) = 0.0_rprec
  wg = wg/SUM(wg)
else
! tmp_p = var_p
endif

! Update VAR with the contribution of VAR_P
var = var + wg*var_p

end subroutine extrap_3d

!------------------------------------------------------------------------------!
!                                 INTERP_2DH                                   !
!------------------------------------------------------------------------------!
subroutine interp_2dh(coord_p, var, var_p)

implicit none

real(rprec), dimension(3), intent(in)      :: coord_p
real(rprec), dimension(:,:), intent(in)    :: var
real(rprec), intent(out)                   :: var_p
integer                                    :: ix,iy,ix1,iy1,nn,jx,jy
integer, dimension(4,2)                    :: cn
real(rprec)                                :: x1,x2
real(rprec), dimension(4)                  :: wg
!------------------------------------------------------------------------------!

! Indexes of the closest surface nodes
ix  = modulo(floor(coord_p(1)/dx),nx)+1
iy  = modulo(floor(coord_p(2)/dy),ny)+1
ix1 = modulo(ix,nx)+1
iy1 = modulo(iy,ny)+1

! Coordinates of the 4 closest surface nodes
cn(1,:) = (/ix,iy/)
cn(2,:) = (/ix1,iy/)
cn(3,:) = (/ix,iy1/)
cn(4,:) = (/ix1,iy1/)

! Fraction of the distance from parcel location to the bottom-left neighboring node
x1 = modulo(coord_p(1),dx)/dx
x2 = modulo(coord_p(2),dy)/dy

! Compute interpolation weights
wg(1)   = (1.-x1)*(1.-x2)
wg(2)   = x1     *(1.-x2)
wg(3)   = (1.-x1)*x2
wg(4)   = x1     *x2

! Interpolate VAR_P from VAR
var_p = 0._rprec
do nn=1,4
    jx = cn(nn,1)
    jy = cn(nn,2)
    var_p = var_p + wg(nn)*var(jx,jy)
end do

end subroutine interp_2dh

!------------------------------------------------------------------------------!
!                                PRINT_3DVAR                                   !
!------------------------------------------------------------------------------!
subroutine print_3dvar(var_in, var_name, time_str)

implicit none

real(rprec), dimension(ld,ny,0:nz), intent(in) :: var_in
character(*), intent(in)                       :: var_name, time_str
real(rprec), dimension(nx,ny,1:nz)               :: var_prt
character(100)                                 :: file_path
integer                                        :: ix, iy, iz
character(len=128) :: fname

write (fname, '(a,a,a,a,a,i0,a,i0,a)') './output/lsm/',var_name,'/',var_name,'_',coord,'_',jt_total,'.out'
fname = trim(fname)
!open(200,file=fname,form='formatted',action='write',position='append')

var_prt(:,:,:) = var_in(1:nx,1:ny,1:nz)

open(200,file=fname,form='unformatted',status='new',access='direct',recl=8*nx*ny*nz)
write(200,rec=1) var_prt
close(200)
!3333 format(3(e11.4,a1,2x),e14.7)

end subroutine print_3dvar

!------------------------------------------------------------------------------!
!                                PRINT_2DVARH                                  !
!------------------------------------------------------------------------------!
subroutine print_2dvarh(var_in, var_name)

implicit none

real(rprec), dimension(ld,ny), intent(in) :: var_in
character(*), intent(in)                :: var_name
character(100)                          :: file_path
integer                                 :: ix, iy, iz
!------------------------------------------------------------------------------!

if (coord .eq. 0) then ! To be printed by bottom rank only

   write (file_path, '(a,a,i0.8,a)') trim(lsm_path)//var_name,'_',jt_total,'.out'
   open(111,file=trim(file_path),status='unknown')
   do ix=1,nx
      do iy=1,ny
         write(111,*) var_in(ix,iy)
      enddo
   enddo
   close(111)

endif

end subroutine print_2dvarh

!------------------------------------------------------------------------------!
!                                PRINT_2DVARV                                  !
!------------------------------------------------------------------------------!
subroutine print_2dvarv(var_in, var_name, time_str)

implicit none

real(rprec), dimension(ny,0:nz), intent(in) :: var_in
character(*), intent(in)                :: var_name, time_str
character(100)                          :: file_path
integer                                 :: ix, iy, iz

file_path = trim(lsm_path)//var_name//'_'//trim(time_str)//'s.csv'

write (file_path, '(a,a,i0.8,a)') trim(lsm_path)//var_name,'_',jt_total,'.out'

open(200,file=trim(file_path),status='unknown')
    write(200,*) "Coord Y,          ", "Coord Z,        ", "Value"
    do iz=1,nz
    do iy=1,ny
        write(200,2222) (iy-1)*dy,',',(iz-1)*dz,',',var_in(iy,iz)
    end do
    end do

close(200)

2222 format(e11.4,2(a1,2x,e11.4))

end subroutine print_2dvarv

!------------------------------------------------------------------------------!
!                               PRINT_PARCELS (old)                            !
!------------------------------------------------------------------------------!
!subroutine print_parcels(time_str)
!
!implicit none
!
!character(*), intent(in) :: time_str
!character(100)           :: file_path
!integer                  :: ip
!character(len=64) :: fname
!character(10) :: coord_str
!
!write(coord_str,'(i0)') coord
!coord_str = adjustl(coord_str)
!
!if(np .gt. 0) then
!
!file_path = lsm_path//'parcels'//'_'//trim(time_str)//'s_coord_'//trim(coord_str)//'p.csv'
!
!open(200,file=trim(file_path),status='unknown')
!write(200,'(8(a15))') "X, ", "Y, ", "Z, ", "U, ", "V, ","W, ","D, ", "M"
!do ip=1,np
!    write(200,1111) x_p(ip),",",y_p(ip),",",z_p(ip),",",    &
!                    u_p(ip),",",v_p(ip),",",w_p(ip),",",    &
!                    d_p(ip),",",mass_p(ip)
!end do
!
!close(200)
!
!endif
!1111 format(e11.4,7(a1,2x,e11.4))
!
!end subroutine print_parcels

!------------------------------------------------------------------------------!
!                                PRINT_PARCELS                                 !
!------------------------------------------------------------------------------!
subroutine print_parcels()

implicit none

integer :: i
character(len=64) :: fname
!------------------------------------------------------------------------------!

if (np .gt. 0) then
   write (fname, '(a,i0,a,i0,a,i0.8,a)') trim(lsm_path)//'/particles_',coord,'_',colour,'_',jt_total,'.out'
   open(200,file=trim(fname),form='formatted',status='replace',action='write',position='append')
   do i=1,np
      IF (sublimation) then
         write(200,*) jt_total,p_ip(i),x_p(i),y_p(i),z_p(i),u_p(i),v_p(i),w_p(i),T_p(i),d_p(i), &
                      re_p(i),nu_p(i),sh_p(i),mass_p(i),delta_m(i),delta_r(i),delta_T(i),timescale_stokes(i)
      else
         write(200,*) jt_total,p_ip(i),x_p(i),y_p(i),z_p(i),u_p(i),v_p(i),w_p(i),ax_p(i),ay_p(i),az_p(i),d_p(i), &
                      re_p(i),mass_p(i),delta_m(i),delta_r(i),delta_T(i),timescale_stokes(i)
      endif
   enddo
   close(200)
endif

end subroutine print_parcels

!------------------------------------------------------------------------------!
!                                PRINT_DISSTKE                                 !
!------------------------------------------------------------------------------!
subroutine print_disstke(time_str)

implicit none

character(*), intent(in) :: time_str
character(100)           :: file_path
integer                  :: ix,iy,iz

file_path = trim(lsm_path)//'/diss_restke_'//trim(time_str)//'.csv'

open(200,file=trim(file_path),status='unknown')

do iz=1,nz
do iy=1,ny
do ix=1,nx
    write(200,'(5(e11.4,3x))') (ix-1)*dx,(iy-1)*dy,(iz-1)*dz,diss(ix,iy,iz),restke(ix,iy,iz)
end do
end do
end do
close(200)

1111 format(e11.4,14(a1,2x,e11.4))

end subroutine print_disstke

!------------------------------------------------------------------------------!
!                             PRINT_MASS_BALANCE                               !
!------------------------------------------------------------------------------!
subroutine print_mass_balance(time)

implicit none

real(rprec), intent(in) :: time
character(100)          :: file_path
real(rprec) :: u_star2
!------------------------------------------------------------------------------!

! Compute entrained mass flux
ve_flux = ve_flux/(l_x*l_y*z_i*z_i)/(ref_dt*z_i*print1_freq)   ![kg/m2/s]

! Compute deposited mass flux
vd_flux = vd_flux/(l_x*l_y*z_i*z_i)/(ref_dt*z_i*print1_freq)   ![kg/m2/s]

! Compute ejected (splashed) mass flux
vs_flux   = vs_flux/(l_x*l_y*z_i*z_i)/(ref_dt*z_i*print1_freq) ![kg/m2/s]

! Compute mean u*^2 along the surface
u_star2  = sum(sum(tau/f_rho,1)/nx)/ny

!call mpi_allreduce(mpi_in_place,ve_flux,1,mpi_rprec,mpi_sum,comm,ierr)
!call mpi_allreduce(mpi_in_place,vd_flux,1,mpi_rprec,mpi_sum,comm,ierr)
!call mpi_allreduce(mpi_in_place,vs_flux,1,mpi_rprec,mpi_sum,comm,ierr)

if (coord .eq. 0) then ! To be printed by bottom rank only
   file_path = trim(lsm_path)//'mass_balance.txt'
   open(300,file=trim(file_path),status='old',position='append')
   write(300,'(f10.3,5(3x,e11.4))') time,u_star2,ve_flux,vd_flux,vs_flux,av_conc/(l_x*l_y*l_z*z_i*z_i*z_i) !last term in [kg/m3]
   close(300)
endif

! Reset mass fluxes
ve_flux = 0._rprec
vd_flux = 0._rprec
vs_flux = 0._rprec

end subroutine print_mass_balance

!------------------------------------------------------------------------------!
!                             PRINT_MOM_BALANCE                                !
!------------------------------------------------------------------------------!
subroutine print_mom_balance(time)

implicit none

real(rprec), intent(in) :: time
real(rprec)             :: px_force,pz_force,u_star2
character(100)          :: file_path
!------------------------------------------------------------------------------!
! Sum in height of mean Fx_lsm and Fz_lsm in (x,y)
px_force = sum(sum(sum(Fx_lsm,1)/nx,1)/ny)*dz
pz_force = sum(sum(sum(Fz_lsm,1)/nx,1)/ny)*dz

! Mean u*^2
u_star2  = sum(sum(tau/f_rho,1)/nx)/ny

write(file_path,'(a,i0,a)') trim(lsm_path)//'mom_balance_c',coord,'.txt'
open(300,file=trim(file_path),status='old',position='append')
write(300,'(f10.3,3(3x,e11.4))') time,px_force,pz_force,u_star2
close(300)

end subroutine print_mom_balance

!------------------------------------------------------------------------------!
!                         PRINT_INDIVIDUAL_PARTICLES                           !
!------------------------------------------------------------------------------!
subroutine print_individual_particles()

implicit none

integer :: i
character(len=64) :: fname

do i=1,np
  if(mod(p_ip(i),1000).eq.0)then
    write (fname, '(a,i0.10,a)') './individual_particle_',p_ip(i) ,'.out'
    fname = trim(fname)
    open(200,file=fname,form='formatted',action='write',position='append')
    if(sublimation) then
    write(200,*) jt_total,p_ip(i),x_p(i),y_p(i),z_p(i),u_p(i),v_p(i),w_p(i),T_p(i),d_p(i)!, &
                 !re_p(i),nu_p(i),sh_p(i),mass_p(i),delta_m(i),delta_r(i),delta_T(i)
    else
    write(200,*) jt_total,p_ip(i),x_p(i),y_p(i),z_p(i),u_p(i),v_p(i),w_p(i),d_p(i)
    endif
    close(200)
  endif
enddo

end subroutine print_individual_particles

!------------------------------------------------------------------------------!
!                                OUTPUT_MANAGER                                !
!------------------------------------------------------------------------------!
subroutine output_manager()

implicit none

real(rprec)   :: time
character(10) :: time_str
character(64) :: filename
!------------------------------------------------------------------------------!

time = lsm_dt*lsm_step*(z_i)
write(time_str,'(f10.3)') time
time_str = adjustl(time_str)

call compute_concentration()
call compute_flux()

! Write on files at PRINT1_FREQ frequency
if ((mod((jt_total-lsm_init),print1_freq).eq.0)) then
   if (lp_stress)   call print_2dvarh(tau, 'stress')
   if (lp_stress)   call print_2dvarh(lsm_surface,'lsm_surface')
   if (lp_stress)   call print_2dvarh(entrainment,'entrainment')
   if (lp_diss)     call print_disstke(time_str)
   if (lp_massbal)  call print_mass_balance(time)
   if (lp_mombal)   call print_mom_balance(time)
endif

! Write total number of parcels aloft at C_CNT4 frequency
if (mod(jt,c_cnt4)==0)then
   if (coord .eq. 0) then
      open(111,file='./np_global.dat',status='replace')
      write(111,*) np_global
      close(111)
   endif
endif

! Write mass and heat output at each 100 time steps
if (mod(jt_total,100) .eq. 0) then
   write(filename,'(a,i0,a)') trim(lsm_path)//'mass_output.c',coord,'.log'
   open(200,file=trim(filename),form='formatted',action='write',position='append')
   write(200,*) jt_total,sum_delta_m
   close(200)

   write(filename,'(a,i0,a)') trim(lsm_path)//'heat_output.c',coord,'.log'
   open(200,file=trim(filename),form='formatted',action='write',position='append')
   write(200,*) jt_total,sum_delta_q
   close(200)
endif

end subroutine output_manager

!------------------------------------------------------------------------------!
!                                PRINT_RESTART                                 !
!------------------------------------------------------------------------------!
subroutine print_restart(time_str)

  implicit none

  character(*), intent(in)    :: time_str
  character(100)              :: file_path
  integer                     :: ip, ix, iy

  file_path = trim(lsm_path)//'rst/parcels_restart.csv'
  open(200,file=trim(file_path),status='unknown')
  write(200,*) np, trim(time_str)
  do ip=1,np
      write(200,1112) x_p(ip),y_p(ip),z_p(ip),u_p(ip),v_p(ip),w_p(ip), &
                      u_sgs(ip),v_sgs(ip),w_sgs(ip),sigma2(ip),        &
                      ax_p(ip),ay_p(ip),az_p(ip),d_p(ip),mass_p(ip)
  end do
  close(200)

  file_path = trim(lsm_path)//'rst/surface_restart.csv'
  open(200,file=trim(file_path),status='unknown')
  do iy=1,ny
      do ix=1,nx
          write(200,2222) (ix-1)*dx,',',(iy-1)*dy,',',lsm_surface(ix,iy)
      end do
  end do
  close(200)

  if (mode.eq.'saltation') then
      file_path = trim(lsm_path)//'rst/entrainment_restart.csv'
      open(200,file=trim(file_path),status='unknown')
      do iy=1,ny
          do ix=1,nx
              write(200,2222) (ix-1)*dx,',',(iy-1)*dy,',',entrainment(ix,iy)
          end do
      end do
      close(200)
  end if

  1112 format(15(e11.4,2x))
  2222 format(e11.4,2(a1,2x,e11.4))

end subroutine print_restart

!------------------------------------------------------------------------------!
!                             SAMPLE_DISTRIBUTION                              !
!------------------------------------------------------------------------------!
subroutine sample_distribution(mean, std, dist, val)

implicit none

real(rprec), intent(in)  :: mean, std ! mean and standard deviation of the distribution
character(*), intent(in) :: dist
real(rprec), intent(out) :: val
real(rprec)              :: m, s2, rand1, rand2
!------------------------------------------------------------------------------!

! Lognormal distribution
if (trim(dist).eq.'lognormal') then
    s2 = log(1._rprec + (std/mean)**2)
    m  = log(mean) - 0.5_rprec*s2
    call random_number(rand1)
    call random_number(rand2)
    val = m + sqrt(s2)*sqrt(-2._rprec*log(rand1))*cos(2._rprec*pi*rand2)
    val = exp(val)

! Exponential distribution
else if (trim(dist).eq.'exponential') then
    call random_number(rand1)
    val = -mean*log(1._rprec-rand1)

! Normal distribution
else if (trim(dist).eq.'normal') then
    call random_number(rand1)
    call random_number(rand2)
    val = mean + std*sqrt(-2._rprec*log(rand1))*cos(2._rprec*pi*rand2)
end if

end subroutine sample_distribution

!------------------------------------------------------------------------------!
!                                INTERP_VEL                                    !
!------------------------------------------------------------------------------!
subroutine interp_vel(coord_p,val_u,val_v,val_w,loc_ip)

implicit none

real(rprec), intent(in), dimension(3)   :: coord_p
real(rprec), intent(out)                :: val_u, val_v, val_w
real(rprec)                             :: phi_p, vel_tp, wall_d
real(rprec), dimension(3)               :: xv, ns, vel, vel_t, x_diff_uvp, x_diff_w
real(rprec), dimension(2,2,2)           :: var_section
integer,intent(in)                      :: loc_ip
integer, dimension(2)                   :: ix_uvp,iy_uvp,iz_uvp, ix_w,iy_w,iz_w ! indexes for interpolation
logical                                 :: at_surf
!------------------------------------------------------------------------------!

if (llbc_special.eq.'IBM_DFA') then
    wall_d = wall_c*dz
else
    wall_d = get_z_local(1) !0.5_rprec*dz
end if

! get indexes of closest W nodes
call get_indexes_3d(coord_p,1._rprec,ix_w,iy_w,iz_w,x_diff_w,at_surf,3,loc_ip)

! Parcel height above surface [m]
if(llbc_special.ne.'IBM_DFA') then
  phi_p = z_p(loc_ip)
else
  ! Get parcel height normal to the surface by interpolating phi between the 8 closest grid nodes
  var_section = phi(ix_w,iy_w,iz_w)
  call interp_var_lsm(var_section, phi_p, x_diff_w(1), x_diff_w(2), x_diff_w(3))
endif

! If parcel is above the 1st UVP node
if (phi_p.gt.wall_d) then
    ! Get indexes of the closest UVP nodes
    call get_indexes_3d(coord_p,0.5_rprec,ix_uvp,iy_uvp,iz_uvp,x_diff_uvp,at_surf,4,loc_ip)

    ! Get flow velocity at parcel location
    ! u
    var_section = u(ix_uvp,iy_uvp,iz_uvp)
    call interp_var_lsm(var_section, val_u, x_diff_uvp(1),x_diff_uvp(2),x_diff_uvp(3))
    ! v
    var_section = v(ix_uvp,iy_uvp,iz_uvp)
    call interp_var_lsm(var_section, val_v, x_diff_uvp(1),x_diff_uvp(2),x_diff_uvp(3))
    ! w
    var_section = w(ix_w,iy_w,iz_w)
    call interp_var_lsm(var_section, val_w, x_diff_w(1),x_diff_w(2),x_diff_w(3))

! If parcel is between zo and the 1st UVP node
else if ((phi_p.le.wall_d).and.(phi_p.gt.lbc_dz0)) then
    ! Get surface normal vector at particle location
    call interp_2dh(coord_p,nis(1,:,:),ns(1))
    call interp_2dh(coord_p,nis(2,:,:),ns(2))
    call interp_2dh(coord_p,nis(3,:,:),ns(3))
    ns = ns/mag(ns)

    ! Compute coordinates of corresponding 1st UVP node
    xv(1) = coord_p(1) + ns(1)*(wall_d-phi_p)
    xv(2) = coord_p(2) + ns(2)*(wall_d-phi_p)
    xv(3) = get_eta_local(wall_d)

    ! get node indexes closest to corresponding 1st UVP node
    call get_indexes_3d(xv,1._rprec,ix_w,iy_w,iz_w,x_diff_w,at_surf,5,loc_ip)
    call get_indexes_3d(xv,0.5_rprec,ix_uvp,iy_uvp,iz_uvp,x_diff_uvp,at_surf,6,loc_ip)

    ! Get flow velocity at the corresponding 1st UVP node
    var_section = u(ix_uvp,iy_uvp,iz_uvp)
    call interp_var_lsm(var_section, vel(1), x_diff_uvp(1),x_diff_uvp(2),x_diff_uvp(3))
    var_section = v(ix_uvp,iy_uvp,iz_uvp)
    call interp_var_lsm(var_section, vel(2), x_diff_uvp(1),x_diff_uvp(2),x_diff_uvp(3))
    var_section = w(ix_w,iy_w,iz_w)
    call interp_var_lsm(var_section, vel(3), x_diff_w(1),x_diff_w(2),x_diff_w(3))

    ! tangential flow velocity vector at 1st UVP node (horizontal component)
    vel_t = vel-dot_product(vel,ns)*ns
    ! Compute horizontal flow velocity from local logarithmic law although it is only valid for averaged speed but at least it approaches the resolved velocity at the 1st UVP node
    vel_tp = mag(vel_t) * log( phi_p/lbc_dz0 ) / log( wall_d/lbc_dz0 )
    val_u = vel_tp *vel_t(1)/mag(vel_t)
    val_v = vel_tp *vel_t(2)/mag(vel_t)
    ! Interpolate vertical flow velocity such that we obtain the same orientation of the velocity vector as at the corresponding 1st UVP node
    val_w = vel(3) * vel_tp / mag(vel_t)

! If parcel is bellow zo
else
    val_u = 0._rprec ! The flow velocity is assumed 0 below zo
    val_v = 0._rprec
    val_w = 0._rprec
end if

end subroutine interp_vel

!------------------------------------------------------------------------------!
!                                GET_INDEXES_3D                                !
!------------------------------------------------------------------------------!
subroutine get_indexes_3d(coord_p,stag,jx_close,jy_close,jz_close,rel_dist,at_surface,from_where,loc_ip)

implicit none

real(rprec), intent(in)               :: stag
real(rprec), dimension(3), intent(in) :: coord_p
real(rprec), dimension(3), intent(out) :: rel_dist
integer, dimension(2), intent(out)    :: jx_close,jy_close,jz_close
logical, intent(out)                  :: at_surface
INTEGER,INTENT(IN)                    :: from_where
integer,intent(in)                    :: loc_ip
!------------------------------------------------------------------------------!

! Index of the closest node with a coordinate lower than or equal to that of the parcel
jx_close(1) = modulo(floor(coord_p(1)/dx),nx)+1
jy_close(1) = modulo(floor(coord_p(2)/dy),ny)+1
jz_close(1) = floor(coord_p(3)/dz+stag)
! Index of closest node with a coordinate higher than that of the parcel
! TRY TO HANDLE BOUNDARIES NICELY FOR THE +1 INDICES
jx_close(2) = modulo(jx_close(1),nx)+1
jy_close(2) = modulo(jy_close(1),ny)+1
jz_close(2) = jz_close(1)+1

! Fraction of the distance from parcel location to the front-bottom-left neighboring node
rel_dist(1) = modulo(coord_p(1),dx)/dx
rel_dist(2) = modulo(coord_p(2),dy)/dy
rel_dist(3) = modulo(coord_p(3)-stag*dz,dz)/dz

! Define indicator for a parcel position between the surface and the first node level above (for later use in extrap_3d)
if( ((stag .eq. 1.0_rprec) .and. (jz_close(1) .eq. 1) .and. (coord .eq. 0)) .OR. &
    ((stag .eq. 0.5_rprec) .and. (jz_close(1) .eq. 0) .and. (coord .eq. 0)) ) then
  at_surface = .true.
else
  at_surface = .false.
endif

!BOUND CHECK
IF((jx_close(1) < 1).OR.(jx_close(1) > NX))THEN
  WRITE(*,*) 'In get_indexes_3d'
  PRINT*,'stag=',stag
  PRINT*,'coord_p=',coord_p
  PRINT*,'OUT OF RANGE: jx_close(1)=',jx_close(1)
  PRINT*,'jy_close(1)=',jy_close(1)
  PRINT*,'jz_close(1)=',jz_close(1)
  PRINT*,'CALLED FROM=',from_where
  PRINT*,'PROC=',coord
  STOP
ENDIF
IF((jy_close(1) < 1).OR.(jy_close(1) > NY))THEN
  WRITE(*,*) 'In get_indexes_3d'
  PRINT*,'stag=',stag
  PRINT*,'coord_p=',coord_p
  PRINT*,'OUT OF RANGE: jy_close(1)=',jy_close(1)
  PRINT*,'jx_close(1)=',jx_close(1)
  PRINT*,'jz_close(1)=',jz_close(1)
  PRINT*,'CALLED FROM=',from_where
  PRINT*,'PROC=',coord
  STOP
ENDIF
IF((jz_close(1) < 0).OR.(jz_close(1) > nz))THEN
  WRITE(*,*) 'In get_indexes_3d'
  PRINT*,'stag=',stag
  PRINT*,'coord_p=',coord_p
  PRINT*,'OUT OF RANGE: jz_close(1)=',jz_close(1)
  PRINT*,'jx_close(1)=',jx_close(1)
  PRINT*,'jy_close(1)=',jy_close(1)
  PRINT*,'CALLED FROM=',from_where
  PRINT*,'PROC=',coord
  write(*,*) 'p_ip:',p_ip(loc_ip)
  write(*,*) 'u,v,w: ',u_p(loc_ip),v_p(loc_ip),w_p(loc_ip)
  write(*,*) 'u_sgs,v_sgs,w_sgs: ', u_sgs(loc_ip),v_sgs(loc_ip),w_sgs(loc_ip)
  write(*,*) 'd_p: ',d_p(loc_ip)
  STOP
ENDIF

end subroutine get_indexes_3d

!------------------------------------------------------------------------------!
!                                INTERP_T                                      !
!------------------------------------------------------------------------------!
subroutine interp_T(coord_p,val_T,loc_ip)

implicit none

integer, intent(in)                     :: loc_ip
real(rprec), intent(in), dimension(3)   :: coord_p
real(rprec), intent(out)                :: val_T
real(rprec)                             :: phi_p,wall_d
real(rprec), dimension(3)               :: xv, ns
real(rprec), dimension(3)               :: tmp_coord, x_diff
real(rprec), dimension(2,2,2)           :: var_section
integer, dimension(2)                   :: ix_close,iy_close,iz_close ! indexes for interpolation
real(rprec)                             :: tmp_T_1
logical                                 :: at_surf

if (llbc_special.eq.'IBM_DFA') then
    wall_d = wall_c*dz
else
    wall_d = get_z_local(1) ! 0.5_rprec*dz
end if

! Parcel height above surface [m]
if (llbc_special.ne.'IBM_DFA') then
  phi_p = z_p(loc_ip)
else
  ! get indexes of closest W nodes
  call get_indexes_3d(coord_p,1._rprec,ix_close,iy_close,iz_close,x_diff,at_surf,7,loc_ip)
  ! Interpolate phi between the 8 closest grid nodes
  var_section = phi(ix_close,iy_close,iz_close)
  call interp_var_lsm(var_section, phi_p, x_diff(1), x_diff(2), x_diff(3))
endif

if (phi_p.gt.wall_d) then
   ! Get indexes of the closest UVP nodes
   call get_indexes_3d(coord_p,0.5_rprec,ix_close,iy_close,iz_close,x_diff,at_surf,8,loc_ip)
   ! Interpolate SC1
   var_section = SC1(ix_close,iy_close,iz_close)
   call interp_var_lsm(var_section, val_T, x_diff(1),x_diff(2),x_diff(3))
else if ((phi_p.le.wall_d).and.(phi_p.gt.LBC_SC1_ZO(1,1))) then
   ! Coordinates for the corresponding UVP node just above the parcel
   tmp_coord(1) = coord_p(1)
   tmp_coord(2) = coord_p(2)
   tmp_coord(3) = get_eta_local(wall_d)
   ! Get indexes of the UVP nodes closest to tmp_coord
   call get_indexes_3d(tmp_coord,0.5_rprec,ix_close,iy_close,iz_close,x_diff,at_surf,9,loc_ip)
   ! Interpolate SC1 at the corresponding 1st UVP node
   var_section = SC1(ix_close,iy_close,iz_close)
   call interp_var_lsm(var_section, tmp_T_1, x_diff(1),x_diff(2),x_diff(3))
  if(LLBC_SC1 .eq. 'DRC') then
   ! Compute air temperature from local logarithmic law although it is only valid for averaged temperature but at least it approaches the resolved temperature at the 1st UVP node
   val_T = DLBC_SC1 + (tmp_T_1 - DLBC_SC1) * log( phi_p/LBC_SC1_ZO(1,1) ) / log( wall_d/LBC_SC1_ZO(1,1) )
  else
   ! Assume vertical gradient of zero and use temperature at first uvp level above the surface
   val_T = tmp_T_1 ! 1.0_rprec
  endif
else
  if(LLBC_SC1 .eq. 'DRC') then
    val_T = DLBC_SC1
  else
   tmp_coord(1) = coord_p(1)
   tmp_coord(2) = coord_p(2)
   tmp_coord(3) = get_eta_local(wall_d)
   ! Get indexes of the UVP nodes closest to tmp_coord
   call get_indexes_3d(tmp_coord,0.5_rprec,ix_close,iy_close,iz_close,x_diff,at_surf,10,loc_ip)
   ! Interpolate SC1
   var_section = SC1(ix_close,iy_close,iz_close)
   call interp_var_lsm(var_section, tmp_T_1, x_diff(1),x_diff(2),x_diff(3))
   val_T = tmp_T_1 !1.0_rprec
  endif
end if

end subroutine interp_T

!------------------------------------------------------------------------------!
!                                INTERP_Q                                      !
!------------------------------------------------------------------------------!
subroutine interp_Q(coord_p,val_Q,loc_ip)

implicit none

integer, intent(in)                     :: loc_ip
real(rprec), intent(in), dimension(3)   :: coord_p
real(rprec), intent(out)                :: val_Q
real(rprec)                             :: phi_p,wall_d
real(rprec), dimension(3)               :: xv, ns
real(rprec) :: tmp_Q_1 !,tmp_Q_2
real(rprec), dimension(3) :: tmp_coord, x_diff
real(rprec), dimension(2,2,2)           :: var_section
integer, dimension(2)                   :: ix_close,iy_close,iz_close ! indexes for interpolation
logical                                 :: at_surf

if (llbc_special.eq.'IBM_DFA') then
    wall_d = wall_c*dz
else
    wall_d = get_z_local(1) !0.5_rprec*dz
end if

! Parcel height above surface [m]
if (llbc_special.ne.'IBM_DFA') then
  phi_p = z_p(loc_ip)
else
  ! get indexes of closest W nodes
  call get_indexes_3d(coord_p,1._rprec,ix_close,iy_close,iz_close,x_diff,at_surf,11,loc_ip)
  ! Interpolate phi between the 8 closest grid nodes
  var_section = phi(ix_close,iy_close,iz_close)
  call interp_var_lsm(var_section, phi_p, x_diff(1), x_diff(2), x_diff(3))
endif

if (phi_p.gt.wall_d) then
    ! Get indexes of the closest UVP nodes
    call get_indexes_3d(coord_p,0.5_rprec,ix_close,iy_close,iz_close,x_diff,at_surf,12,loc_ip)
    ! Interpolate SC2
    var_section = SC2(ix_close,iy_close,iz_close)
    call interp_var_lsm(var_section, tmp_Q_1, x_diff(1),x_diff(2),x_diff(3))
    !call interp_var_lsm(coord_p,SC2,0.5_rprec,tmp_Q_1,15,666,666)
    !call interp_var_lsm(coord_p,SC3,0.5_rprec,tmp_Q_2,16,666,666)
    val_Q = tmp_Q_1 !+tmp_Q_2
else if ((phi_p.le.wall_d).and.(phi_p.gt.LBC_SC2_ZO(1,1))) then
   ! Coordinates for the corresponding UVP node just above the parcel
   tmp_coord(1) = coord_p(1)
   tmp_coord(2) = coord_p(2)
   tmp_coord(3) = get_eta_local(wall_d)
   ! Get indexes of the UVP nodes closest to tmp_coord
   call get_indexes_3d(tmp_coord,0.5_rprec,ix_close,iy_close,iz_close,x_diff,at_surf,13,loc_ip)
   ! Interpolate SC2 at the corresponding 1st UVP node
   var_section = SC2(ix_close,iy_close,iz_close)
   call interp_var_lsm(var_section, tmp_Q_1, x_diff(1),x_diff(2),x_diff(3))
  if(LLBC_SC2 .eq. 'DRC') then
   ! Compute humidity from local logarithmic law although it is only valid for averaged humidity but at least it approaches the resolved humidity at the 1st UVP node
   val_Q = DLBC_SC2 + (tmp_Q_1 - DLBC_SC2) * log( phi_p/LBC_SC2_ZO(1,1) ) / log( wall_d/LBC_SC2_ZO(1,1) )
  else
   ! Assume vertical gradient of zero and use humidity at first uvp level above the surface
   val_Q = tmp_Q_1
  endif
else
   if(LLBC_SC2 .eq. 'DRC') then
     val_Q = DLBC_SC2 !+ DLBC_SC3
   else
    tmp_coord(1) = coord_p(1)
    tmp_coord(2) = coord_p(2)
    tmp_coord(3) = get_eta_local(wall_d)
    ! Get indexes of the UVP nodes closest to tmp_coord
    call get_indexes_3d(tmp_coord,0.5_rprec,ix_close,iy_close,iz_close,x_diff,at_surf,14,loc_ip)
    ! Interpolate SC2
    var_section = SC2(ix_close,iy_close,iz_close)
    call interp_var_lsm(var_section, tmp_Q_1, x_diff(1),x_diff(2),x_diff(3))
    val_Q = tmp_Q_1
   endif
end if

end subroutine interp_Q

!------------------------------------------------------------------------------!
!                                AVERAGE_LES                                   !
!------------------------------------------------------------------------------!
subroutine average_les()

implicit none

real(rprec), parameter             :: c_e = 0.93 ! value from Pope pg. 631
integer                            :: ix,iy,iz
real(rprec), dimension(ld,ny,1:nz) :: LSM_S11, LSM_S22, LSM_S33, LSM_S13, LSM_S23, LSM_S12
real(rprec), dimension(ld,ny,1:nz) :: tmp
integer :: kk
real(rprec) :: fac_1,fac_2,fac_3
! average dissipation (see pope p 587)
call calc_Sij(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,LSM_S11,LSM_S12,LSM_S13,LSM_S23,LSM_S22,LSM_S33)

! assume everything is calculated in w-nodes, even if node=1 is uvp

!$omp parallel do
do iz=1,nz
 !do ix=1,64
 ! do iy=1,64
 !  write(*,*) coord,ix,iy,iz,LSM_S11(ix,iy,iz),LSM_S22(ix,iy,iz),LSM_S33(ix,iy,iz),LSM_S13(ix,iy,iz),LSM_S23(ix,iy,iz),LSM_S12(ix,iy,iz)
 ! enddo
 !enddo
kk = iz+eta_displ(coord)

delta = (dx*dy*(get_z_local(2*kk-1)-get_z_local(2*kk-3)))**(1._rprec/3._rprec)

tmp(:,:,iz) = ((2._rprec*(LSM_S11(:,:,iz)**2+LSM_S22(:,:,iz)**2+LSM_S33(:,:,iz)**2+      &
                     2._rprec*LSM_S13(:,:,iz)**2+2._rprec*LSM_S23(:,:,iz)**2+               &
                     2._rprec*LSM_S12(:,:,iz)**2))**(3./2.))

sum_strain(:,:,iz) = (1.0_rprec - averaging_weight)*sum_strain(:,:,iz) + averaging_weight*tmp(:,:,iz)

sum_w2(:,:,iz) = (1.0_rprec - averaging_weight)*sum_w2(:,:,iz) + (w(:,:,iz)*w(:,:,iz))*averaging_weight
sum_w(:,:,iz) = (1.0_rprec - averaging_weight)*sum_w(:,:,iz) + w(:,:,iz)*averaging_weight
sum_cs(:,:,iz) = (1.0_rprec - averaging_weight)*sum_cs(:,:,iz) + cs_opt2(:,:,iz)*averaging_weight
sum_diss(:,:,iz) = (1.0_rprec - averaging_weight)*sum_diss(:,:,iz) + averaging_weight*(cs_opt2(:,:,iz)*(delta**2.0_rprec)* tmp(:,:,iz))
enddo
!$omp end parallel do

!$omp parallel do
do iz=0,nz

sum_u2(:,:,iz) = (1.0_rprec - averaging_weight)*sum_u2(:,:,iz) + (u(:,:,iz)*u(:,:,iz))*averaging_weight
sum_v2(:,:,iz) = (1.0_rprec - averaging_weight)*sum_v2(:,:,iz) + (v(:,:,iz)*v(:,:,iz))*averaging_weight

sum_u(:,:,iz) = (1.0_rprec - averaging_weight)*sum_u(:,:,iz) + u(:,:,iz) * averaging_weight
sum_v(:,:,iz) = (1.0_rprec - averaging_weight)*sum_v(:,:,iz) + v(:,:,iz) * averaging_weight

enddo
!$omp end parallel do

!if(mod(jt,average_freq).eq.0.0_rprec) then

    tke_u = 0.5_rprec*(sum_u2-(sum_u**2.0_rpreC))
    tke_v = 0.5_rprec*(sum_v2-(sum_v**2.0_rprec))
    tke_w = 0.5_rprec*(sum_w2-(sum_w**2.0_rprec))

    !write(*,*) 'jt,coord,u:',jt,coord,minval(tke_u)
    !write(*,*) 'jt,coord,v:',jt,coord,minval(tke_v)
    !write(*,*) 'jt,coord,w:',jt,coord,minval(tke_w)

    !$omp parallel do
    do iz=1,nz
      kk = iz+eta_displ(coord)

        fac_1=get_z_local(2*kk-1)
        fac_2=get_z_local(2*kk-3)
        fac_3=get_z_local(2*kk-2)

        delta = (dx*dy*(get_z_local(2*kk-1)-get_z_local(2*kk-3)))**(1._rprec/3._rprec)
        !delta = (dx*dy*get_dz_local(real(kk,rprec)))**(1._rprec/3._rprec)
        !diss(:,:,iz)   = sum_cs(:,:,iz)*(delta**2.0_rprec)*sum_strain(:,:,iz)
        diss(:,:,iz) = sum_diss(:,:,iz)
        !restke(:,:,iz) = 0.5_rprec*(tke_u(:,:,iz)+tke_u(:,:,iz-1)) + 0.5_rprec*(tke_v(:,:,iz)+tke_v(:,:,iz-1)) + &
        !                 tke_w(:,:,iz)

        restke(:,:,iz) = tke_u(:,:,iz-1) + ((tke_u(:,:,iz)-tke_u(:,:,iz-1))/(fac_1-fac_2))*(fac_3-fac_2) + &
                         tke_v(:,:,iz-1) + ((tke_v(:,:,iz)-tke_v(:,:,iz-1))/(fac_1-fac_2))*(fac_3-fac_2) + &
                         tke_w(:,:,iz)

        sgstke(:,:,iz) = (diss(:,:,iz)*delta/c_e)**(2._rprec/3._rprec)
    end do
    !$omp end parallel do

call mpi_sendrecv(diss(1,1,1),ld*ny,mpi_rprec,down,510,diss(1,1,nz),ld*ny,mpi_rprec,up,510,comm,status,ierr)
call mpi_sendrecv(diss(1,1,nz-1),ld*ny,mpi_rprec,up,610,diss(1,1,0),ld*ny,mpi_rprec,down,610,comm,status,ierr)

call mpi_sendrecv(restke(1,1,1),ld*ny,mpi_rprec,down,520,restke(1,1,nz),ld*ny,mpi_rprec,up,520,comm,status,ierr)
call mpi_sendrecv(restke(1,1,nz-1),ld*ny,mpi_rprec,up,620,restke(1,1,0),ld*ny,mpi_rprec,down,620,comm,status,ierr)

call mpi_sendrecv(sgstke(1,1,1),ld*ny,mpi_rprec,down,530,sgstke(1,1,nz),ld*ny,mpi_rprec,up,530,comm,status,ierr)
call mpi_sendrecv(sgstke(1,1,nz-1),ld*ny,mpi_rprec,up,630,sgstke(1,1,0),ld*ny,mpi_rprec,down,630,comm,status,ierr)


    if(coord .eq. 0) then
     diss(:,:,0)   = 0.0_rprec !diss(:,:,1)
     restke(:,:,0) = 0.0_rprec !restke(:,:,1)
     sgstke(:,:,0) = 0.0_rprec !sgstke(:,:,1)
    endif
     !sum_u(:,:,:)  = 0._rprec
     !sum_v(:,:,:)  = 0._rprec
     !sum_w(:,:,:)  = 0._rprec
     !sum_u2(:,:,:) = 0.0_rprec
     !sum_v2(:,:,:) = 0.0_rprec
     !sum_w2(:,:,:) = 0.0_rprec
     !tke_u(:,:,:) = 0.0_rprec
     !tke_v(:,:,:) = 0.0_rprec
     !tke_w(:,:,:) = 0.0_rprec
     !sum_strain(:,:,:) = 0._rprec
     !sum_cs(:,:,:) = 0.0_rprec
     !sum_diss(:,:,:) = 0.0_rprec
!end if

end subroutine average_les

!------------------------------------------------------------------------------!
!                               EXCHANGE_PARCELS                               !
!------------------------------------------------------------------------------!
subroutine exchange_parcels()

implicit none

integer :: ip
integer :: n_rank_needed
integer :: np_exchange,global_exchange_p
integer :: i
integer,dimension(nproc) :: incoming_data_size,recvcounts,displs
integer :: ierror
real(rprec),dimension(3) :: coord_p
real(rprec) :: g_coord_uvp
integer :: proc_z
real(rprec),dimension(:),allocatable :: exchange_data
integer(kind=int64),dimension(:),allocatable :: exchange_data_int
integer :: np_old
!------------------------------------------------------------------------------!

! Loop on parcels
!!$omp parallel do private(g_coord_uvp,proc_z), shared(np,z_p,p_exchange,p_exist)
do ip=1,np

   ! Compute rank where parcel belongs
   g_coord_uvp = (get_eta_local(z_p(ip))/dz) + 0.5_rprec
   proc_z = COUNT( REAL(NZ_CUMUL, RPREC) < g_coord_uvp ) !ceiling(g_coord_uvp/real((nz-1),rprec)) - 1

   if (proc_z .ne. coord) then ! Parcel is in the wrong rank
      p_exchange(ip) = .true.
      p_exist(ip) = .false.
   endif

   !if(n_rank_needed .eq. coord) then
   !  p_exchange(ip) = .false.
   !  p_exist(ip) = .true.
   !endif
enddo
!!$omp end parallel do

! Set send/receive parameters
displs = 0
do i=0,nproc-1
   displs(i+1) = i
enddo
recvcounts(:) = 1

! Number of parcels current rank has to exchange
np_exchange = count(p_exchange(1:np))

! Gather number of parcels to be exchanged from all ranks
call mpi_allgatherv(np_exchange,1,MPI_INTEGER,incoming_data_size, &
                    recvcounts,displs,MPI_INTEGER,comm,ierror)
call mpi_barrier(comm,ierr)

! Compute total number of parcels to be exchanged by all ranks
global_exchange_p = sum(incoming_data_size,1)

! Allocate parcel data to be exchanged
allocate(exchange_x_p(global_exchange_p))
allocate(exchange_y_p(global_exchange_p))
allocate(exchange_z_p(global_exchange_p))

allocate(exchange_u_p(global_exchange_p))
allocate(exchange_v_p(global_exchange_p))
allocate(exchange_w_p(global_exchange_p))

allocate(exchange_u_sgs_p(global_exchange_p))
allocate(exchange_v_sgs_p(global_exchange_p))
allocate(exchange_w_sgs_p(global_exchange_p))

allocate(exchange_ax_p(global_exchange_p))
allocate(exchange_ay_p(global_exchange_p))
allocate(exchange_az_p(global_exchange_p))

allocate(exchange_mass_p(global_exchange_p))
allocate(exchange_d_p(global_exchange_p))
allocate(exchange_sigma_p(global_exchange_p))

IF(temperature) then
allocate(exchange_T_p(global_exchange_p))
endif

allocate(exchange_np_global(global_exchange_p))

! Set send/receive parameters
displs=0
do i=2,nproc
 displs(i)=displs(i-1)+incoming_data_size(i-1)
enddo

! Allocate temporary variable to store parcels data to send
allocate(exchange_data(np_exchange))
allocate(exchange_data_int(np_exchange))

! Gather parcel data to be exchanged
exchange_data(1:np_exchange)=pack(x_p(1:np),p_exchange)
call mpi_allgatherv(exchange_data,incoming_data_size(coord+1),MPI_DOUBLE_PRECISION,exchange_x_p,&
                    incoming_data_size,displs,MPI_DOUBLE_PRECISION,COMM,ierror)

exchange_data(1:np_exchange)=pack(y_p(1:np),p_exchange)
call mpi_allgatherv(exchange_data,incoming_data_size(coord+1),MPI_DOUBLE_PRECISION,exchange_y_p,&
                    incoming_data_size,displs,MPI_DOUBLE_PRECISION,COMM,ierror)

exchange_data(1:np_exchange)=pack(z_p(1:np),p_exchange)
call mpi_allgatherv(exchange_data,incoming_data_size(coord+1),MPI_DOUBLE_PRECISION,exchange_z_p,&
                    incoming_data_size,displs,MPI_DOUBLE_PRECISION,COMM,ierror)

exchange_data(1:np_exchange)=pack(u_p(1:np),p_exchange)
call mpi_allgatherv(exchange_data,incoming_data_size(coord+1),MPI_DOUBLE_PRECISION,exchange_u_p,&
                    incoming_data_size,displs,MPI_DOUBLE_PRECISION,COMM,ierror)

exchange_data(1:np_exchange)=pack(v_p(1:np),p_exchange)
call mpi_allgatherv(exchange_data,incoming_data_size(coord+1),MPI_DOUBLE_PRECISION,exchange_v_p,&
                    incoming_data_size,displs,MPI_DOUBLE_PRECISION,COMM,ierror)

exchange_data(1:np_exchange)=pack(w_p(1:np),p_exchange)
call mpi_allgatherv(exchange_data,incoming_data_size(coord+1),MPI_DOUBLE_PRECISION,exchange_w_p,&
                    incoming_data_size,displs,MPI_DOUBLE_PRECISION,COMM,ierror)

exchange_data(1:np_exchange)=pack(u_sgs(1:np),p_exchange)
call mpi_allgatherv(exchange_data,incoming_data_size(coord+1),MPI_DOUBLE_PRECISION,exchange_u_sgs_p,&
                    incoming_data_size,displs,MPI_DOUBLE_PRECISION,COMM,ierror)

exchange_data(1:np_exchange)=pack(v_sgs(1:np),p_exchange)
call mpi_allgatherv(exchange_data,incoming_data_size(coord+1),MPI_DOUBLE_PRECISION,exchange_v_sgs_p,&
                    incoming_data_size,displs,MPI_DOUBLE_PRECISION,COMM,ierror)

exchange_data(1:np_exchange)=pack(w_sgs(1:np),p_exchange)
call mpi_allgatherv(exchange_data,incoming_data_size(coord+1),MPI_DOUBLE_PRECISION,exchange_w_sgs_p,&
                    incoming_data_size,displs,MPI_DOUBLE_PRECISION,COMM,ierror)

exchange_data(1:np_exchange)=pack(ax_p(1:np),p_exchange)
call mpi_allgatherv(exchange_data,incoming_data_size(coord+1),MPI_DOUBLE_PRECISION,exchange_ax_p,&
                    incoming_data_size,displs,MPI_DOUBLE_PRECISION,COMM,ierror)

exchange_data(1:np_exchange)=pack(ay_p(1:np),p_exchange)
call mpi_allgatherv(exchange_data,incoming_data_size(coord+1),MPI_DOUBLE_PRECISION,exchange_ay_p,&
                    incoming_data_size,displs,MPI_DOUBLE_PRECISION,COMM,ierror)

exchange_data(1:np_exchange)=pack(az_p(1:np),p_exchange)
call mpi_allgatherv(exchange_data,incoming_data_size(coord+1),MPI_DOUBLE_PRECISION,exchange_az_p,&
                    incoming_data_size,displs,MPI_DOUBLE_PRECISION,COMM,ierror)

exchange_data(1:np_exchange)=pack(mass_p(1:np),p_exchange)
call mpi_allgatherv(exchange_data,incoming_data_size(coord+1),MPI_DOUBLE_PRECISION,exchange_mass_p,&
                    incoming_data_size,displs,MPI_DOUBLE_PRECISION,COMM,ierror)

exchange_data(1:np_exchange)=pack(d_p(1:np),p_exchange)
call mpi_allgatherv(exchange_data,incoming_data_size(coord+1),MPI_DOUBLE_PRECISION,exchange_d_p,&
                    incoming_data_size,displs,MPI_DOUBLE_PRECISION,COMM,ierror)

exchange_data(1:np_exchange)=pack(sigma2(1:np),p_exchange)
call mpi_allgatherv(exchange_data,incoming_data_size(coord+1),MPI_DOUBLE_PRECISION,exchange_sigma_p,&
                    incoming_data_size,displs,MPI_DOUBLE_PRECISION,COMM,ierror)

if (temperature) then
   exchange_data(1:np_exchange)=pack(T_p(1:np),p_exchange)
   call mpi_allgatherv(exchange_data,incoming_data_size(coord+1),MPI_DOUBLE_PRECISION,exchange_T_p,&
                       incoming_data_size,displs,MPI_DOUBLE_PRECISION,COMM,ierror)
endif

exchange_data_int(1:np_exchange)=pack(p_ip(1:np),p_exchange)
call mpi_allgatherv(exchange_data_int,incoming_data_size(coord+1),MPI_LONG_LONG_INT,exchange_np_global,&
                    incoming_data_size,displs,MPI_LONG_LONG_INT,COMM,ierror)

deallocate(exchange_data)
deallocate(exchange_data_int)

call mpi_barrier(comm,ierr)

! Remove parcels that do not belong to my rank
np_old = np
np = count(p_exist(1:np_old))
x_p(1:np)    = pack(x_p(1:np_old), p_exist)
y_p(1:np)    = pack(y_p(1:np_old), p_exist)
z_p(1:np)    = pack(z_p(1:np_old), p_exist)

u_p(1:np)    = pack(u_p(1:np_old), p_exist)
v_p(1:np)    = pack(v_p(1:np_old), p_exist)
w_p(1:np)    = pack(w_p(1:np_old), p_exist)

u_sgs(1:np)  = pack(u_sgs(1:np_old),p_exist)
v_sgs(1:np)  = pack(v_sgs(1:np_old),p_exist)
w_sgs(1:np)  = pack(w_sgs(1:np_old),p_exist)

ax_p(1:np) = pack(ax_p(1:np_old),p_exist)
ay_p(1:np) = pack(ay_p(1:np_old),p_exist)
az_p(1:np) = pack(az_p(1:np_old),p_exist)

mass_p(1:np) = pack(mass_p(1:np_old),p_exist)
d_p(1:np)    = pack(d_p(1:np_old),p_exist)
sigma2(1:np) = pack(sigma2(1:np_old),p_exist)

if(temperature) then
T_p(1:np) = pack(T_p(1:np_old),p_exist)
endif

p_ip(1:np) = pack(p_ip(1:np_old),p_exist)
p_exist(1:np) = .true.
p_exchange(1:np) = .false.

!!$OMP END WORKSHARE

x_p(np+1:np_old)    = 0._rprec
y_p(np+1:np_old)    = 0._rprec
z_p(np+1:np_old)    = 0._rprec

u_p(np+1:np_old)    = 0._rprec
v_p(np+1:np_old)    = 0._rprec
w_p(np+1:np_old)    = 0._rprec

u_sgs(np+1:np_old)    = 0._rprec
v_sgs(np+1:np_old)    = 0._rprec
w_sgs(np+1:np_old)    = 0._rprec

ax_p(np+1:np_old) = 0.0_rprec
ay_p(np+1:np_old) = 0.0_rprec
az_p(np+1:np_old) = 0.0_rprec

mass_p(np+1:np_old) = 0.0_rprec
d_p(np+1:np_old) = 0.0_rprec
sigma2(np+1:np_old)   = 0._rprec

if (temperature) then
   T_p(np+1:np_old) = 0._rprec
endif

p_ip(np+1:np_old) = 0
p_exist(np+1:np_old)  = .false.
p_exchange(np+1:np_old) = .false.

! Loop on parcels to be exchanged
!!$omp parallel do private(g_coord_uvp,proc_z), shared(p_exist,p_exchange,x_p,y_p,z_p,u_p,v_p,w_p), &
!!$omp& shared(u_sgs,v_sgs,w_sgs,ax_p,ay_p,az_p,mass_p,d_p,sigma2,T_p,p_ip,np)
do i=1,global_exchange_p

   ! Compute rank where parcel belongs
   g_coord_uvp = (get_eta_local(exchange_z_p(i))/dz) + 0.5_rprec
   proc_z = COUNT( REAL(NZ_CUMUL, RPREC) < g_coord_uvp )

   if (coord .eq. proc_z) then ! Parcel belongs to current rank
      !write(*,*) 'Particle being exchanged !',coord,exchange_x_p(i),exchange_y_p(i),exchange_z_p(i)
      np = np + 1 ! Define new parcel index
      if (np.gt.np_max) then
         print*, "INCREASE NUMBER OF PARCELS (D)"
         stop
      end if

      ! Save new parcel data
      p_exist(np) = .true.
      p_exchange(np) = .false.

      x_p(np) = exchange_x_p(i)
      y_p(np) = exchange_y_p(i)
      z_p(np) = exchange_z_p(i)

      u_p(np) = exchange_u_p(i)
      v_p(np) = exchange_v_p(i)
      w_p(np) = exchange_w_p(i)

      u_sgs(np) = exchange_u_sgs_p(i)
      v_sgs(np) = exchange_v_sgs_p(i)
      w_sgs(np) = exchange_w_sgs_p(i)

      ax_p(np) = exchange_ax_p(i)
      ay_p(np) = exchange_ay_p(i)
      az_p(np) = exchange_az_p(i)

      mass_p(np) = exchange_mass_p(i)
      d_p(np) = exchange_d_p(i)
      sigma2(np) = exchange_sigma_p(i)
      if (temperature) then
         T_p(np) = exchange_T_p(i)
      endif
      p_ip(np) = exchange_np_global(i)

   endif
enddo
!!$omp end parallel do

! Deallocate temporary arrays
deallocate(exchange_x_p)
deallocate(exchange_y_p)
deallocate(exchange_z_p)

deallocate(exchange_u_p)
deallocate(exchange_v_p)
deallocate(exchange_w_p)

deallocate(exchange_u_sgs_p)
deallocate(exchange_v_sgs_p)
deallocate(exchange_w_sgs_p)

deallocate(exchange_ax_p)
deallocate(exchange_ay_p)
deallocate(exchange_az_p)

deallocate(exchange_mass_p)
deallocate(exchange_d_p)
deallocate(exchange_sigma_p)

if (temperature) then
   deallocate(exchange_T_p)
endif

deallocate(exchange_np_global)

call mpi_barrier(comm,ierr)

end subroutine exchange_parcels

!------------------------------------------------------------------------------!
!                                 COMPUTE_TAU                                  !
!------------------------------------------------------------------------------!
subroutine compute_tau()

implicit none

!real(rprec), dimension(6) :: tij_surf
!real(rprec), dimension(3) :: t_surf, tn_surf, tt_surf
!------------------------------------------------------------------------------!

!if (llbc_special.eq.'IBM_DFA') then
!    wall_d = wall_c*dz
!else
!    wall_d = 0.5_rprec*dz ! Normal height from the surface where 1st UV velocities are computed
!end if

if (coord .eq. 0) then ! To be compute by bottom rank only
   ! Compute surface shear stress from logarithmic law
   ! assume horizontal surface here and use horizontally smoothed/filtered friction velocity computed in MO_SIMILARITY.f90
   tau(1:nx,1:ny) = f_rho * USTAR_AVG(1:nx,1:ny)**2.0_rprec
endif

end subroutine compute_tau

!------------------------------------------------------------------------------!
!                              INTERP_VAR_LSM                                  !
!------------------------------------------------------------------------------!
subroutine interp_var_lsm(VAR,IVAR,X1,X2,X3)

IMPLICIT NONE

REAL(RPREC),INTENT(IN) :: X1, X2, X3
REAL(RPREC),DIMENSION(2,2,2) :: VAR
REAL(RPREC),INTENT(OUT) :: IVAR
REAL(RPREC),DIMENSION(2,2,2) :: WG
!------------------------------------------------------------------------------!

!CALCULATE INTERPOLATION WEIGHTS
WG(1,1,1) = (1._RPREC - X1)*(1._RPREC - X2)*(1._RPREC - X3)
WG(2,1,1) = (     X1      )*(1._RPREC - X2)*(1._RPREC - X3)
WG(1,2,1) = (1._RPREC - X1)*(     X2      )*(1._RPREC - X3)
WG(2,2,1) = (     X1      )*(     X2      )*(1._RPREC - X3)
WG(1,1,2) = (1._RPREC - X1)*(1._RPREC - X2)*(     X3      )
WG(2,1,2) = (     X1      )*(1._RPREC - X2)*(     X3      )
WG(1,2,2) = (1._RPREC - X1)*(     X2      )*(     X3      )
WG(2,2,2) = (     X1      )*(     X2      )*(     X3      )

IVAR = SUM(WG * VAR)

end subroutine interp_var_lsm

!------------------------------------------------------------------------------!
!                         UPDATE_SUBDOMAIN_BORDERS                             !
!------------------------------------------------------------------------------!
subroutine update_subdomain_borders(VAR)

implicit none

REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(INOUT) :: VAR
REAL(RPREC),DIMENSION(LD,NY) :: tmp_layer

! Update first pair of grid layers at the border of the subdomains: Do not overwrite ghost layers but
! add them to the corresponding grid layer in the other subdomain.
call mpi_sendrecv(VAR(1,1,1),ld*ny,mpi_rprec,down,500, &
                  tmp_layer(1,1),ld*ny,mpi_rprec,up,500,comm,status,ierr)
if (coord /= (nproc-1)) VAR(:,:,nz) = VAR(:,:,nz) + tmp_layer
call mpi_sendrecv(VAR(1,1,nz),ld*ny,mpi_rprec,up,700, &
                  VAR(1,1,1),ld*ny,mpi_rprec,down,700,comm,status,ierr)
! Update second pair of grid layers
call mpi_sendrecv(VAR(1,1,nz-1),ld*ny,mpi_rprec,up,600, &
                  tmp_layer(1,1),ld*ny,mpi_rprec,down,600,comm,status,ierr)
if (coord /= 0) VAR(:,:,0) = VAR(:,:,0) + tmp_layer
call mpi_sendrecv(VAR(1,1,0),ld*ny,mpi_rprec,down,800, &
                  VAR(1,1,nz-1),ld*ny,mpi_rprec,up,800,comm,status,ierr)

end subroutine update_subdomain_borders

!------------------------------------------------------------------------------!
!                                   LSM_RUN                                    !
!------------------------------------------------------------------------------!
subroutine lsm_run()

implicit none
!------------------------------------------------------------------------------!

! Print parcels data to screen
if ((np .ne. 0) .and. (mod((jt_total-lsm_init),300).eq.0) ) then
    print*, "JT,LSM STEP: ",jt_total,lsm_step,": flying parcels",np," in proc: ",coord,"COLOR:",colour,"RE",maxval(re_p(:))
    write(*,*) 'max u,v,w: ',maxval(abs(u_p(:))),maxval(abs(v_p(:))),maxval(abs(w_p(:))),maxval(abs(ax_p(:))),'in proc:',coord
end if

! SALTATION
if (trim(mode).eq.'saltation') then
    ! Compute surface shear stress
    call compute_tau()
    call mpi_barrier(comm,ierr)
    ! Compute aerodynamic entrainment
    call entrain_parcels()
    call mpi_barrier(comm,ierr)

! POINT RELEASE
else if (trim(mode).eq.'point_release') then
    !call point_release()

! SNOWFALL
else if (trim(mode).eq.'snowfall') then
    call snowfall()

! LINE RELEASE
else if (trim(mode).eq.'line_release') then
    !call line_release()
end if

! Compute surface processes
if (trim(mode) .eq. 'saltation') then
   if (coord .eq. 0) then
      call check_lsm_surface()
   endif
endif

! Remove parcels that got deposited or went out the domain
call remove_parcels()
call mpi_barrier(comm,ierr)

!call compute_sgs()

! Compute force terms to be included in LES equations
if (trim(p_type).eq.'inertial') then
   call compute_forces()
   ! Remove parcels that sublimated completely (diameter below 10 microns)
   call remove_parcels()
   call mpi_barrier(comm,ierr)
endif

! Compute particles advection
call compute_advection()

! Remove parcels that got deposited or went out the domain
call remove_parcels()

! Update LSM time step
lsm_step = lsm_step + 1 !.1000000000_rprec !(lsm_dt/ref_dt)

! write parcel properties at PRINT1_FREQ frequency to file before calling exchang_parcels between 
! ranks to avoid the need for exchanging quantities, e.g. delta_m, that are only used for output purposes at 
! this point because they do not influence the next time step
if ((mod((jt_total-lsm_init),print1_freq).eq.0).and.(np.gt.0)) then
   call print_parcels()
end if
!if ((mod((jt_total-lsm_init),1).eq.0).and.(np.gt.0).and.(mod(lsm_step,nint(1.0_rprec/lsm_freq)).eq.0)) then
!     call print_individual_particles()
!end if
call mpi_barrier(comm,ierr)

! Exchange parcels between ranks
call exchange_parcels()
call mpi_barrier(comm,ierr)

! Write remaining output after calling exchange_parcels because parcels should not be beyond the subdomain of the respective coord
! when extrapolating quantities such as parcel mass
call output_manager()

! Extrapolated particle quantities: Update ghost layers and corresponding layers in the other subdomain
IF (HUMIDITY) THEN
   call update_subdomain_borders(Fq_lsm)
ENDIF
IF (TEMPERATURE) THEN
   call update_subdomain_borders(Ft_lsm)
   call update_subdomain_borders(Ft_lsm_1)
   call update_subdomain_borders(Ft_lsm_2)
ENDIF
call update_subdomain_borders(Fx_lsm)
call update_subdomain_borders(Fy_lsm)
call update_subdomain_borders(Fz_lsm)
call update_subdomain_borders(conc)
call update_subdomain_borders(h_flux)
call update_subdomain_borders(v_flux)
call update_subdomain_borders(c_h_flux)
call update_subdomain_borders(c_v_flux)

end subroutine lsm_run

!------------------------------------------------------------------------------!
end module ls_model
!------------------------------------------------------------------------------!
