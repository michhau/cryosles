!!------------------------------------------------------------------------------
!!  OBJECT: subroutine output_loop ()
!!------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 26/10/2011
!!
!!  DESCRIPTION:
!!
!!  This sub calls all the others STD needed for plotting at specified intervals.
!!
!!
!!------------------------------------------------------------------------------
subroutine output_loop

use std_out
use param
use sim_param, only: U,V,W,P,BUOYANCY,div_vel_vs
use ifport
#IF(SCALARS)
use sgsmodule
#ENDIF
use sp_sc1
use sp_sc2
!use sp_sc3
use ls_model
implicit none

integer :: jx,jy,jz,i
real(rprec) :: cfl,cflOld,viscStab,visc_sc_stab,rmsdivvel,maxrmsdivvel,ke, &
               u_vavg,sc1_vavg,sc1_star_havg,      &
               u_star_havg,u_havg,sc1_havg,phi_m_havg,psi_m_havg,obk_l_havg,sc2_star_havg, &
               sc2_havg, tzsc1_havg, tzsc2_havg !,sc3_star_havg
integer,save :: id1=1,id2=1,id3=1,id4=1
integer,dimension(3) :: slc_id
integer :: nterms,nsampl,result
real(rprec) :: tot_mb,step_mb
character(256) :: fname,opath
real(rprec),dimension(ld,ny,0:nz) :: scr
real(rprec) :: g_eddy_cov, w_eddy_cov
integer :: proc_eddy_cov, jz_eddy_cov, n, file_exists

if(verbose) print*,'in output...'

! init scratch
scr=0.0_rprec

!write some statistics to file
if(mod(jt,wbase).eq.0)then

   call Avg_Ke(ke)
   
   !######### write fluxes at the lower boundary to file  ##########
   if((coord==0).and.(LLBC_SPECIAL=='WALL_LAW'))then
      call Havg_Val(ustar_avg(1:nx,1:ny),u_star_havg,'ustar')
      if(TEMPERATURE) call Havg_Val(sc1_star(1:nx,1:ny),sc1_star_havg,'tstar')
      if(HUMIDITY)    call Havg_Val(sc2_star(1:nx,1:ny),sc2_star_havg,'qstar')
      !if(HUMIDITY) call Havg_Val(sc3_star(1:nx,1:ny),sc3_star_havg,'q_part_star')
   endif
   !######### write SGS fluxes of scalars at the upper boundary to file  ##########
   if(coord==nproc-1) then
      if(TEMPERATURE) call Havg_Val(tzsc1(1:nx,1:ny,nz),tzsc1_havg,'t_sgsFluxUpperBoundary')
      if(HUMIDITY)    call Havg_Val(tzsc2(1:nx,1:ny,nz),tzsc2_havg,'q_sgsFluxUpperBoundary')
   endif
endif

call rmsdiv(rmsdivvel,maxrmsdivvel)

!screen prints
if(mod(jt,wbase).eq.0)then

   !returns the divergence of the vel field
   
   !computes cfl etc. (MPI_REDUCE is inside 
   if(sgs)then
      call timestep_conditions(cfl,cflOld,viscStab,(1/re),nu_t)
   else
      call timestep_conditions(cfl,cflOld,viscStab,(1/re))
   endif
   
   call z_mprofile(u,'u',1.0_rprec,10)
   IF(TEMPERATURE) call z_mprofile(sc1,'t',1.0_rprec,10)
   IF(HUMIDITY) call z_mprofile(sc2,'q',1.0_rprec,10)
   !IF(HUMIDITY) call z_mprofile(sc3,'q_part',1.0_rprec,10)
   IF(TEMPERATURE) call z_mprofile(buoyancy*buoyancy,'buoyancy',1.0_rprec,10)

   CALL MPI_BARRIER(COMM,IERR)

   !extra prints just for coord 0 for now
   if(coord==0)then
      call screen_out(jt_total,dt,tt,ke,rmsdivvel,cfl,cflOld,viscStab)
      print*,'--------'
      call sscreen_out(u_star_havg,'uStarHavg')
      if(TEMPERATURE)then
         call sscreen_out(sc1_star_havg,'tStarHavg')
      endif
      If(HUMIDITY) then
         call sscreen_out(sc2_star_havg,'qStarHavg')
      endif
      !If(HUMIDITY) then
      !   call sscreen_out(sc3_star_havg,'qStarHavg')
      !endif
   endif
   !######### write SGS fluxes of scalars at the upper boundary to screen  ##########
   if(coord==nproc-1) then
      if(TEMPERATURE) call sscreen_out(tzsc1_havg,'t_sgsFluxUpperBoundary')
      if(HUMIDITY)    call sscreen_out(tzsc2_havg,'q_sgsFluxUpperBoundary')
   endif

   IF(KE .GT. 1.0E+8)THEN
      PRINT*,"-----------------------------------------------"
      PRINT*,"--- STOPPED SIMULATION AS SOLUTION DIVERGED ---"
      PRINT*,"-----------:q------------------------------------"
      CALL MPI_FINALIZE(IERR)
   ENDIF

ENDIF

!######################################################################
! OUTPUT TO FILE


   !write inst field
   write(opath,'(a)') trim(out_path)//'ts1_field/'
   if((coord==0).and.(jt==1))then
      print*,'Creating folder for ts1_field...'
      result=system('mkdir '//trim(opath))
   endif
   if(jt==1) call mpi_barrier(comm,ierr)
   if(mod(jt,c_cnt1)==0)then
      if(coord==0) print*,'writing ts1_field...'
      call inst_field(u,'u',1,id1,c_cnt1,opath,1,1,1)
      call inst_field(v,'v',1,id1,c_cnt1,opath,1,1,1)
      call inst_field(w,'w',2,id1,c_cnt1,opath,1,1,1)
      call inst_field(p,'p',1,id1,c_cnt1,opath,1,1,1)
      call inst_field(buoyancy,'buoyancy',1,id1,c_cnt1,opath,1,1,1)
      !call inst_field(dudx,'dudx',1,id1,c_cnt1,opath,1,1,1)
      if(TEMPERATURE) call inst_field(sc1,'t',1,id1,c_cnt1,opath,1,1,1)
      if(HUMIDITY) call inst_field(sc2,'q',1,id1,c_cnt1,opath,1,1,1)
      !if(HUMIDITY) call inst_field(sc3,'q_part',1,id1,c_cnt1,opath,1,1,1)
      !call inst_field(div_vel_vs,'divergence',1,id1,c_cnt1,opath,1,1,1) 
      if(LSM) call inst_field(conc,'conc',1,id1,c_cnt1,opath,1,1,1) 
      !call inst_field(weight_mmatrix,'weight',1,id1,c_cnt1,opath,1,1,1)
      id1=id1+1
   endif

   !slices
   write(opath,'(a)') trim(out_path)//'ts1_slices/'
   if((coord==0).and.(jt==1))then
      print*,'Creating folder for ts1_slices...'
      result=system('mkdir '//trim(opath))
   endif
   if(jt==1) call mpi_barrier(comm,ierr)
   if(mod(jt,c_cnt2)==0)then
      if(coord==0) print*,'writing ts1_slice...'
      call inst_slice(u,'u',c_cnt2,opath,id2)
      call inst_slice(v,'v',c_cnt2,opath,id2)
      call inst_slice(w,'w',c_cnt2,opath,id2)
      if(TEMPERATURE) call inst_slice(sc1,'t',c_cnt2,opath,id2)
      if(HUMIDITY) call inst_slice(sc2,'q',c_cnt2,opath,id2)
      !if(HUMIDITY) call inst_slice(sc3,'q_part',c_cnt2,opath,id2)

      id2=id2+1
   endif

   write(opath,'(a)') trim(out_path)//'inlet_slices/'
   if((coord==0).and.(jt==1))then
      print*,'Creating folder for inlet_slices...'
      result=system('mkdir '//trim(opath))
   endif
   if(jt==1) call mpi_barrier(comm,ierr)
   if(jt_total .gt. inlet_count)then
      if(coord==0) print*,'writing inlet_slice #:',id4
      call inlet_slice(u,'u',c_cnt2,opath,id4)
      call inlet_slice(v,'v',c_cnt2,opath,id4)
      call inlet_slice(w,'w',c_cnt2,opath,id4)
      if(TEMPERATURE) call inlet_slice(sc1,'t',c_cnt2,opath,id4)
      if(HUMIDITY) call inlet_slice(sc2,'q',c_cnt2,opath,id4)
      !if(HUMIDITY) call inlet_slice(sc3,'q_part',c_cnt2,opath,id4)
      id4=id4+1
   endif

!save phi function
!if(jt.eq.nsteps)then
!  if(coord==0) result=system('mkdir '//trim(phi_path))
!    call mpi_barrier(comm,ierr)
!    write(fname,'(a,i0)') trim(phi_path)//'phi.c',coord
!    open(11,file=trim(fname),access='direct',recl=8*nx*ny*(nz-1))
!    write(11,rec=1) phi(1:nx,1:ny,1:nz-1)
!    PRINT*,'Phi written for coord = ', coord
!    close(11)
!endif

!ic for other simulations
if(mod(jt,c_cnt4)==0)then
   call w_restart()
endif

!imitate eddy covariance sampling at a grid point at a given height in the center of the domain
g_eddy_cov    = get_eta_local(z_eddy_cov)/dz ! height (in index units) of eddy covariance sampling above surface
proc_eddy_cov = COUNT( REAL(NZ_CUMUL, RPREC) < g_eddy_cov ) ! MPI process that contains the grid point for eddy covariance sampling
! create directory and text file with header
write(opath,'(a)') trim(out_path)//'eddy_cov/'
if((coord == proc_eddy_cov) .and. (jt == 1))then
  print*,'Creating folder and text file for eddy covariance sampling at z=', &
         z_eddy_cov,' m in the center of the domain...'
  result=system('mkdir '//trim(opath))
  inquire(file=trim(opath)//'uvwTq.txt', EXIST=file_exists)
  if(.not. file_exists) then  ! create file for high-frequency records
    open(300,file=trim(opath)//'uvwTq.txt',status='new')
      if (TEMPERATURE .and. HUMIDITY) then
        write(300,'(a60)') "U(m/s)   V(m/s)   W(m/s)   T(K)   q(kg/kg)"
      elseif (TEMPERATURE) then
        write(300,'(a50)') "U(m/s)   V(m/s)   W(m/s)   T(K)"
      elseif (HUMIDITY) then
        write(300,'(a50)') "U(m/s)   V(m/s)   W(m/s)   q(kg/kg)"
      else
        write(300,'(a40)') "U(m/s)   V(m/s)   W(m/s)"
      endif
    close(300)
  endif
  inquire(file=trim(opath)//'vertical_SGS_scalar_fluxes.txt', EXIST=file_exists)
  if(.not. file_exists) then ! create file for time-averaged SGS fluxes
    if(TEMPERATURE .or. HUMIDITY) then
      open(400,file=trim(opath)//'vertical_SGS_scalar_fluxes.txt',status='new')
      if (TEMPERATURE .and. HUMIDITY) then
        write(400,'(a40)') "Time(s)   ttz(Km/s)   tqz(kgm/kgs)"
      elseif (TEMPERATURE) then
        write(400,'(a30)') "Time(s)   ttz(Km/s)"
      elseif (HUMIDITY) then
        write(400,'(a30)') "Time(s)   tqz(kgm/kgs)"
      endif
      close(400)
    endif
  endif
  ! initialize time-averaged SGS fluxes
  if(TEMPERATURE) attz_eddy_cov = 0.0_rprec
  if(HUMIDITY)    atqz_eddy_cov = 0.0_rprec
endif
! synchronize after folder is typically created such that prints are nicely ordered
if(jt_total==1) call mpi_barrier(comm,ierr)

! amend text files
if (coord == proc_eddy_cov) then
  ! high-frequency file with eddy-covariance time series
  if ( mod(jt, c_cnt5) == 0 ) then
    jz_eddy_cov = idnint( g_eddy_cov + 0.5_rprec ) - eta_displ(coord) ! height index for sampling
    ! average vertical velocities at the two w nodes surrounding the uvp node under consideration
    w_eddy_cov  = ( w(nx/2,ny/2,jz_eddy_cov) + w(nx/2,ny/2,jz_eddy_cov+1) ) / 2
    open(300,file=trim(opath)//'uvwTq.txt',status='old',position='append')
    if (TEMPERATURE .and. HUMIDITY) then
      write(300,'(4(1x,f11.6),(1x,f12.9))') u(nx/2,ny/2,jz_eddy_cov), &
            v(nx/2,ny/2,jz_eddy_cov), &
            w_eddy_cov, sc1(nx/2,ny/2,jz_eddy_cov)*t_scale, &
            sc2(nx/2,ny/2,jz_eddy_cov)
    elseif (TEMPERATURE) then
      write(300,'(4(1x,f11.6))') u(nx/2,ny/2,jz_eddy_cov), w_eddy_cov, &
            v(nx/2,ny/2,jz_eddy_cov), &
            sc1(nx/2,ny/2,jz_eddy_cov)*t_scale
    elseif (HUMIDITY) then
      write(300,'(3(1x,f11.6),(1x,f12.9))') u(nx/2,ny/2,jz_eddy_cov), &
            v(nx/2,ny/2,jz_eddy_cov), &
            w_eddy_cov, sc2(nx/2,ny/2,jz_eddy_cov)
    else
      write(300,'(3(1x,f11.6))') u(nx/2,ny/2,jz_eddy_cov), &
            v(nx/2,ny/2,jz_eddy_cov), w_eddy_cov
    endif
    close(300)
  endif
  ! amend text file with time-averaged SGS fluxes
  if ( mod(jt, c_cnt3) == 0 .and. (TEMPERATURE .or. HUMIDITY) ) then
    jz_eddy_cov = idnint( g_eddy_cov + 0.5 ) - eta_displ(coord) ! height index for sampling
    ! average SGS sensible heat flux at specific grid point
    if(TEMPERATURE) then
      attz_eddy_cov = attz_eddy_cov + ( tzsc1(nx/2,ny/2,jz_eddy_cov) + &
                                        tzsc1(nx/2,ny/2,jz_eddy_cov+1) ) / 2
    endif
    ! average SGS moisture flux at specific grid point
    if(HUMIDITY) then
      atqz_eddy_cov = atqz_eddy_cov + ( tzsc2(nx/2,ny/2,jz_eddy_cov) + &
                                        tzsc2(nx/2,ny/2,jz_eddy_cov+1) ) / 2
    endif
    if ( mod(jt, p_cnt3) == 0 ) then
      n = p_cnt3/c_cnt3
      open(400,file=trim(opath)//'vertical_SGS_scalar_fluxes.txt',status='old',position='append')
      if (TEMPERATURE .and. HUMIDITY) then
        write(400,'(f10.3, 2(1x,e11.4))') tt, attz_eddy_cov * t_scale / n, &
              atqz_eddy_cov / n
      elseif (TEMPERATURE) then
        write(400,'(f10.3, (1x,e11.4))') tt, attz_eddy_cov * t_scale / n
      elseif (HUMIDITY) then
        write(400,'(f10.3, 2(1x,e11.4))') tt, atqz_eddy_cov / n
      endif
      close(400)
      ! reinitialize averages
      if(TEMPERATURE) attz_eddy_cov = 0.0_rprec
      if(HUMIDITY)    atqz_eddy_cov = 0.0_rprec
    endif
  endif
endif


if(verbose) print*,'out output!'


END SUBROUTINE OUTPUT_LOOP

