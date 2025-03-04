module avgField

use types,only:rp=>rprec
use param
use sim_param
use sp_sc1
use sp_sc2
!use sp_sc3
use sgsmodule,only:beta,cs_opt2,nu_t
use deriv,only:ddz_uv
use ifport
use ls_model
use stretch

implicit none

!private
character(25) :: o_Path
logical,parameter :: tkeBudget=.false.    
real(rp),dimension(:,:,:),allocatable :: au,av,aw,ap,auu,avv,aww,auv,auw,avw,   &
                                          adudx,adudy,adudz,advdx,advdy,advdz,adwdx,adwdy,adwdz,&
                                          adudx2,adudy2,adudz2,advdx2,advdy2,advdz2,adwdx2,adwdy2,adwdz2,&
                                          atxx,atyy,atzz,atxy,atxz,atyz,         &
                                          acs,abeta,anut,                        &
                                          auuu,avvu,awwu,                        &
                                          auuv,avvv,awwv,                        &
                                          auuw,avvw,awww,                        &
                                          apw,apu,apv,                           &
                                          adpdx,adpdy,adpdz,                     & 
                                          apdudx,apdvdy,apdwdz,                  &
                                          autxx,avtxy,awtxz,                     & 
                                          autxy,avtyy,awtyz,                     & 
                                          autxz,avtyz,awtzz,                     & 
                                          adxx,adyy,adzz,                        &
                                          adxy,adxz,adyz
real(rp),dimension(:,:,:),allocatable :: at,aut,att,avt,awt,attx,atty,attz,     &
                                          adt,abetat
real(rp),dimension(:,:,:),allocatable :: aq,auq,aqq,avq,awq,atqx,atqy,atqz,     &
                                          adq,abetaq
real(rp),dimension(:,:,:),allocatable :: aq_part,auq_part,aqq_part,avq_part,awq_part,atq_partx,atq_party,atq_partz,     &
                                          adq_part,abetaq_part
real(rp),dimension(:,:,:),allocatable :: atfx,atfy,atfz
real(rp),dimension(:,:,:),allocatable :: lsm_conc,lsm_hflux,lsm_vflux
real(rp),dimension(:,:,:),allocatable :: lsm_conc2,lsm_c_hflux,lsm_c_vflux
real(rp),dimension(:,:,:),allocatable :: a_lsm_fx,a_lsm_fy,a_lsm_fz
real(rp),dimension(:,:,:),allocatable :: a_lsm_t_1,a_lsm_t_2,a_lsm_t
real(rp),dimension(:,:,:),allocatable :: a_lsm_q
real(rp),dimension(:,:,:),allocatable :: a_buoy,awbuoy

integer :: ii


!public saveAvgField,initavgfield


!##############################################################################
contains


!##################################################################################



!!------------------------------------------------------------------------------
!!  OBJECT: subroutine ComputeBudget()
!!------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 26/07/2014
!!
!!  DESCRIPTION:
!!
!!  This subroutine extracts the needed statistics for the KAT simulation
!!  all data are interpolated in w nodes
!!
!!
!!  Production is OK
!!  Buoyancy production is also OK
!!  Pressure transport is OK
!!  Transport by molecular diffusion and turbulence is OK
!!  Dissipation is OK
!!  
!!------------------------------------------------------------------------------
subroutine initAvgField()

implicit none

integer :: result
character(256) :: fname

   o_Path=trim(out_path)//'ta1_field/' 

   call allocVar()
   if(TEMPERATURE) call allocTemperatureVar()
   if(HUMIDITY) call allocHumidityVar()
   if(HUMIDITY) call allocHumidity_Part_Var()
   !initialize record index
    ii=floor(real(jt_total,rprec)/real(p_cnt3,rprec))+1
    !if((coord==0).and.(mod(dble(nsteps),dble(p_cnt3)).ne.0.0_rp))then
    !  print*,'Attention, mod(nsteps,p_cnt3).ne.0.0!'
    !  print*,'nsteps=',nsteps
    !  print*,'p_cnt3=',p_cnt3
    !  print*,'c_cnt3=',c_cnt3
    !  stop
    !endif

   !create folder and write parameters file
   if(coord==0)then
   
      !create output folder
      print*,'Creating folder for output...'
      result=system('mkdir '//trim(o_Path))
      if(result.eq.0)then
         print*,'Created output folder!'
      else
         print*,'Output already exists!'
         print*,'Result=',Result
      endif
    
      !write/replace parameters file
      write(fname,'(a)') trim(o_Path)//'parameters.txt'
      open(111,file=fname,status='replace')      !simply replaces in case exists
      print*,'Writing parameters for budget...'
      write(111,*) nx
      write(111,*) ny
      write(111,*) nz_tot-1
      write(111,*) l_x
      write(111,*) l_y
      write(111,*) l_z
      write(111,*) dx
      write(111,*) dy
      write(111,*) dz
      write(111,*) nsteps/p_cnt3
      write(111,*) nproc
      write(111,*) Re
      write(111,*) Ri_sc1
      write(111,*) Pr_sc1
      write(111,*) alpha
      if(SGS)then
         write(111,*) 1
      else
         write(111,*) 0
      endif
      if(TEMPERATURE)then
         write(111,*) 1
      else
         write(111,*) 0
      endif
      if(HUMIDITY)then
         write(111,*) 1
      else
         write(111,*) 0
      endif

      if(LLBC_SPECIAL=='IBM_DFA')then
         write(111,*) 1
      else
         write(111,*) 0
      endif
      write(111,*) phi_rot
      write(111,*) nz_vector
      close(111)
   endif   !coord==0


end subroutine initAvgField


subroutine saveAvgField()

implicit none

integer :: i,j,k

real(rp),dimension(ld,ny,0:nz) :: scr1,scr2,scr3,scr4  !scratch useful for interpolation
character(256) :: fname


!mean vars on uvp nodes
call avg_var(u,au,'u',ii)                                                            !vel
call avg_var(v,av,'v',ii)
!scr1=0.0_rprec ; call int_var(v,scr1); call avg_var(scr1,av,'v',ii)
!scr1=0.0_rprec ; call int_var(u,scr1); call avg_var(scr1,au,'u',ii)
!call avg_var(w,aw,'w',ii)
!call avg_var(scr1*w,auw,'uw',ii)

scr1 = 0.0_rprec ;call int_var(w,scr1); call avg_var(scr1,aw,'w',ii)
if(TEMPERATURE) call avg_var(sc1,at,'t',ii)                                               !vel
if(HUMIDITY) call avg_var(sc2,aq,'q',ii)

if(TEMPERATURE) call avg_var(buoyancy,a_buoy,'buoyancy',ii)

!if(HUMIDITY) call avg_var(sc3,aq,'q_part',ii)

call avg_var(dudx,adudx,'dudx',ii)                                                            !vel
call avg_var(dudy,adudy,'dudy',ii)                                                            !vel
call avg_var(dvdx,advdx,'dvdx',ii)                                                            !vel
call avg_var(dvdy,advdy,'dvdy',ii)                                                            !vel
call avg_var(dwdz,adwdz,'dwdz',ii)                                                            !vel
scr2=0.0_rprec;call int_var(dudz,scr2);call avg_var(scr2,adudz,'dudz',ii)                                                            !vel
 scr2=0.0_rprec;call int_var(dvdz,scr2);call avg_var(scr2,advdz,'dvdz',ii)                                                            !vel
 scr2=0.0_rprec;call int_var(dwdx,scr2);call avg_var(scr2,adwdx,'dwdx',ii)                                                            !vel
 scr2=0.0_rprec;call int_var(dwdy,scr2);call avg_var(scr2,adwdy,'dwdy',ii)                                                            !vel

call avg_var(u*u,auu,'uu',ii)                                                        !res flux
call avg_var(v*v,avv,'vv',ii)
call avg_var(u*v,auv,'uv',ii)
call avg_var(scr1*u,auw,'uw',ii)
call avg_var(scr1*v,avw,'vw',ii)
call avg_var(scr1*scr1,aww,'ww',ii)

if(TEMPERATURE) call avg_var(sc1*sc1,att,'tt',ii)    
if(TEMPERATURE) call avg_var(u*sc1,aut,'tu',ii)    
if(TEMPERATURE) call avg_var(v*sc1,avt,'tv',ii)
if(TEMPERATURE) call avg_var(scr1*sc1,awt,'tw',ii)
if(TEMPERATURE) call avg_var(scr1*buoyancy,awbuoy,'wbuoy',ii)


if(HUMIDITY) call avg_var(sc2*sc2,aqq,'qq',ii)    
if(HUMIDITY) call avg_var(u*sc2,auq,'qu',ii)    
if(HUMIDITY) call avg_var(v*sc2,avq,'qv',ii)
if(HUMIDITY) call avg_var(scr1*sc2,awq,'qw',ii)

!if(HUMIDITY) call avg_var(sc3*sc3,aqq_part,'qq_part',ii)    
!if(HUMIDITY) call avg_var(u*sc3,auq_part,'qu_part',ii)    
!if(HUMIDITY) call avg_var(v*sc3,avq_part,'qv_part',ii)
!if(HUMIDITY) call avg_var(scr1*sc3,awq_part,'qw_part',ii)

call avg_var(txx,atxx,'txx',ii)                     !sgs flux or dns stress
call avg_var(tyy,atyy,'tyy',ii)
call avg_var(tzz,atzz,'tzz',ii)
call avg_var(txy,atxy,'txy',ii)
 scr1=0.0_rprec;call int_var(txz,scr1); call avg_var(scr1,atxz,'txz',ii)
 scr1=0.0_rprec;call int_var(tyz,scr1); call avg_var(scr1,atyz,'tyz',ii)
if(TEMPERATURE) call avg_var(txsc1,attx,'ttx',ii)
if(TEMPERATURE) call avg_var(tysc1,atty,'tty',ii)
if(TEMPERATURE)then
 scr1=0.0_rprec;  call int_var(tzsc1,scr1); call avg_var(scr1,attz,'ttz',ii)
endif

if(HUMIDITY) call avg_var(txsc2,atqx,'tqx',ii)
if(HUMIDITY) call avg_var(tysc2,atqy,'tqy',ii)
if(HUMIDITY)then
 scr1=0.0_rprec;  call int_var(tzsc2,scr1); call avg_var(scr1,atqz,'tqz',ii)
endif

!if(HUMIDITY) call avg_var(txsc3,atq_partx,'tq_partx',ii)
!if(HUMIDITY) call avg_var(tysc3,atq_party,'tq_party',ii)
!if(HUMIDITY)then
!  scr1=0.0_rprec; call int_var(tzsc3,scr1); call avg_var(scr1,atq_partz,'tq_partz',ii)
!endif

if(sgs)then
    scr1(:,:,1:nz)=cs_opt2(:,:,:)
    scr1(:,:,0)=scr1(:,:,1)
    call avg_var(scr1,acs,'cs',ii)                     !les stuff not interpolated
    call avg_var(BETA,abeta,'beta',ii)
    call avg_var(NU_T,anut,'anut',ii)   
    if(TEMPERATURE)then
       scr1(:,:,1:nz)=ds_opt2_sc1(:,:,:)
       scr1(:,:,0)=scr1(:,:,1)   !this way I take care of the fact that Cs is defined in UVP in jz=1
       call avg_var(scr1,adt,'ds_opt2_t',ii)                     !les stuff 
       call avg_var(beta_sc1,abetat,'beta_t',ii)
    endif
   if(HUMIDITY)then
       scr1(:,:,1:nz)=ds_opt2_sc2(:,:,:)
       scr1(:,:,0)=scr1(:,:,1)   !this way I take care of the fact that Cs is defined in UVP in jz=1
       call avg_var(scr1,adq,'ds_opt2_q',ii)                     !les stuff 
       call avg_var(beta_sc2,abetaq,'beta_q',ii)
   endif
   !if(HUMIDITY)then
   !    scr1(:,:,1:nz)=ds_opt2_sc3(:,:,:)
   !    scr1(:,:,0)=scr1(:,:,1)   !this way I take care of the fact that Cs is defined in UVP in jz=1
   !    call avg_var(scr1,adq_part,'ds_opt2_q_part',ii)                     !les stuff 
   !    call avg_var(beta_sc3,abetaq_part,'beta_q_part',ii)
   !endif
endif

scr1=0.0_rprec
if(lsm) then
 !scr1=0.0_rprec ; call int_var(conc,scr1);   
 call avg_var(conc,lsm_conc,'conc',ii)
 call avg_var(conc*conc,lsm_conc2,'conc2',ii)
 scr1=0.0_rprec ; call int_var(h_flux,scr1); call avg_var(scr1,lsm_hflux,'hflux',ii)
 scr1=0.0_rprec ; call int_var(v_flux,scr1); call avg_var(scr1,lsm_vflux,'vflux',ii)
 scr1=0.0_rprec ; call int_var(c_h_flux,scr1); call avg_var(scr1,lsm_c_hflux,'c_hflux',ii)
 scr1=0.0_rprec ; call int_var(c_v_flux,scr1); call avg_var(scr1,lsm_c_vflux,'c_vflux',ii)

 call avg_var(fx_lsm,a_lsm_fx,'force_lsm_x',ii)
 call avg_var(fy_lsm,a_lsm_fy,'force_lsm_y',ii)
 scr1=0.0_rprec; call int_var(fz_lsm,scr1); call avg_var(scr1,a_lsm_fz,'force_lsm_z',ii)

if(temperature) then
 call avg_var(ft_lsm_1,a_lsm_t_1,'lsm_theta_1',ii)
 call avg_var(ft_lsm_2,a_lsm_t_2,'lsm_theta_2',ii)
 call avg_var(ft_lsm,a_lsm_t,'lsm_theta',ii)
endif

if(humidity) then
  call avg_var(fq_lsm,a_lsm_q,'lsm_q',ii)
endif
 
endif


if(tkeBudget)then
   !turb transport
   scr1=0.0_rprec ; call int_var(w,scr1)
   call avg_var(u**2*scr1,auuw,'uuw',ii)
   call avg_var(v**2*scr1,avvw,'vvw',ii)
   call avg_var(scr1**3,awww,'www',ii)
   call avg_var(u**3,auuu,'uuu',ii)
   call avg_var(v**2*u,avvu,'vvu',ii)
   call avg_var(scr1**2*u,awwu,'wwu',ii)
   call avg_var(u**2*v,auuv,'uuv',ii)
   call avg_var(v**3,avvv,'vvv',ii)
   call avg_var(scr1**2*v,awwv,'wwv',ii)
   !sgs transp
   call avg_var(scr1*tzz,awtzz,'wtzz',ii)
   scr2=0.0_rprec ; call int_var(txz,scr2); call avg_var(scr2*u,autxz,'utxz',ii)
   scr3=0.0_rprec ; call int_var(tyz,scr3); call avg_var(scr3*v,avtyz,'vtyz',ii)
   call avg_var(u*txx,autxx,'utxx',ii)
   call avg_var(v*txy,avtxy,'vtxy',ii)
   call avg_var(scr1*scr2,awtxz,'wtxz',ii)
   call avg_var(u*txy,autxy,'utxy',ii)
   call avg_var(v*tyy,avtyy,'vtyy',ii)
   call avg_var(scr1*scr3,awtyz,'wtyz',ii)
   !dissipation
   call avg_var(dudx*txx,adxx,'dxx',ii)
   call avg_var(dvdy*tyy,adyy,'dyy',ii)
   call avg_var(dwdz*tzz,adzz,'dzz',ii)
   call avg_var(0.5_rp*(dudy+dvdx)*txy,adxy,'dxy',ii)
   call int_var(dudz+dwdx,scr1); call avg_var(0.5_rp*scr2*scr1,adxz,'dxz',ii)
   call int_var(dvdz+dwdy,scr1); call avg_var(0.5_rp*scr3*scr1,adyz,'dyz',ii)
   !psdeudo dissipation terms (need the diffusion of tke, but that I can compute after
   scr1=0.0_rprec ; call int_var(dudz,scr1)
   scr2=0.0_rprec ; call int_var(dvdz,scr2)
   scr3=0.0_rprec ; call int_var(dwdx,scr3)
   scr4=0.0_rprec ; call int_var(dwdy,scr4)
   call avg_var(dudx*dudx,adudx2,'dudx2',ii)
   call avg_var(dudy*dudy,adudy2,'dudy2',ii)
   call avg_var(scr1*scr1,adudz2,'dudz2',ii)
   call avg_var(dvdx*dvdx,advdx2,'dvdx2',ii)
   call avg_var(dvdy*dvdy,advdy2,'dvdy2',ii)
   call avg_var(scr2*scr2,advdz2,'dvdz2',ii)
   call avg_var(scr3*scr3,adwdx2,'dwdx2',ii)
   call avg_var(scr4*scr4,adwdy2,'dwdy2',ii)
   call avg_var(dwdz*dwdz,adwdz2,'dwdz2',ii)
endif

scr1=0.0_rprec
scr2=0.0_rprec
scr3=0.0_rprec
!pressure stuff
do k=1,nz-1
  scr1(:,:,k)=p(:,:,k)+0.5_rp*(u(:,:,k)**2+v(:,:,k)**2+(0.5_rp*(w(:,:,k)+w(:,:,k+1)))**2) 
enddo
!derivatives
call filt_da(scr1,scr2,scr3)
call avg_var(scr1,ap,'p',ii)     
call avg_var(scr2,adpdx,'dpdx',ii)
call avg_var(scr3,adpdy,'dpdy',ii)
call ddz_uv(scr2,scr1)
 scr3=0.0_rprec; call int_var(scr2,scr3); call avg_var(scr3,adpdz,'dpdz',ii)    !cleaned pressure
if(tkeBudget)then
   scr2=0.0_rprec ; call int_var(w,scr2)
   call avg_var(scr1*scr2,apw,'pw',ii)        !pressure transp tke
   call avg_var(scr1*u,apu,'pu',ii)           !pressure transp tke
   call avg_var(scr1*v,apv,'pv',ii)           !pressure transp tke
   call avg_var(scr1*dudx,apdudx,'pdudx',ii)  !pressure transp tke
   call avg_var(scr1*dvdy,apdvdy,'pdvdy',ii)  !pressure transp tke
   call avg_var(scr1*dwdz,apdwdz,'pdwdz',ii)  !pressure transp tke
endif

!update record
if(mod(jt,p_cnt3).eq.0)then
   ii=ii+1
endif

if(jt==nsteps) call deAllocVar()
if(TEMPERATURE.and.(jt==nsteps)) call deAllocTemperatureVar()
if(HUMIDITY .and. (jt==nsteps)) call deAllocHumidityVar()
if(HUMIDITY .and. (jt==nsteps)) call deAllocHumidity_Part_Var()
end subroutine saveAvgField


!!------------------------------------------------------------------------------
!!  OBJECT: subroutine avg_var()
!!------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 20/01/2014
!!
!!  DESCRIPTION:
!!
!!  This subroutine updates the value of the current variable
!!
!!------------------------------------------------------------------------------
subroutine avg_var(v_loc,av_loc,str_v,ii)


implicit none

real(rp),dimension(ld,ny,0:nz),intent(in) :: v_loc
real(rp),dimension(ld,ny,0:nz),intent(inout) :: av_loc
character(*),intent(in) :: str_v
integer,intent(in) :: ii
real(rp),dimension(1:nz-1) :: h_av_loc
character(128) :: fname
integer :: n, jz

!counter
n=p_cnt3/c_cnt3

!average
av_loc=av_loc+v_loc

!save if index is reached
if(mod(jt,p_cnt3).eq.0)then
   write(fname,'(a,i0)') trim(o_Path)//trim(str_v)//'.c',coord
   ! Write to file: either horizontal average or whole 3D field
   if (havg_cnt3) then
     do jz=1,nz-1
       h_av_loc(jz) = sum( av_loc(1:nx,1:ny,jz)/(nx*ny) )
     enddo
     open(11,file=trim(fname),status='unknown',access='direct',recl=8*(nz-1))
     write(11,rec=ii) h_av_loc/n
     close(11)
   else
     open(11,file=trim(fname),status='unknown',access='direct',recl=8*nx*ny*(nz-1))
     write(11,rec=ii) av_loc(1:nx,1:ny,1:nz-1)/n
     close(11)
   endif
   av_loc=0.0_rprec      !reinitialize averages
   h_av_loc=0.0_rprec
endif


end subroutine avg_var



!!------------------------------------------------------------------------------
!!  OBJECT: subroutine int_var_2()
!!------------------------------------------------------------------------------
!!
!!  CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 20/01/2014
!!  LAST CHECKED/MODIFIED: Varun Sharma, 2018
!!
!!  DESCRIPTION:
!!
!!  This subroutine interpolates from uvp to w
!!
!!------------------------------------------------------------------------------
subroutine int_var_2(v,iv)

implicit none

real(rp),dimension(ld,ny,0:nz),intent(in) :: v
real(rp),dimension(ld,ny,0:nz),intent(out) :: iv
integer :: k,kk
real(rp) :: fac_1,fac_2,fac_3
integer :: jz_min

if(coord .eq. 0) then
  jz_min = 2
else
  jz_min = 1
endif

!interpolate
do k=jz_min,nz
   kk=k+eta_displ(coord)

   fac_1 = get_z_local(2*kk-1)
   fac_2 = get_z_local(2*kk-3)
   fac_3 = get_z_local(2*kk-2)

   iv(:,:,k)=v(:,:,k-1)+ ((v(:,:,k)-v(:,:,k-1))/(fac_1 - fac_2))*(fac_3 - fac_2)
enddo

end subroutine int_var_2

!!------------------------------------------------------------------------------
!!  OBJECT: subroutine int_var()
!!------------------------------------------------------------------------------
!!
!!  CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 20/01/2014
!!  LAST CHECKED/MODIFIED: Daniela Melo, 25/06/2020
!!
!!  DESCRIPTION:
!!
!!  This subroutine interpolates from w to uvp
!!
!!------------------------------------------------------------------------------
subroutine int_var(v,iv)

implicit none

real(rp),dimension(ld,ny,0:nz),intent(in) :: v
real(rp),dimension(ld,ny,0:nz),intent(out) :: iv
integer :: k,kk
real(rp) :: fac_1,fac_2,fac_3
integer :: jz_min

if(coord .eq. 0) then
  jz_min = 1
else
  jz_min = 0
endif

!interpolate
do k=jz_min,nz-1
   kk=k+eta_displ(coord) ! Global index

   fac_1 = get_z_local(2*kk)
   fac_2 = get_z_local(2*kk-2)
   fac_3 = get_z_local(2*kk-1)

   iv(:,:,k)=v(:,:,k)+ ((v(:,:,k+1)-v(:,:,k))/(fac_1 - fac_2))*(fac_3 - fac_2)
enddo

end subroutine int_var





! SCALAR VARIABLES
subroutine allocTemperatureVar()

implicit none

allocate(at(ld,ny,0:nz)); at=0.0_rp

allocate(a_buoy(ld,ny,0:nz)); a_buoy=0.0_rp
allocate(awbuoy(ld,ny,0:nz)); awbuoy=0.0_rp
allocate(aut(ld,ny,0:nz)); aut=0.0_rp
allocate(att(ld,ny,0:nz)); att=0.0_rp
allocate(avt(ld,ny,0:nz)); avt=0.0_rp
allocate(awt(ld,ny,0:nz)); awt=0.0_rp
allocate(attx(ld,ny,0:nz)); attx=0.0_rp
allocate(atty(ld,ny,0:nz)); atty=0.0_rp
allocate(attz(ld,ny,0:nz)); attz=0.0_rp
allocate(adt(ld,ny,0:nz)); adt=0.0_rp
allocate(abetat(ld,ny,0:nz)); abetat=0.0_rp

end subroutine allocTemperatureVar

subroutine deAllocTemperatureVar()

implicit none

DEallocate(at)
DEallocate(aut)
DEallocate(att)
DEallocate(avt)
DEallocate(awt)
DEallocate(attx)
DEallocate(atty)
DEallocate(attz)
DEallocate(adt)
DEallocate(abetat)
DEALLOCATE(a_buoy)
deallocate(awbuoy)
end subroutine deAllocTemperatureVar

subroutine allocHumidityVar()

implicit none

allocate(aq(ld,ny,0:nz)); aq=0.0_rp
allocate(auq(ld,ny,0:nz)); auq=0.0_rp
allocate(aqq(ld,ny,0:nz)); aqq=0.0_rp
allocate(avq(ld,ny,0:nz)); avq=0.0_rp
allocate(awq(ld,ny,0:nz)); awq=0.0_rp
allocate(atqx(ld,ny,0:nz)); atqx=0.0_rp
allocate(atqy(ld,ny,0:nz)); atqy=0.0_rp
allocate(atqz(ld,ny,0:nz)); atqz=0.0_rp
allocate(adq(ld,ny,0:nz)); adq=0.0_rp
allocate(abetaq(ld,ny,0:nz)); abetaq=0.0_rp

end subroutine allocHumidityVar

subroutine deAllocHumidityVar()

implicit none

DEallocate(aq)
DEallocate(auq)
DEallocate(aqq)
DEallocate(avq)
DEallocate(awq)
DEallocate(atqx)
DEallocate(atqy)
DEallocate(atqz)
DEallocate(adq)
DEallocate(abetaq)

end subroutine deAllocHumidityVar

subroutine allocHumidity_Part_Var()

implicit none

allocate(aq_part(ld,ny,0:nz)); aq_part=0.0_rp
allocate(auq_part(ld,ny,0:nz)); auq_part=0.0_rp
allocate(aqq_part(ld,ny,0:nz)); aqq_part=0.0_rp
allocate(avq_part(ld,ny,0:nz)); avq_part=0.0_rp
allocate(awq_part(ld,ny,0:nz)); awq_part=0.0_rp
allocate(atq_partx(ld,ny,0:nz)); atq_partx=0.0_rp
allocate(atq_party(ld,ny,0:nz)); atq_party=0.0_rp
allocate(atq_partz(ld,ny,0:nz)); atq_partz=0.0_rp
allocate(adq_part(ld,ny,0:nz)); adq_part=0.0_rp
allocate(abetaq_part(ld,ny,0:nz)); abetaq_part=0.0_rp

end subroutine allocHumidity_Part_Var

subroutine deAllocHumidity_Part_Var()

implicit none

DEallocate(aq_part)
DEallocate(auq_part)
DEallocate(aqq_part)
DEallocate(avq_part)
DEallocate(awq_part)
DEallocate(atq_partx)
DEallocate(atq_party)
DEallocate(atq_partz)
DEallocate(adq_part)
DEallocate(abetaq_part)

end subroutine deAllocHumidity_Part_Var




! VARIABLES
subroutine allocVar()

implicit none

allocate(au(ld,ny,0:nz)); au=0.0_rp
allocate(av(ld,ny,0:nz)); av=0.0_rp
allocate(aw(ld,ny,0:nz)); aw=0.0_rp
allocate(ap(ld,ny,0:nz)); ap=0.0_rp
allocate(adudx(ld,ny,0:nz)); adudx=0.0_rp
allocate(adudy(ld,ny,0:nz)); adudy=0.0_rp
allocate(adudz(ld,ny,0:nz)); adudz=0.0_rp
allocate(advdx(ld,ny,0:nz)); advdx=0.0_rp
allocate(advdy(ld,ny,0:nz)); advdy=0.0_rp
allocate(advdz(ld,ny,0:nz)); advdz=0.0_rp
allocate(adwdx(ld,ny,0:nz)); adwdx=0.0_rp
allocate(adwdy(ld,ny,0:nz)); adwdy=0.0_rp
allocate(adwdz(ld,ny,0:nz)); adwdz=0.0_rp
allocate(adpdx(ld,ny,0:nz)); adpdx=0.0_rp
allocate(adpdy(ld,ny,0:nz)); adpdy=0.0_rp
allocate(adpdz(ld,ny,0:nz)); adpdz=0.0_rp
allocate(auu(ld,ny,0:nz)); auu=0.0_rp
allocate(avv(ld,ny,0:nz)); avv=0.0_rp
allocate(aww(ld,ny,0:nz)); aww=0.0_rp
allocate(auv(ld,ny,0:nz)); auv=0.0_rp
allocate(auw(ld,ny,0:nz));  auw=0.0_rp
allocate(avw(ld,ny,0:nz)); avw=0.0_rp
allocate(atxx(ld,ny,0:nz)); atxx=0.0_rp
allocate(atyy(ld,ny,0:nz)); atyy=0.0_rp
allocate(atzz(ld,ny,0:nz)); atzz=0.0_rp
allocate(atxy(ld,ny,0:nz)); atxy=0.0_rp
allocate(atxz(ld,ny,0:nz)); atxz=0.0_rp
allocate(atyz(ld,ny,0:nz)); atyz=0.0_rp
if(sgs) allocate(acs(ld,ny,0:nz)); acs=0.0_rp
if(sgs) allocate(abeta(ld,ny,0:nz)); abeta=0.0_rp
if(sgs) allocate(anut(ld,ny,0:nz)); anut=0.0_rp

if(lsm) then
allocate(lsm_conc(ld,ny,0:nz))
lsm_conc = 0.0_rprec

allocate(lsm_conc2(ld,ny,0:nz))
lsm_conc2 = 0.0_rprec
allocate(lsm_hflux(ld,ny,0:nz))
 lsm_hflux = 0.0_rprec
allocate(lsm_vflux(ld,ny,0:nz))
 lsm_vflux = 0.0_rprec
allocate(lsm_c_hflux(ld,ny,0:nz))
 lsm_c_hflux = 0.0_rprec
allocate(lsm_c_vflux(ld,ny,0:nz))
 lsm_c_vflux = 0.0_rprec

endif

if(lsm) then
   allocate(a_lsm_fx(ld,ny,0:nz)); a_lsm_fx = 0.0_rprec
   allocate(a_lsm_fy(ld,ny,0:nz)); a_lsm_fy = 0.0_rprec
   allocate(a_lsm_fz(ld,ny,0:nz)); a_lsm_fz = 0.0_rprec
   if(temperature) then
       allocate(a_lsm_t_1(ld,ny,0:nz)); a_lsm_t_1 = 0.0_rprec
       allocate(a_lsm_t_2(ld,ny,0:nz)); a_lsm_t_2 = 0.0_rprec
       allocate(a_lsm_t(ld,ny,0:nz)); a_lsm_t = 0.0_rprec
   endif
   if(humidity) then
       allocate(a_lsm_q(ld,ny,0:nz)); a_lsm_q = 0.0_rprec
   endif
endif 
 
if(tkeBudget)then
   allocate(auuw(ld,ny,0:nz)); auuw=0.0_rp
   allocate(avvw(ld,ny,0:nz)); avvw=0.0_rp
   allocate(awww(ld,ny,0:nz)); awww=0.0_rp
   allocate(auuu(ld,ny,0:nz)); auuu=0.0_rp
   allocate(avvu(ld,ny,0:nz)); avvu=0.0_rp
   allocate(awwu(ld,ny,0:nz)); awwu=0.0_rp
   allocate(auuv(ld,ny,0:nz)); auuv=0.0_rp
   allocate(avvv(ld,ny,0:nz)); avvv=0.0_rp
   allocate(awwv(ld,ny,0:nz)); awwv=0.0_rp
   allocate(apw(ld,ny,0:nz)); apw=0.0_rp
   allocate(apu(ld,ny,0:nz)); apu=0.0_rp
   allocate(apv(ld,ny,0:nz)); apv=0.0_rp
   allocate(autxz(ld,ny,0:nz)); autxz=0.0_rp
   allocate(avtyz(ld,ny,0:nz)); avtyz=0.0_rp
   allocate(awtzz(ld,ny,0:nz)); awtzz=0.0_rp
   allocate(autxx(ld,ny,0:nz)); autxx=0.0_rp
   allocate(avtxy(ld,ny,0:nz)); avtxy=0.0_rp
   allocate(awtxz(ld,ny,0:nz)); awtxz=0.0_rp
   allocate(autxy(ld,ny,0:nz)); autxy=0.0_rp
   allocate(avtyy(ld,ny,0:nz)); avtyy=0.0_rp
   allocate(awtyz(ld,ny,0:nz)); awtyz=0.0_rp
   allocate(adxx(ld,ny,0:nz)); adxx=0.0_rp
   allocate(adyy(ld,ny,0:nz)); adyy=0.0_rp
   allocate(adzz(ld,ny,0:nz)); adzz=0.0_rp
   allocate(adxy(ld,ny,0:nz)); adxy=0.0_rp
   allocate(adxz(ld,ny,0:nz)); adxz=0.0_rp
   allocate(adyz(ld,ny,0:nz)); adyz=0.0_rp
   allocate(apdudx(ld,ny,0:nz)); apdudx=0.0_rp
   allocate(apdvdy(ld,ny,0:nz)); apdvdy=0.0_rp
   allocate(apdwdz(ld,ny,0:nz)); apdwdz=0.0_rp
allocate(adudx2(ld,ny,0:nz)); adudx2=0.0_rp
allocate(adudy2(ld,ny,0:nz)); adudy2=0.0_rp
allocate(adudz2(ld,ny,0:nz)); adudz2=0.0_rp
allocate(advdx2(ld,ny,0:nz)); advdx2=0.0_rp
allocate(advdy2(ld,ny,0:nz)); advdy2=0.0_rp
allocate(advdz2(ld,ny,0:nz)); advdz2=0.0_rp
allocate(adwdx2(ld,ny,0:nz)); adwdx2=0.0_rp
allocate(adwdy2(ld,ny,0:nz)); adwdy2=0.0_rp
allocate(adwdz2(ld,ny,0:nz)); adwdz2=0.0_rp
endif 


end subroutine allocVar

subroutine deAllocVar()

implicit none

DEallocate(au)
DEallocate(av)
DEallocate(aw)
DEallocate(ap)
DEallocate(adudx)
DEallocate(adudy)
DEallocate(adudz)
DEallocate(advdx)
DEallocate(advdy)
DEallocate(advdz)
DEallocate(adwdx)
DEallocate(adwdy)
DEallocate(adwdz)
DEallocate(adpdx)
DEallocate(adpdy)
DEallocate(adpdz)
DEallocate(auu)
DEallocate(avv)
DEallocate(aww)
DEallocate(auv)
DEallocate(auw)
DEallocate(avw)
DEallocate(atxx)
DEallocate(atyy)
DEallocate(atzz)
DEallocate(atxy)
DEallocate(atxz)
DEallocate(atyz)
if(sgs) DEallocate(acs)
if(sgs) DEallocate(abeta)
if(sgs) DEallocate(anut)

if(lsm) deallocate(lsm_conc)
if(lsm) deallocate(lsm_conc2)
if(lsm) deallocate(lsm_hflux)
if(lsm) deallocate(lsm_vflux)
if(lsm) deallocate(lsm_c_hflux)
if(lsm) deallocate(lsm_c_vflux)

if(lsm) then
 deallocate(a_lsm_fx)
 deallocate(a_lsm_fy)
 deallocate(a_lsm_fz)
 IF(temperature) then
    deallocate(a_lsm_t_1)
    deallocate(a_lsm_t_2)
    deallocate(a_lsm_t)
 ENDIF
 IF(humidity) then
    deallocate(a_lsm_q)
 ENDIF
endif

if(tkeBudget)then
   DEallocate(auuw)
   DEallocate(avvw)
   DEallocate(awww)
   DEallocate(auuu)
   DEallocate(avvu)
   DEallocate(awwu)
   DEallocate(auuv)
   DEallocate(avvv)
   DEallocate(awwv)
   DEallocate(apw)
   DEallocate(apv)
   DEallocate(apu)
   DEallocate(autxz)
   DEallocate(avtyz)
   DEallocate(awtzz)
   DEallocate(autxx)
   DEallocate(avtxy)
   DEallocate(awtxz)
   DEallocate(autxy)
   DEallocate(avtyy)
   DEallocate(awtyz)
   DEallocate(adxx)
   DEallocate(adyy)
   DEallocate(adzz)
   DEallocate(adxy)
   DEallocate(adxz)
   DEallocate(adyz)
   DEallocate(apdudx)
   DEallocate(apdvdy)
   DEallocate(apdwdz)
DEallocate(adudx2)
DEallocate(adudy2)
DEallocate(adudz2)
DEallocate(advdx2)
DEallocate(advdy2)
DEallocate(advdz2)
DEallocate(adwdx2)
DEallocate(adwdy2)
DEallocate(adwdz2)
endif


end subroutine deAllocVar


end module avgField
