module tkeBudget1D

use ifport
use types,only:rp=>rprec
use param
use sim_param
use sp_sc1
use sgsmodule,only:beta,cs_opt2,nu_t
use deriv,only:ddz_uv

implicit none

private
character(256) :: oPath
real(rp),dimension(:),allocatable :: au,av,aw,ap,auu,avv,aww,auv,auw,avw,                      &
                             adudx,adudy,adudz,advdx,advdy,advdz,adwdx,adwdy,adwdz,           &
                             adudx2,adudy2,adudz2,advdx2,advdy2,advdz2,adwdx2,adwdy2,adwdz2,  &
                             atxx,atyy,atzz,atxy,atxz,atyz,         &
                             acs,abeta,                                    &
                             auuw,avvw,awww,                        &
                             apw,apu,                                      &
                             adpdx,adpdy,adpdz,                     & 
                             apdudx,apdvdy,apdwdz,                  &
                             autxz,avtyz,awtzz,                     & 
                             adxx,adyy,adzz,                        &
                             adxy,adxz,adyz,                        &
                             apdudz,apdwdx,auww,                    &
                             adudxdwdx,adudydwdy,adudzdwdz
real(rp),dimension(:),allocatable :: as,ass,aus,avs,aws,aps,atsx,atsy,atsz,     &
                                 adsdx,adsdy,adsdz,ads,abetas,apdsdz,asww,                &
                                 adsdxdwdx,adsdydwdy,adsdzdwdz
public computeTkeBudget


!##############################################################################
contains


!##################################################################################

!!------------------------------------------------------------------------------
!!  OBJECT: subroutine alloc_tkeBudget1D()
!!------------------------------------------------------------------------------

subroutine alloc_tkeBudget1D()

implicit none

allocate(au(0:nz)); au=0.0_rp
allocate(av(0:nz)); av=0.0_rp
allocate(aw(0:nz)); aw=0.0_rp
allocate(ap(0:nz)); ap=0.0_rp
allocate(auu(0:nz)); auu=0.0_rp
allocate(avv(0:nz)); avv=0.0_rp
allocate(aww(0:nz)); aww=0.0_rp
allocate(auv(0:nz)); auv=0.0_rp
allocate(auw(0:nz)); auw=0.0_rp
allocate(avw(0:nz)); avw=0.0_rp
allocate(adudx(0:nz)); adudx=0.0_rp
allocate(adudy(0:nz)); adudy=0.0_rp
allocate(adudz(0:nz)); adudz=0.0_rp
allocate(advdx(0:nz)); advdx=0.0_rp
allocate(advdy(0:nz)); advdy=0.0_rp
allocate(advdz(0:nz)); advdz=0.0_rp
allocate(adwdx(0:nz)); adwdx=0.0_rp
allocate(adwdy(0:nz)); adwdy=0.0_rp
allocate(adwdz(0:nz)); adwdz=0.0_rp
allocate(adudx2(0:nz)); adudx2=0.0_rp
allocate(adudy2(0:nz)); adudy2=0.0_rp
allocate(adudz2(0:nz)); adudz2=0.0_rp
allocate(advdx2(0:nz)); advdx2=0.0_rp
allocate(advdy2(0:nz)); advdy2=0.0_rp
allocate(advdz2(0:nz)); advdz2=0.0_rp
allocate(adwdx2(0:nz)); adwdx2=0.0_rp
allocate(adwdy2(0:nz)); adwdy2=0.0_rp
allocate(adwdz2(0:nz)); adwdz2=0.0_rp
allocate(atxx(0:nz)); atxx=0.0_rp
allocate(atyy(0:nz)); atyy=0.0_rp
allocate(atzz(0:nz)); atzz=0.0_rp
allocate(atxy(0:nz)); atxy=0.0_rp
allocate(atxz(0:nz)); atxz=0.0_rp
allocate(atyz(0:nz)); atyz=0.0_rp
allocate(acs(0:nz)); acs=0.0_rp
allocate(abeta(0:nz)); abeta=0.0_rp
allocate(auuw(0:nz)); auuw=0.0_rp
allocate(avvw(0:nz)); avvw=0.0_rp
allocate(awww(0:nz)); awww=0.0_rp
allocate(apw(0:nz)); apw=0.0_rp
allocate(apu(0:nz)); apu=0.0_rp
allocate(adpdx(0:nz)); adpdx=0.0_rp
allocate(adpdy(0:nz)); adpdy=0.0_rp
allocate(adpdz(0:nz)); adpdz=0.0_rp
allocate(apdudx(0:nz)); apdudx=0.0_rp
allocate(apdvdy(0:nz)); apdvdy=0.0_rp
allocate(apdwdz(0:nz)); apdwdz=0.0_rp
allocate(autxz(0:nz)); autxz=0.0_rp
allocate(avtyz(0:nz)); avtyz=0.0_rp
allocate(awtzz(0:nz)); awtzz=0.0_rp
allocate(adxx(0:nz)); adxx=0.0_rp
allocate(adyy(0:nz)); adyy=0.0_rp
allocate(adzz(0:nz)); adzz=0.0_rp
allocate(adxy(0:nz)); adxy=0.0_rp
allocate(adxz(0:nz)); adxz=0.0_rp
allocate(adyz(0:nz)); adyz=0.0_rp
allocate(apdudz(0:nz)); apdudz=0.0_rp
allocate(apdwdx(0:nz)); apdwdx=0.0_rp
allocate(auww(0:nz)); auww=0.0_rp
allocate(adudxdwdx(0:nz)); adudxdwdx=0.0_rp
allocate(adudydwdy(0:nz)); adudydwdy=0.0_rp
allocate(adudzdwdz(0:nz)); adudzdwdz=0.0_rp

allocate(as(0:nz)); as=0.0_rp
allocate(ass(0:nz)); ass=0.0_rp
allocate(aus(0:nz)); aus=0.0_rp
allocate(avs(0:nz)); avs=0.0_rp
allocate(aws(0:nz)); aws=0.0_rp
allocate(aps(0:nz)); aps=0.0_rp
allocate(atsx(0:nz)); atsx=0.0_rp
allocate(atsy(0:nz)); atsy=0.0_rp
allocate(atsz(0:nz)); atsz=0.0_rp
allocate(adsdx(0:nz)); adsdx=0.0_rp
allocate(adsdy(0:nz)); adsdy=0.0_rp
allocate(adsdz(0:nz)); adsdz=0.0_rp
allocate(ads(0:nz)); ads=0.0_rp
allocate(abetas(0:nz)); abetas=0.0_rp
allocate(apdsdz(0:nz)); apdsdz=0.0_rp
allocate(asww(0:nz)); asww=0.0_rp
allocate(adsdxdwdx(0:nz)); adsdxdwdx=0.0_rp
allocate(adsdydwdy(0:nz)); adsdydwdy=0.0_rp
allocate(adsdzdwdz(0:nz)); adsdzdwdz=0.0_rp

end subroutine alloc_tkeBudget1D


!!------------------------------------------------------------------------------
!!  OBJECT: subroutine alloc_tkeBudget1D()
!!------------------------------------------------------------------------------

subroutine dealloc_tkeBudget1D()

implicit none

deallocate(au)
deallocate(av)
deallocate(aw)
deallocate(ap)
deallocate(auu)
deallocate(avv)
deallocate(aww)
deallocate(auv)
deallocate(auw)
deallocate(avw)
deallocate(adudx)
deallocate(adudy)
deallocate(adudz)
deallocate(advdx)
deallocate(advdy)
deallocate(advdz)
deallocate(adwdx)
deallocate(adwdy)
deallocate(adwdz)
deallocate(adudx2)
deallocate(adudy2)
deallocate(adudz2)
deallocate(advdx2)
deallocate(advdy2)
deallocate(advdz2)
deallocate(adwdx2)
deallocate(adwdy2)
deallocate(adwdz2)
deallocate(atxx)
deallocate(atyy)
deallocate(atzz)
deallocate(atxy)
deallocate(atxz)
deallocate(atyz)
deallocate(acs)
deallocate(abeta)
deallocate(auuw)
deallocate(avvw)
deallocate(awww)
deallocate(apw)
deallocate(apu)
deallocate(adpdx)
deallocate(adpdy)
deallocate(adpdz)
deallocate(apdudx)
deallocate(apdvdy)
deallocate(apdwdz)
deallocate(autxz)
deallocate(avtyz)
deallocate(awtzz)
deallocate(adxx)
deallocate(adyy)
deallocate(adzz)
deallocate(adxy)
deallocate(adxz)
deallocate(adyz)
deallocate(apdudz)
deallocate(apdwdx)
deallocate(auww)
deallocate(adudxdwdx)
deallocate(adudydwdy)
deallocate(adudzdwdz)

deallocate(as)
deallocate(ass)
deallocate(aus)
deallocate(avs)
deallocate(aws)
deallocate(aps)
deallocate(atsx)
deallocate(atsy)
deallocate(atsz)
deallocate(adsdx)
deallocate(adsdy)
deallocate(adsdz)
deallocate(ads)
deallocate(abetas)
deallocate(apdsdz)
deallocate(asww)
deallocate(adsdxdwdx)
deallocate(adsdydwdy)
deallocate(adsdzdwdz)

end subroutine dealloc_tkeBudget1D


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
subroutine computeTkeBudget()

implicit none

integer :: i,j,k

real(rp),dimension(ld,ny,0:nz) :: scr1,scr2,scr3,scr4  !scratch useful for interpolation
character(256) :: fname
logical,save :: init=.false.
integer,save :: ii=-10000
integer :: result


if(.not.init)then
   oPath=trim(out_path)//'ta1_profiles/' 

   !initialize record index
   ii=floor(real(jt_total,rprec)/real(p_cnt3,rprec)) + 1
   if((coord==0).and.(mod(dble(nsteps),dble(p_cnt3)).ne.0.0_rp))then
      print*,'Attention, mod(nsteps,p_cnt3).ne.0.0!'
      print*,'nsteps=',nsteps
      print*,'p_cnt3=',p_cnt3
      print*,'c_cnt3=',c_cnt3
      stop
   endif

   !create folder and write parameters file
   if(coord==0)then
   
      !create output folder
      print*,'Creating folder for output...'
      result=system('mkdir '//trim(oPath))
      if(result.eq.0)then
         print*,'Created output folder!'
      else
         print*,'Output already exists!'
         print*,'Result=',Result
      endif
    
      !write/replace parameters file
      write(fname,'(a)') trim(oPath)//'parameters.txt'
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
      if(LLBC_SPECIAL=='IBM_DFA')then
         write(111,*) 1
      else
         write(111,*) 0
      endif
      write(111,*)  phi_rot
      close(111)
   endif   !coord==0
endif      !init

!switch
init=.true.

!mean vars on uvp nodes
call avg_var(u,au,'u',ii)                                                            !vel
call avg_var(v,av,'v',ii)
call int_var(w,scr1); call avg_var(scr1,aw,'w',ii)
if(TEMPERATURE) call avg_var(sc1,as,'s',ii)                                               !vel

call avg_var(dudx,adudx,'dudx',ii)                                                            !vel
call avg_var(dudy,adudy,'dudy',ii)                                                            !vel
call int_var(dudz,scr2); call avg_var(scr2,adudz,'dudz',ii)                                   !vel
call avg_var(dvdx,advdx,'dvdx',ii)                                                            !vel
call avg_var(dvdy,advdy,'dvdy',ii)                                                            !vel
call int_var(dvdz,scr2); call avg_var(scr2,advdz,'dvdz',ii)                                   !vel
call int_var(dwdx,scr2);
call int_var(dwdy,scr3);
call avg_var(scr2,adwdx,'dwdx',ii)                                                            !vel
call avg_var(scr3,adwdy,'dwdy',ii)                                                            !vel
call avg_var(dwdz,adwdz,'dwdz',ii)                                                            !vel
if(TEMPERATURE)then
   call avg_var(dsc1dx,adsdx,'dsdx',ii)                                                            !vel
   call avg_var(dsc1dy,adsdy,'dsdy',ii)                                                            !vel
   call int_var(dsc1dz,scr2);
   call avg_var(scr2,adsdz,'dsdz',ii)                                                            !vel
endif

call avg_var(u*u,auu,'uu',ii)                                                        !res flux
call avg_var(v*v,avv,'vv',ii)
call avg_var(scr1*scr1,aww,'ww',ii)
call avg_var(u*v,auv,'uv',ii)
call avg_var(scr1*u,auw,'uw',ii)
call avg_var(scr1*v,avw,'vw',ii)
if(TEMPERATURE)then
   call avg_var(sc1*sc1,ass,'ss',ii)    
   call avg_var(u*sc1,aus,'su',ii)    
   call avg_var(v*sc1,avs,'sv',ii)
   call avg_var(scr1*sc1,aws,'sw',ii)
endif

call avg_var(txx,atxx,'txx',ii)                     !sgs flux or dns stress
call avg_var(tyy,atyy,'tyy',ii)
call avg_var(tzz,atzz,'tzz',ii)
call avg_var(txy,atxy,'txy',ii)
call int_var(txz,scr1); call avg_var(scr1,atxz,'txz',ii)
call int_var(tyz,scr1); call avg_var(scr1,atyz,'tyz',ii)

if(TEMPERATURE) call avg_var(txsc1,atsx,'tsx',ii)
if(TEMPERATURE) call avg_var(tysc1,atsy,'tsy',ii)
if(TEMPERATURE)then
   call int_var(tzsc1,scr1); call avg_var(scr1,atsz,'tsz',ii)
endif

if(sgs)then
    scr1(:,:,1:nz)=cs_opt2(:,:,:)
    scr1(:,:,0)=scr1(:,:,1)
    call avg_var(scr1,acs,'cs',ii)                     !les stuff not interpolated
    call avg_var(BETA,abeta,'beta',ii)
    if(TEMPERATURE)then
       scr1(:,:,1:nz)=ds_opt2_sc1(:,:,:)
       scr1(:,:,0)=scr1(:,:,1)   !this way I take care of the fact that Cs is defined in UVP in jz=1
       call avg_var(scr1,ads,'ds',ii)                     !les stuff 
       call avg_var(beta_sc1,abetas,'betas',ii)
    endif
endif

   !turb transport
   call int_var(w,scr1); call avg_var(u**2*scr1,auuw,'uuw',ii)
   call avg_var(v**2*scr1,avvw,'vvw',ii)
   call avg_var(scr1**3,awww,'www',ii)
   call avg_var(u*scr1**2,auww,'uww',ii)
   if(TEMPERATURE) call avg_var(sc1*scr1**2,asww,'sww',ii)
   !sgs transp
   call avg_var(scr1*tzz,awtzz,'wtzz',ii)
   call int_var(txz,scr3); 
   if((coord==0).and.(llbc_special.eq.'WALL_LAW')) scr3(:,:,1)=txz(:,:,1)
   call avg_var(scr3*u,autxz,'utxz',ii)
   call int_var(tyz,scr2); 
   if((coord==0).and.(llbc_special.eq.'WALL_LAW')) scr2(:,:,1)=tyz(:,:,1)
   call avg_var(scr2*v,avtyz,'vtyz',ii)

   !dissipation (sgs)
   call avg_var(dudx*txx,adxx,'dxx',ii)
   call avg_var(dvdy*tyy,adyy,'dyy',ii)
   call avg_var(dwdz*tzz,adzz,'dzz',ii)
   call avg_var(0.5_rp*(dudy+dvdx)*txy,adxy,'dxy',ii)
   call int_var(dudz+dwdx,scr4); 
   if((coord==0).and.(llbc_special.eq.'WALL_LAW')) scr4(:,:,1)=dudz(:,:,1)+0.5*(dwdx(:,:,1)+dwdx(:,:,2))
   call avg_var(0.5_rp*scr4*scr3,adxz,'dxz',ii)
   call int_var(dvdz+dwdy,scr1); 
   if((coord==0).and.(llbc_special.eq.'WALL_LAW')) scr1(:,:,1)=dvdz(:,:,1)+0.5*(dwdy(:,:,1)+dwdy(:,:,2))
   call avg_var(0.5*scr2*scr1,adyz,'dyz',ii)

   !dissipation
   call int_var(dudz,scr1); 
   call int_var(dvdz,scr2); 
   call int_var(dwdx,scr3); 
   call int_var(dwdy,scr4); 
   call avg_var(dudx*dudx,adudx2,'dudx2',ii)
   call avg_var(dudy*dudy,adudy2,'dudy2',ii)
   call avg_var(scr1*scr1,adudz2,'dudz2',ii)
   call avg_var(dvdx*dvdx,advdx2,'dvdx2',ii)
   call avg_var(dvdy*dvdy,advdy2,'dvdy2',ii)
   call avg_var(scr2*scr2,advdz2,'dvdz2',ii)
   call avg_var(scr3*scr3,adwdx2,'dwdx2',ii)
   call avg_var(scr4*scr4,adwdy2,'dwdy2',ii)
   call avg_var(dwdz*dwdz,adwdz2,'dwdz2',ii)
 
   !dissipation (UW)
   call int_var(dwdx,scr1)
   call int_var(dwdy,scr2)
   call int_var(dudz,scr3)
   call avg_var(dudx*scr1,adudxdwdx,'dudxdwdx',ii)
   call avg_var(dudy*scr2,adudydwdy,'dudydwdy',ii)
   call avg_var(scr3*dwdz,adudzdwdz,'dudzdwdz',ii)

if(TEMPERATURE)then
   !dissipation (SW)
   call int_var(dwdx,scr1)
   call int_var(dwdy,scr2)
   call int_var(dsc1dz,scr3)
   call avg_var(dsc1dx*scr1,adsdxdwdx,'dsdxdwdx',ii)
   call avg_var(dsc1dy*scr2,adsdydwdy,'dsdydwdy',ii)
   call avg_var(scr3*dwdz,adsdzdwdz,'dsdzdwdz',ii)
endif

!pressure stuff
do k=1,nz-1
  scr1(:,:,k)=p(:,:,k)-0.5_rp*(u(:,:,k)**2+v(:,:,k)**2+(0.5_rp*(w(:,:,k)+w(:,:,k+1)))**2) 
enddo
!derivatives
call filt_da(scr1,scr2,scr3)
call avg_var(scr1,ap,'p',ii)    !cleaned pressure
call avg_var(scr2,adpdx,'dpdx',ii)    !cleaned pressure
call avg_var(scr3,adpdy,'dpdy',ii)    !cleaned pressure
call ddz_uv(scr2,scr1)
call int_var(scr2,scr3); call avg_var(scr3,adpdz,'dpdz',ii)    !cleaned pressure
call int_var(w,scr2)
call avg_var(scr1*scr2,apw,'pw',ii)        !pressure transp tke
call avg_var(scr1*u,apu,'pu',ii)
if(TEMPERATURE) call avg_var(scr1*sc1,aps,'ps',ii)
call avg_var(scr1*dudx,apdudx,'pdudx',ii)
call avg_var(scr1*dvdy,apdvdy,'pdvdy',ii)
call avg_var(scr1*dwdz,apdwdz,'pdwdz',ii)
call int_var(dudz,scr2)
call int_var(dwdx,scr3)
call avg_var(scr1*scr3,apdwdx,'pdwdx',ii)  !pressure transp tke
call avg_var(scr1*scr2,apdudz,'pdudz',ii)  !pressure transp tke
if(TEMPERATURE)then
   call int_var(dsc1dz,scr4)
   call avg_var(scr1*scr4,apdsdz,'pdsdz',ii)  !pressure transp tke
endif
!update record
if(mod(jt,p_cnt3).eq.0)then
   ii=ii+1
endif

end subroutine computeTkeBudget


!!------------------------------------------------------------------------------
!!  OBJECT: subroutine avg_var()
!!------------------------------------------------------------------------------
subroutine avg_var(v,av,str_v,ii)

implicit none

real(rp),dimension(ld,ny,0:nz),intent(in) :: v
real(rp),dimension(0:nz),intent(inout) :: av
character(*),intent(in) :: str_v
integer,intent(in) :: ii
character(128) :: fname
integer :: i,j,k
real(rp) :: n
real(rp),dimension(1:nz_tot-1) :: avTot
integer,dimension(nproc) :: recvcounts,displs

   do k=1,nz-1
       av(k)=av(k)+sum(v(1:nx,1:ny,k))
   enddo

!save if index is reached
if(mod(jt,p_cnt3).eq.0)then
   avTot=1234567.0_rp
   recvcounts = nz_vector-1
   displs = eta_displ
   call mpi_gatherv (av(1),nz-1,mpi_rprec,avTot(1),recvcounts,displs,mpi_rprec,rank_of_coord(0),comm,ierr)
   if(coord==0)then
      n=dble(p_cnt3/c_cnt3)*dble(nx*ny)
      write(fname,'(a,i0)') trim(oPath)//trim(str_v)//'.out'
      open(11,file=trim(fname),access='direct',recl=8*(nz_tot-1))
      write(11,rec=ii) avTot/n
      close(11)
   endif
   av=0.0_rp      !reinitialize averages
endif

end subroutine avg_var


!!------------------------------------------------------------------------------
!!  OBJECT: subroutine int_var()
!!------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 20/01/2014
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
integer :: k

!inteprolate
do k=0,nz-1
   iv(:,:,k)=0.5_rp*(v(:,:,k)+v(:,:,k+1))
enddo

end subroutine int_var


end module tkeBudget1D
