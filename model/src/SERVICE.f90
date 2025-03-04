MODULE SERVICE

USE IFPORT
USE TYPES
USE PARAM
USE SIM_PARAM,ONLY:U,V,W,RHSx,RHSy,RHSz,DT
USE SP_SC1,ONLY:SC1
USE SP_SC2,ONLY:SC2
!USE SP_SC3,only:sc3
use stretch
use sgsmodule,only:cs_opt2,F_LM,F_MM,F_QN,F_NN
!USE LS_MODEL,only:np_global=>pip_global
!USE LS_MODEL,ONLY:lsm_diss=>diss,lsm_restke=>restke,lsm_sgstke=>sgstke


!------------------------------------------------------------------------------------------
CONTAINS



!!-----------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE ADAPTIVE_DT(U,V,W,DX,DY,DZ,CFL,DT)
!!-----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 18/12/2012
!!
!!  DESCRIPTION:
!!
!!
!!  It actually calculates something that is not the RMS divergence of the velocity is a norm;
!!  is to check the quality of our divergence free field (how much is div free);
!!
!!
!!
!!
!!-----------------------------------------------------------------------------------
SUBROUTINE ADAPTIVE_DT(U,V,W,DX,DY,DZ,CFL,DT)

IMPLICIT NONE

REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(IN) :: U,V,W
REAL(RPREC),INTENT(IN) :: DX,DY,DZ,CFL
REAL(RPREC),INTENT(OUT) :: DT
REAL(RPREC) :: U_MAX,V_MAX,W_MAX

!VALUES
U_MAX = MAXVAL(U(1:NX,1:NY,1:NZ-1))
V_MAX = MAXVAL(V(1:NX,1:NY,1:NZ-1))
W_MAX = MAXVAL(W(1:NX,1:NY,1:NZ))

!COMPUTE DT
DT = MIN(CFL*DX/U_MAX,CFL*DY/V_MAX,CFL*DZ/W_MAX)

END SUBROUTINE ADAPTIVE_DT




!!-----------------------------------------------------------------------------------
!!  OBJECT: subroutine rmsdiv(rms)
!!-----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 08/10/2011
!!
!!  DESCRIPTION:
!!
!!
!!  It actually calculates something that is not the RMS divergence of the velocity is a norm;
!!  is to check the quality of our divergence free field (how much is div free);
!!
!!
!!
!!
!!-----------------------------------------------------------------------------------
subroutine rmsdiv(div,max_div)

use sim_param, only : dudx,dvdy,dwdz,div_vel_vs

implicit none
integer :: jx, jy, jz
real(rprec),intent(out) :: div,max_div

div=0.0_rprec
div_vel_vs=0.0_rprec

do jz=1,nz-1
do jy=1,ny
do jx=1,nx
   div=div+ABS(dudx(jx,jy,jz)+dvdy(jx,jy,jz)+dwdz(jx,jy,jz))
   div_vel_vs(jx,jy,jz) = dudx(jx,jy,jz)+dvdy(jx,jy,jz)+dwdz(jx,jy,jz)
enddo
enddo
enddo
div=div/(nx*ny*(nz-1))

max_div=(maxval(dudx(1:nx,:,1:nz)+dvdy(1:nx,:,1:nz)+dwdz(1:nx,:,1:nz)))

end subroutine rmsdiv





!!------------------------------------------------------------------------------
!!  OBJECT: subroutine Z_MPROFILE(VAR,STR_VAR,NORM)
!!------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 10/10/2011
!!
!!  DESCRIPTION:
!!
!!  useful subroutine to output results for e.g. "DEBUG" mode.
!!  writes on uvp nodes remember that!
!!
!!  STNC = STENCILS
!!
!------------------------------------------------------------------------------
SUBROUTINE Z_MPROFILE(VAR,STR_VAR,NORM,STNC)

IMPLICIT NONE

REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(IN) :: VAR
REAL(RPREC),DIMENSION(1:NZ-1) :: avgVar
REAL(RPREC),DIMENSION(1:NZ_TOT-1) :: avgVarTot
REAL(RPREC),INTENT(IN) :: NORM
INTEGER,INTENT(IN) :: STNC
CHARACTER(*),INTENT(IN) :: STR_VAR
INTEGER :: JZ
REAL(RPREC) :: Z
INTEGER,DIMENSION(NPROC) :: RECVCOUNTS,DISPLS

!BUILD AVERAGE
DO JZ=1,NZ-1
    avgVar(jz)=sum(var(1:nx,1:ny,jz))/(nx*ny)  
ENDDO

RECVCOUNTS = NZ_VECTOR-1
DISPLS = ETA_DISPL !COORD_OF_RANK*RECVCOUNTS
CALL mpi_gatherv (avgVar(1),NZ-1,MPI_RPREC,avgVarTot(1),RECVCOUNTS,DISPLS,MPI_RPREC,RANK_OF_COORD(0),COMM,IERR)
IF(COORD==0)THEN
   DO JZ=1,NZ_TOT-1,STNC
      Z=get_z_local(2*JZ-1)
      WRITE(*,111) STR_VAR,JZ,Z,avgVarTot(JZ),colour
   ENDDO
ENDIF
 
111 format("VAR, JZ, Z, DIM_VALUE:",(A,3x,I6,3x,F20.10,3x,E20.10,I6))

END SUBROUTINE Z_MPROFILE



!!------------------------------------------------------------------------------
!! OBJECT: SUBROUTINE W_RESTART(U,V,W,SC1,SC2,PATH)
!!------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 18/11/2012
!!
!!  DESCRIPTION:
!!
!!  this subroutine writes in an external file the fields u,v,w,sc1,sc2 at the last step
!!  of the simulation for a given domain and discretization.
!!  Its thought in order to provide IC (euler step) to finer simulations
!!
!!  To be improved: probably worth saving also other values for lagrangian models e.g.
!!
!!------------------------------------------------------------------------------
subroutine w_restart()


implicit none

integer :: jx,jy,jz,result
character(256) :: filename,dir_name

if(verbose) print*,'inside w_restart'

write(dir_name,'(a,i0,a)') './restart_files_',colour,'/'

!ascii mode
if(coord==0)then
      print*,'writing restart file...'
   result=system('mkdir '//trim(dir_name))
endif
call mpi_barrier(comm,ierr)

write(filename,'(a,a,i0,a,i0)') trim(dir_name),'restart_',jt_total,'.c',coord
open(111,file=trim(filename),status='replace')

write(111,*) model
write(111,*) l_x
write(111,*) l_y
write(111,*) l_z
write(111,*) nx
write(111,*) ny
write(111,*) nz
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

do jz=1,nz-1
do jy=1,ny
      write(111,*) u(1:nx,jy,jz)
enddo
enddo
do jz=1,nz-1
do jy=1,ny
      write(111,*) v(1:nx,jy,jz)
enddo
enddo
do jz=1,nz-1
do jy=1,ny
      write(111,*) w(1:nx,jy,jz)
enddo
enddo

do jz=1,nz-1
do jy=1,ny
      write(111,*) RHSx(1:nx,jy,jz)
enddo
enddo

do jz=1,nz-1
do jy=1,ny
      write(111,*) RHSy(1:nx,jy,jz)
enddo
enddo

do jz=1,nz-1
do jy=1,ny
      write(111,*) RHSz(1:nx,jy,jz)
enddo
enddo

do jz=1,nz-1
do jy=1,ny
      write(111,*) cs_opt2(1:nx,jy,jz)
enddo
enddo

do jz=1,nz-1
do jy=1,ny
   write(111,*) F_LM(1:nx,jy,jz)
enddo
enddo

do jz=1,nz-1
do jy=1,ny
   write(111,*) F_MM(1:nx,jy,jz)
enddo
enddo

do jz=1,nz-1
do jy=1,ny
   write(111,*) F_NN(1:nx,jy,jz)
enddo
enddo

do jz=1,nz-1
do jy=1,ny
   write(111,*) F_QN(1:nx,jy,jz)
enddo
enddo

if(TEMPERATURE)then
   do jz=1,nz-1
   do jy=1,ny
         write(111,*) sc1(1:nx,jy,jz)
   enddo
   enddo
endif
if(HUMIDITY)then
   do jz=1,nz-1
   do jy=1,ny
         write(111,*) sc2(1:nx,jy,jz)
   enddo
   enddo
endif
!if(HUMIDITY)then
!   do jz=1,nz-1
!   do jy=1,ny
!         write(111,*) sc3(1:nx,jy,jz)
!   enddo
!   enddo
!endif

if(LSM) then
   do jz=1,nz-1
   do jy=1,ny
      write(111,*) diss(1:nx,jy,jz)
   enddo
   enddo
endif 
if(LSM) then
   do jz=1,nz-1
   do jy=1,ny
      write(111,*) restke(1:nx,jy,jz)
   enddo
   enddo
endif 
if(LSM) then
   do jz=1,nz-1
   do jy=1,ny
      write(111,*) sgstke(1:nx,jy,jz)
   enddo
   enddo
endif 


close(111)

!10 format(nx(f24.12))

print*,'leaving w_restart'

if(verbose) print*,'leaving w_restart'

if(coord .eq. 0) then
  open(111,file='./total_time.dat',status='replace')
  write(111,*) jt_total
  close(111)
endif

!if(coord .eq. 0) then
!  open(111,file='./np_global.dat',status='replace')
!  write(111,*) pip_global
!  close(111)
!endif


end subroutine w_restart



!!------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE R_RESTART()
!!------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 13/07/2014
!!
!!  DESCRIPTION:
!!
!!  this subroutine reads from an file the fields u,v,w,sc1,sc2
!!  Its thought in order to provide IC (euler step) to finer simulations
!!
!!  To be improved: probably worth saving also other values for lagrangian models e.g.
!!
!!---------------------------------------------------------------------------
subroutine r_restart()

implicit none

integer :: jx,jy,jz,i
integer :: nx_old,ny_old,nz_old,model_old,T_flag_old,Q_flag_old
real(rprec) :: lx_old,ly_old,lz_old
real(rprec),dimension(:,:,:),allocatable :: u_old,v_old,w_old,sc1_old,sc2_old !,sc3_old
character(256) :: filename
logical :: file_exists,interp

if(verbose) print*,'inside r_restart'

!check for existence
   write(filename,'(a,i0,a,i0,a,i0)') './restart_files_',colour,'/restart_',jt_total,'.c',coord
   write(*,*) filename
   inquire(file=trim(filename), exist=file_exists) 
   if(coord==0)then
      print*,'filename=',filename
      if(file_exists) print*,'External restart file found, going to read from that.'
   endif

!check
if(file_exists)then

   !ASCII MODE
   open(111,file=trim(filename),status='old')

   !read first lines
   read(111,*) model_old
   read(111,*) lx_old
   read(111,*) ly_old
   read(111,*) lz_old
   read(111,*) nx_old  
   read(111,*) ny_old
   read(111,*) nz_old
   read(111,*) T_flag_old
   read(111,*) Q_flag_old
   !checks
   if((nx_old.ne.nx).or.(ny_old.ne.ny))then
      if(coord==0) print*,'attention, nx_old .ne. nx OR ny_old .ne. ny GOING TO INTERPOLATE!'
      interp=.true.
   else
      interp=.false.
   endif 
   if(nz_old.ne.nz)then
      if(coord==0) print*,'attention, nz_old .ne. nz!'
      if(coord==0) print*,'nz_old=',nz_old
      if(coord==0) print*,'nz=',nz
      !call mpi_finalize(ierr)
   endif
   if(lx_old.ne.l_x)then
      if(coord==0) print*,'attention, lx_old .ne. lx!'
      call mpi_finalize(ierr)
   endif
   if(ly_old.ne.l_y)then
      if(coord==0)  print*,'attention, ly_old .ne. ly!'
      call mpi_finalize(ierr)
   endif   
   if(lz_old.ne.l_z)then
      if(coord==0) print*,'attention, lz_old .ne. lz!'
      !call mpi_finalize(ierr)
   endif
   if((TEMPERATURE .and. (T_flag_old==0)).or.((.not. TEMPERATURE).and.(T_flag_old==1)))then
      if(coord==0) print*,'attention, incompatibilities with s_flags between new and old simulations!'
      call mpi_finalize(ierr)
   endif

   if((HUMIDITY .and. (Q_flag_old==0)).or.((.not. HUMIDITY).and.(Q_flag_old==1)))then
      if(coord==0) print*,'attention, incompatibilities with s_flags between new and old simulations!'
      call mpi_finalize(ierr)
   endif

   !allocate
   allocate(u_old(1:nx_old,1:ny_old,1:nz_old-1))
   allocate(v_old(1:nx_old,1:ny_old,1:nz_old-1))
   allocate(w_old(1:nx_old,1:ny_old,1:nz_old-1))
   if(TEMPERATURE) allocate(sc1_old(1:nx_old,1:ny_old,1:nz_old-1))
   if(HUMIDITY) allocate(sc2_old(1:nx_old,1:ny_old,1:nz_old-1))
   !if(HUMIDITY) allocate(sc3_old(1:nx_old,1:ny_old,1:nz_old-1))

   !read matrices
   do jz=1,nz-1
   do jy=1,ny_old
      read(111,*) u_old(1:nx_old,jy,jz)
   enddo
   enddo
   if(coord==0) print*,'read u'

   do jz=1,nz-1
   do jy=1,ny_old
      read(111,*) v_old(1:nx_old,jy,jz)
   enddo
   enddo
   if(coord==0) print*,'read v'

   do jz=1,nz-1
   do jy=1,ny_old
      read(111,*) w_old(1:nx_old,jy,jz)
   enddo
   enddo
   if(coord==0) print*,'read w'

   do jz=1,nz-1
   do jy=1,ny
      read(111,*) RHSx(1:nx,jy,jz)
   enddo
   enddo
   if(coord==0) print*,'read RHS x'

   do jz=1,nz-1
   do jy=1,ny
      read(111,*) RHSy(1:nx,jy,jz)
   enddo
   enddo
   if(coord==0) print*,'read RHS y'

   do jz=1,nz-1
   do jy=1,ny
      read(111,*) RHSz(1:nx,jy,jz)
   enddo
   enddo
   if(coord==0) print*,'read RHS z'

   do jz=1,nz-1
   do jy=1,ny
      read(111,*) cs_opt2(1:nx,jy,jz)
   enddo
   enddo
   if(coord==0) print*,'read CS_OPT2'

   do jz=1,nz-1
   do jy=1,ny
      read(111,*) F_LM(1:nx,jy,jz)
   enddo
   enddo
   if(coord==0) print*,'read F_LM'

   do jz=1,nz-1
   do jy=1,ny
      read(111,*) F_MM(1:nx,jy,jz)
   enddo
   enddo
   if(coord==0) print*,'read F_MM'

   do jz=1,nz-1
   do jy=1,ny
      read(111,*) F_NN(1:nx,jy,jz)
   enddo
   enddo
   if(coord==0) print*,'read F_NN'

   do jz=1,nz-1
   do jy=1,ny
      read(111,*) F_QN(1:nx,jy,jz)
   enddo
   enddo
   if(coord==0) print*,'read F_QN'

   if(TEMPERATURE)then
      do jz=1,nz-1
      do jy=1,ny_old
         read(111,*) sc1_old(1:nx_old,jy,jz)
      enddo
      enddo
      if(coord==0) print*,'read T'
   endif

   if(HUMIDITY)then
      do jz=1,nz-1
      do jy=1,ny_old
         read(111,*) sc2_old(1:nx_old,jy,jz)
      enddo
      enddo
      if(coord==0) print*,'read Q'
   endif

   !if(HUMIDITY)then
   !   do jz=1,nz-1
   !   do jy=1,ny_old
   !      read(111,*) sc3_old(1:nx_old,jy,jz)
   !   enddo
   !   enddo
   !   if(coord==0) print*,'read Q particles'
   !endif

   if(LSM) then
      do jz=1,nz-1
      do jy=1,ny_old
         read(111,*) diss(1:nx_old,jy,jz)
      enddo
      enddo
      if(coord==0) print*,'read diss'
   endif
   if(LSM) then
      do jz=1,nz-1
      do jy=1,ny_old
         read(111,*) restke(1:nx_old,jy,jz)
      enddo
      enddo
      if(coord==0) print*,'read restke'
    endif
   if(LSM) then
      do jz=1,nz-1
      do jy=1,ny_old
         read(111,*) sgstke(1:nx_old,jy,jz)
      enddo
      enddo
      if(coord==0) print*,'read sgstke'
    endif




   !close file
   close(111)

   !initialize new vars
   if(.not.interp)then
        u(1:nx,1:ny,1:nz-1)=u_old
        v(1:nx,1:ny,1:nz-1)=v_old
        w(1:nx,1:ny,1:nz-1)=w_old
        if(TEMPERATURE) sc1(1:nx,1:ny,1:nz-1)=sc1_old
        if(HUMIDITY) sc2(1:nx,1:ny,1:nz-1)=sc2_old
        !if(HUMIDITY) sc3(1:nx,1:ny,1:nz-1)=sc3_old
        if(coord==0) print*,'Initialized new vars from old ones (DIRECT COPY)!'
   elseif(interp)then
        call bilin_interp(u(1:nx,1:ny,1:nz-1),u_old,nx_old,ny_old)
        call bilin_interp(v(1:nx,1:ny,1:nz-1),v_old,nx_old,ny_old)
        call bilin_interp(w(1:nx,1:ny,1:nz-1),w_old,nx_old,ny_old)
        if(TEMPERATURE) call bilin_interp(sc1(1:nx,1:ny,1:nz-1),sc1_old,nx_old,ny_old)
        if(HUMIDITY) call bilin_interp(sc2(1:nx,1:ny,1:nz-1),sc2_old,nx_old,ny_old)
        !if(HUMIDITY) call bilin_interp(sc3(1:nx,1:ny,1:nz-1),sc3_old,nx_old,ny_old)
        if(coord==0) print*,'Initialized new vars from old ones (INTERPOLATION)!'
   endif

endif

write(*,*) "----------------------------------"
write(*,*) "----------------------------------"
write(*,*) "INITIAL_PROFILES"
write(*,*) "----------------------------------"
write(*,*) "----------------------------------"

call z_mprofile(u,'u',1.0_rprec,1)
IF(TEMPERATURE) call z_mprofile(sc1,'t',1.0_rprec,1)
IF(HUMIDITY) call z_mprofile(sc2,'q',1.0_rprec,1)
!IF(HUMIDITY) call z_mprofile(sc3,'q_part',1.0_rprec,1)


if(verbose) print*,'leaving r_restart'

end subroutine r_restart


!!------------------------------------------------------------------------------
!! OBJECT: subroutine Z_MPROFILE_SUBFL(VAR,STR_VAR,NORM)
!!------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 28/08/2012
!!
!!  DESCRIPTION:
!!
!!  useful subroutine to output results for e.g. "DEBUG" mode for the subfluid grid
!!  writes on W nodes remember that (different from Z_MPROFILE that writes on uvp)
!!
!!------------------------------------------------------------------------------
!SUBROUTINE Z_MPROFILE_SUBFL(VAR,STR_VAR,NORM)

!IMPLICIT NONE

!REAL(RPREC),DIMENSION(NX1,NY1,NZ1),INTENT(IN) :: VAR
!REAL(RPREC),INTENT(IN) :: NORM
!CHARACTER(*),INTENT(IN) :: STR_VAR
!INTEGER :: JX,JY,JZ
!REAL(RPREC) :: Z


!DO JZ=1,NZ1
!   Z=(DBLE(JZ_MIN+JZ)-2.0_RPREC)*DZ         !W NODES
!   WRITE(*,111) STR_VAR,JZ,Z,((SUM(VAR(1:NX1,1:NY1,JZ))*NORM)/DBLE(NX1*NY1))
!ENDDO
!111 format("VAR, JZ, Z, DIM_VALUE:",(A,3x,I6,3x,F10.5,3x,E25.15))

!END SUBROUTINE Z_MPROFILE_SUBFL


!MARCO: SUPPORT SUBROUTINE
SUBROUTINE PRINT_PARAMS()

USE SIM_PARAM,ONLY: DT

IMPLICIT NONE


!PRINT SOME STUFF
print *, 'TIME STEPS', NSTEPS
print *, 'DT = ', DT
print *, 'NX,NY = ', NX,NY
print *, 'NZ_VECTOR = ', NZ_VECTOR
print *, 'L_X,L_Y,L_Z = ', L_X,L_Y,L_Z
print *, 'C_COUNT = ', C_CNT1

!ADD SOME MORE STUFF I WOULD SAY


END SUBROUTINE PRINT_PARAMS



!!------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE TRILIN_INTERP()
!!------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 18/11/2012
!!
!!  DESCRIPTION:
!!
!!  This subroutine performs trilinear interpolation between two meshes.
!!  The variables can be staggered then it has to be specified, this works well
!!  just going from coarse to fine assuming there is a 0 and NZ level for the 
!!  staggered mesh.
!!
!!  IMPORTANT:
!!
!!  Remember to pass the correct values at 0 and NZ levels for the staggered nodes
!!  since they're used in the interp (e.g. fine_uvp(jz=1) needs coarse_uvp(jz=0)
!!
!!------------------------------------------------------------------------------
SUBROUTINE TRILIN_INTERP(CVAR,FVAR,CNX,CNY,CNZ,FNX,FNY,FNZ,LX,LY,LZ,STAGGERED)

IMPLICIT NONE

INTEGER,INTENT(IN) :: CNX,CNY,CNZ,FNX,FNY,FNZ
REAL(RPREC),DIMENSION(CNX,CNY,0:CNZ),INTENT(IN):: CVAR
REAL(RPREC),DIMENSION(FNX,FNY,0:FNZ),INTENT(OUT):: FVAR
REAL(RPREC),INTENT(IN) :: LX,LY,LZ
LOGICAL,INTENT(IN) :: STAGGERED
REAL(RPREC) :: CCOORD_X,CCOORD_Y,CCOORD_Z
INTEGER :: JX,JY,JZ,REF_JX,REF_JY,REF_JZ
REAL(RPREC),DIMENSION(4) :: INTERP_X
REAL(RPREC),DIMENSION(2) :: INTERP_Y
REAL(RPREC) :: CDX,CDY,CDZ,FDX,FDY,FDZ,F_SHIFT,C_SHIFT
REAL(RPREC) :: LAGR_X1,LAGR_X2,LAGR_Y1,LAGR_Y2,LAGR_Z1,LAGR_Z2
REAL(RPREC),DIMENSION(:,:,:),ALLOCATABLE :: CVAR_EXT

IF(VERBOSE) PRINT*,"IN TRILIN_INTERP..."


!DEFINE CDX,CDY,CDZ
CDX=LX/CNX
CDY=LY/CNY
CDZ=LZ/(CNZ-1)
FDX=LX/FNX
FDY=LY/FNY
FDZ=LZ/(FNZ-1)

!STAGGERED TRICK
IF(STAGGERED)THEN
   C_SHIFT=0.5_RPREC*CDZ
   F_SHIFT=0.5_RPREC*FDZ
ELSE
   C_SHIFT=0.0_RPREC
   F_SHIFT=0.0_RPREC
ENDIF

!EXTEND AND COPY EVERYTHING ON LOCAL VECTOR (FOR EXTRA NODE TO CLOSE THE DOMAIN)
ALLOCATE(CVAR_EXT(CNX+1,CNY+1,0:CNZ))
DO JZ=0,CNZ
   CVAR_EXT(1:CNX,1:CNY,JZ)=CVAR(1:CNX,1:CNY,JZ)
ENDDO

!PERIODICITY
CVAR_EXT(CNX+1,1:CNY,0:CNZ)=CVAR(1,1:CNY,0:CNZ)
CVAR_EXT(1:CNX,CNY+1,0:CNZ)=CVAR(1:CNX,1,0:CNZ)
CVAR_EXT(CNX+1,CNY+1,0:CNZ)=CVAR(1,1,0:CNZ)

!MAIN LOOP
DO JZ=0,FNZ
DO JY=1,FNY
DO JX=1,FNX

   !COORD X FOR NODE OF FINE MESH
   CCOORD_X=FDX*(JX-1)
   CCOORD_Y=FDY*(JY-1)
   CCOORD_Z=FDZ*(JZ-1)+F_SHIFT
   
   !REFERENCE JX IN COARSE MESH
   REF_JX=FLOOR(CCOORD_X/CDX)+1
   REF_JY=FLOOR(CCOORD_Y/CDY)+1
   REF_JZ=FLOOR((CCOORD_Z-C_SHIFT)/CDZ)+1
   
   !FIX FOR 0 LEVEL
   IF(CCOORD_Z.LT.0.0_RPREC) REF_JZ=0
   
   !COMPUTE WEIGHT FOR LAGRANGIAN PIECEWISE LINEAR INTERPOLATION
   LAGR_X1 =  CCOORD_X/CDX-(REF_JX-1)
   LAGR_Y1 =  CCOORD_Y/CDY-(REF_JY-1)
   LAGR_Z1 = (CCOORD_Z-C_SHIFT)/CDZ-(REF_JZ-1)
   LAGR_X2 = 1.0_RPREC-LAGR_X1
   LAGR_Y2 = 1.0_RPREC-LAGR_Y1
   LAGR_Z2 = 1.0_RPREC-LAGR_Z1
   
   IF(LAGR_X1.LT.0.0_RPREC) PRINT*,"ATTENTION NEGATIVE LAGRANGIAN BASIS FUNCTION!!"
   IF(LAGR_Y1.LT.0.0_RPREC) PRINT*,"ATTENTION NEGATIVE LAGRANGIAN BASIS FUNCTION!!"
   IF(LAGR_Z1.LT.0.0_RPREC) PRINT*,"ATTENTION NEGATIVE LAGRANGIAN BASIS FUNCTION!!"
   
!    IF(CCOORD_Z .GT. 62.0_RPREC)THEN
!       PRINT*,"CCOORD_Z=",CCOORD_Z
!       PRINT*,"C_SHIFT=",C_SHIFT
!   PRINT*,"REF_JZ=",REF_JZ
!       PRINT*,"LAGR_X1=",LAGR_X1
!       PRINT*,"LAGR_Y1=",LAGR_Y1
!       PRINT*,"LAGR_Z1=",LAGR_Z1
!    ENDIF
   
   !ADD ONE CHECK THAT CAN AVOID PROBLEMS WITH END-NODES (SHOULDN'T ADD ANY COMP. PRICE)
   IF((LAGR_X1==0.0_RPREC) .AND.(LAGR_Y1==0.0_RPREC).AND.(LAGR_Z1==0.0_RPREC))THEN
      INTERP_X(1) = CVAR_EXT(REF_JX,REF_JY,REF_JZ)
      INTERP_X(2) = 0.0_RPREC
      INTERP_X(3) = 0.0_RPREC
      INTERP_X(4) = 0.0_RPREC
   ELSEIF((LAGR_X1==0.0_RPREC) .AND.(LAGR_Y1==0.0_RPREC))THEN
      INTERP_X(1) = CVAR_EXT(REF_JX,REF_JY,REF_JZ)
      INTERP_X(2) = 0.0_RPREC
      INTERP_X(3) = CVAR_EXT(REF_JX,REF_JY,REF_JZ+1)
      INTERP_X(4) = 0.0_RPREC
   ELSEIF((LAGR_X1.EQ.0.0_RPREC).AND.(LAGR_Z1.EQ.0.0_RPREC))THEN
      INTERP_X(1) = CVAR_EXT(REF_JX,REF_JY,REF_JZ)
      INTERP_X(2) = CVAR_EXT(REF_JX,REF_JY+1,REF_JZ)
      INTERP_X(3) = 0.0_RPREC
      INTERP_X(4) = 0.0_RPREC
   ELSEIF((LAGR_Y1.EQ.0.0_RPREC).AND.(LAGR_Z1.EQ.0.0_RPREC))THEN
      INTERP_X(1) = CVAR_EXT(REF_JX,REF_JY,REF_JZ)*LAGR_X2   &
                  + CVAR_EXT(REF_JX+1,REF_JY,REF_JZ)*LAGR_X1
      INTERP_X(2) = 0.0_RPREC
      INTERP_X(3) = 0.0_RPREC
      INTERP_X(4) = 0.0_RPREC
   ELSEIF((LAGR_Z1.EQ.0.0_RPREC))THEN
      INTERP_X(1) = CVAR_EXT(REF_JX,REF_JY,REF_JZ)*LAGR_X2               &
                     + CVAR_EXT(REF_JX+1,REF_JY,REF_JZ)*LAGR_X1
      INTERP_X(2) = CVAR_EXT(REF_JX,REF_JY+1,REF_JZ)*LAGR_X2             &
                     + CVAR_EXT(REF_JX+1,REF_JY+1,REF_JZ)*LAGR_X1
      INTERP_X(3) = 0.0_RPREC
      INTERP_X(4) = 0.0_RPREC
   ELSEIF((LAGR_X1.EQ.0.0_RPREC))THEN
      INTERP_X(1) = CVAR_EXT(REF_JX,REF_JY,REF_JZ)
      INTERP_X(2) = CVAR_EXT(REF_JX,REF_JY+1,REF_JZ)
      INTERP_X(3) = CVAR_EXT(REF_JX,REF_JY,REF_JZ+1)
      INTERP_X(4) = CVAR_EXT(REF_JX,REF_JY+1,REF_JZ+1)
   ELSEIF(LAGR_Y1.EQ.0.0_RPREC)THEN
      INTERP_X(1) = CVAR_EXT(REF_JX,REF_JY,REF_JZ)*LAGR_X2               &
                     + CVAR_EXT(REF_JX+1,REF_JY,REF_JZ)*LAGR_X1
      INTERP_X(2) = 0.0_RPREC
      INTERP_X(3) = CVAR_EXT(REF_JX,REF_JY,REF_JZ+1)*LAGR_X2             &
                     + CVAR_EXT(REF_JX+1,REF_JY,REF_JZ+1)*LAGR_X1
      INTERP_X(4) = 0.0_RPREC
   ELSE
      INTERP_X(1) = CVAR_EXT(REF_JX,REF_JY,REF_JZ)*LAGR_X2               &
                     + CVAR_EXT(REF_JX+1,REF_JY,REF_JZ)*LAGR_X1
      INTERP_X(2) = CVAR_EXT(REF_JX,REF_JY+1,REF_JZ)*LAGR_X2             &
                     + CVAR_EXT(REF_JX+1,REF_JY+1,REF_JZ)*LAGR_X1
      INTERP_X(3) = CVAR_EXT(REF_JX,REF_JY,REF_JZ+1)*LAGR_X2             &
                     + CVAR_EXT(REF_JX+1,REF_JY,REF_JZ+1)*LAGR_X1
      INTERP_X(4) = CVAR_EXT(REF_JX,REF_JY+1,REF_JZ+1)*LAGR_X2           &
                     + CVAR_EXT(REF_JX+1,REF_JY+1,REF_JZ+1)*LAGR_X1   
   ENDIF
   
   !Y INTERPOLATION
   INTERP_Y(1) = INTERP_X(1)*LAGR_Y2 + INTERP_X(2)*LAGR_Y1
   INTERP_Y(2) = INTERP_X(3)*LAGR_Y2 + INTERP_X(4)*LAGR_Y1

   !Z INTERPOLATION
   FVAR(JX,JY,JZ) = INTERP_Y(1)*LAGR_Z2 + INTERP_Y(2)*LAGR_Z1
   
!    !CHECK PROBLEM RELATED TO DOUBLE PRECISION (THIS CASE IS SIMPLE CAUSE FEW VALUES)
!    IF((LAGR_X1.EQ.0.0_RPREC) .AND. (LAGR_Y1.EQ.0.0_RPREC))THEN
!       PRINT*,"LAGR_X1=",LAGR_X1
!       PRINT*,"LAGR_Y1=",LAGR_Y1
!       PRINT*,"LAGR_X2=",LAGR_X2
!       PRINT*,"LAGR_Y2=",LAGR_Y2
!       PRINT*,"LAGR_Z1=",LAGR_Z1
!       PRINT*,"LAGR_Z2=",LAGR_Z2
!       PRINT*,"FLUID_VAR(JX,JY,JZ)",FLUID_VAR(JX,JY,JZ)
!       PRINT*,"FLUID_VAR(JX,JY,JZ+1)=",FLUID_VAR(JX,JY,JZ+1)
!       PRINT*,"FLUID_LVAR(JX,JY,JZ)",FLUID_LVAR(JX,JY,JZ)
!       PRINT*,"FLUID_LVAR(JX,JY,JZ+1)=",FLUID_LVAR(JX,JY,JZ+1)
!       PRINT*,"CASE 1 SOLID_VAR(I) = ",SOLID_VAR(I)
!    ENDIF

ENDDO
ENDDO
ENDDO

!DEALLOCATE MEMORY
DEALLOCATE(CVAR_EXT)


IF(VERBOSE) PRINT*, "OUT TRILIN_INTERP"


END SUBROUTINE TRILIN_INTERP



!!------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE BILIN_INTERP()
!!------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 13/10/2014
!!
!!  DESCRIPTION:
!!
!!  This subroutine performs bilinear interpolation between two meshes, it is 
!!  assumed that nz_old=nz and works perfectly with MPI version of the code
!!
!!------------------------------------------------------------------------------
SUBROUTINE BILIN_INTERP(var,varOld,nxOld,nyOld)

IMPLICIT NONE

REAL(RPREC),DIMENSION(nxOld,nyOld,1:nz-1),INTENT(IN):: varOld
REAL(RPREC),DIMENSION(nx,ny,1:nz-1),INTENT(OUT):: var
INTEGER,INTENT(IN) :: nxOld,nyOld
real(rprec) :: dxOld,dyOld
real(rprec) :: cx,cy
integer :: jx,jy,jz,refJx,refJy,refJz
real(rprec),dimension(2) :: interp_x
real(rprec) :: lagr_x1,lagr_x2,lagr_y1,lagr_y2
real(rprec),dimension(:,:),allocatable :: varExt

if(verbose) print*,"Inside bilin_interp..."

!define old stencil
dxOld = l_x/nxOld
dyOld = l_y/nyOld
   allocate(varext(nxold+1,nyold+1))

!interpolate to new grid
do jz=1,nz-1
   
   !extended for periodicity
   varExt(1:nxOld,1:nyOld)=varOld(1:nxOld,1:nyOld,jz)
   varExt(nxOld+1,1:nyOld)=varOld(1,1:nyOld,jz)
   varExt(1:nxOld,nyOld+1)=varOld(1:nxOld,1,jz)
   varExt(nxOld+1,nyOld+1)=varOld(1,1,jz)

   !interpolate
   do jy=1,ny
   do jx=1,nx

      !coords of node and reference old index
      cx=dx*(jx-1)
      cy=dy*(jy-1)
      refJx=floor(cx/dxOld)+1
      refJy=floor(cy/dyOld)+1
   
      !compute weight for lagrangian piecewise linear interpolation
      lagr_x1 =  cx/dxOld-(refJx-1)
      lagr_y1 =  cy/dyOld-(refJy-1)
      lagr_x2 = 1.0_rprec-lagr_x1
      lagr_y2 = 1.0_rprec-lagr_y1
   
      if(lagr_x1.lt.0.0_rprec) print*,"Attention negative lagrangian basis function!"
      if(lagr_y1.lt.0.0_rprec) print*,"Attention negative lagrangian basis function!"
      if(lagr_x1.gt.1.0_rprec) print*,"Attention wrong lagrangian basis function!"
      if(lagr_y1.gt.1.0_rprec) print*,"Attention wrong lagrangian basis function!"

      !interpolate in x (2 interpolations)
      interp_x(1) = varExt(refJx,refJy)*lagr_x2 + varExt(refJx+1,refJy)*lagr_x1
      interp_x(2) = varExt(refJx,refJy+1)*lagr_x2 + varExt(refJx+1,refJy+1)*lagr_x1
 
      !interpolate in y (1 interpolation)
      var(jx,jy,jz) = interp_x(1)*lagr_y2 + interp_x(2)*lagr_y1
      
   enddo
   enddo

enddo

!deallocate memory
deallocate(varExt)


iF(VERBOSE) PRINT*, "OUT BILIN_INTERP"


END SUBROUTINE BILIN_INTERP



!!-------------------------------------------------------------------------------
!!  OBJECT: subroutine READ_PARAM()
!!-------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 09/05/2013
!!
!!  DESCRIPTION:
!!
!!  Reads from external file the parameters for the actualy simulation
!!
!!-------------------------------------------------------------------------------
SUBROUTINE READ_PARAM(FILENAME)

IMPLICIT NONE

CHARACTER(*),INTENT(IN) :: FILENAME
CHARACTER(1) :: H

IF(COORD==0) PRINT*,'READING PARAMETERS...'

OPEN(112,FILE=TRIM(FILENAME),STATUS='OLD')
READ(112,*) 
READ(112,*) H,L_X
READ(112,*) H,L_Y
READ(112,*) H,L_Z
READ(112,*) H,c_factor
READ(112,*) H,f_factor
READ(112,*) H,noise_pressure
READ(112,*) 
READ(112,*) H,NSTEPS
READ(112,*) H,REF_DT
READ(112,*) H,TR_STEPS
READ(112,*) H,TR_DT_FACTOR
READ(112,*) 
READ(112,*) H,RE
READ(112,*) H,U_STAR
READ(112,*) H,RI_SC1
READ(112,*) H,PR_SC1
READ(112,*) H,PR_SGS_SC1
READ(112,*) H,PR_SC2
READ(112,*) H,PR_SGS_SC2
READ(112,*) 
READ(112,*) H,SGS
READ(112,*) H,TEMPERATURE
READ(112,*) H,HUMIDITY
READ(112,*) H,SLOPE
READ(112,*) H,ALPHA
READ(112,*) H,PFX
READ(112,*) H,PFY
READ(112,*)
READ(112,*) H,coriol
READ(112,*) H,ug
READ(112,*) H,vg
READ(112,*) 
READ(112,*) H,MODEL
READ(112,*) H,IFILTER
READ(112,*) H,CS_COUNT
READ(112,*) H,DYN_INIT
READ(112,*) H,CO
READ(112,*) H,WALL_DMP
READ(112,*)
READ(112,*) 
READ(112,*) 
READ(112,*) H,VI1
READ(112,*) H,VI2
READ(112,*) H,VI3
READ(112,*) H,VI4
READ(112,*) H,VNF
READ(112,*) H,U0
READ(112,*) H,V0
READ(112,*) H,W0
READ(112,*) H,DUDZ0
READ(112,*) H,Z_TURB
READ(112,*)
READ(112,*) H,LLBC
READ(112,*) H,LUBC
READ(112,*) H,LLBC_SPECIAL
READ(112,*) H,LUBC_SPECIAL
READ(112,*) H,LLBC_U
READ(112,*) H,LLBC_V
READ(112,*) H,LLBC_W
READ(112,*) H,LUBC_U
READ(112,*) H,LUBC_V
READ(112,*) H,LUBC_W
READ(112,*) H,DLBC_U
READ(112,*) H,DLBC_V
READ(112,*) H,DLBC_W
READ(112,*) H,DUBC_U
READ(112,*) H,DUBC_V
READ(112,*) H,DUBC_W
READ(112,*) H
READ(112,*) H,LBC_DZ0
READ(112,*) H,DSPL_H
READ(112,*) H,WALL_C
READ(112,*) H
READ(112,*) H,TRIINTERFACE
READ(112,*) H,LPS
READ(112,*) H,PHI_ROT
READ(112,*) H,IT_LPS
READ(112,*) H,NITER_F
READ(112,*)
READ(112,*) H,LDSRM
READ(112,*)
READ(112,*)
READ(112,*)
READ(112,*) H,VERBOSE
READ(112,*) H,PROFILING
READ(112,*) H
READ(112,*) H,WBASE
READ(112,*) H,C_CNT1
READ(112,*) H,C_CNT2
READ(112,*) H,P_CNT3
READ(112,*) H,C_CNT3
READ(112,*) H,HAVG_CNT3
READ(112,*) H,C_CNT4
READ(112,*) H,inlet_count
READ(112,*) H,C_CNT5
READ(112,*) H,Z_EDDY_COV
READ(112,*) 
READ(112,*) 
READ(112,*) H,SC1_INIT
READ(112,*) H,DSC1DZ_INIT
READ(112,*) H,SC1STAR_INIT
READ(112,*) H,NF_SC1
READ(112,*) H,Z_TURB_SC1
READ(112,*) H,Z_ENT_SC1
READ(112,*) H,ENT_SC1_INIT
READ(112,*) H,ENT_DSC1DZ_INIT
READ(112,*) H,LLBC_SC1
READ(112,*) H,LUBC_SC1
READ(112,*) H,DLBC_SC1
READ(112,*) H,DUBC_SC1
READ(112,*) 
READ(112,*) H,SC2_INIT
READ(112,*) H,DSC2DZ_INIT
READ(112,*) H,SC2STAR_INIT
READ(112,*) H,NF_SC2
READ(112,*) H,Z_TURB_SC2
READ(112,*) H,Z_ENT_SC2
READ(112,*) H,ENT_SC2_INIT
READ(112,*) H,ENT_DSC2DZ_INIT
READ(112,*) H,LLBC_SC2
READ(112,*) H,LUBC_SC2
READ(112,*) H,DLBC_SC2
READ(112,*) H,DUBC_SC2
READ(112,*)
READ(112,*)
READ(112,*) H,LSM
READ(112,*) H,NP_MAX
READ(112,*) H,LSM_FREQ
READ(112,*) H,AVERAGE_FREQ
READ(112,*) H,PRINT1_FREQ
READ(112,*) H,LSM_INIT
READ(112,*) H,P_RHO
READ(112,*) H,D_M
READ(112,*) H,D_S
READ(112,*) H,D_MIN
READ(112,*) H,D_MAX
READ(112,*) H,PPP_MIN
READ(112,*) H,PPP_MAX
READ(112,*) H,B_ENE
READ(112,*) H,R_DEP
READ(112,*) H,FILE_DEP
READ(112,*) H,INIT_DEP
READ(112,*) H,MODE
READ(112,*) H,OUTLET
READ(112,*) H,P_TYPE
READ(112,*) H,SETTLING_VEL
READ(112,*) H,molec_weight_water
READ(112,*) H,R_universal
READ(112,*) H,surface_tension
READ(112,*) H,water_rho
READ(112,*) H,c_p_a
READ(112,*) H,c_L
READ(112,*) H,c_p_v
READ(112,*) H,sublimation
READ(112,*) H,Schiller
READ(112,*) H,Mason
READ(112,*) H,thorpe_mason
READ(112,*) H,lp_conc
READ(112,*) H,lp_flux
READ(112,*) H,lp_stress
READ(112,*) H,lp_diss
READ(112,*) H,lp_massbal
READ(112,*) H,lp_restart
READ(112,*) H,lp_mombal
READ(112,*) 
READ(112,*) 
READ(112,*) H,num_of_rows
READ(112,*) H,num_of_columns
READ(112,*) H,spacing_x
READ(112,*) H,spacing_y
READ(112,*) H,wind_farm
READ(112,*) H,Ct
READ(112,*) H,Cp
READ(112,*) H,local_tip_speed_ratio
READ(112,*) H,wt_averaging_time
READ(112,*) H,staggered
READ(112,*) H,start_disk
READ(112,*) H,wind_farm_debug
CLOSE(112)

if ((.not. USE_MPI) .or. (USE_MPI .and. coord == 0)) then
!SCREEN
PRINT*,'VERBOSE=',VERBOSE
PRINT*,'PROFILING=',PROFILING
PRINT*,'----------------'
PRINT*,'L_X=',L_X
PRINT*,'L_X=',L_Y
PRINT*,'L_X=',L_Z
PRINT*,'----------------'
PRINT*,'NSTEPS=',NSTEPS
PRINT*,'REF_DT=',REF_DT
PRINT*,'TR_STEPS=',TR_STEPS
PRINT*,'TR_DT_FACTOR=',TR_DT_FACTOR
PRINT*,'----------------'
PRINT*,'WBASE=',WBASE
PRINT*,'----------------'
PRINT*,'RE=',RE
PRINT*,'U_STAR=',U_STAR
PRINT*,'RI_SC1=',RI_SC1
PRINT*,'PR_SC1=',PR_SC1
PRINT*,'PR_SGS_SC1=',PR_SGS_SC1
PRINT*,'----------------'
PRINT*,'SGS=',SGS
PRINT*,'SC1=',TEMPERATURE
PRINT*,'SC2=',HUMIDITY
PRINT*,'SLOPE=',SLOPE
PRINT*,'ALPHA=',ALPHA
PRINT*,'PFX=',PFX
PRINT*,'PFY=',PFY
PRINT*,'----------------'
PRINT*,'MODEL=',MODEL
PRINT*,'IFILTER=',IFILTER
PRINT*,'CS_COUNT=',CS_COUNT
PRINT*,'DYN_INIT=',DYN_INIT
PRINT*,'CO=',CO
PRINT*,'WALL_DMP=',WALL_DMP
PRINT*,'----------------'
PRINT*,'VI1=',VI1
PRINT*,'VI2=',VI2
PRINT*,'VI3=',VI3
PRINT*,'VI4=',VI4
PRINT*,'VNF=',VNF
PRINT*,'U0=',U0
PRINT*,'V0=',V0
PRINT*,'W0=',W0
PRINT*,'DUDZ0=',DUDZ0
PRINT*,'Z_TURB=',Z_TURB
PRINT*,'DTM_FILE=',TRIM(DTM_FILE)
PRINT*,'---------------'
PRINT*,'NITER_F=',NITER_F
PRINT*,'---------------'
PRINT*,'LDSRM=',LDSRM
PRINT*,'SC1_INIT=',SC1_INIT
PRINT*,'DSC1DZ_INIT=',DSC1DZ_INIT
PRINT*,'NF_SC1=',NF_SC1
PRINT*,'Z_TURB_SC1=',Z_TURB_SC1
PRINT*,'Z_ENT_SC1=',Z_ENT_SC1
PRINT*,'ENT_SC1_INIT=',ENT_SC1_INIT
PRINT*,'ENT_DSC1DZ_INIT=',ENT_DSC1DZ_INIT
PRINT*,'LLBC_SC1=',LLBC_SC1
PRINT*,'LUBC_SC1=',LUBC_SC1
PRINT*,'DLBC_SC1=',DLBC_SC1
PRINT*,'DUBC_SC1=',DUBC_SC1
PRINT*,'---------------------'
end if

END SUBROUTINE READ_PARAM



!!--------------------------------------------------------------------------------
!!  OBJECT: subroutine WRITE_CSCALAR(VALUES,INTERP,LABEL,IDNODE,JT,EPFL_PFORMAT)
!!--------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 14/08/2012
!!
!!  DESCRIPTION:
!!
!!  This subroutine simply interpolates from UVP to W
!!  Values are from 0 to NZ for input and output
!!      
!!--------------------------------------------------------------------------------
SUBROUTINE INTERP_UW(UVAR,WVAR)

IMPLICIT NONE

REAL(RPREC),DIMENSION(:,:,0:),INTENT(IN) :: UVAR
REAL(RPREC),DIMENSION(:,:,0:),INTENT(OUT) :: WVAR
INTEGER :: K


!INTERPOLATE
DO K=1,SIZE(UVAR(1,1,:))-1
   WVAR(:,:,K)=0.5_RPREC*(UVAR(:,:,K)+UVAR(:,:,K-1))
ENDDO


END SUBROUTINE INTERP_UW


!SUPPORT SUBROUTINE TO WRITE VARS TO FILE

! S=1 means W variable, S=0.5 means UVP
!SUBROUTINE INTERP_VAR(X,VAR,S,IVAR)
!IMPLICIT NONE
!
!REAL(RPREC),INTENT(IN) :: X(3)
!REAL(RPREC),INTENT(IN) :: S
!REAL(RPREC),DIMENSION(LD,NY,0:NZ) :: VAR
!REAL(RPREC),INTENT(OUT) :: IVAR
!
!INTEGER :: I,J,K
!INTEGER :: I1,J1,K1
!REAL(RPREC) :: X1,X2,X3,SCR
!REAL(RPREC) :: W1,W2,W3,W4,W5,W6,W7,W8
!REAL(RPREC) :: F1,F2,F3,F4,F5,F6,F7,F8

!---------------------------------------------------------------------

!COMPUTE CELL I,J,K
!I=MODULO(FLOOR(X(1)/DX),NX)+1
!J=MODULO(FLOOR(X(2)/DY),NY)+1
!K=FLOOR(X(3)/DZ+S)

!TRY TO HANDLE BOUNDARIES NICELY FOR THE +1 INDICES
!I1=MODULO(I,NX)+1
!J1=MODULO(J,NY)+1
!K1=K+1

!BOUND CHECK
!IF((I<1).OR.(I>NX))THEN
!  PRINT*,'S=',S
!  PRINT*,'X=',X
!  PRINT*,'I OUT OF RANGE IN INTERP_VEL, I=',I
!  PRINT*,'J=',J
!  PRINT*,'K=',K
!  STOP
!ENDIF
!IF((J<1).OR.(J>NY))THEN
!  PRINT*,'S=',S
!  PRINT*,'X=',X
!  PRINT*,'J OUT OF RANGE IN INTERP_VEL, J=',J
!  PRINT*,'I=',I
!  PRINT*,'K=',K
!  STOP
!ENDIF
!IF((K<1).OR.(K>NZ+1))THEN
!  PRINT*,'S=',S
!  PRINT*,'X=',X
!  PRINT*,'K OUT OF RANGE, K =',K
!  PRINT*,'I=',I
!  PRINT*,'J=',J
!  STOP
!ENDIF

!CALCULATE INTERPOLATION WEIGHTS
!X1=MODULO(X(1),DX)/DX
!X2=MODULO(X(2),DY)/DY
!X3=MODULO((X(3)-S*DZ),DZ)/DZ   

!W1=(1._RPREC - X1)*(1._RPREC - X2)*(1._RPREC - X3)
!W2=(    X1    )*(1._RPREC - X2)*(1._RPREC - X3)
!W3=(1._RPREC - X1)*(    X2    )*(1._RPREC - X3)
!W4=(    X1    )*(    X2    )*(1._RPREC - X3)
!W5=(1._RPREC - X1)*(1._RPREC - X2)*(   X3     )
!W6=(    X1    )*(1._RPREC - X2)*(   X3     )
!W7=(1._RPREC - X1)*(    X2    )*(   X3     )
!W8=(    X1    )*(    X2    )*(   X3     )

!F1=VAR(I ,J ,K)
!F2=VAR(I1,J ,K)
!F3=VAR(I ,J1,K)
!F4=VAR(I1,J1,K)
!F5=VAR(I ,J ,K1)
!F6=VAR(I1,J ,K1)
!F7=VAR(I ,J1,K1)
!F8=VAR(I1,J1,K1)
!
!IVAR=W1*F1+W2*F2+W3*F3+W4*F4+W5*F5+W6*F6+W7*F7+W8*F8
!
!
!END SUBROUTINE INTERP_VAR


FUNCTION MAG (input_loc)
IMPLICIT NONE

REAL(RPREC) :: MAG
REAL(RPREC),INTENT(IN) :: input_loc(:)


MAG = SQRT(DOT_PRODUCT(input_loc, input_loc))

END FUNCTION MAG



END MODULE SERVICE
