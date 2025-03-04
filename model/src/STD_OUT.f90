MODULE STD_OUT

USE TYPES,ONLY:RPREC
USE PARAM
USE SERVICE 
USE DERIV
USE SIM_PARAM
USE SP_SC1
USE SGSMODULE,ONLY: CS_OPT2,NU_T

IMPLICIT NONE



!-------------------------------------------------------------------------------
CONTAINS !----------------------------------------------------------------------
!-------------------------------------------------------------------------------


!!FUNCTION! TO CREATE 1D VECTORS FROM 3D MATRICES
!FUNCTION RSHP(V1,INTP)
!IMPLICIT NONE

!REAL,DIMENSION(NX*NY*NZ) :: RSHP
!REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(IN) :: V1
!INTEGER,INTENT(IN) :: INTP
!INTEGER :: I,J,K,II
!
!!BUILD VAR 1D
!II=1
!IF(INTP.EQ.1)THEN
!   DO K=1,NZ
!   DO J=1,NY
!   DO I=1,NX
!      RSHP(II)=REAL(0.5_RPREC*(V1(I,J,K)+V1(I,J,K-1)))
!      II=II+1
!   ENDDO
!   ENDDO
!   ENDDO
!ELSEIF(INTP.EQ.0)THEN
!   DO K=1,NZ
!   DO J=1,NY
!   DO I=1,NX
!      RSHP(II)=REAL(V1(I,J,K))
!      II=II+1
!   ENDDO
!   ENDDO
!   ENDDO
!ELSE
!   PRINT*,'PLEASE PROVIDE EITHER 1 OR 0 FOR INTERPOLATION'
!   STOP
!ENDIF
!
!END FUNCTION RSHP


! !!--------------------------------------------------------------------------------
! !!  OBJECT: subroutine GID_OUT()
! !!--------------------------------------------------------------------------------
! !!
! !!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 08/11/2012
! !!
! !!  DESCRIPTION:
! !!
! !!  HL subroutine for GiD, makes use of GID_POST
! !!
! !!  CNX,CNY,CNZ  -->   COARSE NX,NY,NZ
! !!  CDX,CDY,CDZ  -->   COARSE DX,DY,DZ
! !!  
! !!--------------------------------------------------------------------------------
! SUBROUTINE GID_OUT(STEP)
! 
! USE GID_POST
! 
! IMPLICIT NONE
! 
! INTEGER,INTENT(IN) :: STEP
! INTEGER,SAVE :: FIRST_ELEM=1,FIRST_NODE=1
! INTEGER :: I,JX,JY,JZ
! INTEGER, DIMENSION(((NX-1)/GNX+1)*((NY-1)/GNY+1)*((NZ-1)/GNZ+1)) :: FLUID_NODE_ID
! INTEGER :: OUTPUT_SIZE
! CHARACTER(256),SAVE :: FILE
! INTEGER :: CNX,CNY,CNZ
! REAL(RPREC) :: CDX,CDY,CDZ
! REAL(RPREC),DIMENSION(LD,NY,0:NZ) :: SCR1,SCR2
! LOGICAL,SAVE :: INIT=.TRUE.
! 
! IF(VERBOSE) PRINT*,'IN GID_OUT'
! 
! !BUILD SUBSIZES
! CNX=(NX-1)/GNX+1
! CNY=(NY-1)/GNY+1
! CNZ=(NZ-1)/GNZ+1
! CDX=GNX*DX
! CDY=GNY*DY
! CDZ=GNZ*DZ
! 
! 
! 
! !###########################################################################
! !SIZE OF OUTPUT
! OUTPUT_SIZE=((NX-1)/GNX+1)*((NY-1)/GNY+1)*((NZ-1)/GNZ+1)
! 
! !CREATES THE FLUID_NODE_ID
! DO I=1,OUTPUT_SIZE
!    FLUID_NODE_ID(I)=I
! ENDDO
! 
! !FILENAME FOR CURRENT STEP
! WRITE(FILE,'(A,I6.6,A)') TRIM(PATH_GID)//TRIM(GID_FILENAME)//".gid/"//        &
!                          TRIM(GID_FILENAME)//"_",STEP,".post.res"
! 
! !OPEN POST
! CALL GID_OPENPOSTRESULTFILE(TRIM(FILE),GIDMODE)
! 
! !WRITE MESH FOR COARSE GRID
! CALL WRITE_CMESH(CNX,CNY,CNZ,CDX,CDY,CDZ,FIRST_ELEM,FIRST_NODE)
! 
! !WRITE POST MOMENTUM
! CALL INTERP_UW(U,SCR1)
! CALL INTERP_UW(V,SCR2)
! CALL WRITE_CVECTOR(SCR1(1:NX:GNX,1:NY:GNY,1:NZ:GNZ),          &
!                    SCR2(1:NX:GNX,1:NY:GNY,1:NZ:GNZ),          &
!                    W(1:NX:GNX,1:NY:GNY,1:NZ:GNZ),'VELOCITY',  &
!                    FLUID_NODE_ID,JT,NX,NY,NZ)
! CALL INTERP_UW(P-0.5_RPREC*(U**2+V**2),SCR1)
! DO I=1,NZ-1
!    SCR1(:,:,I)=SCR1(:,:,I)-0.5_RPREC*(0.5_RPREC*(W(:,:,I)+W(:,:,I+1))**2)
! ENDDO
! CALL WRITE_CSCALAR(SCR1(1:NX:GNX,1:NY:GNY,1:NZ:GNZ),"P_FLUID",    &
!                    FLUID_NODE_ID,JT,NX,NY,NZ)
! CALL WRITE_CSCALAR(TXZ(1:NX:GNX,1:NY:GNY,1:NZ:GNZ),"TXZ_FLUID",   &
!                    FLUID_NODE_ID,JT,NX,NY,NZ)
! CALL WRITE_CSCALAR(TYZ(1:NX:GNX,1:NY:GNY,1:NZ:GNZ),"TYZ_FLUID",   &
!                    FLUID_NODE_ID,JT,NX,NY,NZ)
! CALL WRITE_CSCALAR(DUDZ(1:NX:GNX,1:NY:GNY,1:NZ:GNZ),"DUDZ_FLUID", &
!                    FLUID_NODE_ID,JT,NX,NY,NZ)
! CALL WRITE_CSCALAR(DVDZ(1:NX:GNX,1:NY:GNY,1:NZ:GNZ),"DVDZ_FLUID",   &
!                    FLUID_NODE_ID,JT,NX,NY,NZ)
! 
! 
! 
! !WRITE POST SCALARS
! IF(S_FLAG)THEN
!    CALL INTERP_UW(SC1,SCR1)
!    CALL WRITE_CSCALAR(SCR1(1:NX:GNX,1:NY:GNY,1:NZ:GNZ),"SC1",                  &
!                       FLUID_NODE_ID,JT,NX,NY,NZ)
!    CALL WRITE_CSCALAR(DSC1DZ(1:NX:GNX,1:NY:GNY,1:NZ:GNZ),"DSC1DZ",            &
!                       FLUID_NODE_ID,JT,NX,NY,NZ)
!    CALL WRITE_CSCALAR(TZSC1(1:NX:GNX,1:NY:GNY,1:NZ:GNZ),"TZSC1",              &
!                       FLUID_NODE_ID,JT,NX,NY,NZ)
! ENDIF
! 
! 
! 
! IF(INIT)THEN
!    CALL WRITE_CSCALAR(PHI(1:NX:GNX,1:NY:GNY,1:NZ:GNZ),"LEVEL_SET",   &
!                       FLUID_NODE_ID,JT,NX,NY,NZ)
!    CALL WRITE_CVECTOR(NORM(1,1:NX:GNX,1:NY:GNY,1:NZ:GNZ),            &
!                       NORM(2,1:NX:GNX,1:NY:GNY,1:NZ:GNZ),            &
!                       NORM(3,1:NX:GNX,1:NY:GNY,1:NZ:GNZ),'NORM',  &
!                       FLUID_NODE_ID,JT,NX,NY,NZ)
!    !FILTER PHI (in place so pay attention)
!    SCR1=PHI
!    DO JZ=1,NZ
!       CALL TEST_FILTER(SCR1(:,:,JZ),G_TEST)    
!    ENDDO
!    CALL WRITE_CSCALAR(SCR1(1:NX:GNX,1:NY:GNY,1:NZ:GNZ),"LEVEL_SET_FILTERED",   &
!                       FLUID_NODE_ID,JT,NX,NY,NZ)
!    INIT=.FALSE.
! ENDIF
! 
! !CLOSE POST
! CALL GID_CLOSEPOSTRESULTFILE
! 
! 
! 
! IF(VERBOSE) PRINT*,'OUT GID_OUT'
! 
! END SUBROUTINE GID_OUT
! 
! 
! 
! 
! !!-----------------------------------------------------------------------------------
! !!  OBJECT: subroutine GID_SLC_OUT()
! !!-----------------------------------------------------------------------------------
! !!
! !!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 08/11/2012
! !!
! !!  DESCRIPTION:
! !!
! !!  HL subroutine for GiD, makes use of GID_POST
! !!
! !!
! !!-----------------------------------------------------------------------------------
! SUBROUTINE GID_SLC_OUT(SLCNAME,STEP,SLC,ID_SLC)
! 
! USE LEVEL_SET
! USE GID_POST
! 
! IMPLICIT NONE
! 
! CHARACTER(4),INTENT(IN) :: SLCNAME
! INTEGER,INTENT(IN) :: STEP,SLC
! INTEGER,DIMENSION(:) :: ID_SLC
! INTEGER,DIMENSION(:,:,:),ALLOCATABLE :: IDN
! CHARACTER(256),SAVE :: FILE
! LOGICAL,SAVE :: INIT=.TRUE.
! 
! 
! !ALLOCATE INDEXES
! IF(SLC==1) ALLOCATE(IDN(NX,SIZE(ID_SLC),NZ))
! IF(SLC==2) ALLOCATE(IDN(NX,NY,SIZE(ID_SLC)))
! IF(SLC==3) ALLOCATE(IDN(SIZE(ID_SLC),NY,NZ))
! 
! !FILENAME
! WRITE(FILE,'(A,I6.6,A)') TRIM(PATH_GID)//TRIM(GID_FILENAME)//".gid/"//        &
!                          TRIM(GID_FILENAME)//TRIM(SLCNAME),STEP,".post.res"
! 
! !OPEN FILE AND WRITE MESH
! CALL GID_OPENPOSTRESULTFILE(TRIM(FILE),GIDMODE)
! CALL WSLC_MSH(SLC,ID_SLC,IDN,NX,NY,NZ,DX,DY,DZ,1,1)
! 
! !WRITE SLICE
! CALL WSLC_CSC(TXZ(:,:,1:NZ),"TXZ",SLC,ID_SLC,IDN,JT,NX,NY,NZ)
! CALL WSLC_CSC(TYZ(:,:,1:NZ),"TYZ",SLC,ID_SLC,IDN,JT,NX,NY,NZ)
! CALL WSLC_CSC(CS_OPT2(:,:,1:NZ),"CS",SLC,ID_SLC,IDN,JT,NX,NY,NZ)
! CALL WSLC_CSC(NU_T(:,:,1:NZ),"NU_T",SLC,ID_SLC,IDN,JT,NX,NY,NZ)
! CALL WSLC_CSC(P(:,:,1:NZ),"P",SLC,ID_SLC,IDN,JT,NX,NY,NZ)
! CALL WSLC_CVEC(U(:,:,1:NZ),V(:,:,1:NZ),W(:,:,1:NZ),"VEL",SLC,ID_SLC,IDN,JT,NX,NY,NZ)
! 
! !IF(INIT)THEN
!    CALL WSLC_CSC(PHI(:,:,1:NZ),"LEVEL_SET",SLC,ID_SLC,IDN,JT,NX,NY,NZ)
!    CALL WSLC_CVEC(NORM(1,:,:,1:NZ),NORM(2,:,:,1:NZ),  &
!                   NORM(3,:,:,1:NZ),"NORM",SLC,ID_SLC,IDN,JT,NX,NY,NZ)
!  ! INIT=.FALSE.
! !ENDIF
! 
! !CLOSE FILE
! CALL GID_CLOSEPOSTRESULTFILE
! 
! !DEALLOCATION
! DEALLOCATE(IDN)
! 
! 
! END SUBROUTINE GID_SLC_OUT


!!-----------------------------------------------------------------------------------
!!  OBJECT: subroutine SCREEN_OUT
!!-----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 27/10/2011
!!
!!  DESCRIPTION:
!!
!!  This subroutine writes a file with all the simulation's parameters.
!!
!!-----------------------------------------------------------------------------------
subroutine screen_out(jt_total,dt,tt,ke,rmsdivvel,cfl,cfl_old,visc_stab,visc_sc_stab)

implicit none

integer,intent(in) :: jt_total
real(rprec),intent(in) :: dt,tt,ke,rmsdivvel,cfl,visc_stab,cfl_old
real(rprec),intent(in),optional :: visc_sc_stab

!printings (make this an ext subroutine)
if(present(visc_sc_stab))then

   write(6,7777) jt_total,dt,tt,rmsdivvel,ke
   write(6,7779) cfl,cfl_old,visc_stab,visc_sc_stab

else

   write(6,7776) jt_total,dt,tt,rmsdivvel,ke
   write(6,7778) cfl,cfl_old,visc_stab
   
endif
   
7777 format ('jt,dt,tt,divu,ke,cfl,visc_stab,visc_sc_stab:',1x,i7.7,7(1x,e12.5))
7776 format ('jt,dt,tt,divu,ke:',1x,i7.7,4(1x,e12.5))
7778 format ('cfl,cfl_old,visc_stab:',1x,3(1x,e12.5))
7779 format ('cfl,cfl_old,visc_stab,visc_sc_stab:',1x,4(1x,e12.5))

end subroutine screen_out



!!-----------------------------------------------------------------------------------
!!  OBJECT: subroutine SSCREEN_OUT
!!-----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 27/10/2011
!!
!!  DESCRIPTION:
!!
!!  This subroutine writes on screen a real valued number
!!
!!-----------------------------------------------------------------------------------
subroutine sscreen_out(var,strvar)

implicit none

character(*),intent(in) :: strvar
real(rprec),intent(in) :: var

!printings
write(6,'(a,t25,e12.3)') strvar//' = ',var

end subroutine sscreen_out



!!-------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE INST_FIELD(VAR,STR_VAR,STEP,PATH,NX1,NX2,NY1,NY2,NZ1,NZ2)
!!-------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 20/03/2013
!!
!!  DESCRIPTION:
!!
!!  WRITES THE OUTPUT VAR THAT IS GIVEN
!!
!!  SX,SY,SZ	-->	MEANS A NODE EACH SX, START FROM 1 (e.g. 1,1+SX,1+2SX..)
!!			ITS BETTER IF WE USE EVEN VALUES FOR THESE NODES
!! 
!!
!!  INTP   -->  IF 1 MEANS INTERPOLATE TO W  A VAR WHICH IS FULLY ON UVP
!!         -->  IF 2 MEANS NO INTERPOLATION
!!  NOTE THAT THE CASE OF DUDZ WHERE THE VALUES ARE STORED IN UVP(1),W IS NOT
!!  CONTEMPLATED BECAUSE WE DON'T KNOW HOW MUCH DUDZ(0) IS..
!!
!!  DECIDED NOT TO PUT HERE THE INTERF STUFF CAUSE THAT IS SOMETHING SPECIFIC
!!  WHILE THIS SUBROUTINE IS TOP GENERAL FOR OUTPUTTING, BETTER THEN SAVING
!!  ON AN EXTERNAL FILE THE INTERF_W VALUES AND RETRIEVE THEM WITH MATLAB IN 
!!  PPROCESS PHASE.
!!
!!
!!-------------------------------------------------------------------------------
subroutine Inst_Field(var,strVar,intp,ind,cCnt,path,sx,sy,sz)

implicit none

real(rprec), dimension(ld,ny,0:nz),intent(in) :: var
integer,intent(in) :: intp,ind,sx,sy,sz,cCnt
character(*),intent(in) :: strVar,path
real(rprec),dimension(nx,ny,nz) :: scr1
integer :: i,j,k
character(len=256) :: filename
integer :: nx_fld,ny_fld,nz_fld
logical :: file_exists

!buil recl (brutal but effective way)
nx_fld=0
ny_fld=0
nz_fld=0
do i=1,nx,sx
  nx_fld=nx_fld+1
enddo
do j=1,ny,sy
   ny_fld=ny_fld+1
enddo
do k=1,nz-1,sz
   nz_fld=nz_fld+1
enddo 

!write parameters file
if(coord .eq. 0) then
   write(filename,'(a)') trim(path)//'parameters.txt'
   !INQUIRE(FILE=trim(filename),EXIST=file_exists)
   !if(.not.file_exists) then
     open(111,file=trim(filename),status='replace')
     write(111,*)  nx_fld
     write(111,*)  ny_fld
     write(111,*)  nz_tot-1
     write(111,*)  dx*sx
     write(111,*)  dy*sy
     write(111,*)  dz*sz
     write(111,*)  ind
     write(111,*)  nproc
     write(111,*)  RE
     write(111,*)  alpha
     if(LLBC_SPECIAL=='IBM_DFA')then
       write(111,*) 1
     else
       write(111,*) 0
     endif
     write(111,*)  phi_rot
     close(111)
   !endif
endif
call mpi_barrier(comm,ierr)

!interpolate var or simply copy (loop till nz in this case)
!if(intp.eq.1)then
!   do k=1,nz
!      scr1(1:nx:sx,1:ny:sy,k)=0.5_rprec*(var(1:nx:sx,1:ny:sy,k-1)+var(1:nx:sx,1:ny:sy,k))
!   enddo
!elseif(intp.eq.2)then
!   scr1(1:nx:sx,1:ny:sy,1:nz:sz)=var(1:nx:sx,1:ny:sy,1:nz:sz)
!endif

!write binary (nz-1 cause of MPI)
write(filename,'(a,i0)') trim(path)//trim(strVar)//'.c',coord
open(111,file=trim(filename),status='unknown',access='direct',recl=8*nx_fld*ny_fld*nz_fld)
write(111,rec=ind) var(1:nx,1:ny,1:nz-1)
close(111)

end subroutine Inst_Field



!!-------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE INST_SLICEs(VAR,STR_VAR,STEP,PATH,NX1,NX2,NY1,NY2,NZ1,NZ2)
!!-------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 20/03/2013
!!
!!  DESCRIPTION:
!!
!!
!!
!!-------------------------------------------------------------------------------
subroutine Inst_Slice(var,strVar,c_cnt,opath,ind)

implicit none

real(rprec), dimension(ld,ny,0:nz),intent(in) :: var
integer,intent(in) :: ind,c_cnt
character(*),intent(in) :: strVar,opath
integer :: i,j,k
character(len=256) :: filename
integer,dimension(3) :: nxs,nys
logical :: file_exists

! define slices
nys(1) = floor(ny/3.0)
nys(2) = floor(ny/2.0)
nys(3) = floor(ny*2.0/3.0)
nxs(1) = floor(nx/3.0)
nxs(2) = floor(nx/2.0)
nxs(3) = floor(nx*2.0/3.0)

if(coord .eq. 0) then
   write(filename,'(a)') trim(opath)//'parameters.txt'
!   INQUIRE(file=filename,EXIST=file_exists)
!   if(.not. file_exists) then 
     open(111,file=trim(filename),status='replace')
     write(111,*) nx
     write(111,*) ny
     write(111,*) nz_tot-1
     write(111,*) l_x
     write(111,*) l_y
     write(111,*) l_z
     write(111,*)  dx
     write(111,*)  dy
     write(111,*)  dz
     write(111,*)  ind
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
     write(111,*)  nxs
     write(111,*)  nys
     close(111)
!   endif
endif
call mpi_barrier(comm,ierr)

! write 3 slices xz binary (nz-1 cause of MPI)
write(filename,'(a,i0)') trim(opath)//trim(strVar)//'_XZ.c',coord
open(111,file=trim(filename),status='unknown',access='direct',recl=4*nx*size(nys)*(nz-1))
write(111,rec=ind) real(var(1:nx,nys,1:nz-1))
close(111)

! write 3 slices yz binary (nz-1 cause of MPI)
write(filename,'(a,i0)') trim(opath)//trim(strVar)//'_YZ.c',coord
open(111,file=trim(filename),status='unknown',access='direct',recl=4*ny*size(nxs)*(nz-1))
write(111,rec=ind) real(var(nxs,1:ny,1:nz-1))
close(111)



! write 1 slice xy binary 
write(filename,'(a,i0)') trim(opath)//trim(strVar)//'_XY.c',coord
open(111,file=trim(filename),status='unknown',access='direct',recl=4*nx*ny)
write(111,rec=ind) real(var(1:nx,1:ny,floor((nz-1.0)/2.0)))
close(111)


end subroutine Inst_Slice

subroutine inlet_slice(var,strVar,c_cnt,opath,ind)

implicit none

real(rprec), dimension(ld,ny,0:nz),intent(in) :: var
integer,intent(in) :: ind,c_cnt
character(*),intent(in) :: strVar,opath
integer :: i,j,k
character(len=256) :: filename
logical :: file_exists

if(coord .eq. 0) then
   write(filename,'(a)') trim(opath)//'parameters.txt'
   INQUIRE(file=filename,EXIST=file_exists)
   if(.not. file_exists) then 
     open(111,file=trim(filename),status='new')
     write(111,*) nx
     write(111,*) ny
     write(111,*) nz_tot-1
     write(111,*) l_x
     write(111,*) l_y
     write(111,*) l_z
     write(111,*)  dx
     write(111,*)  dy
     write(111,*)  dz
     write(111,*)  dt*c_cnt
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
     close(111)
   endif
endif
call mpi_barrier(comm,ierr)

! write 3 slices yz binary (nz-1 cause of MPI)
write(filename,'(a,i0)') trim(opath)//trim(strVar)//'_YZ.c',coord
open(111,file=trim(filename),status='unknown',access='direct',recl=8*ny*(nz-1))
write(111,rec=ind) var((nx/2),1:ny,1:nz-1)
close(111)

end subroutine inlet_slice

!!-------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE HAVG_VAL
!!-------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 22/10/2012
!!
!!  DESCRIPTION:
!!
!!  This subroutine writes a ASCII file with the values at dt of a given scalar
!!
!!  TO IMPROVE: adopt binary format
!!
!!-----------------------------------------------------------------------------------
SUBROUTINE Avg_Ke(keGlobOut)

implicit none

real(rprec),intent(out) :: keGlobOut
integer :: jx,jy,jz
character(len=256) :: filename
real(rprec) :: tmp_w,denom,ke

!initialize
ke=0.0_rprec
denom=1.0_rprec/(2.0_rprec*nx*ny*(nz-1))

!compute ke
do jz=1,nz-1
do jy=1,ny
do jx=1,nx
   tmp_w=.5_rprec*(w(jx,jy,jz)+w(jx,jy,jz+1))
   ke=ke+(u(jx,jy,jz)**2+v(jx,jy,jz)**2+tmp_w**2)
enddo
enddo
enddo
ke=ke*denom

!reduce sum to coord0
#IF(MPI) 
   call mpi_reduce(ke,keGlobOut,1,mpi_rprec,mpi_sum,0,comm,ierr)
#ENDIF

!write on file
if(coord==0)then
   keGlobOut=keGlobOut/nproc   !missing part in averaging
   open(111,file='ke.txt',status='unknown',position='append')
   write(111,'(f20.10)') keGlobOut
endif 

END SUBROUTINE Avg_Ke



!!-------------------------------------------------------------------------------
!!  OBJECT: SUBROUTINE HAVG_VAL
!!-------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 02/11/2012
!!
!!  DESCRIPTION:
!!
!!  This subroutine writes a ASCII file with the values at dt of a given scalar
!!
!!  TO IMPROVE: adopt binary format
!!
!!-------------------------------------------------------------------------------
subroutine Havg_Val(var,varOut,strVar)

implicit none

real(rprec),dimension(:,:),intent(in) :: var
real(rprec),intent(out) :: varOut
character(*),intent(in) :: strVar
integer :: nx,ny
character(256) :: file_loc

!define sizes
nx=size(var(:,1))
ny=size(var(1,:))

!compute vol avg
varOut=sum(var(:,:))/(nx*ny)

write(file_loc,'(a,a,i0,a)') trim(strVar),'_',colour,'.txt'      
! write to file (ascii)
open(111,file=trim(file_loc),status='unknown',position='append')
write(111,*) varOut
close(111)
         

end subroutine Havg_Val



subroutine timestep_conditions(cflGlob,cflOldGlob,viscStabGlob,nu_mlc,nu_t)

implicit none

real(kind=rprec) :: delta, u_res_max, nu_max
real(kind=rprec),dimension(1,5) :: max_vels
real(rprec) :: cfl,cflOld,viscStab
real(rprec),intent(out) :: cflGlob,cflOldGlob,viscStabGlob
real(rprec),intent(in) :: nu_mlc
real(rprec),dimension(:,:,0:),optional :: nu_t

max_vels(1,1)=maxval(u(1:nx,1:ny,1:nz-1))
max_vels(1,2)=maxval(v(1:nx,1:ny,1:nz-1))
max_vels(1,3)=maxval(w(1:nx,1:ny,1:nz-1))
u_res_max = sqrt(max_vels(1,1)**2+max_vels(1,2)**2+max_vels(1,3)**2) ! don't bother with interp here

if(present(nu_t))then 
   nu_max=maxval(nu_t(1:nx,1:ny,1:nz-1))
else
   nu_max=nu_mlc
endif

! compute values
cfl = max(max_vels(1,1)*dt/dx,max_vels(1,2)*dt/dy,max_vels(1,3)*dt/dz)
delta = min(dx,dy,dz)
cflOld = u_res_max*dt/delta
viscStab = dt*nu_max/delta**2

! reduce sum to coord0
#IF(MPI)
   call mpi_reduce(cfl,cflGlob,1,mpi_rprec,mpi_max,0,comm,ierr)
   call mpi_reduce(cflOld,cflOldGlob,1,mpi_rprec,mpi_max,0,comm,ierr)
   call mpi_reduce(viscStab,viscStabGlob,1,mpi_rprec,mpi_max,0,comm,ierr)
#ENDIF

end subroutine timestep_conditions



END MODULE STD_OUT
