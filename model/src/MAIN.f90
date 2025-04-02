!!------------------------------------------------------------------------------
!!  OBJECT: Main file of LES-LSM
!!------------------------------------------------------------------------------
!!
!!  MAIN DEVELOPER: Marco Giometto ( mgiometto@gmail.com )
!!  CONTRIBUTORS: Varun Sharma (varun.sharma@epfl.ch)
!!                Francesco Comola (francesco.comola.gmail.com)
!!                Armin Sigmund (armin.sigmund@epfl.ch)
!!                Daniela Brito Melo (daniela.britomelo@epfl.ch)
!!
!!  LAST UPDATE: May 19, 2021
!!
!!  DESCRIPTION: This application uses a Large-Eddy-Simulation (LES) technique
!!               to solve the air flow, temperature and humidity fields on an
!!               atmospheric boundary layer. The code is coupled to a Lagrangian
!!               Stochastic Model for particle dynamics (see LSM_simple.f90).
!!               The code is organized as follows:
!!               - Initialization
!!               - Loop in time:
!!                   - Impose boundary conditions (BC)
!!                   - Compute velocity derivatives, sgs stresses and RHS
!!                   - Compute convective term
!!                   - Compute scalars
!!                   - Add coriolis force
!!                   - Time advancement
!!                   - Particle solver (LSM)  OR  Wind farm model
!!                   - Compute velocities
!!                   - Pressure solver
!!                   - Write outputs
!!
!! MAIN REFERENCES: Albertson, J. D. (1996) Large-Eddy Simulation of Land
!!                  Atmosphere Interaction, PhD thesis, University of California.
!!
!!                  Albertson and Parlange (1999) Natural Integration of Scalar
!!                  Fluxes from Complex Terrain, Advances in Water Resources.
!!
!!                  Giometto M, (2016) Theoretical and numerical studies of
!!                  atmospheric boundary layer flows over complex terrain,
!!                  PhD thesis, EPFL.
!!
!!                  Sharma, V. et al. (2018) On the suitability of the Thorpe-
!!                  -Mason model for calculating sublimation of saltating snow,
!!                  The Cryosphere.
!!
!!                  Comola, F. et al. (2019) Preferential Deposition of Snow
!!                  and Dust Over Hills: Governing processes and relevant scales,
!!                  Journal of Geophysical Research: Atmospheres.
!!------------------------------------------------------------------------------ 

PROGRAM MAIN

USE TYPES
USE PARAM
USE SIM_PARAM
USE DERIV
USE FFT
USE IC
USE BC
USE TEST_FILTERMODULE
USE STD_OUT
USE SPNG
USE SERVICE
USE IFPORT
USE SP_SC1
USE SP_SC2
!USE SP_SC3
USE SGSMODULE
USE SGS_STRESSES
USE SCALARS_MODULE
USE OUTLOOP_SPECIAL
USE MO_SIMILARITY
USE LS_MODEL
USE MWALL_LAW
USE CNVC
USE AVGFIELD
USE STRETCH
USE MMATRIX
USE compute_windfarm
!USE turbine_model_expanded
 
IMPLICIT NONE

INTEGER :: IP_MAIN, MPI_NP, COORDS(1), KKK(1)
INTEGER :: JX,JY,JZ
INTEGER :: ITER,RESULT
CHARACTER(256) :: HEADERS,PHI_FILE,SRF_FILE,sscr1,sscr2,filename, filename_hetero, filename_patches
LOGICAL :: CHECK
REAL(RPREC) :: SCR1,SCR2,P1,P2,P1mom
REAL(RPREC) :: T1,T2,T3,T4,T5,T6
integer :: jz_global
integer :: omp_get_num_procs,omp_get_num_threads,omp_get_thread_num
logical :: file_exists
integer :: kk_main
integer :: mmm
integer :: second_comm,second_colour,second_rank
real(rprec),dimension(3,3,3) :: hehe,tmp_hehe
integer :: mpi_blend
real(rprec) :: rand_num
LOGICAL :: PATCHES(NX, NY)

call mpi_init (ierr)
call mpi_comm_size (MPI_COMM_WORLD, mpi_np, ierr)
call mpi_comm_rank (MPI_COMM_WORLD, global_rank, ierr)

!--check if run-time number of processes agrees with nproc parameter
if (mpi_np /= nproc)then
  write (*, *) 'runtime number of procs = ', np,  &
               ' not equal to nproc = ', nproc
  stop
endif

colour=0
!colour = mod(global_rank,2)
!second_colour = global_rank/2
!call mpi_comm_split(MPI_COMM_WORLD,colour,global_rank,new_comm,ierr)
!call mpi_comm_split(MPI_COMM_WORLD,second_colour,global_rank,second_comm,ierr)

!--set up a 1d cartesian topology 
call mpi_cart_create (mpi_comm_world, 1, (/ nproc /), (/ .false. /),  &
                        .true., comm, ierr)
  !--slight problem here for ghost layers:
  !  u-node info needs to be shifted up to proc w/ rank "up",
  !  w-node info needs to be shifted down to proc w/ rank "down"
call mpi_cart_shift (comm, 0, 1, down, up, ierr)
call mpi_comm_rank (comm, rank, ierr)
call mpi_cart_coords (comm, rank, 1, coords, ierr)
coord = coords(1)  !--use coord (NOT rank) to determine global position
write (chcoord, '(a,i0,a)') '(', coord, ')'  !--() make easier to use
!--rank->coord and coord->rank conversions
do ip_main = 0, nproc-1
  call mpi_cart_rank (comm, (/ ip_main /), rank_of_coord(ip_main), ierr)
  call mpi_cart_coords (comm, ip_main, 1, kkk, ierr)
  coord_of_rank(ip_main) = kkk(1)
end do

!call mpi_comm_rank(second_comm,second_rank,ierr)

!write (*,*) 'Hello! from process with coord,rank,colour = ', coord,rank,colour,second_rank,second_colour,global_rank



!--set the MPI_RPREC variable
MPI_RPREC = MPI_DOUBLE_PRECISION
MPI_CPREC = MPI_DOUBLE_COMPLEX
call MPI_TYPE_CREATE_F90_INTEGER(INTPREC, MPI_INTPREC, ierr)  

! Re distribute the OMP threads
!if(coord .le. 2) then
! call omp_set_num_threads(4)
!else
! call omp_set_num_threads(4)
!endif

!$omp parallel
 write(*,*) 'Hello from:',coord,omp_get_thread_num(),omp_get_num_procs()
!$omp end parallel

! VERTICAL DOMAIN DECOMPOSITION
!--check whether size of NZ_VECTOR is nproc
if ( SIZE(NZ_VECTOR) /= nproc ) then
  write (*, *) 'SIZE(NZ_VECTOR) = ', SIZE(NZ_VECTOR),  &
               ' not equal to nproc = ', nproc
  stop
  !--check whether NZ_VECTOR is consistent with NZ_TOT
elseif ( (SUM(NZ_VECTOR-1) /= NZ_TOT-1) .and. coord == 0 ) then
  write (*, *) 'SUM(NZ_VECTOR-1) = ', SUM(NZ_VECTOR-1),  &
               ' not equal to NZ_TOT-1 = ', NZ_TOT-1
  stop
endif
! Index of uppermost level (ghost layer) in the subdomain array depending on the MPI process
NZ = NZ_VECTOR(coord)
! Number of levels in subdomain and lower subdomains (excluding ghost layers) depending on the MPI process (cumulative)
NZ_CUMUL(0) = NZ_VECTOR(0) - 1
DO ip_main = 1, nproc-1
  NZ_CUMUL(ip_main) = NZ_VECTOR(ip_main) - 1 + NZ_CUMUL(ip_main-1)
END DO
! height of bottom subdomain boundary in index units
ETA_DISPL = (/ 0, NZ_CUMUL(0:nproc-2) /)

write(out_path,'(a,i0,a)') './output_',colour,'/'

!generate main output folder
if(coord==0) result=system('mkdir '//trim(out_path))

write(filename,'(a,i0)') './PARAMETERS.py_',colour
write(filename_patches,'(a)') './patches.in'
write(filename_hetero,'(a)') './hetero_srf_props.in'
!read parames from ext. file
call read_param(filename)
!read patches from file
call read_patches(filename_patches, NX, NY, PATCHES)
!read surface heterogeneity parameters
call read_hetero_srf_props(filename_hetero)
call mpi_barrier(comm,ierr)

!allocation
CALL ALLOC_SIM_PARAM()
CALL ALLOCATE_SGS()

IF(TEMPERATURE)  CALL allocTemperature()

IF(HUMIDITY) CALL allocHumidity()
!IF(HUMIDITY) CALL allocHumidity_particle()
!INITIALIZE JT_TOTAL (ISSUE WITH BC WHEN USING DTM SO PUT HERE)

IF(LSM) CALL ALLOCATE_LSM()


!STENCILS
DX=L_X/NX
DY=L_Y/NY
DZ=1.0_rprec

call init_stretch
!C_STRETCH = 1.0_rprec * LBC_DZ0
!write(*,*) 'c_stretch:',C_STRETCH
!L_eta=dlog(1.0_rprec + (L_Z/C_STRETCH))
!DZ=L_eta/(NZ_tot-1)

!DZ=L_Z/(NZ_TOT-1)
!do jz=0,48
! DZ(jz)=L_Z/(NZ_TOT-1)
!enddo
!do jz=49,65
! DZ(jz)=2.0_rprec * L_Z/real((NZ_TOT-1),rprec)
!enddo


!---------------------------------------------------------------------------------
!INITIAL CONDITIONS
IF(TEMPERATURE)THEN
   !CALL SET_LBC_SC_ZO(LBC_SC1_ZO,1)
   CALL SET_LBC_SC_ZO_HETERO(LBC_SC1_ZO,PATCHES,NX,NY,Z0_0,Z0_1)
   !CALL COMPUTE_LBC_SC(LBC_SC1,1)    !NEED LBC_SC1 IN INIT_SC FOR NOISE NORM.
   CALL COMPUTE_LBC_SC_HETERO(LBC_SC1,PATCHES,NX,NY,T_0,T_1)
   CALL INIT_SC( SC1,SC1_INIT,DSC1DZ_INIT,               &
                 NF_SC1,Z_TURB_SC1,                      &
                 Z_ENT_SC1,ENT_SC1_INIT,ENT_DSC1DZ_INIT,RHS_SC1, &
                 SC1STAR_INIT, LBC_SC1_ZO(1,1) )
ENDIF

IF(HUMIDITY)THEN
   CALL SET_LBC_SC_ZO(LBC_SC2_ZO,2)
   CALL COMPUTE_LBC_SC(LBC_SC2,2)    !NEED LBC_SC1 IN INIT_SC FOR NOISE NORM.
   CALL INIT_SC( SC2,SC2_INIT,DSC2DZ_INIT,               &
                 NF_SC2,Z_TURB_SC2,                      &
                 Z_ENT_SC2,ENT_SC2_INIT,ENT_DSC2DZ_INIT,RHS_SC2, &
                 SC2STAR_INIT, LBC_SC2_ZO(1,1) )
ENDIF
!IF(HUMIDITY)THEN
!   CALL SET_LBC_SC_ZO(LBC_SC3_ZO,2)
!   CALL COMPUTE_LBC_SC(LBC_SC3,2)    !NEED LBC_SC1 IN INIT_SC FOR NOISE NORM.
!   CALL INIT_SC( SC3,SC3_INIT,DSC3DZ_INIT,               &
!                 NF_SC3,Z_TURB_SC3,                      &
!                 Z_ENT_SC3,ENT_SC3_INIT,ENT_DSC3DZ_INIT,RHS_SC3)
!ENDIF


!set momentum ic (important to have set sc1,sc2 before if s_flag)
CALL SET_LBC_ZO(LBC_Z0) !IF WALL_LAW
CALL INIT(U,V,W,U_STAR,PFX,PFY,DSPL_H)
CALL MPI_BARRIER(COMM,IERR)


!initialize time
tt=0
jt_total=0
np_global = 0

INQUIRE(FILE='./total_time.dat',EXIST=file_exists)
if(file_exists) then
    open(111,file='./total_time.dat',action='read')
    read(111,*) jt_total
    close(111)
endif

call mpi_barrier(MPI_COMM_WORLD,ierr)

if(.not. file_exists) then
   if((colour .eq. 0) .and. (coord .eq. 0)) then
     open(111,file='./total_time.dat',action='write')
     write(111,*) jt_total
     close(111)
   endif
endif

INQUIRE(FILE='./np_global.dat',EXIST=file_exists)
if(file_exists) then
    open(111,file='./np_global.dat',action='read')
    read(111,*) np_global
    close(111)
endif

call mpi_barrier(MPI_COMM_WORLD,ierr)

if(.not. file_exists) then
   if((colour .eq. 0) .and. (coord .eq. 0)) then
     open(111,file='./np_global.dat',action='write')
     write(111,*) np_global
     close(111)
   endif
endif

!IF(JT_TOTAL .GT. 4*LSM_INIT) then
!            sublimation = .true.
!ENDIF



call mpi_barrier(MPI_COMM_WORLD,ierr)

IF(LSM) THEN
   CALL INITIALIZE_LSM()
   CALL INITIALIZE_lsm_surface()
   CALL INITIALIZE_PARCELS()
ENDIF
CALL MPI_BARRIER(COMM,IERR)

!check if there's a restart file
CALL R_RESTART()
CALL MPI_BARRIER(COMM,IERR)

!init velocities (put in a service routine)
IF(COORD==0)THEN
   CALL SET_VALUE(W(:,:,:),0,0.0_RPREC)   !IMPORTANT NOT BOGUS FOR RESTART
   CALL SET_VALUE(P(:,:,:),0,0.0_RPREC)
   CALL SET_VALUE(TXX(:,:,:),0,BOGUS)
   CALL SET_VALUE(TXY(:,:,:),0,BOGUS)
   CALL SET_VALUE(TYY(:,:,:),0,BOGUS)
   CALL SET_VALUE(TZZ(:,:,:),0,BOGUS)
ELSEIF(COORD==NPROC-1)THEN
   CALL SET_VALUE(P(:,:,:),NZ,0.0_RPREC)
   CALL SET_VALUE(TXX(:,:,:),NZ,BOGUS)
   CALL SET_VALUE(TXY(:,:,:),NZ,BOGUS)
   CALL SET_VALUE(TYY(:,:,:),NZ,BOGUS)
   CALL SET_VALUE(TZZ(:,:,:),NZ,BOGUS)
ENDIF
   
!sync velocities (put in a service routine)
CALL MPI_SENDRECV (U(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 1,  &
                U(1, 1, NZ), LD*NY, MPI_RPREC, UP, 1,   &
                COMM, STATUS, IERR)
CALL MPI_SENDRECV (V(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 2,  &
                V(1, 1, NZ), LD*NY, MPI_RPREC, UP, 2,   &
                COMM, STATUS, IERR)
CALL MPI_SENDRECV (W(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 3,  &
                W(1, 1, NZ), LD*NY, MPI_RPREC, UP, 3,   &
                COMM, STATUS, IERR)                  
CALL MPI_SENDRECV (U(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 4,  &
                U(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 4,   &
                COMM, STATUS, IERR)
CALL MPI_SENDRECV (V(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 5,  &
                V(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 5,   &
                COMM, STATUS, IERR)
CALL MPI_SENDRECV (W(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 6,  &
                W(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 6,   &
                COMM, STATUS, IERR)

IF(TEMPERATURE) then
CALL MPI_SENDRECV (SC1(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 307,  &
                SC1(1, 1, NZ), LD*NY, MPI_RPREC, UP, 307,   &
                COMM, STATUS, IERR)                  
CALL MPI_SENDRECV (SC1(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 407,  &
                SC1(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 407,   &
                COMM, STATUS, IERR)
ENDIF

IF(HUMIDITY) then
CALL MPI_SENDRECV (SC2(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 303,  &
                SC2(1, 1, NZ), LD*NY, MPI_RPREC, UP, 303,   &
                COMM, STATUS, IERR)                  
CALL MPI_SENDRECV (SC2(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 404,  &
                SC2(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 404,   &
                COMM, STATUS, IERR)
ENDIF
!IF(HUMIDITY) then
!CALL MPI_SENDRECV (SC3(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 311,  &
!                SC3(1, 1, NZ), LD*NY, MPI_RPREC, UP, 311,   &
!                COMM, STATUS, IERR)                  
!CALL MPI_SENDRECV (SC3(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 411,  &
!                SC3(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 411,   &
!                COMM, STATUS, IERR)
!ENDIF

!---------------------------------------------------------------------------------
!VARIOUS INITIALIZATIONS

!initialize kx,ky matrices (with 0 for nyquist)
call init_fft()
if(sgs.or.(llbc_special.eq.'WALL_LAW'))then
   call test_filter_init(2._rprec*filter_size,g_test,ifilter)
   if(sgs.and.(model==3 .or. model==5 .or. model==7))then
      call test_filter_init(4._rprec*filter_size,g_test_test,ifilter)
   endif
   print*,'initialized test filters'
endif
if(TEMPERATURE .and. v_sponge) call init_sponge()

!screen print
if(coord==0) call print_params()

!save time
T1= MPI_WTIME()

    P1 = RI_SC1          ! this is N² in the dimensional equations
    P2 = 1.0_RPREC/RE    ! this is \nu in the dimensional equations
    P1mom = 1.0_RPREC    ! this is =1 in the dimensional equations

call initAvgField()

!call wind_farm_initialize()
if(wind_farm) then
call windfarm_init()
endif
call mpi_barrier(MPI_COMM_WORLD,ierr)

call init_mmatrix()
!call mpi_op_create(compute_mmatrix,.false.,mpi_blend,ierr)
call mpi_barrier(MPI_COMM_WORLD,ierr)

!########################################################
!########################################################
!MAIN LOOP
DO JT=1,NSTEPS

   !VARY TIME STEP
   IF((JT.LE.TR_STEPS))THEN
      DT=REF_DT*TR_DT_FACTOR
   ELSE
      DT=REF_DT
   ENDIF

   !ADVANCE TIME AND STEP
   
   JT_TOTAL=JT_TOTAL+1
   TT = JT_TOTAL*DT*(Z_I/U_STAR)

   !UPDATE OLD RHS
   RHSX_F=RHSX
   RHSY_F=RHSY
   RHSZ_F=RHSZ
  
   !#####################################################  
   !BOUNDARY CONDITIONS

   !INLET AND BUFFER LAYER
   IF(V_INLET) CALL SET_INLET()
   IF(V_BLAYER) CALL BLAYER(U,V,W,8.0_RPREC)
 

   IF (COORD==0) THEN 
         CALL COMPUTE_LBC(LBC_U,LBC_V,LBC_W)
         CALL SET_BC(0,U,LLBC_U,LBC_U,1)
         CALL SET_BC(0,V,LLBC_V,LBC_V,1)
         CALL SET_BC(0,W,LLBC_W,LBC_W,0)
   ENDIF
   
   IF(COORD==NPROC-1) THEN 
      CALL COMPUTE_UBC(UBC_U,UBC_V,UBC_W)
      CALL SET_BC(1,U,LUBC_U,UBC_U,1)
      CALL SET_BC(1,V,LUBC_V,UBC_V,1)
      CALL SET_BC(1,W,LUBC_W,UBC_W,0)
   ENDIF
   
   !###########################################################
   !BOUNDARY CONDITIONS

   !derivatives
   call filt_da(u,dudx,dudy)
   call filt_da(v,dvdx,dvdy)
   call filt_da(w,dwdx,dwdy)
   call ddz_uv(dudz,u)
   call ddz_uv(dvdz,v)
   call ddz_w(dwdz,w)

  call mpi_sendrecv (dwdz(1,1,1),ld*ny,mpi_rprec,down,7,  &
                     dwdz(1,1,nz),ld*ny,mpi_rprec,up,7,   &
                     comm,status,ierr)

  CALL MPI_SENDRECV (dwdz(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 4,  &
                dwdz(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 4,   &
                COMM, STATUS, IERR)


   call mpi_barrier(comm,ierr)
   !#################################################################
   CALL OBUKHOV()   
   !overwrite vel.deriv. (will be used for stress calculations) + stresses at jz=1
   if((coord==0).and.(llbc_special.eq.'WALL_LAW'))then
      call WALL_LAW()
   endif

   !stress calculations
   if(sgs)then
      call sgs_stress(p2)
   else
      call dns_stress(p2)
   endif

   !rhs at uvp (rhsxy) and w (rhsz) nodes
   call convec(rhsx,rhsy,rhsz)  

   !compute divergence of stress
   call divstress(rhsx,txx,txy,txz,.true.)
   call divstress(rhsy,txy,tyy,tyz,.true.)
   call divstress(rhsz,txz,tyz,tzz,.false.)

   !add sponge layer (uvp for rhsx and rhsy need interpolation) (time=0ms)
   if(v_sponge) call set_vsponge(u,v,w,rhsx,rhsy,rhsz)

   IF(TEMPERATURE) then 
    call calc_beta()
    IF(JT_TOTAL .LE. 1000) then
         BUOYANCY=0.0_rprec
    ENDIF
    DO JZ=1,NZ   !RHSZ @JZ=1 COULD AVOID SINCE WE IMPOSE W=0
      RHSz(:,:,JZ) = RHSz(:,:,JZ) + BUOYANCY(:,:,JZ)
    ENDDO
   ENDIF

   !###################################################################
   !###################################################################
   !SCALAR COMPUTATIONS AND COUPLING
   IF(TEMPERATURE)THEN
      IF(SGS)THEN
        IF(TEMPERATURE) THEN 
          IF(COORD==NPROC-1) CALL COMPUTE_UBC_SC( UBC_SC1, 1 )
          CALL SCALAR_SOLVER (SC1,DSC1DX,DSC1DY,DSC1DZ,TXSC1,TYSC1,TZSC1,   &
                             RHS_SC1,RHSF_SC1,                             &
                             LLBC_SC1,DLBC_SC1,LBC_SC1,                    &
                             LUBC_SC1,DUBC_SC1,UBC_SC1,                    &
                             PR_SC1,PR_SGS_SC1,                     &
                             CS_OPT2,NU_T,DS_OPT2_SC1,NU_SC1,SC1_STAR,PSI_H,PHI_H,Ft_lsm)
       ENDIF
       IF(HUMIDITY) THEN 
         IF(COORD==NPROC-1) CALL COMPUTE_UBC_SC( UBC_SC2, 2 )
         CALL SCALAR_SOLVER (SC2,DSC2DX,DSC2DY,DSC2DZ,TXSC2,TYSC2,TZSC2,   &
                            RHS_SC2,RHSF_SC2,                             &
                            LLBC_SC2,DLBC_SC2,LBC_SC2,                    &
                            LUBC_SC2,DUBC_SC2,UBC_SC2,                    &
                            PR_SC2,PR_SGS_SC2,                     & ! NO PARTICLE_FORCING HERE !
                            CS_OPT2,NU_T,DS_OPT2_SC2,NU_SC2,SC2_STAR,PSI_V,PHI_V,Fq_lsm)
       ENDIF
!       IF(HUMIDITY) THEN 
!         CALL SCALAR_SOLVER (SC3,DSC3DX,DSC3DY,DSC3DZ,TXSC3,TYSC3,TZSC3,   &
!                            RHS_SC3,RHSF_SC3,                             &
!                            LLBC_SC3,DLBC_SC3,LBC_SC3,                    &
!                            LUBC_SC3,DUBC_SC3,UBC_SC3,                    &
!                            PR_SC3,PR_SGS_SC3,                     &
!                            CS_OPT2,NU_T,DS_OPT2_SC3,NU_SC3,SC3_STAR,PSI_V,PHI_V)
!       ENDIF
 
  
     ELSE
        IF(COORD==NPROC-1) CALL COMPUTE_UBC_SC( UBC_SC1, 1 )
        CALL SCALAR_SOLVER (SC1,DSC1DX,DSC1DY,DSC1DZ,TXSC1,TYSC1,TZSC1,   &
                            RHS_SC1,RHSF_SC1,                             &
                            LLBC_SC1,DLBC_SC1,LBC_SC1,                    &
                            LUBC_SC1,DUBC_SC1,UBC_SC1,                    &
                            PR_SC1,PR_SGS_SC1,Ft_lsm ) 
     ENDIF
     !ADD COUPLING TERM IN MOMENTUM RHS (BUOYANCY)
  
  ENDIF

   !!!!! CORIOLIS FORCING

   !$omp parallel do
   do jz=1,nz-1
      RHSx(:,:,jz) = RHSx(:,:,jz) + coriol*v(:,:,jz) - coriol*vg
      RHSy(:,:,jz) = RHSy(:,:,jz) - coriol*u(:,:,jz) + coriol*ug
   enddo    
   !$omp end parallel do

   !!!!! END OF CORIOLIS FORCING

   

   !#############################################################
   !TIME ADVANCEMENT

   !choose the time scheme (ab or euler depending if its 1st step or not)
   if (jt_total == 1) then
      tadv1 = 1.0_rprec
      tadv2 = 0.0_rprec
   else
      tadv1 = 1.5_rprec
      tadv2 = 1.0_rprec-tadv1
   end if

   !##############################################################
   ! LAGRANGIAN STOCHASTIC MODEL
    IF (LSM) THEN
       CALL MPI_SENDRECV (U(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 1,  &
                          U(1, 1, NZ), LD*NY, MPI_RPREC, UP, 1,   &
                          COMM, STATUS, IERR)
       CALL MPI_SENDRECV (V(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 2,  &
                          V(1, 1, NZ), LD*NY, MPI_RPREC, UP, 2,   &
                          COMM, STATUS, IERR)
       CALL MPI_SENDRECV (W(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 3,  &
                          W(1, 1, NZ), LD*NY, MPI_RPREC, UP, 3,   &
                          COMM, STATUS, IERR)
       CALL MPI_SENDRECV (U(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 4,  &
                          U(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 4,   &
                          COMM, STATUS, IERR)
       CALL MPI_SENDRECV (V(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 5,  &
                          V(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 5,   &
                          COMM, STATUS, IERR)
       CALL MPI_SENDRECV (W(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 6,  &
                          W(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 6,   &
                          COMM, STATUS, IERR)

       IF(TEMPERATURE) then
       CALL MPI_SENDRECV (SC1(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 300,  &
                          SC1(1, 1, NZ), LD*NY, MPI_RPREC, UP, 300,   &
                          COMM, STATUS, IERR)
       CALL MPI_SENDRECV (SC1(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 400,  &
                          SC1(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 400,   &
                          COMM, STATUS, IERR)
       ENDIF
       IF(HUMIDITY) then
       CALL MPI_SENDRECV (SC2(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 301,  &
                          SC2(1, 1, NZ), LD*NY, MPI_RPREC, UP, 301,   &
                          COMM, STATUS, IERR)
       CALL MPI_SENDRECV (SC2(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 401,  &
                          SC2(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 401,   &
                          COMM, STATUS, IERR)
       ENDIF

!       IF(HUMIDITY) then
!       CALL MPI_SENDRECV (SC3(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 302,  &
!                          SC3(1, 1, NZ), LD*NY, MPI_RPREC, UP, 302,   &
!                          COMM, STATUS, IERR)
!       CALL MPI_SENDRECV (SC3(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 402,  &
!                          SC3(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 402,   &
!                          COMM, STATUS, IERR)
!       ENDIF



       !IF (JT_TOTAL.GE.(LSM_INIT-2*AVERAGE_FREQ)) THEN
       !    CALL AVERAGE_LES()
       !END IF
       IF ((JT_TOTAL.GE.LSM_INIT).AND. (MOD((JT_TOTAL-LSM_INIT), LSM_FREQ) == 0)) THEN

        Fx_lsm(:,:,:)=0.0_rprec
        Fy_lsm(:,:,:)=0.0_rprec
        Fz_lsm(:,:,:)=0.0_rprec
        
        IF(temperature) then
         Ft_lsm(:,:,:) = 0.0_rprec
         Ft_lsm_1(:,:,:) = 0.0_rprec
         Ft_lsm_2(:,:,:) = 0.0_rprec
        ENDIF
        IF(humidity) then
         Fq_lsm(:,:,:) = 0.0_rprec
        ENDIF 
         conc(:,:,:) = 0.0_rprec
         tpxz(:,:) = 0.0_rprec
         tpyz(:,:) = 0.0_rprec

         sum_delta_m=0.0_rprec
         sum_delta_q=0.0_rprec
 
         CALL LSM_RUN()
       END IF
       DO JZ=1,NZ
          RHSX(:,:,JZ) = RHSX(:,:,JZ) + FX_LSM(:,:,JZ)
          RHSY(:,:,JZ) = RHSY(:,:,JZ) + FY_LSM(:,:,JZ)
          RHSZ(:,:,JZ) = RHSZ(:,:,JZ) + FZ_LSM(:,:,JZ)
       END DO
    
    END IF

  if(colour .eq. 0) then
   if (wind_farm) then
    if (jt_total+1 >= start_disk) then !(jt_total+1 = jt).
      call windfarm_main(jt,jt_total)
      CALL MPI_BARRIER(comm,ierr)
    end if
   end if
  endif
   !advection step
   !if(colour .eq. 0) then
   !     if(jt .ge. 1000) then 
   !       pfx = 0.0_rprec
   !       pfy = 0.0_rprec
   !     endif
   !endif
  
   if(noise_pressure) then
    if(jt_total .gt. 5000) then
     call random_number(rand_num)
     if(mod(jt_total,1000) .eq. 0) then
      pfy = 0.2_rprec * pfx * ( 2.0_rprec*rand_num - 1.0_rprec)
      !write(*,*) '****************==== ',pfy,' ====***************'
     endif
    endif
   endif  
 
   do jz=1,nz-1
      u(:,:,jz) = u(:,:,jz)+dt*(tadv1*rhsx(:,:,jz)+tadv2*rhsx_f(:,:,jz)+pfx) 
      v(:,:,jz) = v(:,:,jz)+dt*(tadv1*rhsy(:,:,jz)+tadv2*rhsy_f(:,:,jz)+pfy)
      w(:,:,jz) = w(:,:,jz)+dt*(tadv1*rhsz(:,:,jz)+tadv2*rhsz_f(:,:,jz))
   enddo



   
 do mmm=1,1
   !impose zero "intermediate (in a fractional step sense)" velocity at top (jz=nz)
   !marco: for the same reason impose it at the bottom
   if ((.not. use_mpi) .or. (use_mpi .and. coord == nproc-1)) then 
      w(:,:,nz)=0._rprec   
   elseif ((.not. use_mpi) .or. (use_mpi .and. coord == 0)) then 
      w(:,:,1)=0._rprec   
   endif
   
      !poisson eq.
      call press_stag_array(p,dpdx,dpdy)
  
      CALL MPI_SENDRECV (P(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 13,  &
                      P(1, 1, NZ), LD*NY, MPI_RPREC, UP, 13,   &
                      COMM, STATUS, IERR)                  

     CALL MPI_SENDRECV (P(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 16,  &
                      P(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 16,   &
                      COMM, STATUS, IERR)

      !calculates pressure derivative on z
      !dpdz we don't care at jz=1,jz=nz since we impose bc on w
      !$omp parallel do
      do jz=1,nz-1
         kk_main = jz + eta_displ(coord)
         dpdz(:,:,jz) = (1.0_rprec/(get_dz_local(2*kk_main-2))) *(p(:,:,jz)-p(:,:,jz-1))      !w-nodes
      enddo
      !$omp end parallel do

       

      !$omp parallel do
      do jz=1,nz-1
         u(:,:,jz) = u(:,:,jz) + dt*(-tadv1*dpdx(:,:,jz))
         v(:,:,jz) = v(:,:,jz) + dt*(-tadv1*dpdy(:,:,jz)) 
         w(:,:,jz) = w(:,:,jz) + dt*(-tadv1*dpdz(:,:,jz))
      enddo
      !$omp end parallel do

      

   enddo


   !update rhs with pressure gradient for rhs_f
   !$omp parallel do
   do jz=1,nz-1
      rhsx(:,:,jz) = rhsx(:,:,jz)-dpdx(:,:,jz)
      rhsy(:,:,jz) = rhsy(:,:,jz)-dpdy(:,:,jz)
      rhsz(:,:,jz) = rhsz(:,:,jz)-dpdz(:,:,jz)
   enddo
   !$omp end parallel do


   !###################################################################
   !RE-APPLY BC FOR OUTPUT PURPOSES

   !enforce inlet
   if(v_inlet) call set_inlet()
   
   !blayer
   if(v_blayer) call blayer(u,v,w,8.0_rprec)
   
!SYNC VELOCITIES (to put in service routine)
   CALL MPI_SENDRECV (U(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 1,  &
                      U(1, 1, NZ), LD*NY, MPI_RPREC, UP, 1,   &
                      COMM, STATUS, IERR)
   CALL MPI_SENDRECV (V(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 2,  &
                      V(1, 1, NZ), LD*NY, MPI_RPREC, UP, 2,   &
                      COMM, STATUS, IERR)
   CALL MPI_SENDRECV (W(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 3,  &
                      W(1, 1, NZ), LD*NY, MPI_RPREC, UP, 3,   &
                      COMM, STATUS, IERR)                  
   CALL MPI_SENDRECV (P(1, 1, 1), LD*NY, MPI_RPREC, DOWN, 3,  &
                      P(1, 1, NZ), LD*NY, MPI_RPREC, UP, 3,   &
                      COMM, STATUS, IERR)                  

   CALL MPI_SENDRECV (U(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 4,  &
                      U(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 4,   &
                      COMM, STATUS, IERR)
   CALL MPI_SENDRECV (V(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 5,  &
                      V(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 5,   &
                      COMM, STATUS, IERR)
   CALL MPI_SENDRECV (W(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 6,  &
                      W(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 6,   &
                      COMM, STATUS, IERR)
   CALL MPI_SENDRECV (P(1, 1, NZ-1), LD*NY, MPI_RPREC, UP, 6,  &
                      P(1, 1, 0), LD*NY, MPI_RPREC, DOWN, 6,   &
                      COMM, STATUS, IERR)


   !OUTPUT
   T3= MPI_WTIME(); CALL OUTPUT_LOOP(); T4= MPI_WTIME(); CALL MPI_BARRIER(COMM,IERR)
   T5= MPI_WTIME(); CALL OUTLOOP_SPC(); T6= MPI_WTIME(); CALL MPI_BARRIER(COMM,IERR)


!   if(jt .ge. 1000) then
!   call mpi_barrier(MPI_COMM_WORLD,ierr)
!    if(second_rank .eq. 1) then
!      call mpi_send(u,size(u),mpi_rprec,0,1001,second_comm,ierr)
!    else
 !     call mpi_recv(tmp_for_merge,size(tmp_for_merge),mpi_rprec,1,1001,second_comm,status,ierr)
 !     do xxx=ceiling(0.9_rprec*real(nx,rprec))+1,nx
 !        u(xxx,:,:) = weight_mmatrix(xxx)*tmp_for_merge(xxx,:,:) + (1.0_rprec-weight_mmatrix(xxx))*u(ceiling(0.9_rprec*real(nx,rprec)),:,:)
 !     enddo
 !     !u = tmp_for_merge
  !  endif
  !   call mpi_barrier(MPI_COMM_WORLD,ierr)
  !  tmp_for_merge = 0.0_rprec
  !
  !  if(second_rank .eq. 1) then
 !     call mpi_send(v,size(v),mpi_rprec,0,1002,second_comm,ierr)
 !   else
 !     call mpi_recv(tmp_for_merge,size(tmp_for_merge),mpi_rprec,1,1002,second_comm,status,ierr)
 !     do xxx=ceiling(0.9_rprec*real(nx,rprec))+1,nx
 !        v(xxx,:,:) = weight_mmatrix(xxx)*tmp_for_merge(xxx,:,:) + (1.0_rprec-weight_mmatrix(xxx))*v(ceiling(0.9_rprec*real(nx,rprec)),:,:)
 !     enddo
 !     !u = tmp_for_merge
 !   endif
 !   call mpi_barrier(MPI_COMM_WORLD,ierr)
 !   tmp_for_merge = 0.0_rprec
!
!    if(second_rank .eq. 1) then
!      call mpi_send(w,size(w),mpi_rprec,0,1003,second_comm,ierr)
!    else
!      call mpi_recv(tmp_for_merge,size(tmp_for_merge),mpi_rprec,1,1003,second_comm,status,ierr)
 !     do xxx=ceiling(0.9_rprec*real(nx,rprec))+1,nx
 !        w(xxx,:,:) = weight_mmatrix(xxx)*tmp_for_merge(xxx,:,:) + (1.0_rprec-weight_mmatrix(xxx))*w(ceiling(0.9_rprec*real(nx,rprec)),:,:)
 !     enddo
 !     !u = tmp_for_merge
 !   endif
 !   call mpi_barrier(MPI_COMM_WORLD,ierr)
 !   tmp_for_merge=0.0_rprec
 !
 !   !call mpi_barrier(MPI_COMM_WORLD,ierr)
 !
    !if(temperature) then
    ! if(second_rank .eq. 1) then
    !   call mpi_send(SC1,size(SC1),mpi_rprec,0,1004,second_comm,ierr)
    ! else
    !  call mpi_recv(tmp_for_merge,size(tmp_for_merge),mpi_rprec,1,1004,second_comm,status,ierr)
    !  do xxx=ceiling(0.9_rprec*real(nx,rprec))+1,nx
    !     SC1(xxx,:,:) = weight_mmatrix(xxx)*tmp_for_merge(xxx,:,:) + (1.0_rprec-weight_mmatrix(xxx))*SC1(ceiling(0.9_rprec*real(nx,rprec)),:,:)
    !  enddo
    !  !u = tmp_for_merge
    ! endif
    ! call mpi_barrier(MPI_COMM_WORLD,ierr)
    ! tmp_for_merge=0.0_rprec
    !endif

    !if(humidity) then
    ! if(second_rank .eq. 1) then
    !   call mpi_send(SC2,size(SC2),mpi_rprec,0,1005,second_comm,ierr)
    ! else
    !  call mpi_recv(tmp_for_merge,size(tmp_for_merge),mpi_rprec,1,1005,second_comm,status,ierr)
    !  do xxx=ceiling(0.9_rprec*real(nx,rprec))+1,nx
    !     SC2(xxx,:,:) = weight_mmatrix(xxx)*tmp_for_merge(xxx,:,:) + (1.0_rprec-weight_mmatrix(xxx))*SC2(ceiling(0.9_rprec*real(nx,rprec)),:,:)
    !  enddo
    !  !u = tmp_for_merge
    ! endif
    ! call mpi_barrier(MPI_COMM_WORLD,ierr)
    ! tmp_for_merge=0.0_rprec
    !endif

    !if(humidity) then
    ! if(second_rank .eq. 1) then
    !   call mpi_send(SC3,size(SC3),mpi_rprec,0,1006,second_comm,ierr)
    ! else
    !  call mpi_recv(tmp_for_merge,size(tmp_for_merge),mpi_rprec,1,1006,second_comm,status,ierr)
    !  do xxx=ceiling(0.9_rprec*real(nx,rprec))+1,nx
    !     SC3(xxx,:,:) = weight_mmatrix(xxx)*tmp_for_merge(xxx,:,:) + (1.0_rprec-weight_mmatrix(xxx))*SC3(ceiling(0.9_rprec*real(nx,rprec)),:,:)
    !  enddo
    !  !u = tmp_for_merge
    ! endif
    ! call mpi_barrier(MPI_COMM_WORLD,ierr)
 !   ! tmp_for_merge=0.0_rprec
 !   !endif
 !  
 !  endif
 !
ENDDO

!save time
t2= mpi_wtime()
if(coord==0)then
print*,'total loop time =', t2-t1
print*,'output time = ', t4-t3
print*,'output_spc time = ', t6-t5
endif
call mpi_barrier(comm,ierr)

!DEALLOCATION
CALL DEALLOC_SIM_PARAM()
CALL DEALLOCATE_SGS()
IF(TEMPERATURE) CALL deAllocTemperature()
IF(HUMIDITY) CALL deAllocHumidity()
!if(humidity) call deAllocHumidity_particle()
IF(LSM) CALL DEALLOCATE_LSM()
call mpi_finalize(ierr)


END PROGRAM MAIN

