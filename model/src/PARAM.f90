!---------------------------------------------------------------------------------
!MODULE PARAM
!---------------------------------------------------------------------------------
!
! DESCRIPTION
!
! V_INIT1 = flat velocity profile
! V_INIT2 = logarithmic velocity profile
! V_INIT3 = coriolis velocity
! V_NOISE = add noise to initial velocity profile (WHAT NOISE?)
! T_INIT1 = flat temperature profile
! T_INIT2 = stably stratified (constant gradient of theta)
!
!
!---------------------------------------------------------------------------------
MODULE PARAM

USE TYPES,ONLY:RPREC
#IF(MPI)
  use mpi
#ENDIF

IMPLICIT NONE

!--mpi stuff
#IF(MPI)
  logical,parameter :: USE_MPI = .true.
  integer,parameter :: nproc = 16
#ELSE
  logical,parameter :: USE_MPI = .false.
  integer,parameter :: nproc = 1  !--this must be 1 if no MPI
#ENDIF
#IF(MPI)
  integer :: status(MPI_STATUS_SIZE)
#ENDIF
!--this stuff must be defined, even if not using MPI
character (8) :: chcoord  !--holds character representation of coord
integer :: ierr
integer :: comm
integer :: up, down
integer :: global_rank
integer :: MPI_RPREC, MPI_CPREC, MPI_INTPREC
integer :: rank = 0   !--init to bogus (so its defined, even if no MPI)
integer :: coord = 0  !--same here
integer :: rank_of_coord(0:nproc-1), coord_of_rank(0:nproc-1)
integer :: colour
integer :: new_comm
!--end mpi stuff

!NEEDED
REAL(RPREC),PARAMETER :: PI=3.1415926535897932384626433_RPREC

!######################################################################
!VERBOSE OUTPUT
LOGICAL :: VERBOSE,PROFILING
LOGICAL,PARAMETER :: DEBUG_DDZ = .FALSE.,           &
                     VERBOSE_DTM_DDFA = .FALSE.

!######################################################################
!SIMULATION PARAMETERS
REAL(RPREC) :: RE,U_STAR,                &  !U_STAR IS MAINLY FOR IC
               RI_SC1,PR_SC1,PR_SGS_SC1, &
               RI_SC2,PR_SC2,PR_SGS_SC2, &
               RI_SC3,PR_SC3,PR_SGS_SC3

REAL(RPREC),PARAMETER :: VONK=0.4_RPREC 

!IMPORTANT FLAGS:
LOGICAL :: SGS,TEMPERATURE,HUMIDITY

LOGICAL :: noise_pressure

!SLOPE
LOGICAL ::     SLOPE
REAL(RPREC) :: ALPHA       !(DEG)

!!FORCING
REAL(RPREC) :: PFX,PFY

!IBM
INTEGER ::     NITER_F     !DFA ITERS

!FSI
LOGICAL :: ICM_DFA
LOGICAL,PARAMETER ::     IFEM=.FALSE.

!DYNAMIC ROUGHNESS
LOGICAL :: LDSRM

! Molecular Viscosity
REAL(RPREC),PARAMETER :: NU_MOLEC = 1.24036021044E-5
!######################################################################
!GEOMETRY AND SPATIAL DISCRETIZATION
REAL(RPREC) :: L_X,L_Y,L_Z,L_eta
!INTEGER,PARAMETER ::  NX=64,NY=64,NZ_TOT=129
INTEGER,PARAMETER :: NX=32,NY=32,NZ_TOT=33 !129 ! If SPONGE_TYPE = 1, then NZ_TOT-1 should be multiple of 4
!INTEGER,PARAMETER ::  NX=256,NY=256,NZ_TOT=2049

INTEGER            ::  NZ
INTEGER,DIMENSION(0:nproc-1),PARAMETER :: NZ_VECTOR = 1 + SPREAD(2,1,16) !1 + (/1,2,9,SPREAD(36,1,5)/) !1 + (/1,2,18,SPREAD(19,1,9)/) !1 + SPREAD(16,1,12)
INTEGER,DIMENSION(0:nproc-1)           :: NZ_CUMUL, ETA_DISPL
INTEGER,PARAMETER :: NZ_MAX = MAXVAL(NZ_VECTOR)

REAL(RPREC),PARAMETER :: Z_I=1.0_rprec
REAL(RPREC),PARAMETER :: T_SCALE = 283.15_rprec !263.15_rprec 

!! FOR STRETCHING GRID
!INTEGER :: nz_salt,nz_middle,nz_upper
real(rprec) :: a_grid,b_grid,c_grid,f_grid
real(rprec) :: c_factor,f_factor

real(rprec),dimension(:),allocatable :: get_z_local,get_dz_local
!######################################################################
!OUTPUT PARAMETERS
INTEGER :: WBASE,C_CNT1,C_CNT2,P_CNT3,C_CNT3,C_CNT4,inlet_count,C_CNT5
LOGICAL :: HAVG_CNT3
REAL(RPREC) :: Z_EDDY_COV

!PARAMETERS FOR INST_FIELD
INTEGER :: NXF1,NXF2,NYF1,NYF2,NZF1,NZF2,SX,SY,SZ

!OUTMODES
!CHARACTER(12),PARAMETER :: OUTMODE_FIELD='BINARY'
!CHARACTER(12),PARAMETER :: OUT_MODE="ASCII"

!######################################################################
!INITIAL AND BOUNDARY CONDITIONS ON VELOCITY FIELD    (MAKE LIKE THOSE FOR SCALARS MORE CLEAN!)

!I.C.
LOGICAL :: VI1,VI2,VI3,VI4,VI5,VI6
REAL(RPREC) :: U0,V0,W0,DUDZ0,VNF,Z_TURB

!B.C. 
CHARACTER(25) :: LLBC,LUBC,LLBC_SPECIAL,LUBC_SPECIAL
CHARACTER(25) :: LLBC_U,LLBC_V,LLBC_W,LUBC_U,LUBC_V,LUBC_W
REAL(RPREC)  :: DLBC_U,DLBC_V,DLBC_W,DUBC_U,DUBC_V,DUBC_W

!ROUGHNESS LENGTH
REAL(RPREC) :: LBC_DZ0,DSPL_H,GT_U,GT_V,WALL_C

!SPECIFIC SETTINGS
LOGICAL,PARAMETER      :: V_INLET=.FALSE.,V_BLAYER=.FALSE.,V_SPONGE=.false.
CHARACTER(256),PARAMETER :: BC_INLET='LOG'
!REAL(RPREC),PARAMETER   :: Z_REF = 0.15_RPREC,          &   !POWERLAW STUFF
!                           U_REF = 3.0_RPREC,           &
!                           POWER = 0.16_RPREC,          &   
!                           Z_PWR = 1.0_RPREC

!DTM AS LBC
CHARACTER(256) :: DTM_FILE,triInterface
LOGICAL :: LPS
INTEGER :: PHI_ROT,IT_LPS
REAL(RPREC) :: XSHIFT,YSHIFT,ZSHIFT
CHARACTER(256),PARAMETER :: SURFACE = 'SNOW_PRE-POST'      !'BASEL', GLACIER, VNCV_SUNSET



!######################################################################
!INITIAL AND BOUNDARY CONDITIONS ON SC1 FIELD

REAL(RPREC)   :: SC1_INIT,DSC1DZ_INIT,SC1STAR_INIT,NF_SC1,Z_TURB_SC1,    &
                 Z_ENT_SC1,ENT_SC1_INIT,ENT_DSC1DZ_INIT

CHARACTER(3) :: LLBC_SC1,LUBC_SC1          !DRC,NEU,FLX
REAL(RPREC)  :: DLBC_SC1,DUBC_SC1                   
                           
!######################################################################
!INITIAL AND BOUNDARY CONDITIONS ON SC2 FIELD

REAL(RPREC)   :: SC2_INIT,DSC2DZ_INIT,SC2STAR_INIT,NF_SC2,Z_TURB_SC2,    &
                 Z_ENT_SC2,ENT_SC2_INIT,ENT_DSC2DZ_INIT

CHARACTER(3) :: LLBC_SC2,LUBC_SC2          !DRC,NEU,FLX
REAL(RPREC)  :: DLBC_SC2,DUBC_SC2                   

!######################################################################
!INITIAL AND BOUNDARY CONDITIONS ON SC3 FIELD

REAL(RPREC)   :: SC3_INIT,DSC3DZ_INIT,NF_SC3,Z_TURB_SC3,    &
                 Z_ENT_SC3,ENT_SC3_INIT,ENT_DSC3DZ_INIT

CHARACTER(3) :: LLBC_SC3,LUBC_SC3          !DRC,NEU,FLX
REAL(RPREC)  :: DLBC_SC3,DUBC_SC3                   



!######################################################################
!SPONGE LAYER
LOGICAL,PARAMETER       :: T_SPONGE=.false.
INTEGER,PARAMETER       :: SPONGE_TYPE=2
REAL(RPREC),PARAMETER   :: SPONGE_H=0.33334


!!######################################################################
!SUBGRID SETTINGS

INTEGER :: MODEL,IFILTER,CS_COUNT 

!DYN FOR MOMENTUM AND SCALARS
INTEGER :: DYN_INIT
INTEGER,PARAMETER :: SC_A_INIT=20000   !TO PUT IN PARAMETER WITH SCALARS
!MASON WALL DAMPING
REAL(RPREC) :: CO
INTEGER,PARAMETER ::     NNN=2
LOGICAL ::     WALL_DMP

!BOTTOM BC
LOGICAL,PARAMETER ::       WALL_HOMOG_FLAG=.FALSE.,   &
                           WALL_INST_FLAG=.FALSE.      !AVG RATHER THAN FILT
INTEGER,PARAMETER ::       U_STAR_IT=5             !PICARD FOR USTAR

LOGICAL,PARAMETER ::       CLIP_SGS=.TRUE.        !FOR MODEL 5 WHEN 0 GRADIENTS


!---------------xxxxxxxxx--Initialization--xxxxxxxxx---------------
!--initu = true to read from a file; false to create with random noise
!--initlag = true to initialize cs, FLM & FMM; false to read from vel.out
!--initlag must be true for the first time run, i.e., no vel.out yet.
!------------------------------------------------------------------
!--initu=.FALSE. & S_FLAG=.TRUE. initialize velocity & scalar fields using ran
!--initu=.FALSE. & S_FLAG=.FALSE. initialize velocity fields using ran
!------------------------------------------------------------------
!logical,parameter:: initu=.FALSE.,initsc=.FALSE.,inilag=.TRUE.,interp=.FALSE.
LOGICAL,PARAMETER :: INITU=.FALSE., &
                     INILAG=.TRUE.

!######################################################################
!SLOPE FLOW SETTINGS
LOGICAL,PARAMETER :: sf_norm_1 = .false.


!######################################################################
!WORKING FOLDERS
CHARACTER(*),PARAMETER :: out_avgs        = './output/out_avgs/',       &
                          dsr_path        = './output/dsrm/',           &
                          surf_path       = './output/surf_values/',    &
                          tower_out       = './output/tower_values/',   &
                          rst_path        = './output/rst_files/',      &
                          phi_path        = './output/phi_functions/'

CHARACTER(256) :: out_path
                          
!######################################################################
!EXTRA PARAMETERS FOR SPACE-TIME DISCR

!######################################################################

!STEPS AND DTS
REAL(RPREC) :: REF_DT
INTEGER ::     NSTEPS
  
!TRANSIENT STEPS AND DT
INTEGER::      TR_STEPS
REAL(RPREC) :: TR_DT_FACTOR

INTEGER:: JT,JT_TOTAL    !TO PUT IN SIM_PARAM EVENTUALLY

INTEGER,PARAMETER:: NX2=3*NX/2,NY2=3*NY/2
INTEGER,PARAMETER:: LH=NX/2+1,LD=2*LH,LH_BIG=NX2/2+1,LD_BIG=2*LH_BIG

!SIZE OF THE GRID
REAL(RPREC) :: DX,DY,DZ,C_STRETCH

!REAL(RPREC),DIMENSION(0:NZ_TOT) :: DZ

!TIME ADVANCE PARAMETERS (AB2)
REAL(RPREC):: TADV1,TADV2

!GRAVITY AND DENSITY
REAL(RPREC),PARAMETER::    G=9.81_RPREC
REAL(RPREC),PARAMETER ::   F_RHO=1.18_rprec !1.341305073152195_rprec         !1.0_RPREC

!VARIOUS STUFF (MPI AND CONSTANTS)
INTEGER,PARAMETER ::     IBOGUS=-1234567890
REAL(RPREC),PARAMETER :: BOGUS=-1234567890._RPREC

!HOW WE PERFORM AVERAGES
!INTEGER,PARAMETER ::     SC1_AVGTYPE = 3         !1=PLANAR AVG, 2=VOLUME AVG, 3=REF TEMP
!REAL(RPREC),PARAMETER :: SC1_REF = 9.81_RPREC   !REF FOR BUOYANCY(IF=9.81 WE'RE SOLVING FOR BUOY)
!REAL(RPREC),PARAMETER :: BETA_G1 = 1.0_RPREC     !FILTER SIZE FOR VOLUME AVG (HORIZ WE USE FFT)
!INTEGER,PARAMETER ::     IFILTER_T = 1           !1->cut off 2->Gaussian 3->Top-hat


!######################################################################
!KATABATIC FLOW SIMULATION IC_SCALAR+MOMENTUM   !put in normal bc as INIT_4
!LOGICAL,PARAMETER ::     KAT_INIT=.FALSE.          !PRANDTL ANALITHICAL SOLUTION AS IC
!REAL(RPREC),PARAMETER :: PR_KAT_IC=2.0_RPREC,   &
!                         KB=0.02_RPREC,         &
!                         KM=0.02_RPREC,         &
!                         BS=-0.29_RPREC,        &
!                         MAG_FACT=0.3_RPREC


!######################################################################
! VANCOUVER PARAMETERS
LOGICAL,PARAMETER :: vncv_sim = .FALSE.

!######################################################################
! Geostrophic Wind Parameters
REAL(RPREC):: coriol 
REAL(RPREC):: ug 
REAL(RPREC):: vg
!######################################################################

! LSM VARIABLE
LOGICAL :: LSM
INTEGER :: NP_MAX, AVERAGE_FREQ, PRINT1_FREQ, LSM_INIT, NBINS
REAL(RPREC) :: P_RHO, D_M, D_S, D_MIN, D_MAX, PPP_MIN, PPP_MAX, B_ENE, INIT_DEP, SETTLING_VEL
INTEGER :: LSM_FREQ
CHARACTER(50) :: MODE, FILE_DEP, P_TYPE
LOGICAL :: R_DEP,OUTLET
CHARACTER(256) :: lsm_path ! = './output/lsm/'
REAL(RPREC) :: molec_weight_water,R_universal,surface_tension,water_rho,c_p_a,c_L,c_p_v
LOGICAL :: sublimation,Schiller,Mason,thorpe_mason
LOGICAL :: lp_conc,lp_flux,lp_stress,lp_diss,lp_massbal,lp_restart,lp_mombal

REAL(RPREC),dimension(:,:,:), allocatable :: diss,restke,sgstke

! FOR NEW MPI ADDITIONS
real(rprec),dimension(nx) :: weight_mmatrix
!real(rprec),dimension(ld*ny*(nz+1)) :: flat_weight
!real(rprec),dimension(ld,ny,0:NZ) :: tmp_for_merge
integer :: xxx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#######################################
real(rprec) :: Ct,Cp,local_tip_speed_ratio,wt_averaging_time
logical :: wind_farm,wind_farm_debug
logical :: staggered
integer :: num_of_rows,num_of_columns,spacing_x,spacing_y
integer :: start_disk
!!!!!##################################################

!!#####################################################
!additional parameters for heterogeneous surfaces (MH)
integer :: MODE_0, MODE_1
real(rprec) :: T_0, T_1, Z0_0, Z0_1

!!#####################################################


END MODULE PARAM
