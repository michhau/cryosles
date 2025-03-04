MODULE SGS_STRESSES

USE PARAM
USE TYPES
USE SIM_PARAM
USE SERVICE
use stretch

!-----------------------------------------------------------------------------------
CONTAINS   


!-----------------------------------------------------------------------------------
!  OBJECT: subroutine SGS_STRESSES ()
!!----------------------------------------------------------------------------------
!!
!!  LAST MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 16/05/2013
!!
!!  DESCRIPTION:
!!
!!
!
!!  Through the calc_Sij subroutine: 
!!
!!      - if "wall" then puts all the stresses S11,S12,S13,S23,S22,S33 
!!      are put in the first set of UVP nodes (at jz=1);
!!
!!      - otherwise, and for all the remaining nodes, it puts the stresses: 
!!      S11,S12,S13,S23,S22,S33 to the W nodes;
!!
!!  If dynamic:
!!
!!      - fixes initial cs=0.03 if jt=1 and inilag=true;
!!      - then for jt.GE.DYN_init).OR.(initu)) .AND.  &  (mod(jt,cs_count) == 0
!!      it calls the sgs calculation to update the Cs_coefficient;
!!      - every sgs literally computes the new Cs which is filled in the domain;
!!      - Nu_T is calculated ( c_s^2 l^2 |S| )
!!      - |S| is calculated
!!      - for the wall: computes txx,txy,tyy,tzz as -2 Nu_T S_ij (txz & tyz at wall set by wallstress routine)
!!      - for stress free: computes txx,txy,tyy,tzz as -2 Nu_T(interp) S_ij (interp)(txz & tyz are put = 0 by wallstress routine)
!!      - for 
!!
!!
!!      NOTE: basically we put all stuff on uvp for jz=1 and w for the rest before sending
!!            values to sgs subroutines, then, once the smagorinsky coefficient is computed
!!            we bring back all stuff to their normal position in the staggered fashion.
!!            Namely after this subroutine we have all stresses in their proper position.
!!
!!----------------------------------------------------------------------------------
SUBROUTINE SGS_STRESS (PARAM)
USE SGSMODULE,ONLY: U_LAG,V_LAG,W_LAG,NU_T,CS_OPT2
USE SIM_PARAM
USE TEST_FILTERMODULE

IMPLICIT NONE

REAL(RPREC),INTENT(IN) :: PARAM
                                                 
REAL(RPREC),DIMENSION(NZ) :: CS_TMP
REAL(RPREC),DIMENSION(LD,NY,NZ) :: S11,S12,S22,S33,S13,S23
REAL(RPREC) :: STMP
INTEGER::JX,JY,JZ,JZMIN
REAL(RPREC),DIMENSION(LD,NY) :: CONST,MLC_NU
REAL(RPREC),SAVE :: LAGR_DT=0.0_RPREC
LOGICAL :: CHECK=.TRUE.
REAL(RPREC),DIMENSION(NZ) :: L
REAL(RPREC) :: ZZ,ZZ_LOCAL

IF(VERBOSE) PRINT*,"STARTED SGS_STRESSES"

!###########################################################################
!FILTER SIZE DEFINES MIXING LENGTH
!L=FILTER_SIZE*(DX*DY*DZ)**(1._RPREC/3._RPREC)       !NONDIMENSIONAL

DO JZ=1,NZ
            IF((COORD==0).and.(JZ == 1))THEN
               ZZ=(JZ-0.5_RPREC)
               ZZ_LOCAL = get_z_local(int(2*zz)) ! AS: I assume ZZ_LOCAL should be the height of the first uvp level above surface here because integral scale is on the order of the distance to the wall (Bou-Zeid et al., 2005);
               ! so the real zz=0.5 should first be multiplied by 2 before converting to an integer !c_stretch*(exp(zz)-1.0_rprec)
             ELSE  ! W-NODES
               ZZ= (JZ-1) + ETA_DISPL(COORD) ! ETA_DISPL(COORD) is the height of the bottom of the subdomain in index units
               ZZ_LOCAL = get_z_local(int(2*zz)+1) - get_z_local(int(2*zz)-1) !c_stretch*(exp(zz)-1.0_rprec)
            ENDIF
            L(JZ)=FILTER_SIZE*(DX*DY*(ZZ_LOCAL))**(1._RPREC/3._RPREC)
ENDDO


S11=0.0_rprec
S22=0.0_rprec
S33=0.0_rprec
S12=0.0_rprec
S13=0.0_rprec
     
!SIJ ON UVP NODE AT JZ=1, ON W NODES FOR THE OTHERS (JZ>1)
CALL CALC_SIJ(DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ,           &
              S11,S12,S13,S23,S22,S33)

IF(JT .EQ. 1) then
 U_LAG = 0.0_rprec
 V_LAG = 0.0_rprec
 W_LAG = 0.0_rprec
endif              
!AVERAGE VEL_LAG/THETA_LAG DURING CS
IF(JT_TOTAL.GE.DYN_INIT)THEN
   IF(MODEL==5 .OR. MODEL==7)THEN

      !$OMP PARALLEL DO
      DO JZ=0,NZ
         U_LAG(:,:,JZ) = U_LAG(:,:,JZ)+U(:,:,JZ)
         V_LAG(:,:,JZ) = V_LAG(:,:,JZ)+V(:,:,JZ)
         W_LAG(:,:,JZ) = W_LAG(:,:,JZ)+W(:,:,JZ)
      ENDDO
      !$OMP END PARALLEL DO

      !INCREASE LAGR_DT
      LAGR_DT=LAGR_DT+DT

      !AT DYN_INIT DON'T HAVE YET CS_COUNT VALUES
      IF((JT_TOTAL==DYN_INIT))THEN
         LAGR_DT=LAGR_DT*CS_COUNT
         U_LAG=U_LAG*CS_COUNT
         V_LAG=V_LAG*CS_COUNT
         W_LAG=W_LAG*CS_COUNT
      ENDIF

   ENDIF
ENDIF

!###########################################################################
!###########################################################################
!SMAGORINSKY COEFFICIENT


!SMAGORINSKY FOR JT=1 IF NO RESTART
IF((JT_TOTAL.EQ.1) .OR. MODEL==1)THEN
      
   !SCREEN PRINT
   IF(JT_TOTAL.EQ.DYN_INIT) PRINT*,'RUNNING IN STATIC MODE MODEL = ',MODEL
   IF((MODEL==1).AND.(JT_TOTAL==1)) PRINT*,'RUNNING MODEL = ',MODEL

   !DEFINE SMAG. COEFF
   CS_OPT2(:,:,:)=CO**2      
      
   !DAMPING FUNCTIONS (NOT CHECKED WITH MPI)
   IF(WALL_DMP)THEN
      
      !IF(LLBC_SPECIAL.EQ.'IBM_DFA')THEN 
      !   CS_OPT2=(CO**(-NNN)+(DELTA/VONK/(PHI(:,:,1:NZ)+LBC_DZ0))**NNN)**(-2._RPREC/NNN)
      !ELSE
         DO JZ=1,NZ
            IF((COORD==0).and.(JZ == 1))THEN
               ZZ=(JZ-0.5_RPREC)
               ZZ_LOCAL = get_z_local(int(2*zz)) !c_stretch*(exp(zz)-1.0_rprec)
             ELSE  ! W-NODES
               ZZ=(JZ-1) + ETA_DISPL(COORD)
               ZZ_LOCAL = get_z_local(int(2*zz)+1) - get_z_local(int(2*zz)-1) !get_dz_local(real(zz,rprec)) !c_stretch*(exp(zz)-1.0_rprec)
            ENDIF
            L(JZ)=FILTER_SIZE*(DX*DY*(ZZ_LOCAL))**(1._RPREC/3._RPREC)
            CS_OPT2(:,:,JZ)=(CO**(-NNN)+(L(JZ)/VONK/(get_z_local(int(2*zz))+LBC_DZ0))**NNN)**(-2._RPREC/NNN)
         ENDDO
      !ENDIF
   ENDIF
ENDIF

IF((JT_TOTAL==DYN_INIT).OR.((JT_TOTAL.GT.DYN_INIT).AND.(MOD(JT_TOTAL,CS_COUNT).EQ.0))) THEN
   
   !SCREEN PRINT
   IF((JT_TOTAL.GT.DYN_INIT) .and. (JT .eq. 1)) PRINT*,'RUNNING DYNAMIC MODEL = ',MODEL
      
   !DEFINE SMAG. COEFF
   IF((MODEL==2).OR.(MODEL==6))THEN
      CALL STD_DYNAMIC(CS_TMP,S11,S12,S13,S22,S23,S33)
      DO JZ=1,NZ
         CS_OPT2(:,:,JZ)=CS_TMP(JZ)
      ENDDO
   ELSEIF((MODEL==5).OR.(MODEL==7))THEN
      CALL LAGRANGE_SDEP(S11,S12,S13,S22,S23,S33,LAGR_DT)
      LAGR_DT=0.0_RPREC
   ENDIF
   
ENDIF

 
!###########################################################################
!###########################################################################
!EDDY VISCOSITY
DO JZ=1,NZ
DO JY=1,NY
DO JX=1,NX
   STMP = SQRT(2._RPREC*(S11(JX,JY,JZ)**2.0_rprec + S22(JX,JY,JZ)**2.0_rprec +                &
                         S33(JX,JY,JZ)**2.0_rprec + 2._RPREC*(S12(JX,JY,JZ)**2.0_rprec +      &
                         S13(JX,JY,JZ)**2.0_rprec + S23(JX,JY,JZ)**2.0_rprec)))
   !SMAGORINSKY HYPOTHESIS: EVALUATE NU_T = C_S^2 L^2 |S|
   NU_T(JX,JY,JZ)=STMP*CS_OPT2(JX,JY,JZ)*L(JZ)**2
ENDDO
ENDDO
ENDDO

!###########################################################################
!###########################################################################
!SUBGRID STRESSES
!
!NOTE THAT FROM NOW ON EVERYTHING IS BACK TO ITS ORIGINAL STAGGERED POS
!
!CALL Z_MPROFILE1(CS_OPT2,'CS',1.0_RPREC,1)
!CALL Z_MPROFILE(NU_T,'NU_T',1.0_RPREC,1)
!CALL Z_MPROFILE1(S13,'S13',1.0_RPREC)
!CALL Z_MPROFILE1(S23,'S23',1.0_RPREC)

!PRINT*,'MAXVAL(NU_T)=',MAXVAL(NU_T)

!TEST
!NU_T=0.0_RPREC


!###########################################################################
!JZ=1 SPECIAL TREATMENT 
MLC_NU(:,:)= nu_molec

!COMPUTE STRESSES
IF(COORD==0)THEN
   TXX(:,:,1) = 2._RPREC*(NU_T(:,:,1)+MLC_NU(:,:))*DUDX(:,:,1)    !UVP
   TXY(:,:,1) = 2._RPREC*(NU_T(:,:,1)+MLC_NU(:,:))*0.5_RPREC*(DUDY(:,:,1)+DVDX(:,:,1))    !UVP
   TYY(:,:,1) = 2._RPREC*(NU_T(:,:,1)+MLC_NU(:,:))*DVDY(:,:,1)    !UVP
   TZZ(:,:,1) = 2._RPREC*(NU_T(:,:,1)+MLC_NU(:,:))*DWDZ(:,:,1)    !UVP

   !IF THIS PART HAS NOT BEEN COMPUTED (BY WALL_LAW), PROVIDE DNS VALUES
   if(llbc_special .ne. 'WALL_LAW') then
      TXZ(:,:,1)=PARAM*(DUDZ(:,:,1)+DWDX(:,:,1))       !W
      TYZ(:,:,1)=PARAM*(DVDZ(:,:,1)+DWDY(:,:,1))       !W
   endif

   JZMIN=2
ELSE
   JZMIN=1
ENDIF


!###########################################################################
!STANDARD PROCEDURE

!STUFF ON UVP AND THEN ON W
DO JZ=JZMIN,NZ-1
   
   CONST(:,:) =0.5_RPREC*(NU_T(:,:,JZ)+NU_T(:,:,JZ+1))            !UVP
   TXX(:,:,JZ)=2._RPREC*(CONST+MLC_NU(:,:))*DUDX(:,:,JZ)          !UVP
   TXY(:,:,JZ)=2._RPREC*(CONST+MLC_NU(:,:))   &
                        *0.5_RPREC*(DUDY(:,:,JZ)+DVDX(:,:,JZ))    !UVP
   TYY(:,:,JZ)=2._RPREC*(CONST+MLC_NU(:,:))*DVDY(:,:,JZ)          !UVP
   TZZ(:,:,JZ)=2._RPREC*(CONST+MLC_NU(:,:))*DWDZ(:,:,JZ)          !UVP

ENDDO
DO JZ=JZMIN,NZ

   !STRESSES
   TXZ(:,:,JZ)=2.0_RPREC*(NU_T(:,:,JZ)+MLC_NU(:,:))*S13(:,:,JZ)   !W
   TYZ(:,:,JZ)=2.0_RPREC*(NU_T(:,:,JZ)+MLC_NU(:,:))*S23(:,:,JZ)   !W

ENDDO

!UNA TANTUM CHECK
IF(CHECK.AND.COORD==0)THEN
   txx(:,:,0)=BOGUS        !cause first uvp is outside domain
   txy(:,:,0)=BOGUS        !cause first uvp is outside domain
   tyy(:,:,0)=BOGUS        !cause first uvp is outside domain
   tzz(:,:,0)=BOGUS        !cause first uvp is outside domain
   txz(:,:,0)=BOGUS        !cause first W is outside domain
   tyz(:,:,0)=BOGUS        !cause first W is outside domain
   CHECK=.FALSE.
ELSEIF(CHECK.AND.COORD==NPROC-1)THEN
   txx(:,:,nz)=BOGUS        !cause last uvp is outside domain
   txy(:,:,nz)=BOGUS        !cause last uvp is outside domain
   tyy(:,:,nz)=BOGUS        !cause last uvp is outside domain
   tzz(:,:,nz)=BOGUS        !cause last uvp is outside domain
   CHECK=.FALSE.
ENDIF


#IF(MPI)
  !JF: can avoid if we change nz-1 to nz above
  !call mpi_sendrecv (txz(1, 1, 1), ld*ny, MPI_RPREC, down, 3,  &
  !                   txz(1, 1, nz), ld*ny, MPI_RPREC, up, 3,   &
  !                   comm, status, ierr)
  !call mpi_sendrecv (tyz(1, 1, 1), ld*ny, MPI_RPREC, down, 4,  &
  !                   tyz(1, 1, nz), ld*ny, MPI_RPREC, up, 4,   &
  !                   comm, status, ierr)

  !--also set 0-layer to bogus values
  !txx(:, :, 0) = BOGUS
  !txy(:, :, 0) = BOGUS
  !txz(:, :, 0) = BOGUS
  !tyy(:, :, 0) = BOGUS
  !tyz(:, :, 0) = BOGUS
  !tzz(:, :, 0) = BOGUS

  !--send stuff up to ghost nodes for tzz
  !--moved to here from main
  call mpi_sendrecv (tzz(:, :, nz-1), ld*ny, MPI_RPREC, up, 8,  &
                     tzz(:, :, 0), ld*ny, MPI_RPREC, down, 8, &
                     comm, status, ierr)
#ENDIF


IF(VERBOSE) PRINT*,"FINISHED SGS_STRESSES"

END SUBROUTINE SGS_STRESS


!!------------------------------------------------------------------------------
!!  OBJECT: subroutine CALC_SIJ()
!!------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ) on 27/03/2013
!!
!!  DESCRIPTION:
!!
!!  COMPUTES SIJ ON UVP NODE AT JZ=1, ON W NODES FOR THE OTHERS (JZ>1)
!!
!!------------------------------------------------------------------------------
SUBROUTINE CALC_SIJ(DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ,           &
                    S11,S12,S13,S23,S22,S33)

IMPLICIT NONE

REAL(RPREC) :: UX,UY,UZ,VX,VY,VZ,WX,WY,WZ
REAL(RPREC),DIMENSION(LD,NY,0:NZ),INTENT(IN) :: DUDX,DUDY,DUDZ,DVDX,DVDY, &
                                                DVDZ,DWDX,DWDY,DWDZ
REAL(RPREC),DIMENSION(LD,NY,NZ),INTENT(OUT) ::  S11,S12,S13,S23,S22,S33
INTEGER :: JX,JY,JZ,jzMin
integer :: kk
real(rprec) :: fac_1,fac_2,fac_3

if(coord==0)then  !trick for MPI

   !first level has everything on uvp nodes
   do jy=1,ny
   do jx=1,nx     

      ux=dudx(jx,jy,1)       !uvp-node @ dz/2
      uy=dudy(jx,jy,1)       !uvp-node @ dz/2
      if(llbc_special.eq.'WALL_LAW')then
         uz=dudz(jx,jy,1)
         vz=dvdz(jx,jy,1)
      else
         uz=0.5_rprec*(dudz(jx,jy,1)+dudz(jx,jy,2))
         vz=0.5_rprec*(dvdz(jx,jy,1)+dvdz(jx,jy,2))
      endif
      vx=dvdx(jx,jy,1)                            !uvp-node @ dz/2
      vy=dvdy(jx,jy,1)                            !uvp-node @ dz/2

      fac_1 = get_z_local(2)
      fac_2 = get_z_local(0)
      fac_3 = get_z_local(1)

      wx = dwdx(jx,jy,1) + ((dwdx(jx,jy,2)-dwdx(jx,jy,1))/(fac_1-fac_2))*(fac_3-fac_2)
      wy = dwdy(jx,jy,1) + ((dwdy(jx,jy,2)-dwdy(jx,jy,1))/(fac_1-fac_2))*(fac_3-fac_2)

      !wx=0.5_rprec*(dwdx(jx,jy,1)+dwdx(jx,jy,2))  !uvp-node so i need to interp
      !wy=0.5_rprec*(dwdy(jx,jy,1)+dwdy(jx,jy,2))  !uvp-node so i need to interp
      wz=dwdz(jx,jy,1)                            !uvp-node so i need to interp
     
      s11(jx,jy,1)=ux                             !uvp-node @ dz/2
      s12(jx,jy,1)=0.5_rprec*(uy+vx)              !uvp-node @ dz/2
      s13(jx,jy,1)=0.5_rprec*(uz+wx)              !uvp-node @ dz/2
      s22(jx,jy,1)=vy                             !uvp-node @ dz/2
      s23(jx,jy,1)=0.5_rprec*(vz+wy)              !uvp-node @ dz/2
      s33(jx,jy,1)=wz                             !uvp-node @ dz/2

   enddo
   enddo
  
   jzMin=2

else
   
   jzMin=1

endif

!necessary cause we pass just u,v,w and not dwdz,
#IF(MPI)
  call mpi_sendrecv (dwdz(1, 1, 1), ld*ny, MPI_RPREC, down, 10,  &
                     dwdz(1, 1, nz), ld*ny, MPI_RPREC, up, 10,   &
                     comm, status, ierr)
#ENDIF

!standard nodes (put everything on w, where cs and nu_t are defined)
do jz=jzMin,nz
do jy=1,ny
do jx=1,nx

   kk=jz + ETA_DISPL(COORD)
   
   fac_1=get_z_local(2*kk-1)
   fac_2=get_z_local(2*kk-3)
   fac_3=get_z_local(2*kk-2)

   ux=dudx(jx,jy,jz-1) + ((dudx(jx,jy,jz)-dudx(jx,jy,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)
   uy=dudy(jx,jy,jz-1) + ((dudy(jx,jy,jz)-dudy(jx,jy,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)
   vx=dvdx(jx,jy,jz-1) + ((dvdx(jx,jy,jz)-dvdx(jx,jy,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)
   vy=dvdy(jx,jy,jz-1) + ((dvdy(jx,jy,jz)-dvdy(jx,jy,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)

   wz=dwdz(jx,jy,jz-1) + ((dwdz(jx,jy,jz)-dwdz(jx,jy,jz-1))/(fac_1-fac_2))*(fac_3-fac_2)




   !ux=0.5_rprec*(dudx(jx,jy,jz) + dudx(jx,jy,jz-1))     !w-node
   !uy=0.5_rprec*(dudy(jx,jy,jz) + dudy(jx,jy,jz-1))     !w-node
   uz=dudz(jx,jy,jz)                                    !w-node
   !vx=0.5_rprec*(dvdx(jx,jy,jz) + dvdx(jx,jy,jz-1))     !w-node
   !vy=0.5_rprec*(dvdy(jx,jy,jz) + dvdy(jx,jy,jz-1))     !w-node
   vz=dvdz(jx,jy,jz)                                    !w-node
   wx=dwdx(jx,jy,jz)                                    !w-node
   wy=dwdy(jx,jy,jz)                                    !w-node
   !wz=0.5_rprec*(dwdz(jx,jy,jz) + dwdz(jx,jy,jz-1))     !w-node

   s11(jx,jy,jz)=ux                                     !w-node
   s12(jx,jy,jz)=0.5_rprec*(uy+vx)                      !w-node
   s13(jx,jy,jz)=0.5_rprec*(uz+wx)                      !w-node
   s22(jx,jy,jz)=vy                                     !w-node
   s23(jx,jy,jz)=0.5_rprec*(vz+wy)                      !w-node
   s33(jx,jy,jz)=wz                                     !w-node
  !if(jz .eq. 3) then  
  ! write(*,*) coord,jx,jy,jz,vz,wy
  !endif
enddo
enddo
enddo


end subroutine calc_sij


end module sgs_stresses





