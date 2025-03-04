MODULE MO_SIMILARITY

USE SIM_PARAM
USE TYPES,ONLY:RPREC
USE SP_SC1
USE SP_SC2
!USE SP_SC3
USE PARAM
use stretch

IMPLICIT NONE

REAL(RPREC),DIMENSION(LD,NY) :: U1,V1,TEMP_SC,U_AVG
REAL(RPREC),DIMENSION(LD,NY) :: SC1_AVG,SC2_AVG !,SC3_AVG
CONTAINS

SUBROUTINE OBUKHOV()

IMPLICIT NONE

INTEGER :: JX,JY
INTEGER :: LOC_T
real(rprec) :: z_local

!COMPUTE U1,V1,U_AVG
!FILTERED WAY
U1=U(:,:,1)
V1=V(:,:,1)
CALL TEST_FILTER(U1,G_TEST)
CALL TEST_FILTER(V1,G_TEST)
DO JX=1,NX
  DO JY=1,NY
     U_AVG(JX,JY)=SQRT(U1(JX,JY)**2+V1(JX,JY)**2)
  ENDDO
ENDDO

IF(TEMPERATURE) then
 TEMP_SC = SC1(:,:,1)
 CALL TEST_FILTER(TEMP_SC,G_TEST)
 SC1_AVG = TEMP_SC
ENDIF

IF(HUMIDITY) then
 TEMP_SC = SC2(:,:,1) !+ SC3(:,:,1)
 CALL TEST_FILTER(TEMP_SC,G_TEST)
 SC2_AVG = TEMP_SC
ENDIF

!IF(HUMIDITY) then
! TEMP_SC = SC3(:,:,1)
! CALL TEST_FILTER(TEMP_SC,G_TEST)
! SC3_AVG = TEMP_SC
!ENDIF


!IF(JT .LT. SC_A_INIT) THEN

  PSI_M=0.0_rprec
  PHI_M=1.0_rprec
  z_local = get_z_local(1) !C_STRETCH * (exp(0.5_rprec*dz)-1.0_rprec)
  USTAR_AVG(:,:) = VONK*U_AVG(:,:)/DLOG(z_local/LBC_Z0(:,:))
   
  IF(TEMPERATURE) THEN
    PSI_H=0.0_rprec
    PHI_H=1.0_rprec
    IF(trim(LLBC_SC1) .eq. 'DRC') then
      SC1_STAR(:,:)=-1.0_rprec * VONK*(DLBC_SC1-SC1_AVG(:,:))/DLOG(z_local/LBC_SC1_ZO(:,:))
    ENDIF
 
    IF(trim(LLBC_SC1) .eq. 'NEU') then
      SC1_STAR(:,:) = DLBC_SC1/USTAR_AVG(:,:)
    ENDIF
  ENDIF

  IF(HUMIDITY) THEN
    PSI_V=0.0_rprec
    PHI_V=1.0_rprec
    IF(trim(LLBC_SC2) .eq. 'DRC') then
      SC2_STAR(:,:)=-1.0_rprec * VONK*(DLBC_SC2-SC2_AVG(:,:))/DLOG(z_local/LBC_SC2_ZO(:,:))
    ENDIF

    IF(trim(LLBC_SC2) .eq. 'NEU') then
      SC2_STAR(:,:) = 0.0_rprec ! DLBC_SC2/USTAR_AVG(:,:)
    ENDIF

  ENDIF
  !IF(HUMIDITY) THEN
  !  PSI_V=0.0_rprec
  !  PHI_V=1.0_rprec
  !  IF(trim(LLBC_SC3) .eq. 'DRC') then
  !    SC3_STAR(:,:)=-1.0_rprec * VONK*(DLBC_SC3-SC3_AVG(:,:))/DLOG(z_local/LBC_SC3_ZO(:,:))
  !  ENDIF
  !
  !  IF(trim(LLBC_SC3) .eq. 'NEU') then
  !    SC3_STAR(:,:) = 0.0_rprec ! DLBC_SC2/USTAR_AVG(:,:)
  !  ENDIF

  !ENDIF

!ELSE
!  DO LOC_T=1,U_STAR_IT
!  
!  z_local = get_z_local(1)
!  USTAR_AVG(:,:) = VONK*U_AVG(:,:)/(DLOG(z_local/LBC_Z0(:,:))+PSI_M(:,:))  
! 
!    IF(TEMPERATURE) THEN
!      IF(trim(LLBC_SC1) .eq. 'DRC') then
!        SC1_STAR(:,:) = -1.0_rprec * VONK*(DLBC_SC1-SC1_AVG(:,:))/(DLOG(z_local/LBC_SC1_ZO(:,:))+PSI_H(:,:))
!      ENDIF
!      IF(trim(LLBC_SC1) .eq. 'NEU') then
!        SC1_STAR(:,:) = DLBC_SC1/USTAR_AVG(:,:)
!      ENDIF
!    ENDIF
!    IF(HUMIDITY) THEN
!      IF(trim(LLBC_SC2) .eq. 'DRC') then
!        SC2_STAR(:,:) = -1.0_rprec * VONK*(DLBC_SC2-SC2_AVG(:,:))/(DLOG(z_local/LBC_SC2_ZO(:,:))+PSI_V(:,:))
!      ENDIF
!      IF(trim(LLBC_SC2) .eq. 'NEU') then
!        SC2_STAR(:,:) = 0.0_rprec !DLBC_SC2/USTAR_AVG(:,:)
!      ENDIF
!    ENDIF
!
!    !IF(HUMIDITY) THEN
!    !  IF(trim(LLBC_SC3) .eq. 'DRC') then
!    !    SC3_STAR(:,:) = -1.0_rprec * VONK*(DLBC_SC3-SC3_AVG(:,:))/(DLOG(z_local/LBC_SC3_ZO(:,:))+PSI_V(:,:))
!    !  ENDIF
!    !  IF(trim(LLBC_SC3) .eq. 'NEU') then
!    !    SC3_STAR(:,:) = 0.0_rprec !DLBC_SC2/USTAR_AVG(:,:)
!    !  ENDIF
!    !ENDIF
!
!    IF(TEMPERATURE) THEN
!      IF(HUMIDITY) THEN
!         CALL L_OBK_T_Q()
!      ELSE
!         CALL L_OBK_T()
!      ENDIF
!    ENDIF
!
!    !CALL COMPUTE_PSI_M(DZ/2.0_rprec,LBC_Z0,PSI_M)
!    IF(TEMPERATURE) THEN
!    CALL COMPUTE_PSI_M(z_local,LBC_Z0,PSI_M)
!      IF(HUMIDITY) THEN
!         CALL COMPUTE_PSI_H(z_local,LBC_SC2_ZO,PSI_V) 
!      ENDIF
!         CALL COMPUTE_PSI_H(z_local,LBC_SC1_ZO,PSI_H)
!    ENDIF
!  ENDDO
!
!! CALL COMPUTE_PHI_M(DZ/2.0_rprec,LBC_Z0,PHI_M)
! IF(TEMPERATURE) THEN
! CALL COMPUTE_PHI_M(z_local,LBC_Z0,PHI_M)
!   IF(HUMIDITY) THEN
!        CALL COMPUTE_PHI_H(z_local,LBC_SC2_ZO,PHI_V)
!   ENDIF
!        CALL COMPUTE_PHI_H(z_local,LBC_SC1_ZO,PHI_H)
! ENDIF
!
!ENDIF

END SUBROUTINE OBUKHOV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE L_OBK_T

IMPLICIT NONE

REAL(RPREC) :: non_dim_grav
REAL(RPREC),dimension(ld,ny) :: denom
integer :: i,j

non_dim_grav = G * (Z_I)

denom = (VONK*non_dim_grav*SC1_STAR)

do i=1,ld
  do j=1,ny
   if(denom(i,j) .ne. 0) then
     OBK_L(i,j)=(SC1_AVG(i,j)*(USTAR_AVG(i,j)**2.0_rprec))/denom(i,j)
   else
     OBK_L(i,j) = 1000000.0_rprec
   endif
  enddo
enddo

END SUBROUTINE L_OBK_T

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE L_OBK_T_Q

IMPLICIT NONE

REAL(RPREC) :: non_dim_grav
REAL(RPREC),dimension(ld,ny) :: denom
integer :: i,j

non_dim_grav = G * (Z_I)


!denom =(VONK*non_dim_grav*(SC1_STAR + 0.61_rprec*SC1_AVG*(SC2_STAR+SC3_STAR)))
denom =(VONK*non_dim_grav*(SC1_STAR + 0.61_rprec*SC1_AVG*(SC2_STAR)))

do i=1,ld
   do j=1,ny
    if(denom(i,j) .ne. 0) then
     OBK_L(i,j) = (SC1_AVG(i,j)*(USTAR_AVG(i,j)**2.0_rprec))/denom(i,j)
    else
      OBK_L(i,j) = 10000000.0_rprec
    endif
   enddo
enddo

END SUBROUTINE L_OBK_T_Q

SUBROUTINE COMPUTE_PSI_M(Z,Z0,PSI)

IMPLICIT NONE

REAL(RPREC),INTENT(IN) :: Z
REAL(RPREC),DIMENSION(LD,NY),INTENT(IN) :: Z0
REAL(RPREC),DIMENSION(LD,NY),INTENT(OUT) :: PSI
REAL(RPREC) :: CSI,CSI_0
INTEGER :: JX,JY
REAL(RPREC) :: a_stable,b_stable
REAL(RPREC) :: a_unstable,b_unstable,c_unstable,d_unstable,n_unstable,y_unstable,x_unstable
real(rprec) :: temp_psi_1,temp_psi_2
real(rprec) :: psi_0
!BUILD VALUE OF ST_CSI

a_stable = 6.0_rprec
b_stable = 2.5_rprec

a_unstable = 0.33_rprec
b_unstable = 0.41_rprec
c_unstable = 0.33_rprec
d_unstable = 0.057_rprec
n_unstable = 0.78_rprec


!BUILD PSI_M TOT: -PSI_M(Z)+PSI_M(Z0)
DO JY=1,NY
 DO JX=1,NX

   csi = z/obk_l(jx,jy)
   csi_0 = z0(jx,jy)/obk_l(jx,jy)

   IF(OBK_L(JX,JY) .ge. 0) then
    ! locally stable
    psi(JX,JY) = -1.0_rprec * (-1.0_rprec*a_stable*dlog(csi + (1+csi**b_stable)**(1.0_rprec/b_stable))) &
                 + (-1.0_rprec*a_stable*dlog(csi_0 + (1+(csi_0)**b_stable)**(1.0_rprec/b_stable)))

   ELSE ! locally unstable
    temp_psi_1 = 0.0_rprec
    y_unstable = -1.0_rprec * csi
    x_unstable = (y_unstable/a_unstable)**(1.0_rprec/3.0_rprec)
    psi_0 = -1.0_rprec * dlog(a_unstable) + 3.0_rprec**(0.5_rprec) * b_unstable * (a_unstable**(1.0_rprec/3.0_rprec))*pi*(1.0_rprec/6.0_rprec)
    if(y_unstable .le. b_unstable**(-3.0_rprec)) then
    temp_psi_1 = dlog(a_unstable+y_unstable) - 3.0_rprec * b_unstable * y_unstable**(1.0_rprec/3.0_rprec) + &
           b_unstable*(a_unstable**(1.0_rprec/3.0_rprec))*0.5_rprec * dlog(((1+x_unstable)**2.0_rprec)/(1.0_rprec - x_unstable + x_unstable**2.0_rprec )) + &
           (3.0_rprec**(0.5_rprec))*b_unstable*(a_unstable**(1.0_rprec/3.0_rprec))*datan((2.0_rprec*x_unstable - 1.0_rprec)/(3.0_rprec**(0.5_rprec))) + psi_0
    else
    y_unstable = b_unstable**(-3.0_rprec)
    temp_psi_1 = dlog(a_unstable+y_unstable) - 3.0_rprec * b_unstable * y_unstable**(1.0_rprec/3.0_rprec) + &
           b_unstable*(a_unstable**(1.0_rprec/3.0_rprec))*0.5_rprec * dlog(((1+x_unstable)**2.0_rprec)/(1.0_rprec - x_unstable + x_unstable**2.0_rprec )) + &
           (3.0_rprec**(0.5_rprec))*b_unstable*(a_unstable**(1.0_rprec/3.0_rprec))*datan((2.0_rprec*x_unstable - 1.0_rprec)/(3.0_rprec**(0.5_rprec))) + psi_0
    endif

    temp_psi_2 = 0.0_rprec
    y_unstable = -1.0_rprec * csi_0
    x_unstable = (y_unstable/a_unstable)**(1.0_rprec/3.0_rprec)
    psi_0 = -1.0_rprec * dlog(a_unstable) + 3.0_rprec**(0.5_rprec) * b_unstable * (a_unstable**(1.0_rprec/3.0_rprec))*pi*(1.0_rprec/6.0_rprec)
    if(y_unstable .le. b_unstable**(-3.0_rprec)) then
    temp_psi_2 = dlog(a_unstable+y_unstable) - 3.0_rprec * b_unstable * y_unstable**(1.0_rprec/3.0_rprec) + &
           b_unstable*(a_unstable**(1.0_rprec/3.0_rprec))*0.5_rprec * dlog(((1+x_unstable)**2.0_rprec)/(1.0_rprec - x_unstable + x_unstable**2.0_rprec )) + &
           (3.0_rprec**(0.5_rprec))*b_unstable*(a_unstable**(1.0_rprec/3.0_rprec))*datan((2.0_rprec*x_unstable - 1.0_rprec)/(3.0_rprec**(0.5_rprec))) + psi_0
    else
    y_unstable = b_unstable**(-3.0_rprec)
    temp_psi_2 = dlog(a_unstable+y_unstable) - 3.0_rprec * b_unstable * y_unstable**(1.0_rprec/3.0_rprec) + &
           b_unstable*(a_unstable**(1.0_rprec/3.0_rprec))*0.5_rprec * dlog(((1+x_unstable)**2.0_rprec)/(1.0_rprec - x_unstable + x_unstable**2.0_rprec )) + &
           (3.0_rprec**(0.5_rprec))*b_unstable*(a_unstable**(1.0_rprec/3.0_rprec))*datan((2.0_rprec*x_unstable - 1.0_rprec)/(3.0_rprec**(0.5_rprec))) + psi_0
    endif

    psi(jx,jy) = -1.0_rprec*temp_psi_1 + temp_psi_2

   ENDIF

 ENDDO
ENDDO

END SUBROUTINE COMPUTE_PSI_M

SUBROUTINE COMPUTE_PSI_H(Z,Z0,PSI)

IMPLICIT NONE

REAL(RPREC),INTENT(IN) :: Z
REAL(RPREC),DIMENSION(LD,NY),INTENT(IN) :: Z0
REAL(RPREC),DIMENSION(LD,NY),INTENT(OUT) :: PSI
REAL(RPREC) :: CSI,CSI_0
INTEGER :: JX,JY
REAL(RPREC) :: a_stable,b_stable
REAL(RPREC) :: a_unstable,b_unstable,c_unstable,d_unstable,n_unstable,y_unstable,x_unstable
real(rprec) :: temp_psi_1,temp_psi_2
!BUILD VALUE OF ST_CSI

a_stable = 6.0_rprec
b_stable = 2.5_rprec

a_unstable = 0.33_rprec
b_unstable = 0.41_rprec
c_unstable = 0.33_rprec
d_unstable = 0.057_rprec
n_unstable = 0.78_rprec


!BUILD PSI_M TOT: -PSI_M(Z)+PSI_M(Z0)
DO JY=1,NY
 DO JX=1,NX

   csi = z/obk_l(jx,jy)
   csi_0 = z0(jx,jy)/obk_l(jx,jy)

   IF(OBK_L(JX,JY) .ge. 0) then
    ! locally stable
    psi(JX,JY) = -1.0_rprec * (-1.0_rprec*a_stable*dlog(csi + (1+csi**b_stable)**(1.0_rprec/b_stable))) &
                 + (-1.0_rprec*a_stable*dlog(csi_0 + (1+(csi_0)**b_stable)**(1.0_rprec/b_stable)))

   ELSE ! locally unstable
    temp_psi_1 = 0.0_rprec
    y_unstable = -1.0_rprec * csi
    temp_psi_1 = ((1.0_rprec-d_unstable)/n_unstable)*dlog((c_unstable + y_unstable**(n_unstable))/c_unstable)

    temp_psi_2 = 0.0_rprec
    y_unstable = -1.0_rprec * csi_0
    temp_psi_2 = ((1.0_rprec-d_unstable)/n_unstable)*dlog((c_unstable + y_unstable**(n_unstable))/c_unstable)

    psi(jx,jy) = -1.0_rprec*temp_psi_1 + temp_psi_2

   ENDIF

 ENDDO
ENDDO

END SUBROUTINE COMPUTE_PSI_H

SUBROUTINE COMPUTE_PHI_M(Z,Z0,PHI)

IMPLICIT NONE

REAL(RPREC),INTENT(IN) :: Z
REAL(RPREC),DIMENSION(LD,NY),INTENT(IN) :: Z0
REAL(RPREC),DIMENSION(LD,NY),INTENT(OUT) :: PHI
REAL(RPREC) :: CSI,CSI_0
INTEGER :: JX,JY
REAL(RPREC) :: a_stable,b_stable
REAL(RPREC) :: a_unstable,b_unstable,c_unstable,d_unstable,n_unstable,y_unstable,x_unstable
real(rprec) :: temp_psi_1,temp_psi_2
!BUILD VALUE OF ST_CSI

a_stable = 6.0_rprec
b_stable = 2.5_rprec

a_unstable = 0.33_rprec
b_unstable = 0.41_rprec
c_unstable = 0.33_rprec
d_unstable = 0.057_rprec
n_unstable = 0.78_rprec


!BUILD PSI_M TOT: -PSI_M(Z)+PSI_M(Z0)
DO JY=1,NY
 DO JX=1,NX

   csi = z/obk_l(jx,jy)
   csi_0 = z0(jx,jy)/obk_l(jx,jy)

   IF(OBK_L(JX,JY) .ge. 0) then
    ! locally stable
    phi(jx,jy) = 1.0_rprec + a_stable * ( (csi + (csi**(b_stable))*(1.0_rprec + (csi**(b_stable)))**(-1.0_rprec + 1.0_rprec/b_stable)) /&
                                          (csi + (1.0_rprec + csi**(b_stable))**(1.0_rprec/b_stable)) )
   ELSE ! locally unstable
    y_unstable = -1.0_rprec * csi
    if(y_unstable .le. b_unstable**(-3.0_rprec)) then
      phi(jx,jy) = (a_unstable + b_unstable*y_unstable**(4.0_rprec/3.0_rprec)) / (a_unstable + y_unstable)
    else
      phi(jx,jy) = 1.0_rprec
    endif

   ENDIF

 ENDDO
ENDDO

END SUBROUTINE COMPUTE_PHI_M

SUBROUTINE COMPUTE_PHI_H(Z,Z0,PHI)

IMPLICIT NONE

REAL(RPREC),INTENT(IN) :: Z
REAL(RPREC),DIMENSION(LD,NY),INTENT(IN) :: Z0
REAL(RPREC),DIMENSION(LD,NY),INTENT(OUT) :: PHI
REAL(RPREC) :: CSI,CSI_0
INTEGER :: JX,JY
REAL(RPREC) :: a_stable,b_stable
REAL(RPREC) :: a_unstable,b_unstable,c_unstable,d_unstable,n_unstable,y_unstable,x_unstable
!BUILD VALUE OF ST_CSI

a_stable = 6.0_rprec
b_stable = 2.5_rprec

a_unstable = 0.33_rprec
b_unstable = 0.41_rprec
c_unstable = 0.33_rprec
d_unstable = 0.057_rprec
n_unstable = 0.78_rprec


!BUILD PSI_M TOT: -PSI_M(Z)+PSI_M(Z0)
DO JY=1,NY
 DO JX=1,NX

   csi = z/obk_l(jx,jy)
   csi_0 = z0(jx,jy)/obk_l(jx,jy)

   IF(OBK_L(JX,JY) .ge. 0) then
    ! locally stable
    phi(jx,jy) = 1.0_rprec + a_stable * ( (csi + (csi**(b_stable))*(1.0_rprec + (csi**(b_stable)))**(-1.0_rprec + 1.0_rprec/b_stable)) /&
                                          (csi + (1.0_rprec + csi**(b_stable))**(1.0_rprec/b_stable)) )
   ELSE ! locally unstable
    y_unstable = -1.0_rprec * csi

    phi(jx,jy) = (c_unstable + d_unstable*(y_unstable**(n_unstable)))/(c_unstable + (y_unstable**(n_unstable)))

   ENDIF

 ENDDO
ENDDO

END SUBROUTINE COMPUTE_PHI_H

END MODULE MO_SIMILARITY
