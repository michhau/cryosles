module mmatrix

use param
implicit none


contains

subroutine init_mmatrix
implicit none

integer :: i,jz,jy,jx,counter
real(rprec) :: Lf,Ls,Lp,l

Lf = 0.1_rprec*L_x
Ls = L_x-Lf
Lp = L_x-0.25_rprec*Lf

weight_mmatrix=0.0_rprec

write(*,*) '********************'
write(*,*) Lf,Ls,Lp
write(*,*) '********************'

do i=1,nx
   l = real(i-1,rprec)*dx
   if(l .le. Ls) then
    weight_mmatrix(i) = 0.0_rprec
   elseif((l .gt. Ls) .and. (l .le. Lp)) then
    weight_mmatrix(i) = 0.5_rprec*(1.0_rprec - dcos(pi*(l-Ls)/(Lp-Ls)))
   else
    weight_mmatrix(i) = 1.0_rprec
   endif
   !write(*,*) weight_mmatrix(i),i
enddo

!counter=0
!do jz=0,nz
!  do jy=1,ny
!   do jx=1,ld
!   counter=counter+1
!   flat_weight(counter) = weight_mmatrix(jx,jy,jz)
!   enddo
!enddo
!enddo

!flat_weight = pack(weight_mmatrix,.true.)

!write(*,*) flat_weight

end subroutine init_mmatrix

!subroutine compute_mmatrix(vec_in,vec_inout,length,datatype)
!
!implicit none
!
!integer,intent(in) :: length,datatype
!real(rprec) :: vec_in(length)
!real(rprec) :: vec_inout(length)
!integer :: i
!
!do i=1,length
!  vec_inout(i) = (1.0_rprec-flat_weight(i))*vec_inout(i) + flat_weight(i)*vec_in(i)
!enddo
!
!end subroutine compute_mmatrix


end module mmatrix

