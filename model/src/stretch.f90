module stretch

use param
implicit none


contains

subroutine init_stretch
implicit none


a_grid = 0.1_rprec
b_grid = -1.0_rprec * (1.0_rprec / 2.0_rprec) * (Nz_tot-1)

c_grid = c_factor*dx
f_grid = f_factor*dx

! get c_grid and f_grid from 
allocate(get_z_local(-2:2*Nz_tot))
allocate(get_dz_local(-2:2*Nz_tot))

call compute_z_local
call compute_dz_local
end subroutine init_stretch

!function get_z_local(i)
!real(rprec),intent(in) :: i
!real(rprec) :: get_z_local
!
!get_z_local = (A_grid/alpha_grid)*dlog( dcosh(alpha_grid*i-gamma_grid)/dcosh(alpha_grid-gamma_grid)) + (A_grid+B_grid)*i - A_grid
!
!end function get_z_local

!function get_dz_local(i)
!real(rprec), intent(in) :: i
!real(rprec) :: get_dz_local
!
!get_dz_local = A_grid*(dtanh(alpha_grid*i-gamma_grid)+1.0_rprec) + B_grid
!
!end function get_dz_local

subroutine compute_z_local
implicit none

integer :: i
character(len=256) :: filename1, filename2

do i=-2,2*(nz_tot)
!get_z_local(i) = (A_grid/alpha_grid)*dlog( dcosh(alpha_grid*(real(i,rprec)*0.5_rprec)-gamma_grid)/dcosh(alpha_grid-gamma_grid)) + (A_grid+B_grid)*(real(i,rprec)*0.5_rprec) - A_grid
get_z_local(i) = ((f_grid-c_grid)/(2.0_rprec*a_grid)) * &
                 (sqrt(1.0_rprec+(a_grid*(real(i,rprec)*0.5_rprec+b_grid))**2.0_rprec) &
                 - sqrt(1.0_rprec+(a_grid*b_grid)**2.0_rprec)) + &
                 ((f_grid+c_grid)/2.0_rprec) * real(i,rprec)*0.5_rprec  

enddo

if(coord .eq. 0) then
  write(filename1,'(a)') './wheight.log'
  open(111,file=trim(filename1),status='replace')
  do i=0,2*(nz_tot)-2,2
    write(111,*) real(i,rprec)/2.0_rprec+1.0_rprec,get_z_local(i)
  enddo

  write(filename2,'(a)') './uvpheight.log'
  open(112,file=trim(filename2),status='replace')
  do i=1,2*(nz_tot)-1,2
    write(112,*) (real(i,rprec)+1.0_rprec)/2.0_rprec,get_z_local(i)
  enddo
endif

end subroutine compute_z_local

subroutine compute_dz_local
implicit none

integer :: i

do i=-2,2*(nz_tot)
  !get_dz_local(i) = A_grid*(dtanh(alpha_grid*(real(i,rprec)*0.5_rprec)-gamma_grid)+1.0_rprec) + B_grid
  get_dz_local(i) = (f_grid-c_grid)*0.5_rprec*((a_grid*(real(i,rprec)*0.5_rprec+b_grid))/ &
                    (sqrt(1.0_rprec+(a_grid*(real(i,rprec)*0.5_rprec+b_grid))**2.0_rprec)) + 1.0_rprec) + c_grid
  !write(*,*) 'i,z',i,get_dz_local(i)
enddo

end subroutine compute_dz_local

function get_eta_local(i)
real(rprec), intent(in) :: i
real(rprec) :: get_eta_local

!get_eta_local = 5.0102375513521844_rprec*(1.5716860586567082_rprec + 2.217681496581414_rprec*(i) - sqrt(2.470197066975863_rprec - 4.603623765100343_rprec*(i) +3.1475911809789314_rprec*(i)*(i)));
!get_eta_local = 1.1562086656966579_rprec*(9.636779865654959_rprec + 88.70725986325656_rprec*(i) - sqrt(92.86752617909338_rprec - 672.7383457534048_rprec*(i) + 7101.752602083714_rprec*(i)*(i)))
!get_eta_local = 5.637555837147649_rprec*(3.367831612972066_rprec + 17.977888694869334_rprec*(i) - sqrt(11.342289773333995_rprec - 78.63221637707105_rprec*(i) + 306.19674385402926_rprec*(i)*(i)))
! for domain height of 18.317655 m and NZ_tot = 193:
!get_eta_local = 162823600_rprec/4379319_rprec*(i) + 45.5963426738067_rprec - sqrt(992.650706509249_rprec*(i)*(i) - 2691.39679378331_rprec*(i) + 2079.02646522722_rprec)
! for domain height of 6.4 m and NZ_tot = 129:
!get_eta_local = 4.174150902928559_rprec*(7.22783801679223_rprec + 13.30942671875569_rprec*(i) - sqrt(52.241642396986944_rprec - 138.07086947692574_rprec*(i) + 113.37013733243383_rprec*(i)*(i)))
!get_eta_local = 1.9649756542092672_rprec*(16.29776532602988_rprec + 5.219611873785569_rprec*(i) - 1.000000_rprec*sqrt(265.61715462234242_rprec - 153.31078467741929_rprec*(i) + 24.588024171949374_rprec*(i*i)))
!get_eta_local = 20.371837_rprec*(i)
get_eta_local = 1018589_rprec*(i)
end function get_eta_local

end module stretch

