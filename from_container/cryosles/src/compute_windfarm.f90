module compute_windfarm

use sim_param
use param

!use mpi
!use parameters_IO , only:pi,dt_dim,dx,dy,dz,num_of_rows,num_of_columns, &
!                         spacing_x,spacing_y,staggered,rprec, &
!                         ny,nz,p_on_y,p_on_z,wind_farm_debug,&
!                         Ct,Cp,local_tip_speed_ratio,wt_averaging_time,&
!                         running_period
!
!use system_variables , only: u,rhsx,rhsy,rhsz !,array_wt_points
!use decomp_2d , only: nrank
!use decomp_2d_custom , only: coord

implicit none

!include 'mpif.h'

real(rprec),allocatable,dimension(:,:) :: array_wt_points
real(rprec),dimension(:),allocatable  :: time_disk_avg_vel
real(rprec),dimension(:),allocatable :: turbine_flag
real(rprec),dimension(:,:),allocatable :: default_turbine
real(rprec) :: weight
integer :: no_of_turbines

integer :: nrank
integer :: p_on_y
integer :: p_on_z

real(rprec) :: dt_dim
real(rprec) :: dz_local
!private 
!public :: windfarm_init,windfarm_main

contains
   
subroutine windfarm_init()

   integer :: i,turb_no,loc_array_serial_no,j,k
   integer :: counter_points,proc_y,proc_z,n_rank_needed
   integer,dimension(:),allocatable :: turbine_loc_x,turbine_loc_y
   logical :: file_exists
   character(len=64) :: fname

   nrank = coord
   p_on_y = 1
   p_on_z = nproc
   dt_dim = ref_dt
   dz_local = 1000.0_rprec/128.0_rprec

   allocate(turbine_loc_x(num_of_rows))
   allocate(turbine_loc_y(num_of_columns))

   do i=1,num_of_rows
    turbine_loc_x(i) = 12 + spacing_x*i-spacing_x/2
   enddo
  
   do i=1,num_of_columns
    turbine_loc_y(i) = spacing_y*i-spacing_y/2
   enddo
 
   allocate(default_turbine(39,3))
   open(100,file='./default_turbine.txt',action='read',form='formatted')
     do i=1,size(default_turbine,1)
      read(100,*) (default_turbine(i,j),j=1,3)
     enddo
   close(100)
   
k=13
turb_no=0
loc_array_serial_no=0
   do i=1,size(turbine_loc_x)
      do j=1,size(turbine_loc_y)
         turb_no=turb_no+1
         do counter_points=1,size(default_turbine,1)
           if(staggered) then
            if(mod(i,2) .eq. 0) then
             proc_y=floor((real(turbine_loc_y(j)+spacing_y/4,rprec) + default_turbine(counter_points,1)-1)/(ny/p_on_y))
            else
             proc_y=floor((real(turbine_loc_y(j)-spacing_y/4,rprec) + default_turbine(counter_points,1)-1)/(ny/p_on_y))
            endif
           else
            proc_y=floor((real(turbine_loc_y(j),rprec) + default_turbine(counter_points,1)-1)/(ny/p_on_y))
           endif  
           proc_z=floor((real(k,rprec)+default_turbine(counter_points,2)-1)/((nz_tot-1)/p_on_z))
           n_rank_needed = proc_y*p_on_y+proc_z
            if(nrank .eq. n_rank_needed) then
              loc_array_serial_no=loc_array_serial_no+1
            endif
         enddo
      enddo
   enddo   
 
   if(loc_array_serial_no .ne. 0) then
     allocate(array_wt_points(loc_array_serial_no,9))
   else
     allocate(array_wt_points(1,9))
     array_wt_points(1,1:3)=1.0_rprec
     array_wt_points(1,4:9)=0.0_rprec
     num_of_rows=1
     num_of_columns=1
   endif
   
   turb_no=0
   loc_array_serial_no=0
   do i=1,size(turbine_loc_x)
      do j=1,size(turbine_loc_y)
         turb_no=turb_no+1
         do counter_points=1,size(default_turbine,1)
           if(staggered) then
            if(mod(i,2) .eq. 0) then
             proc_y=floor((real(turbine_loc_y(j)+spacing_y/4,rprec) + default_turbine(counter_points,1)-1)/(ny/p_on_y))
            else
             proc_y=floor((real(turbine_loc_y(j)-spacing_y/4,rprec) + default_turbine(counter_points,1)-1)/(ny/p_on_y))
            endif
           else
            proc_y=floor((real(turbine_loc_y(j),rprec) + default_turbine(counter_points,1)-1)/(ny/p_on_y))
           endif 
         proc_z=floor((real(k,rprec)+default_turbine(counter_points,2)-1)/((nz_tot-1)/p_on_z))
         n_rank_needed = proc_y*p_on_y+proc_z
         if(nrank .eq. n_rank_needed) then
           loc_array_serial_no=loc_array_serial_no+1
           array_wt_points(loc_array_serial_no,1) = turbine_loc_x(i)
           if(staggered) then
              if(mod(i,2) .eq. 0) then
           array_wt_points(loc_array_serial_no,2) = turbine_loc_y(j)+spacing_y/4+default_turbine(counter_points,1) - proc_y*(ny/p_on_y)
              else
           array_wt_points(loc_array_serial_no,2) = turbine_loc_y(j)-spacing_y/4+default_turbine(counter_points,1) - proc_y*(ny/p_on_y)
              endif
           else
            array_wt_points(loc_array_serial_no,2) = turbine_loc_y(j) + default_turbine(counter_points,1) - proc_y*(ny/p_on_y)
           endif 
           array_wt_points(loc_array_serial_no,3) = k+default_turbine(counter_points,2) - proc_z*((nz_tot-1)/p_on_z)
           array_wt_points(loc_array_serial_no,4) = default_turbine(counter_points,3)
           array_wt_points(loc_array_serial_no,5) = turb_no
           array_wt_points(loc_array_serial_no,6) = counter_points
           array_wt_points(loc_array_serial_no,7:9) = 0.0_rprec
         endif
         enddo
      enddo
   enddo   

no_of_turbines=1
   if(size(array_wt_points,1) .gt. 1) then
      do i=2,size(array_wt_points,1)
         if(array_wt_points(i,5) .ne. array_wt_points(i-1,5)) then
           no_of_turbines=no_of_turbines+1
         endif
      enddo
   endif  
   write(*,*) 'nrank,coord,no_of_turbines:',nrank,coord,no_of_turbines
  
   !if(wind_farm_debug) then
   write(fname,'("./wind_farm_setup/array_wt_points.",i0)') nrank
   open(UNIT=200,FILE=fname,FORM="FORMATTED",ACTION="WRITE",STATUS="REPLACE")
   do i=1,size(array_wt_points,1)
        write(200,*) array_wt_points(i,1:3)
   enddo
   close(200)
   !endif
   allocate(time_disk_avg_vel(num_of_rows*num_of_columns))
   time_disk_avg_vel(:)=0.0_rprec

   weight=dt_dim/(wt_averaging_time*60.0_rprec + dt_dim)
   write(*,*) 'WEIGHT:', weight,dt_dim,wt_averaging_time,dt
   
  write(fname,'("./turbine_data/restart_vel.",i0)') nrank
  INQUIRE(FILE=fname,EXIST=file_exists)
  if(file_exists) then
     open(UNIT=200,FILE=fname,FORM="FORMATTED",ACTION="READ")
     do i=1,size(time_disk_avg_vel)
     read(200,*) time_disk_avg_vel(i)
     enddo
     close(200)
  endif

  allocate(turbine_flag(num_of_rows*num_of_columns))
  turbine_flag(:)=1.0_rprec
  write(fname,'("./wind_farm_setup/turb_flag.txt")')
  inquire(file=fname,exist=file_exists)
  if(file_exists) then
     open(unit=200,file=fname,form="formatted",action="read") 
     do i=1,size(turbine_flag)
     read(200,*) turbine_flag(i)
     enddo
     close(200)
   endif

   end subroutine windfarm_init

   subroutine windfarm_main(jt,jt_total)
 
   integer,intent(in) :: jt,jt_total
   real(rprec),dimension(no_of_turbines,2) :: u_sum,power_sum,torque_sum
   integer :: i,j
   integer :: loc_vec
   character(len=64) :: fname,temp
   integer,dimension(p_on_y*p_on_z) :: incoming_data_size,recvcounts,displs
   integer :: ierror
   real(rprec),dimension(:),allocatable :: all_turbine_vel,all_turbine_numbers
   real(rprec),dimension(:),allocatable :: all_turbine_power,all_turbine_torque
   real(rprec),dimension(num_of_rows*num_of_columns) :: disk_vel_avg
   real(rprec),dimension(num_of_rows*num_of_columns) :: disk_power_avg
   real(rprec),dimension(num_of_rows*num_of_columns) :: disk_torque_avg
   real(rprec) :: angle,radial_distance
   real(rprec),dimension(num_of_rows*num_of_columns) :: angular_velocity
   integer :: reclen
   ! compute u_avg locally
   loc_vec=size(array_wt_points,1)/no_of_turbines
   do i=1,no_of_turbines
      u_sum(i,1)=0.0_rprec
      do j=1,loc_vec
         u_sum(i,1)= u_sum(i,1) + u(int(array_wt_points((i-1)*loc_vec+j,1)),& 
                                    int(array_wt_points((i-1)*loc_vec+j,2)),& 
                                    int(array_wt_points((i-1)*loc_vec+j,3)) ) 

      enddo
      u_sum(i,2)=array_wt_points((i-1)*loc_vec+1,5)
   enddo  

   !write(*,*) 'nrank,size(u_sum):',nrank,size(u_sum,1),size(u_sum,2)  
   ! transfer all the velocities - MPI_ALLGATHERV
   !write(*,*) 'nrank,u_sum:',nrank,u_sum 

   ! compute forces for each point locally
   if(wind_farm_debug) then
   write(fname,'("./u_sum.",i0)') nrank
   open(UNIT=200,FILE=fname,FORM="FORMATTED",ACTION="WRITE",POSITION="APPEND")
   do i=1,size(u_sum,1)
        write(200,*) u_sum(i,:)
   enddo
   close(200)
   endif

   do i=0,p_on_y*p_on_z-1
      displs(i+1)=i
   enddo 
   recvcounts(:)=1
   call mpi_allgatherv(size(u_sum,1),1,MPI_INTEGER,incoming_data_size, &
                       recvcounts,displs,MPI_INTEGER,comm,ierror)
   
   if(wind_farm_debug) then
      if(nrank .eq. 0) then
       write(*,*) 'jt,incoming data size',jt,incoming_data_size
      endif
   endif
 
   allocate(all_turbine_vel(sum(incoming_data_size,1)))
   allocate(all_turbine_numbers(sum(incoming_data_size,1)))
   allocate(all_turbine_power(sum(incoming_data_size,1)))
   allocate(all_turbine_torque(sum(incoming_data_size,1)))

   displs=0
 
   do i=2,p_on_y*p_on_z  
     displs(i)=displs(i-1)+incoming_data_size(i-1)
   enddo

   call mpi_allgatherv(u_sum(:,1),incoming_data_size(nrank+1),MPI_DOUBLE_PRECISION,all_turbine_vel,&
               incoming_data_size,displs,MPI_DOUBLE_PRECISION,comm,ierror)

   call mpi_allgatherv(u_sum(:,2),incoming_data_size(nrank+1),MPI_DOUBLE_PRECISION,all_turbine_numbers,&
               incoming_data_size,displs,MPI_DOUBLE_PRECISION,comm,ierror)

   if(wind_farm_debug) then
    if(nrank .eq. 0) then
    open(UNIT=200,FILE='incoming_vel_data.txt',FORM="FORMATTED",ACTION="WRITE",POSITION="APPEND")
      do i=1,sum(incoming_data_size)
        write(*,*) all_turbine_vel(i),all_turbine_numbers(i)
      enddo
    close(200)
    endif
   endif

   disk_vel_avg(:) = 0.0_rprec
   do i=1,num_of_rows*num_of_columns
      do j=1,sum(incoming_data_size)
         if (int(all_turbine_numbers(j)) .eq. i) then
               !write(*,*) 'i am here ! - computing disk_vel_avg'
             disk_vel_avg(i) = disk_vel_avg(i)+all_turbine_vel(j)
         endif 
      enddo
   enddo     
   
   disk_vel_avg = disk_vel_avg/size(default_turbine,1)
   time_disk_avg_vel = (1.0_rprec-weight)*time_disk_avg_vel + weight*(disk_vel_avg)
   angular_velocity = local_tip_speed_ratio*disk_vel_avg/(50.0_rprec/1.0_rprec)
  
   angular_velocity = max(real(1E-24),real(angular_velocity)) 

   if(mod(jt_total,p_cnt3) .eq. 0) then
      write(fname,'("./turbine_data/restart_vel.",i0)') nrank
      open(UNIT=200,FILE=fname,FORM="FORMATTED",ACTION="WRITE",STATUS="REPLACE")
      do i=1,size(time_disk_avg_vel)
         write(200,*) time_disk_avg_vel(i)
      enddo
      close(200) 
   endif

   if(wind_farm_debug) then
     if(nrank .eq. 0) then
       do i=1,size(time_disk_avg_vel,1)
         write(*,*) 'time_disk_avg_vel',jt,time_disk_avg_vel(i),disk_vel_avg(i),weight
       enddo
     endif
   endif
  !!$OMP PARALLEL DO
   do i=1,size(array_wt_points,1)

        if(array_wt_points(i,5) .ne. 0) then
      radial_distance = sqrt((default_turbine(array_wt_points(i,6),1)*dy)**2.0_rprec + &
                             (default_turbine(array_wt_points(i,6),2)*dz_local)**2.0_rprec)
      
      angle = atan2(default_turbine(array_wt_points(i,6),2)*dz_local, &
                    default_turbine(array_wt_points(i,6),1)*dy) + pi/2.0_rprec

      rhsx(array_wt_points(i,1),array_wt_points(i,2),array_wt_points(i,3)) = &
      rhsx(array_wt_points(i,1),array_wt_points(i,2),array_wt_points(i,3)) - & 
      0.5_rprec*Ct*(time_disk_avg_vel(array_wt_points(i,5))**2.0_rprec)*(array_wt_points(i,4)/dx)* &
      turbine_flag(array_wt_points(i,5)) 
 
      array_wt_points(i,7) = array_wt_points(i,7) + 0.5_rprec*Ct* &
                             (time_disk_avg_vel(array_wt_points(i,5))**2.0_rprec)*(array_wt_points(i,4)/dx)* &
      turbine_flag(array_wt_points(i,5)) 
 
 
     if ( radial_distance .ne. 0.0_rprec) then 
      rhsy(array_wt_points(i,1),array_wt_points(i,2),array_wt_points(i,3)) = &
      rhsy(array_wt_points(i,1),array_wt_points(i,2),array_wt_points(i,3)) - & 
      0.5_rprec*Cp*(time_disk_avg_vel(array_wt_points(i,5))**2.0_rprec) * &
      (u(array_wt_points(i,1),array_wt_points(i,2),array_wt_points(i,3))/ &
      (angular_velocity(array_wt_points(i,5)) * radial_distance)) &
      * cos(angle) * (array_wt_points(i,4)/dx) * &
      turbine_flag(array_wt_points(i,5)) 

      !write(*,*) angular_velocity(array_wt_points(i,5)), radial_distance, cos(angle),radial_distance

      array_wt_points(i,8) = array_wt_points(i,8) + 0.5_rprec*Cp*((time_disk_avg_vel(array_wt_points(i,5))**2.0_rprec) * &
                                                    (u(array_wt_points(i,1),array_wt_points(i,2),array_wt_points(i,3))/ &
                                                    (angular_velocity(array_wt_points(i,5)) * radial_distance))) &
                                                    * cos(angle) * (array_wt_points(i,4)/dx) * &
      turbine_flag(array_wt_points(i,5)) 



      rhsz(array_wt_points(i,1),array_wt_points(i,2),array_wt_points(i,3)) = &
      rhsz(array_wt_points(i,1),array_wt_points(i,2),array_wt_points(i,3)) - & 
      0.5_rprec*Cp*(time_disk_avg_vel(array_wt_points(i,5))**2.0_rprec) * &
      (u(array_wt_points(i,1),array_wt_points(i,2),array_wt_points(i,3))/ &
      (angular_velocity(array_wt_points(i,5)) * radial_distance)) &
      * sin(angle) * (array_wt_points(i,4)/dx)* &
      turbine_flag(array_wt_points(i,5)) 


       array_wt_points(i,9) = array_wt_points(i,9) + 0.5_rprec*Cp*((time_disk_avg_vel(array_wt_points(i,5))**2.0_rprec) * &
                                                    (u(array_wt_points(i,1),array_wt_points(i,2),array_wt_points(i,3))/ &
                                                    (angular_velocity(array_wt_points(i,5)) * radial_distance))) &
                                                    * sin(angle) * (array_wt_points(i,4)/dx) * &
      turbine_flag(array_wt_points(i,5)) 

  
     endif
       endif
   enddo 
   !!$OMP END PARALLEL DO

   loc_vec=size(array_wt_points,1)/no_of_turbines
   do i=1,no_of_turbines
      power_sum(i,1)=0.0_rprec
      torque_sum(i,1)=0.0_rprec
      do j=1,loc_vec
        if(array_wt_points((i-1)*loc_vec+j,5) .ne. 0) then
         power_sum(i,1)= power_sum(i,1) + u(int(array_wt_points((i-1)*loc_vec+j,1)),     & 
                                            int(array_wt_points((i-1)*loc_vec+j,2)),     & 
                                            int(array_wt_points((i-1)*loc_vec+j,3)) )    &
        *0.5_rprec*Ct*(time_disk_avg_vel(array_wt_points((i-1)*loc_vec+j,5))**2.0_rprec) &
        *(array_wt_points((i-1)*loc_vec+j,4)/dx) * &
      turbine_flag(array_wt_points((i-1)*loc_vec+j,5)) 


         radial_distance = sqrt((default_turbine(array_wt_points((i-1)*loc_vec+j,6),1)*dy)**2.0_rprec + &
                                (default_turbine(array_wt_points((i-1)*loc_vec+j,6),2)*dz_local)**2.0_rprec)

         if ( radial_distance .ne. 0.0_rprec) then 
         torque_sum(i,1) = torque_sum(i,1) + 0.5_rprec*Cp*(time_disk_avg_vel(array_wt_points((i-1)*loc_vec+j,5))**2.0_rprec) * &
                          (u(array_wt_points((i-1)*loc_vec+j,1),array_wt_points((i-1)*loc_vec+j,2),array_wt_points((i-1)*loc_vec+j,3))) &
                           * (array_wt_points((i-1)*loc_vec+j,4)/dx) * &
      turbine_flag(array_wt_points((i-1)*loc_vec+j,5)) 


         endif
       endif
      enddo
      power_sum(i,2)=array_wt_points((i-1)*loc_vec+1,5)
      torque_sum(i,2)=array_wt_points((i-1)*loc_vec+1,5)
   enddo  

   call mpi_allgatherv(power_sum(:,1),incoming_data_size(nrank+1),MPI_DOUBLE_PRECISION,all_turbine_power,&
               incoming_data_size,displs,MPI_DOUBLE_PRECISION,comm,ierror)

   call mpi_allgatherv(torque_sum(:,1),incoming_data_size(nrank+1),MPI_DOUBLE_PRECISION,all_turbine_torque,&
               incoming_data_size,displs,MPI_DOUBLE_PRECISION,comm,ierror)

   disk_power_avg(:)  = 0.0_rprec
   disk_torque_avg(:) = 0.0_rprec
   do i=1,num_of_rows*num_of_columns
      do j=1,sum(incoming_data_size)
         if (int(all_turbine_numbers(j)) .eq. i) then
             disk_power_avg(i)  = disk_power_avg(i)  +  all_turbine_power(j)
             disk_torque_avg(i) = disk_torque_avg(i) +  all_turbine_torque(j)
         endif 
      enddo
   enddo     


   !if(nrank .eq. 0) then
   ! open(200,file='power.txt',form="formatted",action="write",position="append")
   !   write(200,*) jt_total,disk_power_avg
   ! close(200)
   ! open(200,file='torque.txt',form="formatted",action="write",position="append")
   !   write(200,*) jt_total,disk_torque_avg
   ! close(200)
   !endif

   if(mod(jt_total,p_cnt3) .eq. 0) then

   write (fname, '(a,i7.7,a)') './wind_farm_forces/wt_x',jt_total,'.out'
   write (temp, '(".c",i0)') nrank
   fname = trim (fname) // temp
   open(unit=200,file=fname,form="formatted",action='write')
      do i=1,size(array_wt_points,1)
        write(200,*) array_wt_points(i,7)
      enddo
   close(200)
   array_wt_points(:,7) = 0.0_rprec

   write (fname, '(a,i7.7,a)') './wind_farm_forces/wt_y',jt_total,'.out'
   write (temp, '(".c",i0)') nrank
   fname = trim (fname) // temp
   open(unit=200,file=fname,form="formatted",action='write')
       do i=1,size(array_wt_points,1)
         write(200,*) array_wt_points(i,8)
       enddo
   close(200)
   array_wt_points(:,8) = 0.0_rprec

   write (fname, '(a,i7.7,a)') './wind_farm_forces/wt_z',jt_total,'.out'
   write (temp, '(".c",i0)') nrank
   fname = trim (fname) // temp
   open(unit=200,file=fname,form="formatted",action='write')
       do i=1,size(array_wt_points,1) 
        write(200,*) array_wt_points(i,9)
       enddo
   close(200)
   array_wt_points(:,9) = 0.0_rprec

   endif

   deallocate(all_turbine_vel)
   deallocate(all_turbine_numbers)
   deallocate(all_turbine_power)
   deallocate(all_turbine_torque)

   end subroutine windfarm_main

end module
