!!------------------------------------------------------------------------------
!!  DESCRIPTION
!!
!!  This module calls the specific output subroutines and contains the output
!!  parameters (e.g. p_cnt_xx,c_cnt_xx)
!!
!!  NOTE: make use to USE just the specific output and not the full modules
!!        (e.g. IBM_OUT and not IBM or IBM_DFA_INIT).
!!
!!------------------------------------------------------------------------------
module outloop_special

use param
use avgField
use tkeBudget1D
use tkeBudget1DW
implicit none

logical,parameter :: surf_data=.true.
integer,parameter :: p_cnts=400000,c_cnts=40

contains

!!------------------------------------------------------------------------------
!!  OBJECT: subroutine OUTPUT_LOOP_SPECIFIC ()
!!------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 26/10/2011
!!
!!  DESCRIPTION:
!!
!!  This sub calls all the others SPECIFIC subroutines needed for plotting 
!!  at specified intervals.
!!
!!
!!------------------------------------------------------------------------------
subroutine outloop_spc


implicit none

integer :: i,j,k
real(rprec),parameter :: eps=1000._rprec*epsilon(0._rprec)


if(mod(jt,c_cnt3)==0) call saveAvgField()

!if(jt==0)      call_alloc_tkeBudget1D()
!if(mod(jt,c_cnt3)==0) call computeTkeBudget()
!if(jt==nsteps) call_dealloc_tkeBudget1D()

!if(jt==0)      call alloc_tkeBudget1Dw()
!if(mod(jt,c_cnt3)==0) call computeTkeBudgetWnodes()
!if(jt==nsteps) call dealloc_tkeBudget1Dw()

end subroutine outloop_spc


end module outloop_special
