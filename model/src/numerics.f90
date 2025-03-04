!!----------------------------------------------------------------------------------
!!  OBJECT: num_lib
!!----------------------------------------------------------------------------------
!!
!!  LAST MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 30/12/2011
!!
!!  DESCRIPTION:
!!
!!  This module provides a set of computational subroutines
!!
!!----------------------------------------------------------------------------------

module numerics

use types,only:rprec

implicit none


!-----------------------------------------------------------------------------------
CONTAINS   !------------------------------------------------------------------------
!-----------------------------------------------------------------------------------



!!----------------------------------------------------------------------------------
!!  OBJECT: subroutine LU_solver
!!----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 11/12/2011
!!
!!  DESCRIPTION:
!!
!!  This subroutine solves the linear system Ax = b composed of a general matrix 
!!  through an LU decomposition + backsubstitution
!!
!!----------------------------------------------------------------------------------
subroutine LU_solver(A,B)

implicit none

external DGESV

real(rprec), dimension(:,:), intent(inout) :: A
real(rprec), dimension(:,:), intent(inout) :: B
integer, dimension(:), allocatable :: IPIV
integer :: N, NRHS, LDA, LDB, INFO


!allocate memory
allocate(IPIV(size(B(:,1))))

!initialize variables
N = size(B(:,1))
NRHS = size(B(1,:))
LDA = N
LDB = N

!solves the system through lapack subroutine
call DGESV( N,NRHS,A,LDA,IPIV,B,LDB,INFO )

!raise warning if info different than 0
if (INFO .ne. 0) then
  print*, "<<< attention INFO from subroutine LU_solver is .ne. 0 >>>"
  print*,"info = ", info
  stop
endif

!deallocate memory
deallocate(IPIV)


end subroutine LU_solver



!!----------------------------------------------------------------------------------
!!  OBJECT: subroutine LU_solver
!!----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 11/12/2011
!!
!!  DESCRIPTION:
!!
!!  This subroutine solves the linear system Ax = b composed of a general matrix 
!!  through an LU decomposition + backsubstitution
!!
!!----------------------------------------------------------------------------------
subroutine LU_sym_solver(A,B)

implicit none

real(rprec), dimension(:,:), intent(in) :: A
real(rprec), dimension(:,:), intent(inout) :: B
integer, dimension(:), allocatable :: IPIV
real(rprec), dimension(:), allocatable :: WORK
integer :: N,NRHS,LDA,LDB,INFO,LWORK
character(1) :: UPLO


!initialize variables
UPLO = "U"
N=size(A(:,1))
NRHS=size(B(1,:))
LDA=N
LDB=N
LWORK=N

!allocate memory
allocate(IPIV(N))
allocate(WORK(LWORK))

!solves the system through lapack subroutine
call DSYSV( UPLO,N,NRHS,A,LDA,IPIV,B,LDB,WORK,LWORK,INFO )

!raise warning if info different than 0
if (INFO .ne. 0) then
  print*, "<<< attention INFO from subroutine LU_solver is .ne. 0 >>>"
  stop
endif

! print*, "Actual LWORK SIZE = ", LWORK
! print*, "Optimal LWORK SIZE = ", WORK(1)

!deallocate memory
deallocate(IPIV)
deallocate(WORK)

end subroutine LU_sym_solver



!!----------------------------------------------------------------------------------
!!  OBJECT: subroutine CHOL_solver
!!----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 11/12/2011
!!
!!  DESCRIPTION:
!!
!!  This subroutine solves the linear system Ax = b composed of a general matrix 
!!  through an LU decomposition + backsubstitution
!!
!!----------------------------------------------------------------------------------
subroutine CHOL_solver(A,B)

implicit none

external DPOSV

real(rprec), dimension(:,:), intent(in) :: A
real(rprec), dimension(:,:), intent(inout) :: B
integer :: N, NRHS, LDA, LDB, INFO
character(1) :: UPLO


!initialize variables
UPLO = "U"
N=size(A(:,1))
NRHS=size(B(1,:))
LDA=size(A(:,1))
LDB=size(A(:,1))

!solves the system through lapack subroutine
call DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )

!raise warning if info different than 0
if (INFO .ne. 0) then
  print*, "<<< Attention INFO from subroutine CHOL_solver is .ne. 0:", INFO, " >>>"
  stop
endif


end subroutine CHOL_solver



!!----------------------------------------------------------------------------------
!!  OBJECT: subroutine eigv
!!----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 29/12/2011
!!
!!  DESCRIPTION:
!!
!!  Computes all eigenvalues, and optionally, eigenvectors of a real
!!  symmetric matrix.
!!
!!----------------------------------------------------------------------------------
subroutine eigv(A)

implicit none

external DSYEV

real(rprec), dimension(:), allocatable :: W
real(rprec), dimension(:), allocatable :: WORK
real(rprec), dimension(:,:), intent(in) :: A
real(rprec), dimension(:,:), allocatable :: A_TEMP
integer :: N, LDA, LWORK, INFO
integer :: i,j
character :: JOBZ, UPLO


!allocate memory
allocate(WORK(3*size(A(1,:)-1)))
allocate(W(size(A(1,:))))
allocate(A_TEMP(size(A(1,:)),size(A(1,:))))

!initialize A_TEMP
do i=1,size(A(:,1))
do j=1,size(A(1,:))
  A_TEMP(i,j) = A(i,j)
enddo
enddo

!initialize other variables
JOBZ = "N"
UPLO = "U"
N = size(A(1,:))
LDA = size(A(1,:))
LWORK = 3*size(A(1,:))-1

!lapack subroutine
call DSYEV( JOBZ, UPLO, N, A_TEMP, LDA, W, WORK, LWORK, INFO )

!check for successful exit of DSYEV subroutine
if (INFO .ne. 0) then
  print*, "<<< Attention INFO from subroutine eigv is .ne. 0 >>>"
  stop
endif

!raise warning if a 0.0 is found
if (any(W .eq. 0.0)) then
  print*, "<<< Attention the matrix Amat is not invertible >>>"
  stop
endif

!raise warning if a negative is found
if (minval(W) .lt. 0.0) then
  print*, "<<< Attention the matrix Amat is not positive-definite >>>"
  print*, "min_lambda(matrix) =", minval(W)
  stop
endif

!deallocate memory
deallocate(WORK)
deallocate(W)
deallocate(A_TEMP)


end subroutine eigv



!!----------------------------------------------------------------------------------
!!  OBJECT: subroutine checksym
!!----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 29/12/2011
!!
!!  DESCRIPTION:
!!
!!  This subroutine checks if the given matrix is symmetric
!!
!!----------------------------------------------------------------------------------
subroutine checksym(A)

implicit none


real(rprec), dimension(:,:), intent(in) :: A
integer :: i,j

do i=1,size(A(:,1))
do j=i,size(A(1,:))

  if(A(i,j) .ne. A(j,i)) then
    print*, "<<< Attention the given matrix is not symmetric >>>"
    stop
  endif

enddo
enddo


end subroutine checksym



!!----------------------------------------------------------------------------------
!!  OBJECT: subroutine cond_num
!!----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 03/12/2011
!!
!!  DESCRIPTION:
!!
!!  DPOCON estimates the the condition number of a real symmetric positive definite !!  matrix using the Cholesky factorization A = U**T*U or A = L*L**T computed by
!!  DPOTRF.
!!
!!----------------------------------------------------------------------------------
subroutine cond_num(A)

implicit none

external DGECON, DPOTRF

character :: UPLO
real(rprec), dimension(:,:), intent(in):: A
real(rprec), dimension(:,:), allocatable :: A_TEMP
integer :: INFO, N, LDA
integer, dimension(:), allocatable :: IWORK
real(rprec), dimension(:), allocatable :: WORK
real(rprec) :: ANORM
real(rprec) :: RCOND

!parameters for lapack subroutine
UPLO = "U"
ANORM = 1.0
N = size(A(1,:))
LDA = size(A(1,:))

!allocate memory
allocate(WORK(size(A(1,:))*3))
allocate(IWORK(size(A(1,:))))
allocate(A_TEMP(size(A(:,1)),size(A(1,:))))

A_TEMP = A

!cholesky factorization A = U**T*U
call DPOTRF( UPLO, N, A_TEMP, LDA, INFO )

!raise warning if info different than 0
if (INFO .ne. 0) then
  print*, "<<< attention INFO from DPOTRF (CHOL dec) is .ne. 0:", INFO, ">>> "
  stop
endif

!calculates reciprocal of condition number
CALL DPOCON( UPLO, N, A_TEMP, LDA, ANORM, RCOND, WORK, IWORK, INFO )

!raise warning if info different than 0
if (INFO .ne. 0) then
  print*, "<<< attention INFO from DPOCON is .ne. 0:", INFO, ">>> "
  stop
endif

!print condition number as return value
print*, "CONDITION NUMBER OF MATRIX A = ", int(1.0/RCOND)

!deallocate memory
deallocate(WORK)
deallocate(IWORK)
deallocate(A_TEMP)

end subroutine cond_num



!!----------------------------------------------------------------------------------
!!  OBJECT: subroutine invert_matrix_LU
!!----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 03/12/2011
!!
!!  DESCRIPTION:
!!
!!  This subroutine inverts the given matrix through LU decomposition and BackSub
!!  using the lapack library (when compiling add -llapack -lblas)
!!
!!----------------------------------------------------------------------------------
subroutine invert_matrix_LU(A,A_INV)

implicit none

external DGETRF
external DGETRI
external DGECON

real(rprec), dimension(:,:), intent(in):: A
real(rprec), dimension(:,:), intent(out):: A_INV
integer :: INFO, LWORK, M, N, LDA
integer, dimension(:), allocatable :: IPIV, IWORK
real(rprec), dimension(:), allocatable :: WORK
real(rprec), dimension(:), allocatable :: WORK1
character :: NORM


!allocate memory
allocate(IPIV(size(A(1,:))))
allocate(IWORK(size(A(1,:))))
allocate(WORK(size(A(1,:))))
allocate(WORK1(size(A(1,:))*4))

!parameters for lapack subroutine
NORM = "1"
N = size(A(1,:))
M = size(A(1,:))
LDA = size(A(1,:))
LWORK = size(A(1,:))

!calculates the inverse of Amat (posed to 0 if A is all 0s)
call DGETRF(M, N, A, LDA, IPIV, INFO)

!raise warning if info different than 0
if (INFO .ne. 0) then
  print*, "<<< attention INFO from DGETRF (LU dec) is .ne. 0:", INFO, ">>> "
  stop
endif

! call DGECON( NORM, N, A, LDA, ANORM, RCOND, WORK1, IWORK,INFO )
call DGETRI(N, A, LDA, IPIV, WORK, LWORK, INFO)

!raise warning if info different than 0
if (INFO .ne. 0) then
  print*, "<<< attention INFO from DGETRI (sys solve) is .ne. 0:", INFO, ">>> "
  stop
endif

!assign output value
A_INV = A

deallocate(IPIV)
deallocate(WORK)
deallocate(IWORK)
deallocate(WORK1)

end subroutine invert_matrix_LU



!!----------------------------------------------------------------------------------
!!  OBJECT: subroutine invert_matrix_CHOL
!!----------------------------------------------------------------------------------
!!
!!  LAST CHECKED/MODIFIED: Marco Giometto ( mgiometto@gmail.com ), day 03/12/2011
!!
!!  DESCRIPTION:
!!
!!  This subroutine inverts the given matrix through Cholevsky decomposition 
!!  and BackSub using the lapack library (when compiling add -llapack -lblas)
!!
!!----------------------------------------------------------------------------------
subroutine invert_matrix_CHOL(A,A_INV)

implicit none

external DPOTRF
external DPOTRI
external DPOCON

real(rprec), dimension(:,:), intent(in):: A
real(rprec), dimension(:,:), intent(out):: A_INV
integer :: INFO, N, LDA
! real(rprec) :: ANORM, RCOND
character(1) :: UPLO

!parameters for lapack subroutine
UPLO = "U"

N = size(A(1,:))
LDA = size(A(1,:))

!cholesky factorization A = U**T*U
call DPOTRF( UPLO, N, A, LDA, INFO )

!raise warning if info different than 0
if (INFO .ne. 0) then
  print*, "<<< attention INFO from DPOTRF (CHOL dec) is .ne. 0:", INFO, ">>> "
  stop
endif

!call DPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, IWORK, INFO )
call DPOTRI( UPLO, N, A, LDA, INFO )

!raise warning if info different than 0
if (INFO .ne. 0) then
  print*, "<<< attention INFO from DPOTRI (sys solve) is .ne. 0:", INFO, ">>> "
  stop
endif

!assign output value
A_INV = A


end subroutine invert_matrix_CHOL





end module numerics
