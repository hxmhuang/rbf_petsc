#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>
#include "mat_math_type.h"

module	rbf 
    use matrix	
	type MyStruct
	sequence
	PetscScalar :: a,b,c
	end type MyStruct

contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Create the points in rbf application.
! A: Matrix with m*n rows and 2 columns
! m: there are m points in x direction 
! n: there are n points in y direction 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine rbf_createpoints(A,m,n,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	PetscInt,		intent(in)		::	m,n	
	Mat,			intent(out)		::	A 
	PetscErrorCode,	intent(out)		::	ierr
		
	PetscReal	                    ::  xmin,xmax,ymin,ymax,dx,dy,xcord,ycord
	PetscInt	                    ::  ista,iend
	PetscInt	                    ::  nrow,ncol,nlocal
    PetscInt,allocatable            ::  idxm(:),idxn(:) 
    PetscScalar,allocatable         ::  rows(:)
	integer		                    ::  i,j
	
    xmin=0.0
	xmax=1.0
	ymin=0.0
	ymax=1.0
	dx= (xmax-xmin)/(m-1)
	dy= (ymax-ymin)/(n-1)
	
	call MatGetSize(A,nrow,ncol,ierr)
    if(ncol/=2) then
        print *, "Error in rbf_createpoints: The column size of matrix A should be 2."
        stop
    endif
    
    call MatGetOwnershipRange(A,ista,iend,ierr)
    nlocal=iend-ista
    allocate(idxm(nlocal),idxn(ncol),rows(2*nlocal))
	
    do i=ista,iend-1
		xcord = xmin+(i/n)*dx
		ycord = ymin+mod(i,n)*dy
        rows((i-ista)*ncol+1)=xcord
        rows((i-ista)*ncol+2)=ycord
        idxm(i-ista+1)=i
	enddo
    do j=1,ncol
        idxn(j)=j-1
    enddo
    !print *, "nlocal=",nlocal,"idxm=",idxm,"idxn=",idxn,"rows=",rows
    call MatSetValues(A,nlocal,idxm,ncol,idxn,rows,INSERT_VALUES,ierr) 
    
	call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    deallocate(idxm,idxn,rows)
end subroutine


subroutine rbf_testfunctionD(A,v,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)		::	A 
	Vec,			intent(out)	    ::	v 
	PetscErrorCode,	intent(out)		::	ierr
	
	PetscInt						::	nrow,ncol
	PetscInt						::  ista,iend
	PetscScalar,allocatable			::	row(:)
	PetscReal						::	xcord,ycord,res
	integer							:: 	i

	call MatGetSize(A,nrow,ncol,ierr)
	if(ncol/=2)then
		print *, "Error: the column size of Matrix A in testfunctionD must equal to 2"
		stop	
	endif	
	
	call MatGetOwnershipRange(A,ista,iend,ierr)
	!print *,">istat=",ista,"iend=",iend,">ncol=",ncol
	allocate(row(ncol))
	
	do i=ista,iend-1
		call MatGetRow(A,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row,ierr)
		xcord=row(1)
		ycord=row(2)
		call testfunction(xcord,ycord,res)
		call VecSetValue(v,i,res,INSERT_VALUES,ierr)
		call MatRestoreRow(A,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row,ierr)
	enddo
    
    call VecAssemblyBegin(v,ierr)
    call VecAssemblyEnd(v,ierr)
	deallocate(row)
end subroutine


subroutine testfunction(xcord,ycord,res)
	implicit none
#include <petsc/finclude/petscsys.h>
	PetscReal,intent(in):: xcord,ycord
	PetscReal,intent(out):: res 
	res= 0.75*exp(-((9*xcord-2)*(9*xcord-2)+(9*ycord-2)*(9*ycord-2))/4) &
		+0.75*exp(-((9*xcord+1)*(9*xcord+1))/49-((9*ycord+1)*(9*ycord+1))/10) &
		+0.5*exp(-((9*xcord-7)*(9*xcord-7)+(9*ycord-3)*(9*ycord-3))/4) &
		-0.2*exp(-(9*xcord-4)*(9*xcord-4)-(9*ycord-7)*(9*ycord-7))
end subroutine


subroutine rbf_distancematrix(A1,A2,B,ierr)
	use matrix 
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	   	:: A1,A2
	Mat,			intent(out)     :: B 
	PetscErrorCode,	intent(out)	    :: ierr
	Mat                     		:: W1,W2,W3 
	Mat                     		:: P1,P2,P3 
	PetscInt	            		:: nrow1,ncol1,nrow2,ncol2
	PetscScalar						:: alpha
	
    PetscMPIInt		myrank,mysize
	PetscLogEvent	ievent(13)
	call PetscLogEventRegister("1_dm_eprod",0, ievent(1), ierr)
	call PetscLogEventRegister("2_dm_sum",0, ievent(2), ierr)
	call PetscLogEventRegister("3_dm_rep",0, ievent(3), ierr)
	call PetscLogEventRegister("4_dm_xyt",0, ievent(4), ierr)
	call PetscLogEventRegister("5_dm_axpy",0, ievent(5), ierr)
	call PetscLogEventRegister("6_dm_destroy",0, ievent(6), ierr)
	call PetscLogEventRegister("7_dm_eprod",0, ievent(7), ierr)
	call PetscLogEventRegister("8_dm_sum",0, ievent(8), ierr)
	call PetscLogEventRegister("9_dm_trans",0, ievent(9), ierr)
	call PetscLogEventRegister("10_dm_rep",0, ievent(10), ierr)
	call PetscLogEventRegister("11_dm_axpy",0, ievent(11), ierr)
	call PetscLogEventRegister("12_dm_math_sqrt",0, ievent(12), ierr)
	call PetscLogEventRegister("13_dm_destroy",0, ievent(13), ierr)

    call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr)
	call MPI_Comm_rank(PETSC_COMM_WORLD,mysize,ierr)
 
	call MatGetSize(A1,nrow1,ncol1,ierr)
	call MatGetSize(A2,nrow2,ncol2,ierr)
    !print *, "nrow1=",nrow1,"nrow2=",nrow2
	if(ncol1/=ncol2)then
		print *, "Error in rbf_distancematrix: the matrix A1 and A2 should have the same column size."
		stop	
	endif
    
	call PetscLogEventBegin(ievent(1),ierr)
    call mat_eprod(A1,A1,W1,ierr)
	call PetscLogEventEnd(ievent(1),ierr)
    
	call PetscLogEventBegin(ievent(2),ierr)
    call mat_sum(W1,2,W2,ierr)
	call PetscLogEventEnd(ievent(2),ierr)
    
	call PetscLogEventBegin(ievent(3),ierr)
    call mat_rep(W2,1,nrow2,P1,ierr)
	call PetscLogEventEnd(ievent(3),ierr)
	!if(myrank==0) print *, ">P1="
	!call mat_view(P1,ierr)
	
	call PetscLogEventBegin(ievent(4),ierr)
	call mat_xyt(A1,A2,P2,ierr)
	call PetscLogEventEnd(ievent(4),ierr)
	!if(myrank==0) print *, ">P2="
	!call mat_view(P2,ierr)
	call PetscLogEventBegin(ievent(5),ierr)
	alpha=-2.0
	call mat_axpy(P1,alpha,P2,ierr)
	call PetscLogEventEnd(ievent(5),ierr)
	
	call PetscLogEventBegin(ievent(6),ierr)
    call mat_destroy(W1,ierr)
    call mat_destroy(W2,ierr)
	call PetscLogEventEnd(ievent(6),ierr)
	
	call PetscLogEventBegin(ievent(7),ierr)
	call mat_eprod(A2,A2,W1,ierr)
	call PetscLogEventEnd(ievent(7),ierr)
	call PetscLogEventBegin(ievent(8),ierr)
	call mat_sum(W1,2,W2,ierr)
	call PetscLogEventEnd(ievent(8),ierr)
	call PetscLogEventBegin(ievent(9),ierr)
	call mat_trans(W2,W3,ierr)
	call PetscLogEventEnd(ievent(9),ierr)
	call PetscLogEventBegin(ievent(10),ierr)
	call mat_rep(W3,nrow1,1,P3,ierr)
	call PetscLogEventEnd(ievent(10),ierr)
	!if(myrank==0) print *, ">P3="
	!call mat_view(P3,ierr)
	
	call PetscLogEventBegin(ievent(11),ierr)
	alpha=1.0
	call mat_axpy(P1,alpha,P3,ierr)	
	call PetscLogEventEnd(ievent(11),ierr)
    
    
	call PetscLogEventBegin(ievent(12),ierr)
    call mat_math(P1,MAT_MATH_SQRT,B,ierr)
	call PetscLogEventEnd(ievent(12),ierr)
	!if(myrank==0) print *, ">B="
	!call mat_view(B,ierr)

	call PetscLogEventBegin(ievent(13),ierr)
	call mat_destroy(W1,ierr)
    call mat_destroy(W2,ierr)
    call mat_destroy(W3,ierr)
	call mat_destroy(P1,ierr)
    call mat_destroy(P2,ierr)
    call mat_destroy(P3,ierr)
	call PetscLogEventEnd(ievent(13),ierr)
end subroutine

subroutine rbf_guassian(ep,A,B,ierr)
!    use matrix
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include "mat_math_type.h"
	PetscReal,		intent(in)		::  ep		
	Mat,			intent(in)		::	A 
	Mat,			intent(out)		::	B 
	PetscErrorCode,	intent(out)		::	ierr
    PetscScalar                     ::  alpha	
	Mat                 			::	W
    !rbf = @(e,r) exp(-(e*r).^2)
    call mat_eprod(A,A,W,ierr)
    alpha=-ep*ep
    call mat_scale(W,alpha,ierr)
    call mat_math(W,MAT_MATH_EXP,B,ierr)
	call mat_destroy(W,ierr)
end subroutine

end module
