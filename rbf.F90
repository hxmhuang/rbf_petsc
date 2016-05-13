#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>

module	rbf 
	
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


subroutine rbf_distancematrix(dsites,ctrs,dm,ierr)
	use matrix 
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	   	:: dsites,ctrs
	Mat,			intent(out)     :: dm
	PetscErrorCode,	intent(out)	    :: ierr
	Mat                     		:: W1,W2,W3,W4 
	Mat                     		:: P1,P2,P3 
	PetscInt	            		:: nrow1,ncol1,nrow2,ncol2,ista,iend	
	PetscScalar						:: alpha

	PetscMPIInt		myrank,mysize
	call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr)
	call MPI_Comm_rank(PETSC_COMM_WORLD,mysize,ierr)
 
	call MatGetSize(dsites,nrow1,ncol1,ierr)
	call MatGetSize(dsites,nrow2,ncol2,ierr)
	if(ncol1/=ncol2)then
		print *, "Error in rbf_distancematrix: the matrix A1 and A2 should have the same column size."
		stop	
	endif
    
    call mat_eprod(dsites,dsites,W1,ierr)
    call mat_sum(W1,2,W2,ierr)
    call mat_rep(W2,1,nrow2,P1,ierr)
	if(myrank==0) print *, ">P1="
	call mat_view(P1,ierr)
	
	call mat_xyt(dsites,ctrs,P2,ierr)
	if(myrank==0) print *, ">P2="
	call mat_view(P2,ierr)
	alpha=-2.0
	call mat_axpy(P1,alpha,P2,ierr)
	
    call mat_destroy(W1,ierr)
    call mat_destroy(W2,ierr)
	
	call mat_eprod(ctrs,ctrs,W1,ierr)
	call mat_sum(W1,2,W2,ierr)
	call mat_trans(W2,W3,ierr)
	call mat_rep(W3,nrow1,1,P3,ierr)
	if(myrank==0) print *, ">P3="
	call mat_view(P3,ierr)
	
	alpha=1.0
	call mat_axpy(P1,alpha,P3,ierr)	
   
	call mat_copy(P1,dm,ierr) 
	if(myrank==0) print *, ">dm="
	call mat_view(dm,ierr)

	call mat_destroy(W1,ierr)
    call mat_destroy(W2,ierr)
    call mat_destroy(W3,ierr)
    
	call mat_destroy(P1,ierr)
    call mat_destroy(P2,ierr)
    call mat_destroy(P3,ierr)

end subroutine


end module
