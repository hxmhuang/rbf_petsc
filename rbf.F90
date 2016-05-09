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
		
	PetscReal	xmin,xmax,ymin,ymax,dx,dy,xcord,ycord
	PetscInt	ista,iend
	integer		i

	xmin=0.0
	xmax=1.0
	ymin=0.0
	ymax=1.0
	dx= (xmax-xmin)/(m-1)
	dy= (ymax-ymin)/(n-1)
	
	! generate matrix dsites and ctrs with size M*2 
	call MatGetOwnershipRange(A,ista,iend,ierr)
	
	do i=ista,iend-1
		xcord = xmin+(i/n)*dx
		ycord = ymin+mod(i,n)*dy
		call MatSetValue(A,i,0,xcord,INSERT_VALUES,ierr) 
		call MatSetValue(A,i,1,ycord,INSERT_VALUES,ierr) 
	enddo
	
	call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
end subroutine


subroutine rbf_testfunctionD(A,v,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)		::	A 
	Vec,			intent(inout)	::	v 
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


subroutine rbf_distancematrix(dsites,ctrs,dm)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,intent(in)	:: dsites,ctrs
	Mat,intent(out):: dm
	PetscInt	nrow,ncol,rsta,rend	
	PetscErrorCode	ierr
	integer		i
	call MatGetSize(dsites,nrow,ncol,ierr)

	call MatGetOwnershipRange(dsites,rsta,rend,ierr)
	print *,rsta,rend
	do i=rsta,rend-1,1
	
	enddo

end subroutine


end module
