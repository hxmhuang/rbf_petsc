#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>

module particle 
	type MyStruct
	sequence
	PetscScalar :: a,b,c
	end type MyStruct

contains

subroutine createpoints(npoints,xd,yd,rhsd,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>

	PetscInt,		intent(in)		::	npoints
	Vec,			intent(inout)	::	xd,yd,rhsd 
	PetscErrorCode,	intent(out)		::	ierr
		
	integer		i,ista,iend	
	PetscReal xmin,xmax,ymin,ymax,delta,xcord,ycord,res
	xmin=0.0	
	xmax=1.0	
	ymin=0.0	
	ymax=1.0	

	delta= (xmax-xmin)/(npoints-1)

	call VecGetOwnershipRange(xd,ista,iend,ierr)
	print *,ista,iend
	do i=ista,iend-1,1
		xcord = xmin+(i/npoints)*delta;
		ycord = ymin+mod(i,npoints)*delta;
		call testfunction(xcord,ycord,res);
		call VecSetValues(xd,1,i,xcord,INSERT_VALUES,ierr)
		call VecSetValues(yd,1,i,ycord,INSERT_VALUES,ierr)
		call VecSetValues(rhsd,1,i,res,INSERT_VALUES,ierr)
	enddo

	call VecAssemblyBegin(xd,ierr)
	call VecAssemblyEnd(xd,ierr)
	call VecAssemblyBegin(yd,ierr)
	call VecAssemblyEnd(yd,ierr)
	call VecAssemblyBegin(rhsd,ierr)
	call VecAssemblyEnd(rhsd,ierr)
end subroutine


subroutine generatepoints(npoints,dsites,rhs,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>

	PetscInt,		intent(in)		::	npoints
	Mat,			intent(out)		::  dsites	
	Vec,			intent(out)		::  rhs 
	PetscErrorCode,	intent(out)		::	ierr
		
	PetscReal	xmin,xmax,ymin,ymax,delta,xcord,ycord,res
	PetscInt	ista,iend
	integer		i

	xmin=0.0	
	xmax=1.0	
	ymin=0.0	
	ymax=1.0	
	
	delta= (xmax-xmin)/(npoints-1)
	
	! generate matrix dsites and ctrs with size M*2 
	call MatCreate(PETSC_COMM_WORLD,dsites,ierr);
	call MatSetSizes(dsites,PETSC_DECIDE,PETSC_DECIDE,npoints*npoints,2,ierr)
	call MatSetFromOptions(dsites,ierr)
    call MatSetUp(dsites,ierr)
		
	call MatGetOwnershipRange(dsites,ista,iend,ierr)
	
	
	do i=ista,iend-1
		xcord = xmin+(i/npoints)*delta;
		ycord = ymin+mod(i,npoints)*delta;
		call testfunction(xcord,ycord,res);
		call MatSetValue(dsites,i,0,xcord,INSERT_VALUES,ierr) 
		call MatSetValue(dsites,i,1,ycord,INSERT_VALUES,ierr) 
		call VecSetValue(rhs,i,res,INSERT_VALUES,ierr)
	enddo

	call MatAssemblyBegin(dsites,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(dsites,MAT_FINAL_ASSEMBLY,ierr)

	call VecAssemblyBegin(rhs,ierr)
	call VecAssemblyEnd(rhs,ierr)

	call MatView(dsites,PETSC_VIEWER_STDOUT_WORLD,ierr)
	call VecView(rhs,PETSC_VIEWER_STDOUT_WORLD,ierr)
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


subroutine distancematrix(dsites,ctrs,dm)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,intent(in)	:: dsites,ctrs
	Mat,intent(out):: dm
	PetscInt	nrow,ncol,rsta,rend	
	Mat			work1,work2
	PetscErrorCode	ierr
	integer		i
	call MatGetSize(dsites,nrow,ncol,ierr)

	call MatGetOwnershipRange(dsites,rsta,rend,ierr)
	print *,rsta,rend
	do i=rsta,rend-1,1
	
	enddo

end subroutine


end module
