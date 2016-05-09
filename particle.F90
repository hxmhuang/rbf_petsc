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

! -----------------------------------------------------------------------
! The repmat function in matrix algebra library is used to replicate matrix
! with m times in row and n times in column. The name of this function is
! refered from MATLAB. 
! -----------------------------------------------------------------------
subroutine ma_repmat(orgM,m,n,newM,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::  orgM	
	PetscInt,		intent(in)	::	m,n
	Mat,			intent(out)	::	newM 
	PetscErrorCode,	intent(out)	::	ierr
	PetscInt					::	nrow,ncol
	PetscInt					::	num
	PetscInt,allocatable		::	idxm(:),idxn(:)
	PetscScalar,allocatable		::	row(:),rows(:),allrows(:)
	PetscInt					::  ista,iend,ilocal
	integer						::	i,j,k

	call MatGetSize(orgM,nrow,ncol,ierr)
	call MatGetOwnershipRange(orgM,ista,iend,ierr)
	ilocal= iend-ista	
	!print *,">istat=",ista,"iend=",iend
	
	allocate(idxm(nrow),idxn(ncol),row(ncol),rows(ilocal*ncol),allrows(nrow*ncol))
	
	do i=ista,iend-1
		call MatGetRow(orgM,i,num,idxn,row,ierr)
		j=mod(i,ilocal)		
		rows((j*ncol+1):) = row 
		call MatRestoreRow(orgM,i,num,idxn,row,ierr)
	enddo
	call mpi_allgather(rows,ilocal*ncol,MPIU_SCALAR,allrows,ilocal*ncol,MPIU_SCALAR,PETSC_COMM_WORLD,ierr)
	!print *,">>allrows=",allrows
	
	call MatCreate(PETSC_COMM_WORLD,newM,ierr)
	call MatSetSizes(newM,PETSC_DECIDE,PETSC_DECIDE,nrow*m,ncol*n,ierr)
	call MatSetFromOptions(newM,ierr)
	call MatSetUp(newM,ierr)

	call MatAssemblyBegin(newM,MAT_FLUSH_ASSEMBLY,ierr)
	call MatAssemblyEnd(newM,MAT_FLUSH_ASSEMBLY,ierr)
	
	do j=1,n
		do i=1,m
			do k=1,nrow
				idxm(k)=(i-1)*nrow+k-1
			enddo
			do k=1,ncol
				idxn(k)=(j-1)*ncol+k-1
			enddo
			call MatSetValues(newM,nrow,idxm,ncol,idxn,allrows,INSERT_VALUES,ierr)
		enddo
	enddo

	call MatAssemblyBegin(newM,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(newM,MAT_FINAL_ASSEMBLY,ierr)
	call MatView(newM,PETSC_VIEWER_STDOUT_WORLD,ierr)
	
	deallocate(idxm,idxn,row,rows,allrows)

end subroutine

! -----------------------------------------------------------------------
! The eprod function in matrix algebra library is used to implement the
! puoduct of elements with the same positions in two matrixs. 
! -----------------------------------------------------------------------
subroutine ma_eprod(A1,A2,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::  A1,A2 
	Mat,			intent(out)	::	B
	PetscErrorCode,	intent(out)	::	ierr
	PetscInt					::	nrow,ncol
	PetscInt					::	num
	PetscInt,allocatable		::	idxm(:),idxn(:)
	PetscScalar,allocatable		::	row1(:),row2(:),results(:)
	PetscInt					::  ista,iend,ilocal
	integer						::	i,j,k

	call MatGetSize(A1,nrow,ncol,ierr)
	call MatGetOwnershipRange(A1,ista,iend,ierr)
	ilocal= iend-ista	
	print *,">istat=",ista,"iend=",iend
	
	allocate(idxm(nrow),idxn(ncol),row1(ncol),row2(ncol),results(ncol))
	
	call MatDuplicate(A1,MAT_COPY_VALUES,B,ierr)
	
	do i=ista,iend-1
		call MatGetRow(A1,i,num,idxn,row1,ierr)
		call MatRestoreRow(A1,i,num,idxn,row1,ierr)
		call MatGetRow(A2,i,num,idxn,row2,ierr)
		call MatRestoreRow(A2,i,num,idxn,row2,ierr)
		results = row1*row2	
		!print *,">i=",i,"row1=",row1
		!print *,">i=",i,"row2=",row2
		print *,">i=",i,"results=",results
	
		!call MatSetValues(B,1,ista,ncol,idxn,results,INSERT_VALUES,ierr)
		call MatSetValues(B,1,i,ncol,idxn,results,INSERT_VALUES,ierr)
	enddo

	call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatView(B,PETSC_VIEWER_STDOUT_WORLD,ierr)
	
	deallocate(idxm,idxn,row1,row2,results)

end subroutine


end module
