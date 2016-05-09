! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Matrix Algebra Library 
! -----------------------------------------------------------------------
#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>

module matrixalgebra

contains

! -----------------------------------------------------------------------
! The repmat function in matrix algebra library is used to replicate matrix
! with m times in row and n times in column. The name of this function is
! refered from MATLAB. 
! -----------------------------------------------------------------------
subroutine ma_repmat1(orgM,m,n,newM,ierr)
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

	call MatAssemblyBegin(newM,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(newM,MAT_FINAL_ASSEMBLY,ierr)
	
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

	call MatAssemblyBegin(newM,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(newM,MAT_FINAL_ASSEMBLY,ierr)
	
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
	integer						::	i

	call MatGetSize(A1,nrow,ncol,ierr)
	call MatGetOwnershipRange(A1,ista,iend,ierr)
	ilocal= iend-ista	
	!print *,">istat=",ista,"iend=",iend
	
	allocate(idxm(nrow),idxn(ncol),row1(ncol),row2(ncol),results(ncol))
	
	call MatDuplicate(A1,MAT_COPY_VALUES,B,ierr)
	
	do i=ista,iend-1
		call MatGetRow(A1,i,num,idxn,row1,ierr)
		call MatRestoreRow(A1,i,num,idxn,row1,ierr)
		call MatGetRow(A2,i,num,idxn,row2,ierr)
		call MatRestoreRow(A2,i,num,idxn,row2,ierr)
		results = row1*row2	
		!print *,">i=",i,"results=",results
		call MatSetValues(B,1,i,ncol,idxn,results,INSERT_VALUES,ierr)
	enddo

	call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatView(B,PETSC_VIEWER_STDOUT_WORLD,ierr)
	
	deallocate(idxm,idxn,row1,row2,results)

end subroutine

! -----------------------------------------------------------------------
! Sum the elements in a matrix along with the row or column 
! -----------------------------------------------------------------------

subroutine ma_sum(A,dims,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	
	Mat,			intent(in)	::	A
	PetscInt,		intent(in)	::	dims
	Mat,			intent(out)	::	B
	PetscErrorCode,	intent(out)	::	ierr
	
	Mat							::  Tmp
	PetscInt					::	nrow,ncol
	PetscInt					::	num
	PetscInt					::  ista,iend

	call MatGetSize(A,nrow,ncol,ierr)
	call MatGetOwnershipRange(A,ista,iend,ierr)
	print *,">istat=",ista,"iend=",iend
	
	
	call MatDuplicate(A,MAT_COPY_VALUES,B,ierr)
	
end subroutine


subroutine ma_ones(A,m,n,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	
	Mat,			intent(out)	::	A
	PetscInt,		intent(in)	::	m,n	
	PetscErrorCode,	intent(out)	::	ierr
	
	PetscInt					::	nrow,ncol
	PetscInt					::	num
	PetscInt					::  ista,iend
	PetscInt,allocatable		::	idxn(:)
	PetscScalar,allocatable		::	row(:),results(:)
	integer 					:: 	i,j
	
	! generate matrix A with size M*M
	call MatCreate(PETSC_COMM_WORLD,A,ierr);
	call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n,ierr)
	call MatSetFromOptions(A,ierr)
	call MatSetUp(A,ierr)

	call MatGetSize(A,nrow,ncol,ierr)
	call MatGetOwnershipRange(A,ista,iend,ierr)
	!print *,">istat=",ista,"iend=",iend,">ncol=",ncol
	allocate(idxn(ncol),row(ncol),results(ncol))

	!call MatAssemblyBegin(T,MAT_FINAL_ASSEMBLY,ierr)
	!call MatAssemblyEnd(T,MAT_FINAL_ASSEMBLY,ierr)
	do i=ista,iend-1
		do j=1,ncol
			idxn(j)=j-1
			row(j)=1.0
		enddo
		call MatSetValues(A,1,i,ncol,idxn,row,INSERT_VALUES,ierr)
	enddo
	call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)

	deallocate(idxn,row,results)
end subroutine

subroutine ma_zeros(A,m,n,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	
	Mat,			intent(out)	::	A
	PetscInt,		intent(in)	::	m,n	
	PetscErrorCode,	intent(out)	::	ierr
	
	PetscInt					::	nrow,ncol
	PetscInt					::	num
	PetscInt					::  ista,iend
	PetscInt,allocatable		::	idxn(:)
	PetscScalar,allocatable		::	row(:),results(:)
	integer 					:: 	i,j
	
	! generate matrix A with size M*M
	call MatCreate(PETSC_COMM_WORLD,A,ierr);
	call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n,ierr)
	call MatSetFromOptions(A,ierr)
	call MatSetUp(A,ierr)

	call MatGetSize(A,nrow,ncol,ierr)
	call MatGetOwnershipRange(A,ista,iend,ierr)
	!print *,">istat=",ista,"iend=",iend,">ncol=",ncol
	allocate(idxn(ncol),row(ncol),results(ncol))

	!call MatAssemblyBegin(T,MAT_FINAL_ASSEMBLY,ierr)
	!call MatAssemblyEnd(T,MAT_FINAL_ASSEMBLY,ierr)
	do i=ista,iend-1
		do j=1,ncol
			idxn(j)=j-1
			row(j)=0.0
		enddo
		call MatSetValues(A,1,i,ncol,idxn,row,INSERT_VALUES,ierr)
	enddo
	call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)

	deallocate(idxn,row,results)
end subroutine


end module
