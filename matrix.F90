! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Matrix Algebra Library 
! -----------------------------------------------------------------------
#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>

module matrix

contains

subroutine mat_rep(A,m,n,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	
	Mat,			intent(in)	::  A 
	PetscInt,		intent(in)	::	m,n
	Mat,			intent(out)	::	B	
	PetscErrorCode,	intent(out)	::	ierr
	PetscInt					::	nrow1,ncol1,nrow2,ncol2
	PetscInt					::  ista1,iend1,ilocal1
	!PetscInt					::  ista2,iend2,ilocal2
	!PetscScalar,allocatable		::	row(:),rows(:)
	!integer						::	i,j,k

	call MatGetSize(A,nrow1,ncol1,ierr)
	nrow2=m*nrow1
	ncol2=n*ncol2

	! generate matrix B with size M*M
	call MatCreate(PETSC_COMM_WORLD,B,ierr);
	call MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,nrow2,ncol2,ierr)
	call MatSetFromOptions(B,ierr)
	call MatSetUp(B,ierr)

	call MatGetOwnershipRange(A,ista1,iend1,ierr)
	ilocal1= iend1-ista1	
	!print *,">istat=",ista,"iend=",iend
end subroutine


! -----------------------------------------------------------------------
! The repmat function in matrix algebra library is used to replicate matrix
! with m times in row and n times in column. The name of this function is
! refered from MATLAB. 
! -----------------------------------------------------------------------
subroutine mat_repmat1(orgM,m,n,newM,ierr)
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
	print *,">istat=",ista,"iend=",iend
	
	allocate(idxm(nrow),idxn(ncol),row(ncol),rows(ilocal*ncol),allrows(nrow*ncol))
	
	do i=ista,iend-1
		call MatGetRow(orgM,i,num,idxn,row,ierr)
		j=mod(i,ilocal)		
		rows((j*ncol+1):) = row 
		call MatRestoreRow(orgM,i,num,idxn,row,ierr)
	enddo
	call mpi_allgather(rows,ilocal*ncol,MPIU_SCALAR,allrows,ilocal*ncol,MPIU_SCALAR,PETSC_COMM_WORLD,ierr)
	print *,">>allrows=",allrows
	
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
	
	deallocate(idxm,idxn,row,rows,allrows)

end subroutine

! -----------------------------------------------------------------------
! The eprod function in matrix algebra library is used to implement the
! puoduct of elements with the same positions in two matrixs. 
! -----------------------------------------------------------------------
subroutine mat_eprod(A1,A2,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::  A1,A2 
	Mat,			intent(out)	::	B
	PetscErrorCode,	intent(out)	::	ierr
	PetscInt					::	nrow1,ncol1,nrow2,ncol2
	PetscInt					::	num
	PetscInt,allocatable		::	idxm(:),idxn(:)
	PetscScalar,allocatable		::	row1(:),row2(:),results(:)
	PetscInt					::  ista,iend,ilocal
	integer						::	i

	call MatGetSize(A1,nrow1,ncol1,ierr)
	call MatGetSize(A2,nrow2,ncol2,ierr)
	if(nrow1/=nrow2 .or. ncol1/=ncol2)then
		print *, "Error: Matrix A1 and Matrix A2 should have the same sizes"
		stop	
	endif

	call MatGetOwnershipRange(A1,ista,iend,ierr)
	ilocal= iend-ista	
	!print *,">istat=",ista,"iend=",iend
	
	allocate(idxm(nrow1),idxn(ncol1),row1(ncol1),row2(ncol1),results(ncol1))
	
	do i=ista,iend-1
		call MatGetRow(A1,i,num,idxn,row1,ierr)
		call MatRestoreRow(A1,i,num,idxn,row1,ierr)
		call MatGetRow(A2,i,num,idxn,row2,ierr)
		call MatRestoreRow(A2,i,num,idxn,row2,ierr)
		results = row1*row2	
		!print *,">i=",i,"results=",results
		call MatSetValues(B,1,i,ncol1,idxn,results,INSERT_VALUES,ierr)
	enddo

	call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
	
	deallocate(idxm,idxn,row1,row2,results)

end subroutine

! -----------------------------------------------------------------------
! Sum the elements in a matrix along with the row or column 
! -----------------------------------------------------------------------

subroutine mat_sum(A,dims,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	
	Mat,			intent(in)	::	A
	PetscInt,		intent(in)	::	dims
	Mat,			intent(out)	::	B
	PetscErrorCode,	intent(out)	::	ierr
	
	PetscInt					::	nrow,ncol
	PetscInt					::  ista,iend

	call MatGetSize(A,nrow,ncol,ierr)
	call MatGetOwnershipRange(A,ista,iend,ierr)
	print *,">istat=",ista,"iend=",iend
	
	
end subroutine


subroutine mat_ones(A,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	
	Mat,			intent(out)	::	A
	PetscErrorCode,	intent(out)	::	ierr
	
	PetscInt					::	nrow,ncol
	PetscInt					::  ista,iend
	PetscInt,allocatable		::	idxn(:)
	PetscScalar,allocatable		::	row(:),results(:)
	integer 					:: 	i,j
	
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

	deallocate(idxn,row,results)
end subroutine

subroutine mat_zeros(A,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::	A
	PetscErrorCode,	intent(out)	::	ierr
	PetscInt					::	nrow,ncol
	PetscInt					::  ista,iend
	PetscInt,allocatable		::	idxn(:)
	PetscScalar,allocatable		::	row(:)
	integer 					:: 	i,j
	
	call MatGetSize(A,nrow,ncol,ierr)
	call MatGetOwnershipRange(A,ista,iend,ierr)
	!print *,">istat=",ista,"iend=",iend,">ncol=",ncol
	allocate(idxn(ncol),row(ncol))

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

	deallocate(idxn,row)
end subroutine


subroutine mat_diag(A,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::	A
	PetscErrorCode,	intent(out)	::	ierr
	PetscInt					::	nrow,ncol	
	PetscInt					::  ista,iend
	PetscInt,allocatable		::	idxn(:)
	PetscScalar,allocatable		::	row(:)
	integer 					:: 	i,j
	
	! generate matrix A with size m*n
	call MatGetSize(A,nrow,ncol,ierr)
	if(nrow /= ncol) then
		print *, "Error: Matrix should be square"
		stop	
	endif
	call MatGetOwnershipRange(A,ista,iend,ierr)
	allocate(idxn(ncol),row(ncol))

	do i=ista,iend-1
		do j=1,ncol
			idxn(j)=j-1
			if (i==(j-1)) then
				row(j)=1.0
			else
				row(j)=0.0
			endif
		enddo
		call MatSetValues(A,1,i,ncol,idxn,row,INSERT_VALUES,ierr)
	enddo
	call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

	deallocate(idxn,row)
end subroutine


subroutine mat_create(A,m,n,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	PetscInt,		intent(in)	::	m,n	
	Mat,			intent(out)	::	A
	PetscErrorCode,	intent(out)	::	ierr
	! generate matrix A with size m*n
	call MatCreate(PETSC_COMM_WORLD,A,ierr);
	call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n,ierr)
	call MatSetFromOptions(A,ierr)
	call MatSetUp(A,ierr)
end subroutine


subroutine mat_copy(A,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::	A
	Mat,			intent(out)	::	B
	PetscErrorCode,	intent(out)	::	ierr
	! veiw matrix A
	call MatDuplicate(A,MAT_COPY_VALUES,B,ierr)
end subroutine


subroutine mat_destroy(A,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::	A
	PetscErrorCode,	intent(out)	::	ierr
	! destroy matrix A
	call MatDestroy(A,ierr)
end subroutine


subroutine mat_view(A,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::	A
	PetscErrorCode,	intent(out)	::	ierr
	! veiw matrix A
	call MatView(A,PETSC_VIEWER_STDOUT_WORLD,ierr)
end subroutine


! -----------------------------------------------------------------------
! The mat_hjoin function is used to combine two matrixs into one matrix along with x direction 
! -----------------------------------------------------------------------
subroutine mat_hjoin(A1,A2,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::  A1,A2 
	Mat,			intent(out)	::	B
	PetscErrorCode,	intent(out)	::	ierr
	PetscInt					::	nrow1,ncol1,nrow2,ncol2
	PetscInt,allocatable		::	idxn(:)
	PetscScalar,allocatable		::	row1(:),row2(:),results(:)
	PetscInt					::  ista,iend
	integer						::	i,j

	call MatGetSize(A1,nrow1,ncol1,ierr)
	call MatGetSize(A2,nrow2,ncol2,ierr)
	if(nrow1/=nrow2)then
		print *, "Error: Matrix A1 and Matrix A2 should have the same row size"
		stop	
	endif

	call MatGetOwnershipRange(A1,ista,iend,ierr)
	!print *,">istat=",ista,"iend=",iend

	allocate(idxn(ncol1+ncol2),row1(ncol1),row2(ncol2),results(ncol1+ncol2))

	do i=ista,iend-1
		call MatGetRow(A1,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row1,ierr)
		call MatRestoreRow(A1,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row1,ierr)
		call MatGetRow(A2,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row2,ierr)
		call MatRestoreRow(A2,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row2,ierr)
		do j=1,(ncol1+ncol2)
			idxn(j)=j-1
		enddo
		results(1:ncol1) = row1
		results((ncol1+1):(ncol1+ncol2)) = row2
		!print *,">i=",i,"results=",results
		call MatSetValues(B,1,i,(ncol1+ncol2),idxn,results,INSERT_VALUES,ierr)
	enddo

	call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)

	deallocate(idxn,row1,row2,results)
end subroutine



end module
