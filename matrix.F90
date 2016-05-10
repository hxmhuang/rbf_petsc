! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Matrix Algebra Library 
! -----------------------------------------------------------------------
#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>

module matrix

contains

! -----------------------------------------------------------------------
! The repmat function in matrix algebra library is used to replicate matrix
! with m times in row and n times in column. The name of this function is
! refered from MATLAB. P is an extended unit matrix (mM*m) and T is another 
! extended unit matrix (N*nN), we compute using this equation B= P*A*T. 
! For example, if the size of A is M*N=3*2, suppose m=3 and n=2, we have
! P= [1 0 0]
!	 [0 1 0]
!	 [0 0 1]
!    [1 0 0]
!    [0 1 0]
!    [0 0 1]
!    [1 0 0]
!    [0 1 0]
!    [0 0 1]
!
! T= [1 0 1 0]
!	 [0 1 0 1]
! -----------------------------------------------------------------------
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
	PetscInt					::	nrow,ncol

	Mat							:: I1,I2,W

	Mat							:: P, T
	call MatGetSize(A,nrow,ncol,ierr)
	
	call mat_create(I1,nrow,nrow,ierr)
	call mat_create(W,nrow,m*nrow,ierr)
	call mat_create(P,m*nrow,nrow,ierr)
	call mat_create(I2,ncol,ncol,ierr)
	call mat_create(T,ncol,n*ncol,ierr)
	
	call mat_diag(I1,ierr)
    call mat_mhjoin(I1,m,I1,0,W,ierr)
    call MatTranspose(W,MAT_INITIAL_MATRIX,P,ierr)

	call mat_diag(I2,ierr)
	call mat_mhjoin(I2,n,I2,0,T,ierr)
  
	!print *, "=============P==========="
	!call mat_view(P,ierr)
	!print *, "=============A==========="
	!call mat_view(A,ierr)
	!print *, "=============T==========="
	!call mat_view(T,ierr)

	call MatMatMatMult(P,A,T,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,B,ierr)	

	call mat_destroy(I1,ierr)
	call mat_destroy(I2,ierr)
	call mat_destroy(W,ierr)
	call mat_destroy(P,ierr)
	call mat_destroy(T,ierr)
	
	call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)

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
	!print *,">ista=",ista,"iend=",iend
	
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
	print *,">ista=",ista,"iend=",iend
	
	
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
	!print *,">ista=",ista,"iend=",iend,">ncol=",ncol
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
	!print *,">ista=",ista,"iend=",iend,">ncol=",ncol
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
subroutine bk_mat_hjoin(A1,A2,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::  A1,A2 
	Mat,			intent(out)	::	B
	PetscErrorCode,	intent(out)	::	ierr
	PetscInt					::	nrow1,ncol1,nrow2,ncol2,nrow3,ncol3
	PetscInt,allocatable		::	idxm1(:),idxn1(:)
	PetscInt,allocatable		::	idxm2(:),idxn2(:)
	PetscInt,allocatable		::	idxm3(:),idxn3(:)
	PetscScalar,allocatable		::	row1(:),row2(:),row3(:)
	PetscInt					::  ista,iend,ilocal,indx
	integer						::	i,ib1,ie1,ib2,ie2,ib3,ie3

	call MatGetSize(A1,nrow1,ncol1,ierr)
	call MatGetSize(A2,nrow2,ncol2,ierr)
	if(nrow1/=nrow2)then
		print *, "Error: Matrix A1 and Matrix A2 should have the same row size"
		stop	
	endif
	
	nrow3=nrow1
	ncol3=ncol1+ncol2

	call MatGetOwnershipRange(A1,ista,iend,ierr)
	ilocal=iend-ista
	!print *,">ista=",ista,"iend=",iend,"ilocal=",ilocal,"nrow1=",nrow1,"ncol1=",ncol1

	allocate(idxm1(ilocal),idxn1(ncol1),row1(ilocal*ncol1))
	allocate(idxm2(ilocal),idxn2(ncol2),row2(ilocal*ncol2))
	allocate(idxm3(ilocal),idxn3(ncol3),row3(ilocal*ncol3))

	do i=ista,iend-1
		indx=i-ista+1
		idxm1(indx)=i
		idxm2(indx)=i
		idxm3(indx)=i
	enddo
	do i=1,ncol1
		idxn1(i)=i-1
	enddo
	do i=1,ncol2
		idxn2(i)=i-1
	enddo
	do i=1,ncol3
		idxn3(i)=i-1
	enddo
	!print *,">idxm1=",idxm1,"idxn1=",idxn1

	call MatGetValues(A1,ilocal,idxm1,ncol1,idxn1,row1,ierr)
	call MatGetValues(A2,ilocal,idxm2,ncol2,idxn2,row2,ierr)

	do i=1,ilocal
		ib1=(i-1)*ncol1+1	
		ie1=i*ncol1	
		ib2=(i-1)*ncol2+1	
		ie2=i*ncol2

		ib3=(i-1)*ncol3+1
		ie3=ib3+ncol1-1
		row3(ib3:ie3)=row1(ib1:ie1) 
		!print *,">row1(",ib1,":",ie1,")=",row1(ib1:ie1)
		!print *,">row3(",ib3,":",ie3,")=",row3(ib3:ie3)

		ib3=ie3+1
		ie3=ib3+ncol2-1	
		row3(ib3:ie3)=row2(ib2:ie2) 
		!print *,">row2(",ib2,":",ie2,")=",row2(ib2:ie2)
		!print *,">row3(",ib3,":",ie3,")=",row3(ib3:ie3)
	enddo

	call MatSetValues(B,ilocal,idxm3,ncol3,idxn3,row3,INSERT_VALUES,ierr)

	call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)

	deallocate(idxm1,idxn1,row1)
	deallocate(idxm2,idxn2,row2)
	deallocate(idxm3,idxn3,row3)
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
	!print *,">ista=",ista,"iend=",iend

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

! -----------------------------------------------------------------------
! The mat_mhjoin function will extend A1,A2 with a factor of m,n folds respectively, 
! and then combine two matrixs into one matrix along with x direction.
! For example, using mat_mhjion(A1,3,A2,2,B,ierr), we can get
!	B= [A1,A1,A1,A2,A2] 
! -----------------------------------------------------------------------
subroutine mat_mhjoin(A1,m,A2,n,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,				 intent(in)	::  A1 
	PetscInt,			 intent(in)	::  m,n
	Mat,				 intent(in)	::  A2 
	Mat,				 intent(out)::	B
	PetscErrorCode,		 intent(out)::	ierr
	PetscInt						::	nrow1,ncol1,nrow2,ncol2,nrow3,ncol3
	PetscInt,allocatable			::	idxn(:)
	PetscScalar,allocatable			::	row1(:),row2(:),row3(:)
	PetscInt						::  ista,iend
	PetscInt						::  ncycle1,ncycle2
	integer							::	i,j,k
	integer							::	ib,ie

	ncycle1=m
	ncycle2=n

	call MatGetSize(A1,nrow1,ncol1,ierr)
	call MatGetSize(A2,nrow2,ncol2,ierr)
	if(nrow1/=nrow2)then
		print *, "Error: Matrix A1 and Matrix A2 should have the same row size"
		stop	
	endif
	nrow3=nrow1
	ncol3=ncol1*ncycle1+ncol2*ncycle2

	call MatGetOwnershipRange(A1,ista,iend,ierr)
	!print *,">ista=",ista,"iend=",iend

	allocate(idxn(ncol3),row1(ncol1),row2(ncol2),row3(ncol3))

	do i=ista,iend-1
		call MatGetRow(A1,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row1,ierr)
		call MatRestoreRow(A1,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row1,ierr)
		call MatGetRow(A2,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row2,ierr)
		call MatRestoreRow(A2,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row2,ierr)
		
		do j=1,ncol3
			idxn(j)=j-1
		enddo

		do k=1,ncycle1
			ib=(k-1)*ncol1+1
			ie=ib+ncol1-1
			row3(ib:ie)=row1
			!print *,">row3(",ib,":",ie,")=",row1
		enddo

		do k=1,ncycle2
			ib=ie+1
			ie=ib+ncol2-1
			row3(ib:ie)=row2
			!print *,">row3(",ib,":",ie,")=",row2
		enddo

		call MatSetValues(B,1,i,ncol3,idxn,row3,INSERT_VALUES,ierr)
	enddo

	call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)

	deallocate(idxn,row1,row2,row3)
end subroutine


subroutine mat_seq(A,ierr)
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
	!print *,">ista=",ista,"iend=",iend,">ncol=",ncol
	allocate(idxn(ncol),row(ncol),results(ncol))

	!call MatAssemblyBegin(T,MAT_FINAL_ASSEMBLY,ierr)
	!call MatAssemblyEnd(T,MAT_FINAL_ASSEMBLY,ierr)
	do i=ista,iend-1
		do j=1,ncol
			idxn(j)=j-1
			row(j)=i*ncol+j
		enddo
		call MatSetValues(A,1,i,ncol,idxn,row,INSERT_VALUES,ierr)
	enddo
	call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

	deallocate(idxn,row,results)
end subroutine


end module
