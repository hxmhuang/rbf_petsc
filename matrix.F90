! -----------------------------------------------------------------------
! Matrix Algebra Library 
! -----------------------------------------------------------------------
#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>
#include "mat_math_type.h"

module matrix

contains


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


! -----------------------------------------------------------------------
! A=0 
! -----------------------------------------------------------------------
subroutine mat_zeros(A,m,n,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(out)	::	A
	PetscInt,	    intent(in)	::	m,n
	PetscErrorCode,	intent(out)	::	ierr
	
    call mat_create(A,m,n,ierr)
    call MatZeroEntries(A,ierr)
	call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
end subroutine


! -----------------------------------------------------------------------
! A=1. Since it is not sparse, we will remove this interface later 
! -----------------------------------------------------------------------
subroutine mat_ones(A,m,n,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(out)	::	A
	PetscInt,	    intent(in)	::	m,n
	PetscErrorCode,	intent(out)	::	ierr
	
	PetscInt					::  ista,iend
	PetscInt,allocatable		::	idxn(:)
	PetscScalar,allocatable		::	row(:),results(:)
	integer 					:: 	i,j
	
    call mat_create(A,m,n,ierr)
	call MatGetOwnershipRange(A,ista,iend,ierr)
	!print *,">ista=",ista,"iend=",iend,">ncol=",ncol
	allocate(idxn(n),row(n),results(n))

	!call MatAssemblyBegin(T,MAT_FINAL_ASSEMBLY,ierr)
	!call MatAssemblyEnd(T,MAT_FINAL_ASSEMBLY,ierr)
	do i=ista,iend-1
		do j=1,n
			idxn(j)=j-1
			row(j)=1.0
		enddo
		call MatSetValues(A,1,i,n,idxn,row,INSERT_VALUES,ierr)
	enddo
	call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

	deallocate(idxn,row,results)
end subroutine


! -----------------------------------------------------------------------
! A=[1 2 3], This function is only used to generate the test data.
!   [4 5 6]
!   [7 8 9]
! -----------------------------------------------------------------------
subroutine mat_seq(A,m,n,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	
	Mat,			intent(out)	::	A
	PetscInt,	    intent(in)	::	m,n
	PetscErrorCode,	intent(out)	::	ierr
	
	PetscInt					::  ista,iend
	PetscInt,allocatable		::	idxn(:)
	PetscScalar,allocatable		::	row(:)
	integer 					:: 	i,j
	
    call mat_create(A,m,n,ierr)
	call MatGetOwnershipRange(A,ista,iend,ierr)
	!print *,">ista=",ista,"iend=",iend,">ncol=",ncol
	allocate(idxn(n),row(n))

	do i=ista,iend-1
		do j=1,n
			idxn(j)=j-1
			row(j)=i*n+j
		enddo
		call MatSetValues(A,1,i,n,idxn,row,INSERT_VALUES,ierr)
	enddo
	call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

	deallocate(idxn,row)
end subroutine


! -----------------------------------------------------------------------
! The eye function is used to generate the simple and complex identity matrixs. 
! For example, if A is a 2*6 matrix, we can use mat_eye(A,ierr) to obtain 
! A= [1 0 1 0 1 0]
!	 [0 1 0 1 1 0]
! if A is a 6*2 matrix, then mat_eye(A,ierr) will generate
! A= [1 0]
!	 [0 1]
!	 [1 0]
!    [0 1]
!	 [1 0]
!    [0 1]
! -----------------------------------------------------------------------
subroutine mat_eye(A,m,n,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(out)	::	A
	PetscErrorCode,	intent(out)	::	ierr
	PetscInt,       intent(in)	::	m,n	
	PetscInt					::	nmax, nmin 
	PetscInt					::  ista,iend,xpos,ypos
	PetscInt,allocatable		::	idxn(:)
	PetscScalar,allocatable		::	row(:)
	integer 					:: 	i,j,counter
	
	call mat_create(A,m,n,ierr)
    nmin=min(m,n)
	nmax=max(m,n)
 	if(mod(nmax,nmin) /= 0) then
		print *, "Error in mat_eye: the size of matrix A should be (NM)*M or M*(NM)"
		stop	
	endif
	call MatGetOwnershipRange(A,ista,iend,ierr)
	allocate(idxn(nmax/nmin),row(nmax/nmin))

    do i=ista,iend-1
    	xpos=mod(i,nmin)

        counter=0	
		do j=1,n
    		ypos=mod(j-1,nmin)
    		if(ypos==xpos) then
                counter=counter+1
    			row(counter)=1.0
			    idxn(counter)=j-1
    		    !print *,"i=",i,"j=",j,"xpos=",xpos,"ypos=",ypos,"idxn(",counter,")=",idxn(counter),"row(",j,")=",row(counter) 
    		endif
		enddo
	   	call MatSetValues(A,1,i,counter,idxn,row,INSERT_VALUES,ierr)
	enddo
		
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
	deallocate(idxn,row)
end subroutine


! -----------------------------------------------------------------------
! B=A 
! -----------------------------------------------------------------------
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


! -----------------------------------------------------------------------
! B=A1+A2
! -----------------------------------------------------------------------
subroutine mat_add(A1,A2,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::  A1,A2 
	Mat,			intent(out)	::	B
	PetscErrorCode,	intent(out)	::	ierr
	PetscInt					::	nrow1,ncol1,nrow2,ncol2
    
	call MatGetSize(A1,nrow1,ncol1,ierr)
	call MatGetSize(A2,nrow2,ncol2,ierr)
	if(nrow1/=nrow2 .or. ncol1/=ncol2)then
		print *, "Error in mat_add: the matrix A1 and A2 should have the same size."
		stop	
	endif

    call MatCreate(PETSC_COMM_WORLD,B,ierr)
    call MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,nrow1,ncol1,ierr)
    call MatSetType(B,MATCOMPOSITE,ierr)
	call MatSetUp(B,ierr)
    !call MatCompositeSetType(B,MAT_COMPOSITE_MULTIPLICATIVE,ierr)
    
    call MatCompositeAddMat(B,A2,ierr)
    call MatCompositeAddMat(B,A1,ierr)

    call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
end subroutine


! -----------------------------------------------------------------------
! B=[A1 A2] 
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
	PetscInt					::	col1,col2,m,n
	PetscInt,allocatable		::	idxn1(:),idxn2(:),idxn3(:)
	PetscScalar,allocatable		::	row1(:),row2(:),row3(:)
	PetscInt					::  ista,iend
	integer						::	i
	call MatGetSize(A1,nrow1,ncol1,ierr)
	call MatGetSize(A2,nrow2,ncol2,ierr)
	if(nrow1/=nrow2)then
		print *, "Error: Matrix A1 and Matrix A2 should have the same row size"
		stop	
	endif

    call mat_create(B,nrow1,(ncol1+ncol2),ierr)
	call MatGetOwnershipRange(A1,ista,iend,ierr)

	do i=ista,iend-1
	    call MatGetRow(A1,i,col1,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
        m=col1
	    call MatRestoreRow(A1,i,col1,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
	    
        call MatGetRow(A2,i,col2,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
		n=col2
	    call MatRestoreRow(A2,i,col2,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
	    
	    allocate(idxn1(m),row1(m))
	    allocate(idxn2(n),row2(n))
	    allocate(idxn3(m+n),row3(m+n))
        
        call MatGetRow(A1,i,col1,idxn1,row1,ierr)
        idxn3(1:m)=idxn1
        row3(1:m)=row1
		call MatRestoreRow(A1,i,col1,idxn1,row1,ierr)
        
        call MatGetRow(A2,i,col2,idxn2,row2,ierr)
        idxn3((m+1):(m+n))=ncol1+idxn2
        row3((m+1):(m+n))=row2
		call MatRestoreRow(A2,i,col2,idxn2,row2,ierr)
		
		!print *,">i=",i,"idxn3=",idxn3,"row3=",row3
		call MatSetValues(B,1,i,(m+n),idxn3,row3,INSERT_VALUES,ierr)
	    deallocate(idxn1,idxn2,idxn3,row1,row2,row3)
	enddo

	call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
end subroutine


! -----------------------------------------------------------------------
! Compute B=A1*A2
! -----------------------------------------------------------------------
subroutine mat_mult(A1,A2,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::  A1,A2 
	Mat,			intent(out)	::	B	
	PetscErrorCode,	intent(out)	::	ierr
	PetscInt					::	nrow1,ncol1,nrow2,ncol2
	
	call MatGetSize(A1,nrow1,ncol1,ierr)
    call MatGetSize(A2,nrow2,ncol2,ierr)
 	if(ncol1/=nrow2)then
 		print *, "Error in mat_mult: the column of A1 matrix should equal to the row of A2 matrix."
 		stop	
 	endif
    call MatMatMult(A1,A2,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,B,ierr) 
end subroutine


! -----------------------------------------------------------------------
! B=A1.*A2
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
	PetscInt					::	col1,col2,m,n
	PetscInt,allocatable		::	idxn1(:),idxn2(:),idxn3(:),idxtmp(:)
	PetscScalar,allocatable		::	row1(:),row2(:),row3(:),rowtmp(:)
	PetscInt					::  ista,iend
    PetscInt                    ::  pos1,pos2,counter
	integer						::	i
	call MatGetSize(A1,nrow1,ncol1,ierr)
	call MatGetSize(A2,nrow2,ncol2,ierr)
	if(nrow1/=nrow2 .or. ncol1/=ncol2)then
		print *, "Error: Matrix A1 and Matrix A2 should have the same sizes"
		stop	
	endif
    
    call mat_create(B,nrow1,ncol1,ierr)
	call MatGetOwnershipRange(A1,ista,iend,ierr)
	    
	
    do i=ista,iend-1
	    call MatGetRow(A1,i,col1,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
        m=col1
	    call MatRestoreRow(A1,i,col1,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
	    
        call MatGetRow(A2,i,col2,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
		n=col2
	    call MatRestoreRow(A2,i,col2,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
	    
        allocate(idxn1(m),row1(m))
        allocate(idxtmp(m),rowtmp(m))
	    allocate(idxn2(n),row2(n))
	    allocate(idxn3(min(m,n)),row3(min(m,n)))
    	
        call MatGetRow(A1,i,col1,idxn1,row1,ierr)
        m=col1
        idxtmp=idxn1
        rowtmp=row1
	    call MatRestoreRow(A1,i,col1,idxn1,PETSC_NULL_SCALAR,ierr)
 
        call MatGetRow(A2,i,col2,idxn2,row2,ierr)
        counter=0
        pos1=1
        pos2=1
        do while(pos1<=m .and. pos2<=n)
            if(idxtmp(pos1)<idxn2(pos2)) then
                pos1=pos1+1
            elseif(idxtmp(pos1)==idxn2(pos2))then
                counter=counter+1
                idxn3(counter)=idxn2(pos2)
                row3(counter)=rowtmp(pos1)*row2(pos2)
                pos1=pos1+1
                pos2=pos2+1
            else
                pos2=pos2+1
            endif
        enddo
        call MatRestoreRow(A2,i,col2,idxn2,row2,ierr)
           
		!print *,">i=",i,"idxn3=",idxn3,"row3=",row3
		call MatSetValues(B,1,i,counter,idxn3(1:counter),row3(1:counter),INSERT_VALUES,ierr)
	    deallocate(idxn1,idxn2,idxn3,idxtmp,row1,row2,row3,rowtmp)
	enddo

	call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
end subroutine


! -----------------------------------------------------------------------
! The mat_rep function is used to replicate matrix with m times in row and
! n times in column. It's name refers to the repmat function in MATLAB. 
! Suppose P is an extended identity mM*M matrix and T is another 
! extended N*nN identity matrix, we can compute W=P*A firstly and then compute
! B=W*T. These two stpes are faster than computing B=P*A*T directly.
! If the size of A is M*N=3*2, suppose m=3 and n=2, we have
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

	Mat							:: P,T,W

	call MatGetSize(A,nrow,ncol,ierr)
	
    call mat_eye(P,m*nrow,nrow,ierr)
    call mat_mult(P,A,W,ierr) 

	call mat_eye(T,ncol,n*ncol,ierr)
    call mat_mult(W,T,B,ierr) 

	call mat_destroy(P,ierr)
	call mat_destroy(W,ierr)
	call mat_destroy(T,ierr)
	
    call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
end subroutine


! -----------------------------------------------------------------------
! Sum of elements along with the row or column.
! Suppose A=[1,2,3]
!           [4,5,6],
! then mat_sum(A,1,B) will make B=[5,7,9],
!      mat_sum(A,2,B) will make B=[6 ]
!                                 [15]
! -----------------------------------------------------------------------
subroutine mat_sum(A,ndim,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::	A
	PetscInt,       intent(in)  ::  ndim	
	Mat,			intent(out)	::	B
	PetscErrorCode,	intent(out)	::	ierr
	Mat				            ::	W
	PetscInt					::	nrow,ncol	
	
 	if((ndim<1) .or. (ndim>2)) then
		print *, "Error in mat_sum: the dim should be 1 or 2"
		stop	
	endif

    call MatGetSize(A,nrow,ncol,ierr)
    if(ndim==1) then
        call mat_ones(W,1,nrow,ierr)
        call mat_mult(W,A,B,ierr)
        call mat_destroy(W,ierr)
    elseif(ndim==2) then
        call mat_ones(W,ncol,1,ierr)
        call mat_mult(A,W,B,ierr)
        call mat_destroy(W,ierr)
    endif

	call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
end subroutine


! -----------------------------------------------------------------------
! Compute Y = a*X + Y.
! -----------------------------------------------------------------------
subroutine mat_axpy(Y,a,X,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	    ::  X 
	PetscScalar,    intent(in)	    ::	a
	Mat,			intent(inout)   ::	Y	
	PetscErrorCode,	intent(out)	    ::	ierr

	call MatAXPY(Y,a,X,DIFFERENT_NONZERO_PATTERN,ierr)
end subroutine

! -----------------------------------------------------------------------
! Compute Y = a*Y + X.
! -----------------------------------------------------------------------
subroutine mat_aypx(Y,a,X,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	    ::  X 
	PetscScalar,    intent(in)	    ::	a
	Mat,			intent(inout)   ::	Y	
	PetscErrorCode,	intent(out)	    ::	ierr

    call MatAYPX(Y,a,X,DIFFERENT_NONZERO_PATTERN,ierr)
end subroutine

! -----------------------------------------------------------------------
! B = A^T.
! -----------------------------------------------------------------------
subroutine mat_trans(A,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	    ::  A 
	Mat,			intent(out)     ::	B	
	PetscErrorCode,	intent(out)	    ::	ierr

    call MatTranspose(A,MAT_INITIAL_MATRIX,B,ierr)
end subroutine


! -----------------------------------------------------------------------
! B = X*Y^T
! -----------------------------------------------------------------------
subroutine mat_xyt(X,Y,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	    ::  X,Y 
	Mat,			intent(out)     ::	B	
	Mat 			                ::	W	
	PetscErrorCode,	intent(out)	    ::	ierr
    !call MatMatTransposeMult(X,Y,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,B,ierr) 	
    call mat_trans(Y,W,ierr)
    call mat_mult(X,W,B,ierr)
    call mat_destroy(W,ierr)
end subroutine


! -----------------------------------------------------------------------
! B = X^T*Y
! -----------------------------------------------------------------------
subroutine mat_xty(X,Y,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	    ::  X,Y 
	Mat,			intent(out)     ::	B	
	PetscErrorCode,	intent(out)	    ::	ierr
    call MatTransposeMatMult(X,Y,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,B,ierr)
end subroutine


! -----------------------------------------------------------------------
! X = a*X
! -----------------------------------------------------------------------
subroutine mat_scale(X,a,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(inout)	::  X 
	PetscScalar,    intent(in)     ::	a	
	PetscErrorCode,	intent(out)	    ::	ierr
    call MatScale(X,a,ierr)
end subroutine


! -----------------------------------------------------------------------
! B=fun(A,opt) 
! opt options: exp,log,sin,cos,tan
! -----------------------------------------------------------------------
subroutine mat_math(A,opt,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include "mat_math_type.h"
	Mat,			intent(in)	::  A
	Integer,        intent(in)  ::  opt
	Mat,			intent(out)	::	B
	PetscErrorCode,	intent(out)	::	ierr
    
	PetscInt					::	nrow,ncol
	PetscInt					::	col,m
	PetscInt,allocatable        ::	idxn(:),idxtmp(:)
	PetscScalar,allocatable     ::	row(:),rowtmp(:)
	PetscInt					::  ista,iend
	integer						::	i,j
    PetscScalar                 ::  newval 
    
    call MatGetSize(A,nrow,ncol,ierr)

    call mat_create(B,nrow,ncol,ierr)
	
    call MatGetOwnershipRange(A,ista,iend,ierr)
	    
    do i=ista,iend-1
        call MatGetRow(A,i,col,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
        m=col
        call MatRestoreRow(A,i,col,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
        
        allocate(idxn(m),idxtmp(m),row(m),rowtmp(m))

        call MatGetRow(A,i,col,idxtmp,rowtmp,ierr)
        m=col
        idxn=idxtmp
        do j=1,m
            select case(opt)
                case (MAT_MATH_EXP)
                    row=exp(rowtmp)
                case (MAT_MATH_SQRT)
                    row=sqrt(rowtmp)
                case (MAT_MATH_LOG)
                    row=log(rowtmp)
                case (MAT_MATH_SIN)
                    row=sin(rowtmp)
                case (MAT_MATH_COS)
                    row=cos(rowtmp)
                case (MAT_MATH_TAN)
                    row=tan(rowtmp)
                case default
                    newval=0.0    
            end select
        enddo
        call MatRestoreRow(A,i,col,idxtmp,rowtmp,ierr)

    	!print *,">i=",i,"idxn=",idxn,"row=",row
    	call MatSetValues(B,1,i,m,idxn,row,INSERT_VALUES,ierr)
        deallocate(idxn,idxtmp,row,rowtmp)
	enddo

	call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
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


subroutine mat_solve(A,b,x,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>
	Mat,			intent(in)	::	A
	Vec,			intent(in)	::	b
	Vec,			intent(out)	::	x
	PetscErrorCode,	intent(out)	::	ierr
    KSP                         ::  ksp
    PC                          ::  pc
    PetscReal                   ::  norm,tol
    PetscInt                    ::  its
    
    call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
    
    call KSPSetOperators(ksp,A,A,ierr)
    call KSPGetPC(ksp,pc,ierr)
    call PCSetType(pc,PCJACOBI,ierr)
    tol = 1.0e-15
    call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
    call KSPSetFromOptions(ksp,ierr)
    call KSPSolve(ksp,b,x,ierr)
    !call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
    !call KSPGetIterationNumber(ksp,its,ierr)
    !print *, ">Iterations number=",its 
end subroutine



end module
