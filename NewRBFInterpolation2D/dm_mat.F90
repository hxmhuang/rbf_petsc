! -----------------------------------------------------------------------
! Distributed Matrix Computing Library
! -----------------------------------------------------------------------
module dm_mat 
   implicit none

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
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_create",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
	! generate matrix A with size m*n
	call MatCreate(PETSC_COMM_WORLD,A,ierr);
	call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,n,ierr)
	call MatSetFromOptions(A,ierr)
	call MatSetUp(A,ierr)
    call PetscLogEventEnd(ievent,ierr)
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
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_zeros",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call mat_create(A,m,n,ierr)
    call MatZeroEntries(A,ierr)
	call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    
    call PetscLogEventEnd(ievent,ierr)
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
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_ones",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    
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
    
    call PetscLogEventEnd(ievent,ierr)
end subroutine



! -----------------------------------------------------------------------
! A=[1 2 3], This function is only used to generate the test data.
!   [4 5 6]
!   [7 8 9]
! -----------------------------------------------------------------------
subroutine mat_seqs(A,m,n,ierr)
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
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_seqs",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
	
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
    call PetscLogEventEnd(ievent,ierr)
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
subroutine mat_eyes(A,m,n,ierr)
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
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_eyes",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
	
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
    call PetscLogEventEnd(ievent,ierr)
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
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_copy",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    
    call MatDuplicate(A,MAT_COPY_VALUES,B,ierr)
    
    call PetscLogEventEnd(ievent,ierr)
end subroutine


! -----------------------------------------------------------------------
! C=A+B
! -----------------------------------------------------------------------
subroutine bk_mat_add(A,B,C,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::  A,B 
	Mat,			intent(out)	::	C
	PetscErrorCode,	intent(out)	::	ierr
	PetscInt					::	nrow1,ncol1,nrow2,ncol2
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_add",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    
	call MatGetSize(A,nrow1,ncol1,ierr)
	call MatGetSize(B,nrow2,ncol2,ierr)
	if(nrow1/=nrow2 .or. ncol1/=ncol2)then
		print *, "Error in mat_add: the matrix A and B should have the same size."
		stop	
	endif

    call MatCreate(PETSC_COMM_WORLD,C,ierr)
    call MatSetSizes(C,PETSC_DECIDE,PETSC_DECIDE,nrow1,ncol1,ierr)
    call MatSetType(C,MATCOMPOSITE,ierr)
	call MatSetUp(C,ierr)
    !call MatCompositeSetType(B,MAT_COMPOSITE_MULTIPLICATIVE,ierr)
    
    call MatCompositeAddMat(C,B,ierr)
    call MatCompositeAddMat(C,A,ierr)

    call MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY,ierr)
    
    call PetscLogEventEnd(ievent,ierr)
end subroutine


! -----------------------------------------------------------------------
! C=A+B
! -----------------------------------------------------------------------
subroutine mat_add(A,B,C,ierr) 
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,	       intent(in)	::  A 
	Mat,	       intent(in)	::  B 
	Mat,           intent(out)  ::	C
	PetscErrorCode,intent(out)  ::	ierr
	PetscInt					::	nrow1,ncol1,nrow2,ncol2
    PetscScalar                 ::  alpha	
    call MatGetSize(A,nrow1,ncol1,ierr)
	call MatGetSize(B,nrow2,ncol2,ierr)
	if(nrow1/=nrow2 .or. ncol1/=ncol2)then
		print *, "Error: matrix A and matrix B should have the same size"
		stop	
	endif
    
    alpha=1.0
    call MatDuplicate(A,MAT_COPY_VALUES,C,ierr)
	call MatAXPY(C,alpha,B,DIFFERENT_NONZERO_PATTERN,ierr)
end subroutine 


! -----------------------------------------------------------------------
! C=A-B
! -----------------------------------------------------------------------
subroutine mat_minus(A,B,C,ierr) 
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,	       intent(in)	::  A 
	Mat,	       intent(in)	::  B 
	Mat,           intent(out)  ::	C
	PetscErrorCode,intent(out)  ::	ierr
	PetscInt					::	nrow1,ncol1,nrow2,ncol2
    PetscScalar                 ::  alpha	
    call MatGetSize(A,nrow1,ncol1,ierr)
	call MatGetSize(B,nrow2,ncol2,ierr)
	if(nrow1/=nrow2 .or. ncol1/=ncol2)then
		print *, "Error: matrix A and matrix B should have the same size"
		stop	
	endif
    
    alpha=-1.0
    call MatDuplicate(A,MAT_COPY_VALUES,C,ierr)
	call MatAXPY(C,alpha,B,DIFFERENT_NONZERO_PATTERN,ierr)
end subroutine 



! -----------------------------------------------------------------------
! C=[A B] 
! -----------------------------------------------------------------------
subroutine mat_hjoin(A,B,C,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,	        intent(in)	::  A,B 
	Mat,            intent(out)	::	C
	PetscErrorCode,	intent(out)	::	ierr
	PetscInt					::	nrow1,ncol1,nrow2,ncol2
	PetscInt					::	col1,col2,m,n
	PetscInt,allocatable		::	idxn1(:),idxn2(:),idxn3(:)
	PetscScalar,allocatable		::	row1(:),row2(:),row3(:)
	PetscInt					::  ista,iend
	integer						::	i
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_hjoin",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
	
    call MatGetSize(A,nrow1,ncol1,ierr)
	call MatGetSize(B,nrow2,ncol2,ierr)
	if(nrow1/=nrow2)then
		print *, "Error: Matrix A and Matrix B should have the same row size"
		stop	
	endif

    call mat_create(C,nrow1,(ncol1+ncol2),ierr)
	call MatGetOwnershipRange(A,ista,iend,ierr)

	do i=ista,iend-1
	    call MatGetRow(A,i,col1,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
        m=col1
	    call MatRestoreRow(A,i,col1,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
	    
        call MatGetRow(B,i,col2,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
		n=col2
	    call MatRestoreRow(B,i,col2,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
	    
	    allocate(idxn1(m),row1(m))
	    allocate(idxn2(n),row2(n))
	    allocate(idxn3(m+n),row3(m+n))
        
        call MatGetRow(A,i,col1,idxn1,row1,ierr)
        idxn3(1:m)=idxn1
        row3(1:m)=row1
		call MatRestoreRow(A,i,col1,idxn1,row1,ierr)
        
        call MatGetRow(B,i,col2,idxn2,row2,ierr)
        idxn3((m+1):(m+n))=ncol1+idxn2
        row3((m+1):(m+n))=row2
		call MatRestoreRow(B,i,col2,idxn2,row2,ierr)
		
		!print *,">i=",i,"idxn3=",idxn3,"row3=",row3
		call MatSetValues(C,1,i,(m+n),idxn3,row3,INSERT_VALUES,ierr)
	    deallocate(idxn1,idxn2,idxn3,row1,row2,row3)
	enddo

	call MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY,ierr)
    
    call PetscLogEventEnd(ievent,ierr)
end subroutine


! -----------------------------------------------------------------------
! C=A*B
! -----------------------------------------------------------------------
subroutine mat_mult(A,B,C,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::  A,B 
	Mat,			intent(out)	::	C	
	PetscErrorCode,	intent(out)	::	ierr
	PetscInt					::	nrow1,ncol1,nrow2,ncol2
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_mult",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
	
	call MatGetSize(A,nrow1,ncol1,ierr)
    call MatGetSize(B,nrow2,ncol2,ierr)
 	if(ncol1/=nrow2)then
 		print *, "Error in mat_mult: the column of A matrix should equal to the row of B matrix."
 		stop	
 	endif

    call MatMatMult(A,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,C,ierr) 
	call MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY,ierr)
    call PetscLogEventEnd(ievent,ierr)
end subroutine


! -----------------------------------------------------------------------
! C=A.*B
! -----------------------------------------------------------------------
subroutine mat_emult(A,B,C,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::  A,B 
	Mat,			intent(out)	::	C
	PetscErrorCode,	intent(out)	::	ierr
	PetscInt					::	nrow1,ncol1,nrow2,ncol2
	PetscInt					::	col1,col2,m,n
	PetscInt,allocatable		::	idxn1(:),idxn2(:),idxn3(:),idxtmp(:)
	PetscScalar,allocatable		::	row1(:),row2(:),row3(:),rowtmp(:)
	PetscInt					::  ista,iend
    PetscInt                    ::  pos1,pos2,counter
	integer						::	i
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_emult",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
	
    call MatGetSize(A,nrow1,ncol1,ierr)
	call MatGetSize(B,nrow2,ncol2,ierr)
	if(nrow1/=nrow2 .or. ncol1/=ncol2)then
		print *, "Error: Matrix A1 and Matrix A2 should have the same sizes"
		stop	
	endif
    
    call mat_create(C,nrow1,ncol1,ierr)
	call MatGetOwnershipRange(A,ista,iend,ierr)
	    
	
    do i=ista,iend-1
	    call MatGetRow(A,i,col1,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
        m=col1
	    call MatRestoreRow(A,i,col1,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
	    
        call MatGetRow(B,i,col2,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
		n=col2
	    call MatRestoreRow(B,i,col2,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr)
	    
        allocate(idxn1(m),row1(m))
        allocate(idxtmp(m),rowtmp(m))
	    allocate(idxn2(n),row2(n))
	    allocate(idxn3(min(m,n)),row3(min(m,n)))
    	
        call MatGetRow(A,i,col1,idxn1,row1,ierr)
        m=col1
        idxtmp=idxn1
        rowtmp=row1
	    call MatRestoreRow(A,i,col1,idxn1,PETSC_NULL_SCALAR,ierr)
 
        call MatGetRow(B,i,col2,idxn2,row2,ierr)
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
        call MatRestoreRow(B,i,col2,idxn2,row2,ierr)
           
		!print *,">i=",i,"idxn3=",idxn3,"row3=",row3
		call MatSetValues(C,1,i,counter,idxn3(1:counter),row3(1:counter),INSERT_VALUES,ierr)
	    deallocate(idxn1,idxn2,idxn3,idxtmp,row1,row2,row3,rowtmp)
	enddo

	call MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY,ierr)
    call PetscLogEventEnd(ievent,ierr)
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
	Mat							::  P,T,W
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_rep",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

	call MatGetSize(A,nrow,ncol,ierr)
	
    call mat_eyes(P,m*nrow,nrow,ierr)
    call mat_mult(P,A,W,ierr) 

	call mat_eyes(T,ncol,n*ncol,ierr)
    call mat_mult(W,T,B,ierr) 

	call mat_destroy(P,ierr)
	call mat_destroy(W,ierr)
	call mat_destroy(T,ierr)
	
    call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
    call PetscLogEventEnd(ievent,ierr)
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
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_sum",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
	
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
    call PetscLogEventEnd(ievent,ierr)
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
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_axpy",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
	
    call MatAXPY(Y,a,X,DIFFERENT_NONZERO_PATTERN,ierr)
    call PetscLogEventEnd(ievent,ierr)
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
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_aypx",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)

    call MatAYPX(Y,a,X,DIFFERENT_NONZERO_PATTERN,ierr)
    call PetscLogEventEnd(ievent,ierr)
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
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_trans",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    
    call MatTranspose(A,MAT_INITIAL_MATRIX,B,ierr)
    
    call PetscLogEventEnd(ievent,ierr)
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
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_xyt",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    !call MatMatTransposeMult(X,Y,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,B,ierr) 	
    call mat_trans(Y,W,ierr)
    call mat_mult(X,W,B,ierr)
    call mat_destroy(W,ierr)
    
    call PetscLogEventEnd(ievent,ierr)
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
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_xty",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    
    call MatTransposeMatMult(X,Y,MAT_INITIAL_MATRIX,PETSC_DEFAULT_REAL,B,ierr)
    
    call PetscLogEventEnd(ievent,ierr)
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
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_scale",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    
    call MatScale(X,a,ierr)
    
    call PetscLogEventEnd(ievent,ierr)
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
#include "mat_type.h"
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
	PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_math",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    
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
                    row=0.0    
            end select
        enddo
        do j=1,m
            !Maybe there is a potenial bug here. When rowtmp < 1E-16, exp(rowtmp) will be NaN. 
            !To avoid this, we have to set row=0.0 at his situation.
            if(isnan(row(j))) row(j)=0.0
    	    !print *,"------rowtmp(",j,")=",rowtmp(j),"row(",j,")=",row(j)
        enddo
        
        call MatRestoreRow(A,i,col,idxtmp,rowtmp,ierr)
    	
        !print *,">i=",i,"idxn=",idxn,"row=",row
    	call MatSetValues(B,1,i,m,idxn,row,INSERT_VALUES,ierr)
        deallocate(idxn,idxtmp,row,rowtmp)
	enddo

	call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
    
    call PetscLogEventEnd(ievent,ierr)
end subroutine


! -----------------------------------------------------------------------
! Solve Ax=b 
! -----------------------------------------------------------------------
subroutine mat_solve(A,b,x,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscpc.h>
	Mat,			intent(in)	::	A
	Mat,			intent(in)	::	b
	Mat,			intent(out)	::	x
	PetscErrorCode,	intent(out)	::	ierr
    
    Vec                         ::  vec_b
    Vec                         ::  vec_x
    KSP                         ::  ksp
    PC                          ::  pc
    PetscReal                   ::  tol
    PetscLogEvent	            ::  ievent
	call PetscLogEventRegister("mat_solve",0, ievent, ierr)
    call PetscLogEventBegin(ievent,ierr)
    !PetscInt                    ::  its
	
    call mat_mat2vec(b,vec_b,ierr)
    call VecDuplicate(vec_b,vec_x,ierr)

    
    call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
    
    call KSPSetOperators(ksp,A,A,ierr)
    call KSPGetPC(ksp,pc,ierr)
    !call PCSetType(pc,PCJACOBI,ierr)
    tol = 1.0e-10
    call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
    call KSPSetFromOptions(ksp,ierr)
    call KSPSolve(ksp,vec_b,vec_x,ierr)
    !call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
    !call KSPGetIterationNumber(ksp,its,ierr)
    !print *, ">Iterations number=",its 
    
    call mat_vec2mat(vec_x,x,ierr) 
    call KSPDestroy(ksp,ierr)
    call VecDestroy(vec_b,ierr)
    call VecDestroy(vec_x,ierr)

    call PetscLogEventEnd(ievent,ierr)
end subroutine



! -----------------------------------------------------------------------
! convert one m*1 matrix into a vector 
! -----------------------------------------------------------------------
subroutine mat_mat2vec(A,v,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,    intent(in)      :: A 
    Vec,    intent(out)     :: v
    PetscErrorCode      	::ierr

	PetscInt	            	::  nrow,ncol	
	PetscInt					::  ista,iend
	PetscInt					::  ni
	PetscInt,allocatable		::  idx(:) 
	PetscScalar,allocatable		::  y(:) 
	integer 					:: 	i
   
    call MatGetSize(A,nrow,ncol,ierr)
    if(ncol /=1) then
        print *, "Error in mat_mat2vec: the column of A should be 1"
        stop 
    endif

    call VecCreate(PETSC_COMM_WORLD,v,ierr)
    call VecSetSizes(v,PETSC_DECIDE,nrow,ierr)
    call VecSetFromOptions(v,ierr)

    call MatGetOwnershipRange(A,ista,iend,ierr)
    ni=iend-ista

    !print *, "1=ista=",ista,"iend=",iend,"ni=",ni,"==idx=",idx
    allocate(idx(ni),y(ni))
    do i=ista,iend-1
        idx(i-ista+1)=i 
    enddo
    !print *, "2=ista=",ista,"iend=",iend,"ni=",ni,"==idx=",idx

    call MatGetValues(A,ni,idx,1,0,y,ierr)
    call VecSetValues(v,ni,idx,y,INSERT_VALUES,ierr)  
    call VecAssemblyBegin(v,ierr)
    call VecAssemblyEnd(v,ierr)

    deallocate(idx,y)

end subroutine


! -----------------------------------------------------------------------
! convert one m*1 vector into a matrix 
! -----------------------------------------------------------------------
subroutine mat_vec2mat(v,A,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Vec,    intent(in)          ::  v
    Mat,    intent(out)         ::  A 
	PetscErrorCode              ::  ierr

	PetscInt	            	::  nrow	
	PetscInt					::  ista,iend
	PetscInt					::  ni
	PetscInt,allocatable		::  ix(:) 
	PetscScalar,allocatable		::  y(:) 
	integer 					:: 	i
     
    call VecGetSize(v,nrow,ierr)
    call mat_create(A,nrow,1,ierr)

    call MatGetOwnershipRange(A,ista,iend,ierr)
    ni=iend-ista

    allocate(ix(ni),y(ni))
	
	do i=ista,iend-1
        ix(i-ista+1)=i 
    enddo
    call VecGetValues(v,ni,ix,y,ierr)  
    call MatSetValues(A,ni,ix,1,0,y,INSERT_VALUES,ierr)
	call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

	deallocate(ix,y)

end subroutine



! -----------------------------------------------------------------------
! Load a standard row-cloumn file into a matrix 
! -----------------------------------------------------------------------
subroutine mat_load(filename,A,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    character(len=*),   intent(in)  ::  filename 
	Mat,			    intent(out)	::	A
	PetscErrorCode,	    intent(out)	::	ierr
    
    !PetscReal,allocatable       :: x(:,:)
    PetscScalar,allocatable         :: x(:,:)
 	PetscScalar,allocatable         :: row(:) 
	PetscInt,allocatable		    :: idxn(:)
 	PetscInt					    :: nrow,ncol
 	PetscInt					    :: ista,iend
    integer                         :: i,j,fid
    PetscLogEvent               ::  ievent
    call PetscLogEventRegister("mat_load",0, ievent, ierr) 
    call PetscLogEventBegin(ievent,ierr) 
    
    fid=1000
    open(fid,FILE=filename)
    call getfilerowcol(fid,nrow,ncol,ierr)
    !print *, "nrow=",nrow, "ncol=",ncol
    
    call mat_create(A,nrow,ncol,ierr)
    call MatGetOwnershipRange(A,ista,iend,ierr)
    allocate(x(ncol,nrow),row(ncol),idxn(ncol))

    do i=1,nrow
       read(fid,*), x(:,i)
    enddo
    close(fid)

    do i=ista,iend-1
		do j=1,ncol
			idxn(j)=j-1
		enddo
		call MatSetValues(A,1,i,ncol,idxn,x(:,i+1),INSERT_VALUES,ierr)
	enddo
	call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    deallocate(x,row)
    call PetscLogEventEnd(ievent,ierr) 
end subroutine


subroutine getfilerowcol(fid,nrow,ncol,ierr)
    implicit none
    integer, intent(in)     :: fid 
    integer, intent(out)    :: nrow,ncol 
    integer, intent(out)    :: ierr
    character(len=1)        :: cdummy
    character(len=1000)     :: string 
    real(kind=8)            :: value
    integer                 :: i,j
    nrow=0
    ncol=0
    rewind(fid)

    do
        read(fid, * , iostat = ierr ) cdummy
        if(ierr /= 0) exit
        nrow = nrow + 1
    end do
    rewind(fid)

    read(fid, '(a)' ) string 
    !print *, "string=",string 

    do i =1,1000 ! 
        read(string,*, iostat=ierr ) (value, j=1,i)
        !print *,"value=",value,"ierr=",ierr
        if ( ierr /= 0 ) then
        ncol = i - 1
        exit
        endif
    enddo
    rewind(fid)
end subroutine 


! -----------------------------------------------------------------------
! A(m,n)=value. Note that the starting point of m and n is 1. 
! -----------------------------------------------------------------------
subroutine mat_setvalue(A,m,n,value,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(out)	::	A
	PetscInt,	    intent(in)	::	m,n
	PetscScalar,    intent(in)	::	value
	PetscErrorCode,	intent(out)	::	ierr
    PetscLogEvent               ::  ievent
    call PetscLogEventRegister("mat_setvalue",0, ievent, ierr) 
    call PetscLogEventBegin(ievent,ierr) 

    call MatSetValue(A,m-1,n-1,value,INSERT_VALUES,ierr)
    
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    call PetscLogEventEnd(ievent,ierr) 
end subroutine


! -----------------------------------------------------------------------
! B=A(rows,cols). Get sub matrix.
! -----------------------------------------------------------------------
subroutine mat_submatrix(A,Rows,Cols,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	::	A,Rows,Cols
	PetscErrorCode,	intent(out)	::	ierr
	Mat,			intent(out)	:: 	B
	IS							:: 	ISRows,ISCols
    PetscLogEvent               ::  ievent
    call PetscLogEventRegister("mat_submatrix",0, ievent, ierr) 
    call PetscLogEventBegin(ievent,ierr) 
	
	call mat_mat2is(Rows,ISRows,ierr)
	call mat_mat2is(Cols,ISCols,ierr)
	!call ISView(ISRows,PETSC_VIEWER_STDOUT_WORLD,ierr)
	!call ISView(ISCols,PETSC_VIEWER_STDOUT_WORLD,ierr)
	call MatGetSubMatrix(A,ISRows,ISCols,MAT_INITIAL_MATRIX,B,ierr)	

	call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
    call PetscLogEventEnd(ievent,ierr) 
end subroutine


! -----------------------------------------------------------------------
! convert one m*1 matrix into a IS object. 
! -----------------------------------------------------------------------
subroutine mat_mat2is(A,is,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,    intent(in)      :: A 
    IS,	    intent(out)     :: is 
    PetscErrorCode      	::ierr

	PetscInt	            	::  nrow,ncol	
	PetscInt					::  ista,iend
	PetscInt					::  ni
	PetscInt,allocatable		::  idx(:) 
	PetscScalar,allocatable		::  y(:) 
	integer 					:: 	i
   
    call MatGetSize(A,nrow,ncol,ierr)
    if(ncol /=1) then
        print *, "Error in mat_mat2is: the column of A should be 1"
        stop
    endif

    call MatGetOwnershipRange(A,ista,iend,ierr)
	ni=iend-ista
 
    allocate(idx(ni),y(ni))
    do i=ista,iend-1
        idx(i-ista+1)=i 
    enddo
   
	call MatGetValues(A,ni,idx,1,0,y,ierr) 
	
	call ISCreateGeneral(PETSC_COMM_WORLD,ni,int(y-1),PETSC_COPY_VALUES,is,ierr)
		
    deallocate(idx,y)
end subroutine


! -----------------------------------------------------------------------
! B=A(:,i). Get a column from A.
! -----------------------------------------------------------------------
subroutine mat_getcol(A,n,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,        intent(in)      ::  A 
    PetscInt,	intent(in)      ::  n 
    Mat,        intent(in)      ::  B 
    PetscErrorCode      	    ::  ierr

	PetscInt	            	::  nrow,ncol	
	PetscInt					::  ista,iend
	PetscInt					::  ni
	PetscInt,allocatable		::  row(:) 
	PetscInt                    ::  col 
	integer 					:: 	i
    IS                          ::  ISRows, ISCols 
    PetscLogEvent               ::  ievent
    call PetscLogEventRegister("mat_getcol",0, ievent, ierr) 
    call PetscLogEventBegin(ievent,ierr) 
  
    call MatGetSize(A,nrow,ncol,ierr)
    if(n-1 > ncol) then
        print *, "Error in mat_getcol: the second paramenter should not greater than the column size of A"
        stop
    endif
    
    call MatGetOwnershipRange(A,ista,iend,ierr)
	ni=iend-ista
 
    allocate(row(ni))
    do i=ista,iend-1
        row(i-ista+1)=i 
    enddo
    col=n-1

	call ISCreateGeneral(PETSC_COMM_WORLD,ni,row,PETSC_COPY_VALUES,ISRows,ierr)
    call ISCreateGeneral(PETSC_COMM_WORLD,1,col,PETSC_COPY_VALUES,ISCols,ierr)
	call MatGetSubMatrix(A,ISRows,ISCols,MAT_INITIAL_MATRIX,B,ierr)	
    
    deallocate(row)
    call PetscLogEventEnd(ievent,ierr) 
end subroutine



! -----------------------------------------------------------------------
! B=A(m,:). Get a row from A.
! -----------------------------------------------------------------------
subroutine mat_getrow(A,m,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    Mat,        intent(in)      :: A 
    PetscInt,	intent(in)      :: m 
    Mat,        intent(in)      :: B 
    Mat                         :: W 
    PetscErrorCode      	    ::ierr

	PetscInt	            	::  nrow,ncol	
    PetscLogEvent               ::  ievent
    call PetscLogEventRegister("mat_getrow",0, ievent, ierr) 
    call PetscLogEventBegin(ievent,ierr) 
   
    call MatGetSize(A,nrow,ncol,ierr)
    if(m-1 > nrow) then
        print *, "Error in mat_getrow: the second paramenter should not greater than the row size of A"
        stop
    endif
    call mat_trans(A,W,ierr)
    call mat_getcol(W,m,B,ierr)
    call mat_destroy(W,ierr)

    call PetscLogEventEnd(ievent,ierr) 
end subroutine





end module
