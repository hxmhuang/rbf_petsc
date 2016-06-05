! -----------------------------------------------------------------------
! Distributed Matrix Wrapper 
! -----------------------------------------------------------------------
module dm 
    use dm_type
    use dm_mat
    implicit none
#include "mat_type.h"

    interface operator (+)
        module procedure dm_add
    end interface

    interface operator (-)
        module procedure dm_minus
    end interface

    interface operator (*)
        module procedure dm_mult1
        module procedure dm_mult2
        module procedure dm_mult3
        module procedure dm_mult4
        module procedure dm_mult5
    end interface
    
    ! element multiple
    interface operator (.em.)
        module procedure dm_emult
    end interface
    
    ! join horizontally
    interface operator (.hj.)
        module procedure dm_hjoin
    end interface

    interface dm_axpy
        module procedure dm_axpy1
        module procedure dm_axpy2
    end interface

    interface dm_aypx
        module procedure dm_aypx1
        module procedure dm_aypx2
    end interface

    interface operator (.tr.)
        module procedure dm_trans
    end interface

	interface dm_setvalue
        module procedure dm_setvalue1
        module procedure dm_setvalue2
    end interface

    interface assignment(=)
        module procedure dm_copy
    end interface


contains

! -----------------------------------------------------------------------
! Initialize the distributed matrix environment 
! -----------------------------------------------------------------------
function dm_init() result(ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscErrorCode  ::  ierr 
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
end function

! -----------------------------------------------------------------------
! Get the rank number of the current process in the commmunicator 
! -----------------------------------------------------------------------
function dm_comm_rank() result(myrank)
	implicit none
#include <petsc/finclude/petscsys.h>
    integer         ::  myrank
    PetscErrorCode  ::  ierr 
	call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr)
end function


! -----------------------------------------------------------------------
! Get the size of processes in the commmunicator 
! -----------------------------------------------------------------------
function dm_comm_size() result(mysize)
	implicit none
#include <petsc/finclude/petscsys.h>
    integer         ::  mysize
    PetscErrorCode  ::  ierr 
	call MPI_Comm_rank(PETSC_COMM_WORLD,mysize,ierr)
end function


! -----------------------------------------------------------------------
! Get the input paramenters 
! -----------------------------------------------------------------------
function dm_get_int(str) result(input)
	implicit none
#include <petsc/finclude/petscsys.h>
    character(len=*)::  str
    PetscInt        ::  input 
    PetscErrorCode  ::  ierr 
    call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,str,input,PETSC_NULL_BOOL,ierr)
end function

! -----------------------------------------------------------------------
! Get the input paramenters 
! -----------------------------------------------------------------------
function dm_get_bool(str) result(input)
	implicit none
#include <petsc/finclude/petscsys.h>
    character(len=*)::  str
    logical         ::  input 
	PetscErrorCode  ::  ierr
    !call PetscOptionsGetBool(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,str,flag,PETSC_NULL_BOOL,ierr)
    input=.false.
	call PetscOptionsGetBool(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,'-debug',input,PETSC_NULL_BOOL,ierr)
    !print *,"str=",str," flag=",input
end function


! -----------------------------------------------------------------------
! Get the input paramenters 
! -----------------------------------------------------------------------
function dm_get_real(str) result(input)
	implicit none
#include <petsc/finclude/petscsys.h>
    character(len=*)::  str
    PetscReal       ::  input 
    PetscErrorCode  ::  ierr 
    call PetscOptionsGetReal(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,str,input,PETSC_NULL_BOOL,ierr)
end function


! -----------------------------------------------------------------------
! Finalize the distributed matrix environment 
! -----------------------------------------------------------------------
function dm_finalize() result(ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscErrorCode  ::  ierr 
    call PetscFinalize(ierr)
end function


function dm_create(m,n) result(A)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscInt,       intent(in)  ::  m,n 
    type(Matrix)              ::  A
    PetscErrorCode              ::  ierr
    ! generate matrix A with size m*n
    call mat_create(A%x,m,n,ierr)
    A%xtype=MAT_XTYPE_IMPLICIT
end function 

! -----------------------------------------------------------------------
!Destroy a matrix to free the memory
! -----------------------------------------------------------------------
function dm_destroy(A) result(ierr)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    type(Matrix),   intent(in)  ::  A
    PetscErrorCode              ::  ierr
    ! destroy matrix A
    call mat_destroy(A%x,ierr)
end function 


! -----------------------------------------------------------------------
! Print a matrix on screen
! -----------------------------------------------------------------------
function dm_view(A) result(ierr) 
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    type(Matrix),  intent(in)   ::  A
    PetscErrorCode              ::  ierr
	call MatAssemblyBegin(A%x,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A%x,MAT_FINAL_ASSEMBLY,ierr)
    call MatView(A%x,PETSC_VIEWER_STDOUT_WORLD, ierr)
end function 

! -----------------------------------------------------------------------
! A=0 
! -----------------------------------------------------------------------
function dm_zeros(m,n) result(A)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscInt,       intent(in)  ::  m,n 
    type(Matrix)              ::  A
    PetscErrorCode              ::  ierr
    
    call mat_zeros(A%x,m,n,ierr)
    A%xtype=MAT_XTYPE_IMPLICIT 
end function



function dm_ones(m,n) result(A)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscInt,       intent(in)  ::  m,n
    type(Matrix)              ::  A
    PetscErrorCode              ::  ierr
    
    call mat_ones(A%x,m,n,ierr)
    A%xtype=MAT_XTYPE_IMPLICIT 
end function


! -----------------------------------------------------------------------
! A=[1 2 3], This function is only used to generate the test data.
!   [4 5 6]
!   [7 8 9]
! -----------------------------------------------------------------------
function dm_seqs(m,n) result(A)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    PetscInt,       intent(in)  ::  m,n
    type(Matrix)              ::  A
    PetscErrorCode              ::  ierr

    call mat_seqs(A%x,m,n,ierr)
    A%xtype=MAT_XTYPE_IMPLICIT 
end function


! -----------------------------------------------------------------------
! The eyes function is used to generate the simple and complex identity matrixs. 
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
function dm_eyes(m,n) result(A)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	PetscInt,       intent(in)	::	m,n	
	type(Matrix)			    ::	A
	PetscErrorCode	            ::	ierr
    
    call mat_eyes(A%x,m,n,ierr)
    A%xtype=MAT_XTYPE_IMPLICIT 
end function 


! -----------------------------------------------------------------------
! B=A. This function uses the implicit matrix A directly because A is not need to free. 
! -----------------------------------------------------------------------
subroutine dm_copy(B,A)
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    type(Matrix),  intent(in)    ::  A
    type(Matrix),  intent(out)   ::  B
    type(Matrix)                 ::  W
    PetscErrorCode               ::  ierr
    !Free the space of B matrix 
    !call mat_destroy(B%x,ierr)
    if(B%xtype==MAT_XTYPE_EXPLICIT) then
        W%x=B%x
    endif
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        B%x=A%x
    else
        call mat_copy(A%x,B%x,ierr)
    endif

    if(B%xtype==MAT_XTYPE_EXPLICIT) then
        call mat_destroy(W%x,ierr)
    endif
 
    B%xtype=MAT_XTYPE_EXPLICIT
end subroutine



! -----------------------------------------------------------------------
! C=A+B
! -----------------------------------------------------------------------
function dm_add(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)                ::	C
	PetscErrorCode      		::	ierr

    call mat_add(A%x,B%x,C%x,ierr)
    C%xtype=MAT_XTYPE_IMPLICIT 
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! C=A-B
! -----------------------------------------------------------------------
function dm_minus(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              ::	C
	PetscErrorCode      		::	ierr
    call mat_minus(A%x,B%x,C%x,ierr)
    C%xtype=MAT_XTYPE_IMPLICIT 
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! C=[A B] 
! -----------------------------------------------------------------------
function dm_hjoin(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              ::	C
	PetscErrorCode      		::	ierr
    call mat_hjoin(A%x,B%x,C%x,ierr)
    C%xtype=MAT_XTYPE_IMPLICIT 

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! C=A*B
! -----------------------------------------------------------------------
function dm_mult1(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              ::	C
	PetscErrorCode      		::	ierr
    call mat_mult(A%x,B%x,C%x,ierr)
    C%xtype=MAT_XTYPE_IMPLICIT 

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 


function dm_mult2(alpha,A) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	PetscScalar,	intent(in)	::  alpha 
	type(Matrix),	intent(in)	::  A 
	type(Matrix)              ::	B
	PetscErrorCode      		::	ierr
    call mat_copy(A%x,B%x,ierr) 
    call mat_scale(B%x,alpha,ierr)
    B%xtype=MAT_XTYPE_IMPLICIT 

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


function dm_mult3(A,alpha) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	PetscScalar,	intent(in)	::  alpha 
	type(Matrix),	intent(in)	::  A 
	type(Matrix)              ::	B
	PetscErrorCode      		::	ierr
    call mat_copy(A%x,B%x,ierr) 
    call mat_scale(B%x,alpha,ierr)
    B%xtype=MAT_XTYPE_IMPLICIT 

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


function dm_mult4(alpha,A) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	real,           intent(in)	::  alpha 
	type(Matrix),	intent(in)	::  A 
	type(Matrix)                ::	B
	PetscErrorCode      		::	ierr
    call mat_copy(A%x,B%x,ierr) 
    call mat_scale(B%x,real(alpha,8),ierr)
    B%xtype=MAT_XTYPE_IMPLICIT 

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


function dm_mult5(A,alpha) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	real,           intent(in)	::  alpha 
	type(Matrix),	intent(in)	::  A 
	type(Matrix)                ::	B
	PetscErrorCode      		::	ierr
    call mat_copy(A%x,B%x,ierr) 
    call mat_scale(B%x,real(alpha,8),ierr)
    B%xtype=MAT_XTYPE_IMPLICIT 

    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! B=A1.*A2
! -----------------------------------------------------------------------
function dm_emult(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              ::	C
	PetscErrorCode      		::	ierr
    call mat_emult(A%x,B%x,C%x,ierr)
    C%xtype=MAT_XTYPE_IMPLICIT 
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! B=repmat(A,m,n)
! -----------------------------------------------------------------------
function dm_rep(A,m,n) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  m,n 
	type(Matrix)              ::	B
	PetscErrorCode      		::	ierr
    call mat_rep(A%x,m,n,B%x,ierr)
    B%xtype=MAT_XTYPE_IMPLICIT 
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! Sum of elements along with the row or column.
! Suppose A=[1,2,3]
!           [4,5,6],
! then mat_sum(A,1,B) will make B=[5,7,9],
!      mat_sum(A,2,B) will make B=[6 ]
!                                 [15]
! -----------------------------------------------------------------------
function dm_sum(A,ndim) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	integer,	    intent(in)	::  ndim 
	type(Matrix)              ::	B
	PetscErrorCode      		::	ierr
    call mat_sum(A%x,ndim,B%x,ierr)
    B%xtype=MAT_XTYPE_IMPLICIT 
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! Compute Y = a*X + Y.
! -----------------------------------------------------------------------
function dm_axpy1(Y,a,X) result(ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  X 
	PetscScalar,    intent(in)	::	a
	type(Matrix), intent(inout) ::  Y 
	PetscErrorCode      		::	ierr
    call mat_axpy(Y%x,a,X%x,ierr)
    Y%xtype=MAT_XTYPE_IMPLICIT 
    
    if (X%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(X%x,ierr)
    endif
end function 


function dm_axpy2(Y,a,X) result(ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  X 
	real,           intent(in)	::	a
	type(Matrix), intent(inout) ::  Y 
	PetscErrorCode      		::	ierr
    call mat_axpy(Y%x,real(a,kind=8),X%x,ierr)
    Y%xtype=MAT_XTYPE_IMPLICIT 
    
    if (X%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(X%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! Compute Y = a*Y + X.
! -----------------------------------------------------------------------
function dm_aypx1(Y,a,X) result(ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  X 
	PetscScalar,    intent(in)	::	a
	type(Matrix), intent(inout) ::  Y 
	PetscErrorCode      		::	ierr
    call mat_aypx(Y%x,a,X%x,ierr)
    Y%xtype=MAT_XTYPE_IMPLICIT 
    
    if (X%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(X%x,ierr)
    endif
end function 


function dm_aypx2(Y,a,X) result(ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  X 
	real        ,   intent(in)	::	a
	type(Matrix), intent(inout) ::  Y 
	PetscErrorCode      		::	ierr
    call mat_aypx(Y%x,real(a,kind=8),X%x,ierr)
end function 



! -----------------------------------------------------------------------
! B = A^T.
! -----------------------------------------------------------------------
function dm_trans(A) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(Matrix)              ::	B
	PetscErrorCode      		::	ierr
    call mat_trans(A%x,B%x,ierr)
    B%xtype=MAT_XTYPE_IMPLICIT 
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! C = A*B^T
! -----------------------------------------------------------------------
function dm_xyt(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              ::	C
	PetscErrorCode      		::	ierr
    call mat_xyt(A%x,B%x,C%x,ierr)
    C%xtype=MAT_XTYPE_IMPLICIT 
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! C = A^T*B
! -----------------------------------------------------------------------
function dm_xty(A,B) result(C)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)	::  B 
	type(Matrix)              ::	C
	PetscErrorCode      		::	ierr
    call mat_xty(A%x,B%x,C%x,ierr)
    C%xtype=MAT_XTYPE_IMPLICIT 
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (B%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(B%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! B=exp(A) 
! -----------------------------------------------------------------------
function dm_exp(A) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include "mat_type.h"
	type(Matrix),	intent(in)	::  A 
	type(Matrix)              ::	B
	PetscErrorCode      		::	ierr
    call mat_math(A%x,MAT_MATH_EXP,B%x,ierr)
    B%xtype=MAT_XTYPE_IMPLICIT 
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! B=log(A) 
! -----------------------------------------------------------------------
function dm_log(A) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include "mat_type.h"
	type(Matrix),	intent(in)	::  A 
	type(Matrix)              ::	B
	PetscErrorCode      		::	ierr
    call mat_math(A%x,MAT_MATH_LOG,B%x,ierr)
    B%xtype=MAT_XTYPE_IMPLICIT 
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! B=sqrt(A) 
! -----------------------------------------------------------------------
function dm_sqrt(A) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include "mat_type.h"
	type(Matrix),	intent(in)	::  A 
	type(Matrix)              ::	B
	PetscErrorCode      		::	ierr
    call mat_math(A%x,MAT_MATH_SQRT,B%x,ierr)
    B%xtype=MAT_XTYPE_IMPLICIT 
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! Solve Ax=b 
! -----------------------------------------------------------------------
function dm_solve(A,b) result(x)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)	::  A 
	type(Matrix),	intent(in)  ::	b
	type(Matrix)            	::	x
	PetscErrorCode      		::	ierr
    
    call mat_solve(A%x,b%x,x%x,ierr)
    x%xtype=MAT_XTYPE_IMPLICIT 
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (b%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(b%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! Load a standard row-cloumn file into a matrix 
! -----------------------------------------------------------------------
function dm_load(filename) result(A)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
    character(len=*),   intent(in)  ::  filename 
	type(Matrix)	            	::  A 
	PetscErrorCode      		    ::	ierr
    
    call mat_load(filename,A%x,ierr)
    A%xtype=MAT_XTYPE_IMPLICIT 
end function 


! -----------------------------------------------------------------------
! A(m,n)=value 
! -----------------------------------------------------------------------
function dm_setvalue1(A,m,n,value) result(ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>    
	type(Matrix)	            ::  A 
	PetscInt,	    intent(in)	::	m,n
	PetscScalar,    intent(in)	::	value
	PetscErrorCode	        	::	ierr
	
    call mat_setvalue(A%x,m,n,value,ierr)
end function 

function dm_setvalue2(A,m,n,value) result(ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>    
	type(Matrix)	            ::  A 
	PetscInt,	    intent(in)	::	m,n
	PetscInt,    	intent(in)	::	value
	PetscErrorCode	        	::	ierr
	
    call mat_setvalue(A%x,m,n,real(value,kind=8),ierr)
end function 


! -----------------------------------------------------------------------
! B=A(rows,cols). Get sub matrix.
! -----------------------------------------------------------------------
function dm_submatrix(A,Rows,Cols) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>    
	type(Matrix),	intent(in)	::	A
	type(Matrix),	intent(in)	::	Rows
	type(Matrix),	intent(in)	::	Cols
	type(Matrix)				::	B
	PetscErrorCode				::	ierr
    
	call mat_submatrix(A%x,Rows%x,Cols%x,B%x,ierr)
    B%xtype=MAT_XTYPE_IMPLICIT 
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
    if (Rows%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(Rows%x,ierr)
    endif
    if (Cols%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(Cols%x,ierr)
    endif
end function 


! -----------------------------------------------------------------------
! B=A(:,i). Get a column from A.
! -----------------------------------------------------------------------
function dm_getcol(A,n) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>    
	type(Matrix),	intent(in)	::	A
    integer,        intent(in)  ::  n
	type(Matrix)				::	B
	PetscErrorCode				::	ierr
    
	call mat_getcol(A%x,n,B%x,ierr)
    B%xtype=MAT_XTYPE_IMPLICIT 
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 
 

! -----------------------------------------------------------------------
! B=A(m,:). Get a row from A.
! -----------------------------------------------------------------------
function dm_getrow(A,n) result(B)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>    
	type(Matrix),	intent(in)	::	A
    integer,        intent(in)  ::  n
	type(Matrix)				::	B
	PetscErrorCode				::	ierr
    
	call mat_getrow(A%x,n,B%x,ierr)
    B%xtype=MAT_XTYPE_IMPLICIT 
    
    if (A%xtype==MAT_XTYPE_IMPLICIT) then
        call mat_destroy(A%x,ierr)
    endif
end function 


end module 


