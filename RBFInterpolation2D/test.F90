! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Implementing the rbf_interpolation2D with the ASM preconditioner
! -----------------------------------------------------------------------
program test 
	use matrix
	use vector 
	use rbf 
	implicit none

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscviewer.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#include "mat_math_type.h"
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	!                   Variable declarations
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	!KSP				ksp
	!PC				pc
	Mat				A,A1,A2,B
	Vec             rhs,x
    !PetscReal		error
	PetscMPIInt		myrank,mysize
	PetscErrorCode	ierr
	PetscInt		m,n
	PetscBool		debug
    PetscScalar     alpha	
	PetscLogEvent	ievent(20)
 	m=2
 	n=3
    debug = .false.
    alpha=1.0
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	!                 Beginning of program
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

	call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr)
	call MPI_Comm_rank(PETSC_COMM_WORLD,mysize,ierr)
    
    call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,'-m',m,PETSC_NULL_BOOL,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,'-n',n,PETSC_NULL_BOOL,ierr)
    call PetscOptionsGetBool(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,'-debug',debug,PETSC_NULL_BOOL,ierr)
    if(myrank==0) then 
        print *, "============Input paramenters============"
        print *, "m=",m,",n=",n,",debug=",debug
     endif 
	
	call PetscLogEventRegister("mat_zeros",0, ievent(1), ierr)
	call PetscLogEventRegister("mat_ones",0, ievent(2), ierr)
	call PetscLogEventRegister("mat_seq",0, ievent(3), ierr)
	call PetscLogEventRegister("mat_eye",0, ievent(4), ierr)
	call PetscLogEventRegister("mat_copy",0, ievent(5), ierr)
	call PetscLogEventRegister("mat_add",0, ievent(6), ierr)
	call PetscLogEventRegister("mat_hjoin",0, ievent(7), ierr)
	call PetscLogEventRegister("mat_mult",0, ievent(8), ierr)
	call PetscLogEventRegister("mat_eprod",0, ievent(9), ierr)
	call PetscLogEventRegister("mat_rep",0, ievent(10), ierr)
	call PetscLogEventRegister("mat_sum",0, ievent(11), ierr)
	call PetscLogEventRegister("mat_axpy",0, ievent(12), ierr)
	call PetscLogEventRegister("mat_aypx",0, ievent(13), ierr)
	call PetscLogEventRegister("mat_trans",0, ievent(14), ierr)
	call PetscLogEventRegister("mat_xyt",0, ievent(15), ierr)
	call PetscLogEventRegister("mat_xty",0, ievent(16), ierr)
	call PetscLogEventRegister("mat_scale",0, ievent(17), ierr)
	call PetscLogEventRegister("mat_solve",0, ievent(18), ierr)
	call PetscLogEventRegister("mat_math",0, ievent(19), ierr)
	call PetscLogEventRegister("mat_20",0, ievent(20), ierr)


 	if(myrank==0) print *, "==============Test mat_zeros==============="
	call PetscLogEventBegin(ievent(1),ierr)
    call mat_zeros(A,m,n,ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
 	endif
 	call mat_destroy(A,ierr)	
	call PetscLogEventEnd(ievent(1),ierr)

 	if(myrank==0) print *, "==============Test mat_ones==============="
	call PetscLogEventBegin(ievent(2),ierr)
 	call mat_ones(A,m,n,ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
 	endif
 	call mat_destroy(A,ierr)	
	call PetscLogEventEnd(ievent(2),ierr)

 	if(myrank==0) print *, "==============Test mat_seq==============="
	call PetscLogEventBegin(ievent(3),ierr)
 	call mat_seq(A,m,n,ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
 	endif
 	call mat_destroy(A,ierr)	
	call PetscLogEventEnd(ievent(3),ierr)
  
    if(myrank==0) print *, "==============Test mat_eye==============="
	call PetscLogEventBegin(ievent(4),ierr)
    call mat_eye(A,m,m,ierr)	
    if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
 	endif
    call mat_destroy(A,ierr)	

 	call mat_eye(A,m,2*m,ierr)	
    if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
 	endif
 	call mat_destroy(A,ierr)	

   	call mat_eye(A,2*m,m,ierr)	
    if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
 	endif
   	call mat_destroy(A,ierr)	
	call PetscLogEventEnd(ievent(4),ierr)


 	if(myrank==0) print *, "==============Test mat_copy==============="
	call PetscLogEventBegin(ievent(5),ierr)
 	call mat_eye(A,m,m,ierr)
 	call mat_copy(A,B,ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
        if(myrank==0) print *, ">B="
        call mat_view(B,ierr)
 	endif
 	call mat_destroy(A,ierr)	
 	call mat_destroy(B,ierr)
	call PetscLogEventEnd(ievent(5),ierr)

    if(myrank==0) print *, "==============Test mat_add==============="
	call PetscLogEventBegin(ievent(6),ierr)
    call mat_eye(A1,m,m,ierr)
    call mat_ones(A2,m,m,ierr)
    call mat_add(A1,A2,B,ierr)
    if(debug) then
        if(myrank==0) print *, ">A1="
        call mat_view(A1,ierr)
        if(myrank==0) print *, ">A2="
        call mat_view(A2,ierr)
        if(myrank==0) print *, ">B="
        call mat_view(B,ierr)
 	endif
    call mat_destroy(A1,ierr)	
    call mat_destroy(A2,ierr)	
    call mat_destroy(B,ierr)	
	call PetscLogEventEnd(ievent(6),ierr)

    if(myrank==0) print *, "==============Test mat_hjoin==============="
    call PetscLogEventBegin(ievent(7),ierr)
    call mat_zeros(A1,m,n,ierr)
    call mat_seq(A2,m,m,ierr)
    call mat_hjoin(A1,A2,B,ierr)
    if(debug) then
        if(myrank==0) print *, ">A1="
        call mat_view(A1,ierr)
        if(myrank==0) print *, ">A2="
        call mat_view(A2,ierr)
        if(myrank==0) print *, ">B="
        call mat_view(B,ierr)
 	endif
    call mat_destroy(A1,ierr)	
    call mat_destroy(A2,ierr)	
    call mat_destroy(B,ierr)	
    call PetscLogEventEnd(ievent(7),ierr)


    if(myrank==0) print *, "==============Test mat_mult==============="
    call PetscLogEventBegin(ievent(8),ierr)
    call mat_ones(A1,m,m,ierr)
    call mat_eye(A2,m,2*m,ierr)
    call mat_mult(A1,A2,B,ierr)
    if(debug) then
        if(myrank==0) print *, ">A1="
        call mat_view(A1,ierr)
        if(myrank==0) print *, ">A2="
        call mat_view(A2,ierr)
        if(myrank==0) print *, ">B="
        call mat_view(B,ierr)
 	endif
    call mat_destroy(A1,ierr)	
    call mat_destroy(A2,ierr)	
    call mat_destroy(B,ierr)	
	call PetscLogEventEnd(ievent(8),ierr)


 	if(myrank==0) print *, "==============Test mat_eprod==============="
    call PetscLogEventBegin(ievent(9),ierr)
 	call mat_seq(A1,m,m,ierr)
 	call mat_eye(A2,m,m,ierr)
 	call mat_eprod(A1,A2,B,ierr)
    if(debug) then
        if(myrank==0) print *, ">A1="
        call mat_view(A1,ierr)
        if(myrank==0) print *, ">A2="
        call mat_view(A2,ierr)
        if(myrank==0) print *, ">B="
        call mat_view(B,ierr)
 	endif
 	call mat_destroy(A1,ierr)
 	call mat_destroy(A2,ierr)
 	call mat_destroy(B,ierr)	

    call mat_seq(A,m,n,ierr)
 	call mat_eprod(A,A,B,ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
        if(myrank==0) print *, ">B="
        call mat_view(B,ierr)
 	endif
 	call mat_destroy(A,ierr)
 	call mat_destroy(B,ierr)	
 

call PetscLogEventEnd(ievent(9),ierr)


    if(myrank==0) print *, "==============Test mat_rep==============="
    call PetscLogEventBegin(ievent(10),ierr)
    call mat_seq(A,m,m,ierr)
    call mat_rep(A,3,2,B,ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
        if(myrank==0) print *, ">B="
        call mat_view(B,ierr)
 	endif
    call mat_destroy(A,ierr)	
    call mat_destroy(B,ierr)	
	call PetscLogEventEnd(ievent(10),ierr)


    if(myrank==0) print *, "==============Test mat_sum==============="
    call PetscLogEventBegin(ievent(11),ierr)
    !call mat_seq(A,m,m,ierr)
    call mat_seq(A,m,m,ierr)
    call mat_sum(A,1,B,ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
        if(myrank==0) print *, ">B="
        call mat_view(B,ierr)
 	endif
 	call mat_destroy(A,ierr)	
 	call mat_destroy(B,ierr)	
 
    !call mat_seq(A,m,m,ierr)
    call mat_seq(A,m,m,ierr)
    call mat_sum(A,2,B,ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
        if(myrank==0) print *, ">B="
        call mat_view(B,ierr)
 	endif
 	call mat_destroy(A,ierr)	
 	call mat_destroy(B,ierr)	
    call PetscLogEventEnd(ievent(11),ierr)


    if(myrank==0) print *, "==============Test mat_axpy=============="
    call PetscLogEventBegin(ievent(12),ierr)
 	call mat_seq(A,m,n,ierr)
    call mat_ones(B,m,n,ierr)
    alpha=1.0    
    call mat_axpy(B,alpha,A,ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
        if(myrank==0) print *, ">B="
        call mat_view(B,ierr)
 	endif
 	call mat_destroy(A,ierr)	
 	call mat_destroy(B,ierr)	
	call PetscLogEventEnd(ievent(12),ierr)


    if(myrank==0) print *, "==============Test mat_aypx=============="
    call PetscLogEventBegin(ievent(13),ierr)
 	call mat_seq(A,m,n,ierr)
 	call mat_ones(B,m,n,ierr)
    alpha=10.0    
    call mat_aypx(B,alpha,A,ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
        if(myrank==0) print *, ">B="
        call mat_view(B,ierr)
 	endif
 	call mat_destroy(A,ierr)	
 	call mat_destroy(B,ierr)	
	call PetscLogEventEnd(ievent(13),ierr)

    
    if(myrank==0) print *, "==============Test mat_trans=============="
	call PetscLogEventBegin(ievent(14),ierr)
    call mat_seq(A,m,2*m,ierr)
    call mat_trans(A,B,ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
        if(myrank==0) print *, ">B="
        call mat_view(B,ierr)
 	endif
    call mat_destroy(A,ierr)	
 	call mat_destroy(B,ierr)	
	call PetscLogEventEnd(ievent(14),ierr)


    if(myrank==0) print *, "==============Test mat_xyt=============="
    call PetscLogEventBegin(ievent(15),ierr)
    call mat_seq(A1,m,m,ierr)
    call mat_ones(A2,m,m,ierr)
    call mat_xyt(A1,A2,B,ierr)
    if(debug) then
        if(myrank==0) print *, ">A1="
        call mat_view(A1,ierr)
        if(myrank==0) print *, ">A2="
        call mat_view(A2,ierr)
        if(myrank==0) print *, ">B=A1*A2^T="
        call mat_view(B,ierr)
 	endif
 	call mat_destroy(A1,ierr)	
 	call mat_destroy(A2,ierr)	
 	call mat_destroy(B,ierr)	
	call PetscLogEventEnd(ievent(15),ierr)


    if(myrank==0) print *, "==============Test mat_xty=============="
    call PetscLogEventBegin(ievent(16),ierr)
 	call mat_seq(A1,m,m,ierr)
 	call mat_ones(A2,m,m,ierr)
    call mat_xty(A1,A2,B,ierr)
    if(debug) then
        if(myrank==0) print *, ">A1="
        call mat_view(A1,ierr)
        if(myrank==0) print *, ">A2="
        call mat_view(A2,ierr)
        if(myrank==0) print *, ">B=A1^T*A2="
        call mat_view(B,ierr)
 	endif
 	call mat_destroy(A1,ierr)	
 	call mat_destroy(A2,ierr)	
 	call mat_destroy(B,ierr)	
	call PetscLogEventEnd(ievent(16),ierr)

    if(myrank==0) print *, "==============Test mat_scale=============="
    call PetscLogEventBegin(ievent(17),ierr)
 	call mat_eye(A,m,m,ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
 	endif
    alpha=2.0
    call mat_scale(A,alpha,ierr)
    if(debug) then
        if(myrank==0) print *, ">ScaleA="
        call mat_view(A,ierr)
 	endif
 	call mat_destroy(A,ierr)	
	call PetscLogEventEnd(ievent(17),ierr)

    if(myrank==0) print *, "==============Test mat_solve=============="
    call PetscLogEventBegin(ievent(18),ierr)
 	call mat_seq(A,m,m,ierr)
    call vec_create(rhs,m,ierr)
    call vec_create(x,m,ierr)
    alpha=1.0
    call VecSet(rhs,alpha,ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
        if(myrank==0) print *, ">rhs="
        call vec_view(rhs,ierr)
 	endif

    call mat_solve(A,rhs,x,ierr)
    if(debug) then
        if(myrank==0) print *, ">x="
        call vec_view(x,ierr)
 	endif
 	call mat_destroy(A,ierr)	
 	call vec_destroy(rhs,ierr)	
 	call vec_destroy(x,ierr)	
	call PetscLogEventEnd(ievent(18),ierr)

	if(myrank==0) print *, "==============Test mat_math=============="
    call PetscLogEventBegin(ievent(19),ierr)
	!m=10
 	call mat_seq(A,m,m,ierr)
    call mat_math(A,MAT_MATH_SQRT,B,ierr)
	if(debug) then
        if(myrank==0) print *, ">A="
        call mat_view(A,ierr)
        if(myrank==0) print *, ">B="
        call mat_view(B,ierr)
 	endif
 	call mat_destroy(A,ierr)	
 	call mat_destroy(B,ierr)	
	call PetscLogEventEnd(ievent(19),ierr)


    if(myrank==0) print *, "===================End==================="
	call PetscFinalize(ierr)
end program
