! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Implementing the rbf_interpolation2D with the ASM preconditioner
! -----------------------------------------------------------------------
program main
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
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !                   Variable declarations
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !KSP				ksp
    !PC				pc
    Mat				A,dsites,ctrs
	Vec				u,x,b,rhs
	!PetscReal		error
	PetscMPIInt		myrank,mysize
	PetscErrorCode	ierr
	PetscReal		ep
	PetscInt		meval,neval,m,n
	PetscBool		debug
    PetscScalar     alpha	
	PetscLogEvent	ievent(20)
    
    debug = .false.
    alpha=1.0
	ep=4.1
	m=3
	n=3
	meval=30
	neval=30
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
		
	call PetscLogEventRegister("rbf_zeros",0, ievent(1), ierr)
	call PetscLogEventRegister("rbf_ones",0, ievent(2), ierr)
	call PetscLogEventRegister("rbf_seq",0, ievent(3), ierr)
	call PetscLogEventRegister("rbf_eye",0, ievent(4), ierr)
	call PetscLogEventRegister("rbf_copy",0, ievent(5), ierr)
	call PetscLogEventRegister("rbf_add",0, ievent(6), ierr)
	call PetscLogEventRegister("rbf_hjoin",0, ievent(7), ierr)
	call PetscLogEventRegister("rbf_mult",0, ievent(8), ierr)
	call PetscLogEventRegister("rbf_eprod",0, ievent(9), ierr)
	call PetscLogEventRegister("rbf_rep",0, ievent(10), ierr)
	call PetscLogEventRegister("rbf_sum",0, ievent(11), ierr)
	call PetscLogEventRegister("rbf_axpy",0, ievent(12), ierr)
	call PetscLogEventRegister("rbf_aypx",0, ievent(13), ierr)
	call PetscLogEventRegister("rbf_trans",0, ievent(14), ierr)
	call PetscLogEventRegister("rbf_xyt",0, ievent(15), ierr)
	call PetscLogEventRegister("rbf_xty",0, ievent(16), ierr)
	call PetscLogEventRegister("rbf_17",0, ievent(17), ierr)
	call PetscLogEventRegister("rbf_18",0, ievent(18), ierr)
	call PetscLogEventRegister("rbf_19",0, ievent(19), ierr)
	call PetscLogEventRegister("rbf_20",0, ievent(20), ierr)

	! generate some vectors: x,b,u
	call vec_create(x,m*n,ierr)
	call vec_duplicate(x,b,ierr)
	call vec_duplicate(x,u,ierr)
	call vec_duplicate(x,rhs,ierr)

	! generate matrix A with size M*N
	call mat_create(A,m*n,m*n,ierr)
	
	print *, "==============createpoints & testfunctionD==============="
	call mat_create(dsites,m*n,2,ierr)
	call rbf_createpoints(dsites,m,n,ierr)
	call rbf_testfunctionD(dsites,rhs,ierr)
    if(debug) then
        if(myrank==0) print *, ">dsites="
        call mat_view(dsites,ierr)
        if(myrank==0) print *, ">rhs="
	    call vec_view(rhs,ierr)
 	endif
	
    call mat_copy(dsites,ctrs,ierr)

	call vec_destroy(x,ierr)
	call vec_destroy(b,ierr)
	call vec_destroy(u,ierr)
	call vec_destroy(rhs,ierr)
	call mat_destroy(a,ierr)
	call mat_destroy(dsites,ierr)
	call mat_destroy(ctrs,ierr)
	call PetscFinalize(ierr)

end program
