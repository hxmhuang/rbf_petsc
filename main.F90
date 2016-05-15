! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Implementing the rbf_interpolation2D with the ASM preconditioner
! -----------------------------------------------------------------------
program main
    use matrix
    use vector 
    use rbf 
    implicit none

#include "mat_math_type.h"
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
    Mat				dsites,ctrs,epoints
    !Mat				W1,W2,W3
    Mat				DM_data,DM_eval,IM,EM
	Vec				u,x,rhs,exact,s,norm
	!PetscReal		error
	PetscMPIInt		myrank,mysize
	PetscErrorCode	ierr
	PetscReal		ep
	PetscInt		meval,neval,m,n
	PetscBool		debug
    PetscScalar     alpha,rmse	
	PetscLogEvent	ievent(20)
    
    debug = .false.
    alpha=1.0
	ep=6.1
	m=3
	n=3
	meval=4
	neval=4
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	!                 Beginning of program
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

	call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr)
	call MPI_Comm_rank(PETSC_COMM_WORLD,mysize,ierr)
 
    call PetscOptionsGetReal(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,'-ep',ep,PETSC_NULL_BOOL,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,'-m',m,PETSC_NULL_BOOL,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,'-n',n,PETSC_NULL_BOOL,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,'-meval',meval,PETSC_NULL_BOOL,ierr)
    call PetscOptionsGetInt(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,'-neval',neval,PETSC_NULL_BOOL,ierr)
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
	call vec_duplicate(x,u,ierr)
	call vec_duplicate(x,rhs,ierr)
	
    call vec_create(exact,meval*neval,ierr)
	call vec_duplicate(exact,s,ierr)
	call vec_duplicate(exact,norm,ierr)

	! generate matrix A with size M*N
	!call mat_create(A,m*n,m*n,ierr)
	
    call mat_create(dsites,m*n,2,ierr)
    call mat_create(epoints,meval*neval,2,ierr)
    
    call rbf_createpoints(dsites,m,n,ierr)
    
    call rbf_testfunctionD(dsites,rhs,ierr)
    
    call mat_copy(dsites,ctrs,ierr)
    
    call rbf_distancematrix(dsites,ctrs,DM_data,ierr)
    
    call rbf_guassian(ep,DM_data,IM,ierr)
    
    call mat_solve(IM,rhs,x,ierr) 
    
    call rbf_createpoints(epoints,meval,neval,ierr)
    
    call rbf_testfunctionD(epoints,exact,ierr)
    
    call rbf_distancematrix(epoints,ctrs,DM_eval,ierr)
    
    call rbf_guassian(ep,DM_eval,EM,ierr)
    
    call MatMult(EM,x,s,ierr)
   
    call VecCopy(s,norm,ierr)
    alpha=-1.0
    call VecAXPY(norm,alpha,exact,ierr)
    call VecNorm(norm,NORM_2,rmse,ierr) 
        
    rmse=rmse/neval
    
    !if(myrank==0) print *, ">DM_eval="
    !call mat_view(DM_eval,ierr)
    !call vec_view(exact,ierr)

    if(debug) then
        if(myrank==0) print *, ">dsites="
        call mat_view(dsites,ierr)
        if(myrank==0) print *, ">rhs="
        call vec_view(rhs,ierr)
        if(myrank==0) print *, ">DM_data="
        call mat_view(DM_data,ierr)
        if(myrank==0) print *, ">IM="
        call mat_view(IM,ierr)
        
        if(myrank==0) print *, ">epoints="
        call mat_view(epoints,ierr)
        if(myrank==0) print *, ">exact="
        call vec_view(exact,ierr)
        if(myrank==0) print *, ">DM_eval="
        call mat_view(DM_eval,ierr)
        if(myrank==0) print *, ">EM="
        call mat_view(EM,ierr)
        
        if(myrank==0) print *, ">s="
        call vec_view(s,ierr)
        if(myrank==0) print *, ">norm="
        call vec_view(norm,ierr)
    endif
    
    if(myrank==0) then 
        print *, "================Output================="
        if(myrank==0) print *, ">rmse=",rmse
     endif 
    
    call vec_destroy(x,ierr)
    call vec_destroy(u,ierr)
    call vec_destroy(rhs,ierr)
    call vec_destroy(exact,ierr)
    call vec_destroy(s,ierr)
    call vec_destroy(norm,ierr)
    call mat_destroy(DM_data,ierr)
    call mat_destroy(DM_eval,ierr)
    call mat_destroy(IM,ierr)
    call mat_destroy(EM,ierr)
    call mat_destroy(dsites,ierr)
    call mat_destroy(ctrs,ierr)
    call mat_destroy(epoints,ierr)
    call PetscFinalize(ierr)

end program
