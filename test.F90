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
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	!                   Variable declarations
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	!KSP				ksp
	!PC				pc
	Mat				A,A1,A2,B
	!PetscReal		error
	PetscMPIInt		myrank,mysize
	PetscErrorCode	ierr
	PetscInt		m,n
	PetscBool		debug
	
!	debug = .false.
	debug = .true.
	
	m=3
	n=2
! 	m=900
!  	n=900
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	!                 Beginning of program
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

	call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr)
	call MPI_Comm_rank(PETSC_COMM_WORLD,mysize,ierr)

 	if(myrank==0) print *, "==============Test mat_eprod==============="
 	call mat_create(A,m,n,ierr)
 	call mat_create(B,m,n,ierr)
 	call mat_seq(A,ierr)
 	call mat_eprod(A,A,B,ierr)
 	if(debug) call mat_view(B,ierr)
 	call mat_destroy(A,ierr)
 	call mat_destroy(B,ierr)	

 	if(myrank==0) print *, "==============Test mat_ones==============="
 	call mat_create(A,m,n,ierr)
 	call mat_ones(A,ierr)
 	if(debug) call mat_view(A,ierr)
 	call mat_destroy(A,ierr)	

 	if(myrank==0) print *, "==============Test mat_zeros==============="
 	call mat_create(A,m,n,ierr)
 	call mat_zeros(A,ierr)
 	if(debug) call mat_view(A,ierr)
 	call mat_destroy(A,ierr)	

 	if(myrank==0) print *, "==============Test mat_seq==============="
 	call mat_create(A,m,n,ierr)
 	call mat_seq(A,ierr)
 	if(debug) call mat_view(A,ierr)
 	call mat_destroy(A,ierr)	
    
 	if(myrank==0) print *, "==============Test mat_eye==============="
   	call mat_create(A,m,m,ierr)
    call mat_eye(A,ierr)	
    if(debug) call mat_view(A,ierr)
    call mat_destroy(A,ierr)	

 	call mat_create(A,m,2*m,ierr)
 	call mat_eye(A,ierr)	
 	if(debug) call mat_view(A,ierr)
 	call mat_destroy(A,ierr)	

  	call mat_create(A,2*m,m,ierr)
   	call mat_eye(A,ierr)	
  	if(debug) call mat_view(A,ierr)
   	call mat_destroy(A,ierr)	


 	if(myrank==0) print *, "==============Test mat_copy==============="
 	call mat_create(A,m,n,ierr)
 	call mat_seq(A,ierr)
 	call mat_copy(A,B,ierr)
 	if(debug) call mat_view(B,ierr)
 	call mat_destroy(A,ierr)	
 	call mat_destroy(B,ierr)
 
    if(myrank==0) print *, "==============Test mat_hjoin==============="
    call mat_create(A,m,n,ierr)
    call mat_create(B,m,2*n,ierr)
    call mat_seq(A,ierr)
    call mat_hjoin(A,A,B,ierr)
    if(debug) call mat_view(B,ierr)
    call mat_destroy(A,ierr)
    call mat_destroy(B,ierr)	
    
    call mat_create(A1,m,n,ierr)
    call mat_create(A2,m,m,ierr)
    call mat_create(B,m,m+n,ierr)
    call mat_zeros(A1,ierr)
    call mat_eye(A2,ierr)
    call mat_hjoin(A1,A2,B,ierr)
    if(debug) call mat_view(B,ierr)
    call mat_destroy(A1,ierr)	
    call mat_destroy(A2,ierr)	
    call mat_destroy(B,ierr)	
   

	if(myrank==0) print *, "==============Test mat_rep==============="
	call mat_create(A,m,n,ierr)
	call mat_create(B,3*m,2*n,ierr)
	call mat_seq(A,ierr)
	call mat_rep(A,3,2,B,ierr)
	if(debug) call mat_view(B,ierr)
	call mat_destroy(A,ierr)	
	call mat_destroy(B,ierr)	
    
	if(myrank==0) print *, "==============End==============="
	call PetscFinalize(ierr)
end program
