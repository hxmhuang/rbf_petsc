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
	Mat				A,A1,A2,B,dsites,ctrs
	!PetscReal		error
	PetscMPIInt		myrank,mysize
	PetscErrorCode	ierr
	PetscReal		ep
	PetscInt		meval,neval,m,n
	!PetscInt rstart,rend,r;
	!integer i

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

	call mat_create(dsites,m*n,2,ierr)
	call rbf_createpoints(dsites,m,n,ierr)
	
	print *, "==============mat_eprod==============="
	call mat_create(A,m*n,2,ierr)
	call mat_eprod(dsites,dsites,A,ierr)
	!call mat_view(A,ierr)
	call mat_destroy(dsites,ierr)
	call mat_destroy(A,ierr)	

	print *, "==============mat_ones==============="
	call mat_create(A,m,n,ierr)
	call mat_ones(A,ierr)
	!call mat_view(A,ierr)
	call mat_destroy(A,ierr)	

	print *, "==============mat_zeros==============="
	call mat_create(A,m,n,ierr)
	call mat_zeros(A,ierr)
	!call mat_view(A,ierr)
	call mat_destroy(A,ierr)	

	print *, "==============mat_seq==============="
	call mat_create(A,m,n,ierr)
	call mat_seq(A,ierr)
	!call mat_view(A,ierr)
	call mat_destroy(A,ierr)	

	print *, "==============mat_diag==============="
	call mat_create(A,4,4,ierr)
	call mat_diag(A,ierr)	
	!call mat_view(A,ierr)
	call mat_destroy(A,ierr)	
	
	print *, "==============mat_copy==============="
	call mat_create(A,m,n,ierr)
	call mat_seq(A,ierr)
	call mat_copy(A,B,ierr)
	!call mat_view(B,ierr)
	call mat_destroy(A,ierr)	
	call mat_destroy(B,ierr)
 
!   print *, "==============mat_hjoin1==============="
!   call mat_create(A,m,n,ierr)
!   call mat_create(B,m,2*n,ierr)
!   call mat_seq(A,ierr)
!   call mat_hjoin(A,A,B,ierr)
!   !call mat_view(B,ierr)
!   call mat_destroy(A,ierr)
!   call mat_destroy(B,ierr)	
!   
!   print *, "==============mat_hjoin2==============="
!   call mat_create(A1,m,n,ierr)
!   call mat_create(A2,m,m,ierr)
!   call mat_create(B,m,m+n,ierr)
!   call mat_zeros(A1,ierr)
!   call mat_diag(A2,ierr)
!   call mat_hjoin(A1,A2,B,ierr)
!   !call mat_view(B,ierr)
!   call mat_destroy(A1,ierr)	
!   call mat_destroy(A2,ierr)	
!   call mat_destroy(B,ierr)	
!   
!   print *, "==============mat_mhjoin==============="
!   call mat_create(A1,m,m,ierr)
!   call mat_create(A2,m,n,ierr)
!   call mat_diag(A1,ierr)
!   call mat_zeros(A2,ierr)
!   
!   call mat_create(B,2*m,2*m+1*n,ierr)
!   call mat_mhjoin(A1,2,A2,1,B,ierr)
!   !call mat_view(B,ierr)
!   call mat_destroy(A1,ierr)	
!   call mat_destroy(A2,ierr)	
!   call mat_destroy(B,ierr)	
    
    print *, "==============mat_rep==============="
    call mat_create(A,m,n,ierr)
    call mat_create(B,3*m,2*n,ierr)
    call mat_seq(A,ierr)
    call mat_rep(A,3,2,B,ierr)
    call mat_view(B,ierr)
    call mat_destroy(A,ierr)	
    call mat_destroy(B,ierr)	
    
	print *, "============================="
	call PetscFinalize(ierr)
end program
