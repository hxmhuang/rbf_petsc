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
	Mat				A,dsites,ctrs
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
	call mat_view(A,ierr)
	call mat_destroy(A,ierr)	
	
	call mat_destroy(dsites,ierr)	
	call mat_destroy(ctrs,ierr)	

	print *, "==============mat_ones==============="
	!call mat_ones(ctrs,6,6,ierr)
	call mat_create(A,6,5,ierr)
	call mat_ones(A,ierr)
	call mat_view(A,ierr)
	call mat_destroy(A,ierr)	

	!call MatView(ctrs,PETSC_VIEWER_STDOUT_WORLD,ierr)
	print *, "==============mat_zeros==============="
	call mat_create(A,5,6,ierr)
	call mat_zeros(A,ierr)
	call mat_view(A,ierr)
	call mat_destroy(A,ierr)	

	print *, "==============mat_diag==============="
	call mat_create(A,4,4,ierr)
	call mat_diag(A,ierr)	
	call mat_view(A,ierr)
	call mat_destroy(A,ierr)	

	print *, "============================="
	!call mat_repmat(dsites,2,2,ctrs,ierr)
	!call mat_repmat(dsites,2,2,ctrs,ierr)

	call mat_destroy(A,ierr)
	call mat_destroy(dsites,ierr)
	call mat_destroy(ctrs,ierr)
	call PetscFinalize(ierr)

end program
