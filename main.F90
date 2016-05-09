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

	! generate some vectors: x,b,u
	call vec_create(x,m*n,ierr)
	call vec_duplicate(x,b,ierr)
	call vec_duplicate(x,u,ierr)
	call vec_duplicate(x,rhs,ierr)

	! generate matrix A with size M*M
	call mat_create(A,m*n,m*n,ierr)

	
	print *, "==============createpoints & testfunctionD==============="
	call mat_create(dsites,m*n,2,ierr)
	call rbf_createpoints(dsites,m,n,ierr)
	call rbf_testfunctionD(dsites,rhs,ierr)
	call mat_view(dsites,ierr)
	call vec_view(rhs,ierr)
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
