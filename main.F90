! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Implementing the rbf_interpolation2D with the ASM preconditioner
! -----------------------------------------------------------------------
program main
	use matrixalgebra
	use particle 
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
	Vec				u,x,b,xdim,ydim,rhs
	!PetscReal		error
	PetscMPIInt		myrank,mysize
	PetscErrorCode	ierr
	PetscReal		ep
	PetscInt		npoints,neval,M,N
	!PetscInt rstart,rend,r;
	!integer i

	ep=4.1
	npoints=3
	neval=30
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	!                 Beginning of program
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

	call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr)
	call MPI_Comm_rank(PETSC_COMM_WORLD,mysize,ierr)

	M=npoints*npoints
	N=neval*neval

	! generate some vectors: x,b,u
	call VecCreate(PETSC_COMM_WORLD,x,ierr)
	call VecSetSizes(x,PETSC_DECIDE,M,ierr)
	call VecSetFromOptions(x,ierr)
	call VecDuplicate(x,b,ierr)
	call VecDuplicate(x,u,ierr)
	call VecDuplicate(x,xdim,ierr)
	call VecDuplicate(x,ydim,ierr)
	call VecDuplicate(x,rhs,ierr)

	! generate matrix A with size M*M
	call MatCreate(PETSC_COMM_WORLD,A,ierr);
	call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,M,M,ierr)
	call MatSetFromOptions(A,ierr)
	call MatSetUp(A,ierr)

	! generate the sample points
	! call createpoints(npoints,xdim,ydim,rhs,ierr)
	call generatepoints(npoints,dsites,rhs,ierr)

	call MatDuplicate(dsites,MAT_COPY_VALUES,ctrs,ierr)

	call ma_repmat(dsites,2,2,ctrs,ierr)

	call ma_eprod(dsites,dsites,ctrs,ierr)

	call VecDestroy(x,ierr)
	call VecDestroy(b,ierr)
	call VecDestroy(u,ierr)
	call VecDestroy(rhs,ierr)
	call MatDestroy(dsites,ierr)
	call MatDestroy(ctrs,ierr)
	call PetscFinalize(ierr)

end program
