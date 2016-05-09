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

	call ma_matcreate(dsites,m*n,2,ierr)
	call createpoints(dsites,m,n,ierr)
	
	print *, "==============ma_eprod==============="
	call ma_matcreate(A,m*n,2,ierr)
	call ma_eprod(dsites,dsites,A,ierr)
	call ma_matview(A,ierr)
	call ma_matdestroy(A,ierr)	
	
	call ma_matdestroy(dsites,ierr)	
	call ma_matdestroy(ctrs,ierr)	

	print *, "==============ma_ones==============="
	!call ma_ones(ctrs,6,6,ierr)
	call ma_matcreate(A,6,5,ierr)
	call ma_ones(A,ierr)
	call ma_matview(A,ierr)
	call ma_matdestroy(A,ierr)	

	!call MatView(ctrs,PETSC_VIEWER_STDOUT_WORLD,ierr)
	print *, "==============ma_zeros==============="
	call ma_matcreate(A,5,6,ierr)
	call ma_zeros(A,ierr)
	call ma_matview(A,ierr)
	call ma_matdestroy(A,ierr)	

	print *, "==============ma_diag==============="
	call ma_matcreate(A,4,4,ierr)
	call ma_diag(A,ierr)	
	call ma_matview(A,ierr)
	call ma_matdestroy(A,ierr)	

	print *, "============================="
	!call ma_repmat(dsites,2,2,ctrs,ierr)
	!call ma_repmat(dsites,2,2,ctrs,ierr)

	call ma_matdestroy(A,ierr)
	call ma_matdestroy(dsites,ierr)
	call ma_matdestroy(ctrs,ierr)
	call PetscFinalize(ierr)

end program
