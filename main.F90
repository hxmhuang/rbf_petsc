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
	Mat				A,Tmp,dsites,ctrs
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
	call ma_veccreate(x,m*n,ierr)
	call ma_vecduplicate(x,b,ierr)
	call ma_vecduplicate(x,u,ierr)
	call ma_vecduplicate(x,rhs,ierr)

	! generate matrix A with size M*M
	call ma_matcreate(A,m*n,m*n,ierr)

	! generate the sample points
	!call ma_matcreate(dsites,M,2,ierr)
	!call createpoints(dsites,npoints,npoints,ierr)
	
	print *, "==============createpoints & testfunctionD==============="
	call ma_matcreate(dsites,m*n,2,ierr)
	call createpoints(dsites,m,n,ierr)
	call testfunctionD(dsites,rhs,ierr)
	call ma_matview(dsites,ierr)
	call ma_vecview(rhs,ierr)
	call ma_matcopy(dsites,ctrs,ierr)
	
	print *, "==============ma_eprod==============="
	call ma_matcreate(Tmp,m*n,2,ierr)
	call ma_eprod(dsites,dsites,Tmp,ierr)
	call ma_matview(Tmp,ierr)
	call ma_matdestroy(Tmp,ierr)	
	
	call ma_matdestroy(dsites,ierr)	
	call ma_matdestroy(ctrs,ierr)	

	print *, "==============ma_ones==============="
	!call ma_ones(ctrs,6,6,ierr)
	call ma_matcreate(Tmp,6,5,ierr)
	call ma_ones(Tmp,ierr)
	call ma_matview(Tmp,ierr)
	call ma_matdestroy(Tmp,ierr)	

	!call MatView(ctrs,PETSC_VIEWER_STDOUT_WORLD,ierr)
	print *, "==============ma_zeros==============="
	call ma_matcreate(Tmp,5,6,ierr)
	call ma_zeros(Tmp,ierr)
	call ma_matview(Tmp,ierr)
	call ma_matdestroy(Tmp,ierr)	

	print *, "==============ma_diag==============="
	call ma_matcreate(Tmp,4,4,ierr)
	call ma_diag(Tmp,ierr)	
	call ma_matview(Tmp,ierr)
	call ma_matdestroy(Tmp,ierr)	

	print *, "============================="
	!call ma_repmat(dsites,2,2,ctrs,ierr)
	!call ma_repmat(dsites,2,2,ctrs,ierr)

	call ma_vecdestroy(x,ierr)
	call ma_vecdestroy(b,ierr)
	call ma_vecdestroy(u,ierr)
	call ma_vecdestroy(rhs,ierr)
	call ma_matdestroy(a,ierr)
	call ma_matdestroy(dsites,ierr)
	call ma_matdestroy(ctrs,ierr)
	call ma_matdestroy(tmp,ierr)
	call PetscFinalize(ierr)

end program
