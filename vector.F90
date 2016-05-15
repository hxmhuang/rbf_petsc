! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Vector Library 
! -----------------------------------------------------------------------
#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>

module vector 

contains
subroutine vec_create(v,m,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
	PetscInt,		intent(in)	::	m	
	Vec,			intent(out)	::	v
	PetscErrorCode,	intent(out)	::	ierr
	! generate vector v with size m
	call VecCreate(PETSC_COMM_WORLD,v,ierr)
	call VecSetSizes(v,PETSC_DECIDE,m,ierr)
	call VecSetFromOptions(v,ierr)
end subroutine

subroutine vec_duplicate(v1,v2,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
	Vec,		intent(in)	::	v1	
	Vec,			intent(out)	::	v2
	PetscErrorCode,	intent(out)	::	ierr
	! duplicate vector v2 using v1
	call VecDuplicate(v1,v2,ierr)
end subroutine

subroutine vec_destroy(v,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
	Vec,			intent(in)	::	v
	PetscErrorCode,	intent(out)	::	ierr
	! destroy vector v
	call VecDestroy(v,ierr)
end subroutine


subroutine vec_view(v,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Vec,			intent(in)	::	v
	PetscErrorCode,	intent(out)	::	ierr
	! veiw vec A
	call VecView(v,PETSC_VIEWER_STDOUT_WORLD,ierr)
end subroutine



end module
