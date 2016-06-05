program main 

    use dm 
    implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscviewer.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
    type(Matrix)    :: A,B,C,D,E,F,G,H 
    type(Matrix)    :: X,Y,Z,U 
    integer         :: myrank, mysize 
    integer         :: m,n 
    real(kind=8)    :: ep,alpha
    logical         :: debug 
    integer         :: ierr
    character(len=50):: filename
    debug=.false.

    ierr=dm_init()
    
    myrank=dm_comm_rank()
    
    mysize=dm_comm_size()
    
    m=dm_get_int('-m')
    n=dm_get_int('-n')
    ep=dm_get_real('-ep')
    !debug=dm_get_bool('-debug')

    call PetscOptionsGetBool(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,'-debug',debug,PETSC_NULL_BOOL,ierr)
    
    if(myrank==0) then 
       print *, "==============Input paramenters==========="
        print *, "m=",m,",n=",n,"ep=",ep,"debug=",debug
     endif 
	
    
 	if(myrank==0) print *, "==============Test rbf_createpoints============"
!   A=dm_zeros(m,n)
!   if(debug) then
!       if(myrank==0) print *, ">A="
!       ierr=dm_view(A)
!	endif
!   ierr=dm_destroy(A)




    call PetscFinalize(ierr)
end program
