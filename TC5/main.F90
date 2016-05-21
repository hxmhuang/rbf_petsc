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
    Mat			    nodes,pts,PX,PY,PZ,F
    !Mat				W1,W2,W3
    Mat				DM_data,DM_eval,IM,EM
	Vec				u,x,rhs,exact,s,norm
	PetscMPIInt		myrank,mysize
	PetscErrorCode	ierr
	
    PetscReal		ep,tend,gamma,pi,lam_c,thm_c,mR,hm0
	PetscInt	    fd, order, dim, dt, N	
	PetscBool		debug
    PetscScalar     alpha,angle,radius,gh0,u0,g
    character*100   filename 
    
    debug = .false.
    alpha=1.0

    fd=31
    tend=15
    order=4
    dim=2
	ep=2.7
    dt=900
    gamma=-2.98e-17
    N=4
    filename="md001.00004"
    
    
    angle=0         ! Angle of rotation measured from the equator.
    u0=20           ! Speed of rotation in meters/second
    radius=6.37122e6     ! Mean radius of the earth (meters).
    g = 9.80616       ! Gravitational constant (m/s^2).
    gh0 = g*5960      ! Initial condition for the geopotential field (m^2/s^2).
  
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	!                 Beginning of program
	! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

	call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr)
	call MPI_Comm_rank(PETSC_COMM_WORLD,mysize,ierr)
    
    call PetscOptionsGetBool(PETSC_NULL_OBJECT,PETSC_NULL_CHARACTER,'-debug',debug,PETSC_NULL_BOOL,ierr)
     if(myrank==0) then 
        print *, "============Input paramenters============"
        print *, "debug=",debug
     endif 
	 
    call mat_create(nodes,N,4,ierr)
    
    call mat_load(filename,nodes,ierr)
    
    call rbf_cart2sph(nodes,pts,ierr)
    
    call rbf_project_cart2sph(nodes,PX,PY,PZ,ierr)
   
    call rbf_coriolis_force(nodes,angle,F,ierr)



    if(debug) then
        if(myrank==0) print *, ">nodes="
        call mat_view(nodes,ierr)
        if(myrank==0) print *, ">pts="
        call mat_view(pts,ierr)
        if(myrank==0) print *, ">PX="
        call mat_view(PX,ierr)
        if(myrank==0) print *, ">PY="
        call mat_view(PY,ierr)
        if(myrank==0) print *, ">PZ="
        call mat_view(PZ,ierr)
        if(myrank==0) print *, ">F="
        call mat_view(F,ierr)
 	endif

 	call mat_destroy(nodes,ierr)	
 	call mat_destroy(PX,ierr)	
 	call mat_destroy(PY,ierr)	
 	call mat_destroy(PZ,ierr)	

    call PetscFinalize(ierr)

end program
