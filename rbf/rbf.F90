#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>
#include "mat_type.h"

module	rbf 
    use dm 
	implicit none

contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Create the points in rbf application.
! A: Matrix with m*n rows and 2 columns
! m: there are m points in x direction 
! n: there are n points in y direction 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine rbf_createpoints(A,m,n,ierr)
	implicit none
	integer,		intent(in)	::	m,n	
	type(Matrix),   intent(out)	::	A 
	integer,		intent(out)	::	ierr
		
	real				::  xmin,xmax,ymin,ymax
	real				::  dx,dy,xcord,ycord
	integer				::  ista,iend
	integer				::  ncol,nlocal
    integer,allocatable ::  idxm(:),idxn(:) 
    real,allocatable 	::  rows(:)
	integer		        ::  i,j
	
    xmin=0.0
	xmax=1.0
	ymin=0.0
	ymax=1.0
	dx= (xmax-xmin)/(m-1)
	dy= (ymax-ymin)/(n-1)
 	
   	A=dm_zeros(m*n,2)
 	
    if(A%ncol/=2) then
        print *, "Error in rbf_createpoints: The column size of matrix A should be 2."
        stop
    endif
   	
    ncol=A%ncol 
    ista=A%ista
    iend=A%iend
    nlocal=iend-ista
    allocate(idxm(nlocal),idxn(ncol),rows(2*nlocal))
    
    do i=ista,iend-1
    	xcord = xmin+(i/n)*dx
    	ycord = ymin+mod(i,n)*dy
        rows((i-ista)*ncol+1)=xcord
        rows((i-ista)*ncol+2)=ycord
        idxm(i-ista+1)=i
    enddo
    do j=1,ncol
        idxn(j)=j-1
    enddo
    !print *, "nlocal=",nlocal,"idxm=",idxm,"idxn=",idxn,"rows=",rows
    call dm_setvalues(A,nlocal,idxm,ncol,idxn,rows,ierr) 
    
    deallocate(idxm,idxn,rows)
end subroutine 


subroutine rbf_testfunction(A,B,ierr)
	implicit none
	type(Matrix),intent(in)		::	A 
	type(Matrix),intent(out)	::	B 
	integer,intent(out)			::	ierr
	type(Matrix)				::	A1,A2,A3 
	type(Matrix)				::	W1,W2,W3,W4
	integer						::  nrow,ncol
	
	W1=dm_zeros(1,2)
   	W2=dm_zeros(1,2)
  	W3=dm_zeros(1,2)
    !W1=[7 sqrt(10)] 
  	call dm_setvalue(W1,0,0,7,ierr)
  	call dm_setvalue(W1,0,1,sqrt(10.0),ierr)
	!W2=[7 3]
 	call dm_setvalue(W2,0,0,7,ierr)
    call dm_setvalue(W2,0,1,3,ierr)
	!W3=[4 7]
    call dm_setvalue(W3,0,0,4,ierr)
    call dm_setvalue(W3,0,1,7,ierr)
     	
    nrow=A%nrow
    ncol=A%ncol
    A1=dm_rep(W1,nrow,floor((ncol+1.0)/2.0))
    A2=dm_rep(W2,nrow,floor((ncol+1.0)/2.0))
    A3=dm_rep(W3,nrow,floor((ncol+1.0)/2.0))

    W1=dm_exp( (-1.0/4.0)*dm_sum(dm_squ(9*A-2),2) )	
 	W2=dm_exp( (-1.0)*dm_sum(dm_squ((9*A+1).ed.(A1)),2) )	
 	W3=dm_exp( (-1.0/4.0)*dm_sum(dm_squ(9*A-A2),2) )	
 	W4=dm_exp( (-1.0)*dm_sum(dm_squ(9*A-A3),2) )	
 
 	B=0.75*W1+0.75*W2+0.5*W3-0.2*W4
	
 	call dm_destroy(A1,ierr)
    call dm_destroy(A2,ierr)
    call dm_destroy(A3,ierr)
    call dm_destroy(W1,ierr)
    call dm_destroy(W2,ierr)
    call dm_destroy(W3,ierr)
    call dm_destroy(W4,ierr)
end subroutine


subroutine rbf_distancematrix(A,B,C,ierr)
	implicit none
	type(Matrix),intent(in)		::	A,B 
	type(Matrix),intent(out)	::	C 
	integer,intent(out)			::	ierr
	type(Matrix)				::	W1,W2,W3
	integer						::  m,n

	m=A%nrow
	n=B%nrow	

	W1=dm_rep(dm_sum((A .em. A),2),1,n)
	W2=2*dm_xyt(A,B)
	W3=dm_rep( dm_trans( dm_sum((B .em. B),2) ),m,1)
   	C=W1-W2+W3
	
	call dm_destroy(W1,ierr)
    call dm_destroy(W2,ierr)
    call dm_destroy(W3,ierr)
end subroutine


!B= exp(-ep^2*A). Note that A equals to r^2.
subroutine rbf_guassian(ep,A,B,ierr)
	implicit none
	real(kind=8),	intent(in)	::  ep		
	type(Matrix),	intent(in)	::	A 
	type(Matrix),	intent(out)	::	B 
	integer,		intent(out)	::	ierr
   
	B=dm_exp((-ep**2)*A) 
end subroutine


! -----------------------------------------------------------------------
! Variables for projecting an arbitrary Cartesian vector onto the surface of the sphere.
! -----------------------------------------------------------------------
subroutine rbf_projection(A,PX,PY,PZ,ierr)
	implicit none
	type(Matrix),	intent(in)	::  A
	type(Matrix),	intent(out)	::  PX,PY,PZ	
	integer,		intent(out)	::	ierr
    type(Matrix)				::  x,y,z,x2,xy,y2,xz,z2,yz
	
 	x=dm_getcol(A,0)
 	y=dm_getcol(A,1)
 	z=dm_getcol(A,2)
 	x2= 1-dm_squ(x)
    y2= 1-dm_squ(y)
    z2= 1-dm_squ(z)
    xy=	(-1.0)*x .em. y
    xz= (-1.0)*x .em. z
    yz= (-1.0)*y .em. z 
    
    PX= x2 	.hj. xy .hj. xz	
    PY= xy	.hj. y2 .hj. yz	
    PZ= xz 	.hj. yz .hj. z2	

	call dm_destroy(x,ierr)	
	call dm_destroy(y,ierr)	
	call dm_destroy(z,ierr)	
	call dm_destroy(x2,ierr)	
	call dm_destroy(y2,ierr)	
	call dm_destroy(z2,ierr)	
	call dm_destroy(xy,ierr)	
	call dm_destroy(xz,ierr)	
	call dm_destroy(yz,ierr)	
end subroutine


! -----------------------------------------------------------------------
! Compute the coriolis force of each point in A. The angle is the rotation measured from the equator.
! -----------------------------------------------------------------------
subroutine rbf_coriolis_force(A,angle,B,ierr)
	implicit none
	type(Matrix),	intent(in)	::  A
	real(kind=8),   intent(in)  ::  angle	
	type(Matrix),	intent(out)	::  B
	integer,		intent(out)	::	ierr
    type(Matrix)				::  x,z
	real(kind=8)     			::  omega 
    omega=7.292e-5  ! Rotation rate of the earth (1/seconds).
	
 	x=dm_getcol(A,0)
 	z=dm_getcol(A,2)
   	B=2*omega*(x*sin(angle)+z*cos(angle)) 
	ierr=0		
	call dm_destroy(x,ierr)	
	call dm_destroy(z,ierr)	
end subroutine



end module
