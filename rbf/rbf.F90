#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>
#include "mat_type.h"

module	rbf 
    use dm 

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
 	
   	A=dm_create(m*n,2)
 	
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
    ierr=dm_setvalues(A,nlocal,idxm,ncol,idxn,rows) 
    
    deallocate(idxm,idxn,rows)
    !    ierr=dm_view(A)
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
  	ierr=dm_setvalue(W1,0,0,7)
  	ierr=dm_setvalue(W1,0,1,sqrt(10.0))
	!W2=[7 3]
 	ierr=dm_setvalue(W2,0,0,7)
    ierr=dm_setvalue(W2,0,1,3)
	!W3=[4 7]
    ierr=dm_setvalue(W3,0,0,4)
    ierr=dm_setvalue(W3,0,1,7)
     	
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
	
 	ierr=dm_destroy(A1)
    ierr=dm_destroy(A2)
    ierr=dm_destroy(A3)
    ierr=dm_destroy(W1)
    ierr=dm_destroy(W2)
    ierr=dm_destroy(W3)
    ierr=dm_destroy(W4)
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
	
	ierr=dm_destroy(W1)
    ierr=dm_destroy(W2)
    ierr=dm_destroy(W3)
end subroutine

!B= exp(-ep^2*A). Note that A equals to r^2.
subroutine rbf_guassian(ep,A,B,ierr)
	implicit none
	real*8,			intent(in)	::  ep		
	type(Matrix),	intent(in)	::	A 
	type(Matrix),	intent(out)	::	B 
	integer,		intent(out)	::	ierr
   
	B=dm_exp((-ep**2)*A) 
end subroutine



end module
