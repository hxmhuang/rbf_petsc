#include <petsc/finclude/petscsysdef.h>
#include <petsc/finclude/petscvecdef.h>
#include "mat_math_type.h"

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
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	PetscInt,		intent(in)		::	m,n	
	type(Matrix),	intent(out)		::	A 
	PetscErrorCode,	intent(out)		::	ierr
		
	PetscReal	                    ::  xmin,xmax,ymin,ymax,dx,dy,xcord,ycord
	PetscInt	                    ::  ista,iend
	PetscInt	                    ::  nrow,ncol,nlocal
    PetscInt,allocatable            ::  idxm(:),idxn(:) 
    PetscScalar,allocatable         ::  rows(:)
	integer		                    ::  i,j
	
    xmin=0.0
	xmax=1.0
	ymin=0.0
	ymax=1.0
	dx= (xmax-xmin)/(m-1)
	dy= (ymax-ymin)/(n-1)
	
	call MatGetSize(A%x,nrow,ncol,ierr)
    if(ncol/=2) then
        print *, "Error in rbf_createpoints: The column size of matrix A should be 2."
        stop
    endif
    
    call MatGetOwnershipRange(A%x,ista,iend,ierr)
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
    call MatSetValues(A,nlocal,idxm,ncol,idxn,rows,INSERT_VALUES,ierr) 
    
	call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
    deallocate(idxm,idxn,rows)
end subroutine


subroutine rbf_testfunctionD(A,v,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	type(Matrix),	intent(in)		::	A 
	type(Matrix),	intent(out)	    ::	v 
	PetscErrorCode,	intent(out)		::	ierr
	
	PetscInt						::	nrow,ncol
	PetscInt						::  ista,iend
	PetscScalar,allocatable			::	row(:)
	PetscReal						::	xcord,ycord,res
	integer							:: 	i

	call MatGetSize(A%x,nrow,ncol,ierr)
	if(ncol/=2)then
		print *, "Error: the column size of Matrix A in testfunctionD must equal to 2"
		stop	
	endif	
	
	call MatGetOwnershipRange(A%x,ista,iend,ierr)
	!print *,">istat=",ista,"iend=",iend,">ncol=",ncol
	allocate(row(ncol))
	
	do i=ista,iend-1
		call MatGetRow(A%x,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row,ierr)
		xcord=row(1)
		ycord=row(2)
		call testfunction(xcord,ycord,res)
		call VecSetValue(v,i,res,INSERT_VALUES,ierr)
		call MatRestoreRow(A,i,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,row,ierr)
	enddo
    
    call VecAssemblyBegin(v,ierr)
    call VecAssemblyEnd(v,ierr)
	deallocate(row)
end subroutine


subroutine testfunction(xcord,ycord,res)
	implicit none
#include <petsc/finclude/petscsys.h>
	PetscReal,intent(in):: xcord,ycord
	PetscReal,intent(out):: res 
	res= 0.75*exp(-((9*xcord-2)*(9*xcord-2)+(9*ycord-2)*(9*ycord-2))/4) &
		+0.75*exp(-((9*xcord+1)*(9*xcord+1))/49-((9*ycord+1)*(9*ycord+1))/10) &
		+0.5*exp(-((9*xcord-7)*(9*xcord-7)+(9*ycord-3)*(9*ycord-3))/4) &
		-0.2*exp(-(9*xcord-4)*(9*xcord-4)-(9*ycord-7)*(9*ycord-7))
end subroutine


subroutine rbf_distancematrix(A1,A2,B,ierr)
	use matrix 
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
	Mat,			intent(in)	   	:: A1,A2
	Mat,			intent(out)     :: B 
	PetscErrorCode,	intent(out)	    :: ierr
	Mat                     		:: W1,W2,W3 
	Mat                     		:: P1,P2,P3 
	PetscInt	            		:: nrow1,ncol1,nrow2,ncol2
	PetscScalar						:: alpha
	
    PetscMPIInt		myrank,mysize
	PetscLogEvent	ievent(13)
	call PetscLogEventRegister("1_dm_eprod",0, ievent(1), ierr)
	call PetscLogEventRegister("2_dm_sum",0, ievent(2), ierr)
	call PetscLogEventRegister("3_dm_rep",0, ievent(3), ierr)
	call PetscLogEventRegister("4_dm_xyt",0, ievent(4), ierr)
	call PetscLogEventRegister("5_dm_axpy",0, ievent(5), ierr)
	call PetscLogEventRegister("6_dm_destroy",0, ievent(6), ierr)
	call PetscLogEventRegister("7_dm_eprod",0, ievent(7), ierr)
	call PetscLogEventRegister("8_dm_sum",0, ievent(8), ierr)
	call PetscLogEventRegister("9_dm_trans",0, ievent(9), ierr)
	call PetscLogEventRegister("10_dm_rep",0, ievent(10), ierr)
	call PetscLogEventRegister("11_dm_axpy",0, ievent(11), ierr)
	call PetscLogEventRegister("12_dm_copy",0, ievent(12), ierr)
	call PetscLogEventRegister("13_dm_destroy",0, ievent(13), ierr)

    call MPI_Comm_rank(PETSC_COMM_WORLD,myrank,ierr)
	call MPI_Comm_rank(PETSC_COMM_WORLD,mysize,ierr)
 
	call MatGetSize(A1,nrow1,ncol1,ierr)
	call MatGetSize(A2,nrow2,ncol2,ierr)
    !print *, "nrow1=",nrow1,"nrow2=",nrow2
	if(ncol1/=ncol2)then
		print *, "Error in rbf_distancematrix: the matrix A1 and A2 should have the same column size."
		stop	
	endif
    
	call PetscLogEventBegin(ievent(1),ierr)
    call mat_eprod(A1,A1,W1,ierr)
	call PetscLogEventEnd(ievent(1),ierr)
    
	call PetscLogEventBegin(ievent(2),ierr)
    call mat_sum(W1,2,W2,ierr)
	call PetscLogEventEnd(ievent(2),ierr)
    
	call PetscLogEventBegin(ievent(3),ierr)
    call mat_rep(W2,1,nrow2,P1,ierr)
	call PetscLogEventEnd(ievent(3),ierr)
	!if(myrank==0) print *, ">P1="
	!call mat_view(P1,ierr)
	
	call PetscLogEventBegin(ievent(4),ierr)
	call mat_xyt(A1,A2,P2,ierr)
	call PetscLogEventEnd(ievent(4),ierr)
	!if(myrank==0) print *, ">P2="
	!call mat_view(P2,ierr)
	call PetscLogEventBegin(ievent(5),ierr)
	alpha=-2.0
	call mat_axpy(P1,alpha,P2,ierr)
	call PetscLogEventEnd(ievent(5),ierr)
	
	call PetscLogEventBegin(ievent(6),ierr)
    call mat_destroy(W1,ierr)
    call mat_destroy(W2,ierr)
	call PetscLogEventEnd(ievent(6),ierr)
	
	call PetscLogEventBegin(ievent(7),ierr)
	call mat_eprod(A2,A2,W1,ierr)
	call PetscLogEventEnd(ievent(7),ierr)
	call PetscLogEventBegin(ievent(8),ierr)
	call mat_sum(W1,2,W2,ierr)
	call PetscLogEventEnd(ievent(8),ierr)
	call PetscLogEventBegin(ievent(9),ierr)
	call mat_trans(W2,W3,ierr)
	call PetscLogEventEnd(ievent(9),ierr)
	call PetscLogEventBegin(ievent(10),ierr)
	call mat_rep(W3,nrow1,1,P3,ierr)
	call PetscLogEventEnd(ievent(10),ierr)
	!if(myrank==0) print *, ">P3="
	!call mat_view(P3,ierr)
	
	call PetscLogEventBegin(ievent(11),ierr)
	alpha=1.0
	call mat_axpy(P1,alpha,P3,ierr)	
	call PetscLogEventEnd(ievent(11),ierr)
    
	call PetscLogEventBegin(ievent(12),ierr)
    call mat_copy(P1,B,ierr) 
	call PetscLogEventEnd(ievent(12),ierr)
	!call PetscLogEventBegin(ievent(12),ierr)
    !call mat_math(P1,MAT_MATH_SQRT,B,ierr)
	!call PetscLogEventEnd(ievent(12),ierr)
	!if(myrank==0) print *, ">B="
	!call mat_view(B,ierr)

	call PetscLogEventBegin(ievent(13),ierr)
	call mat_destroy(W1,ierr)
    call mat_destroy(W2,ierr)
    call mat_destroy(W3,ierr)
	call mat_destroy(P1,ierr)
    call mat_destroy(P2,ierr)
    call mat_destroy(P3,ierr)
	call PetscLogEventEnd(ievent(13),ierr)
end subroutine

!B= exp(-ep^2*A). Note that A equals to r^2.
subroutine rbf_guassian(ep,A,B,ierr)
!    use matrix
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include "mat_math_type.h"
	PetscReal,		intent(in)		::  ep		
	Mat,			intent(in)		::	A 
	Mat,			intent(out)		::	B 
	PetscErrorCode,	intent(out)		::	ierr
    PetscScalar                     ::  alpha	
	Mat                 			::	W
    !rbf = @(e,r) exp(-(e*r).^2)
    call mat_copy(A,W,ierr)
    alpha=-ep*ep
    call mat_scale(W,alpha,ierr)
    call mat_math(W,MAT_MATH_EXP,B,ierr)
	call mat_destroy(W,ierr)
end subroutine


!B= -2*ep^2*exp(-ep^2*A). Note that A equals to r^2.
subroutine rbf_diff_guassian(ep,A,B,ierr)
!    use matrix
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include "mat_math_type.h"
	PetscReal,		intent(in)		::  ep		
	Mat,			intent(in)		::	A 
	Mat,			intent(out)		::	B 
	PetscErrorCode,	intent(out)		::	ierr
    PetscScalar                     ::  alpha	
	Mat                 			::	W
    !rbf = @(e,r) exp(-(e*r).^2)
    call mat_copy(A,W,ierr)
    alpha=-ep**2
    call mat_scale(W,alpha,ierr)
    call mat_math(W,MAT_MATH_EXP,B,ierr)
    alpha=2*alpha
    call mat_scale(B,alpha,ierr)
	call mat_destroy(W,ierr)
end subroutine


! -----------------------------------------------------------------------
! Transform Cartesian to spherical coordinates.
! [TH,PHI,R] = cart2sph(X,Y,Z) transforms corresponding elements of
!    data stored in Cartesian coordinates X,Y,Z to spherical
!    coordinates (azimuth TH, elevation PHI, and radius R).
!    TH and PHI are returned in radians.
! where,
!   azimuth = atan2(y,x)
!   elevation = atan2(z,sqrt(x.^2 + y.^2))
!   r = sqrt(x.^2 + y.^2 + z.^2)
! -----------------------------------------------------------------------
subroutine rbf_cart2sph(A,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include "mat_math_type.h"
	Mat,			intent(in)	::  A
	Mat,			intent(out)	::	B
	PetscErrorCode,	intent(out)	::	ierr
    
	PetscInt					::	nrow1,ncol1,nrow2,ncol2
	PetscInt,allocatable        ::	idx1(:),idx2(:)
	PetscScalar,allocatable     ::	row1(:),row2(:)
	PetscInt					::  ista,iend
	integer						::	i,j
    PetscScalar                 ::  newval 
    
    call MatGetSize(A,nrow1,ncol1,ierr)
    call MatGetOwnershipRange(A,ista,iend,ierr)
   
    nrow2=nrow1
    ncol2=3
    call mat_create(B,nrow2,ncol2,ierr)
    allocate(idx1(ncol1),idx2(ncol2),row1(ncol1),row2(ncol2))
    
    do i=ista,iend-1
        do j=1,ncol1
            idx1(j)=j-1
        enddo
        do j=1,ncol2
            idx2(j)=j-1
        enddo
        
        call MatGetRow(A,i,ncol1,idx1,row1,ierr)
        row2(1) = atan2(row1(2),row1(1)) 
        row2(2) = atan2(row1(3),sqrt(row1(1)**2+row1(2)**2))
        row2(3) = sqrt(row1(1)**2+row1(2)**2+row1(3)**2) 
        call MatRestoreRow(A,i,ncol1,idx1,row1,ierr)
        
    	call MatSetValues(B,1,i,ncol2,idx2,row2,INSERT_VALUES,ierr)
	enddo
    
    deallocate(idx1,idx2,row1,row2)
    call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
end subroutine




! -----------------------------------------------------------------------
! Variables for projecting an arbitrary Cartesian vector onto the surface of the sphere.
! -----------------------------------------------------------------------
subroutine rbf_project_cart2sph(A,PX,PY,PZ,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include "mat_math_type.h"
	Mat,			intent(in)	::  A
	Mat,			intent(out)	::  PX,PY,PZ	
	PetscErrorCode,	intent(out)	::	ierr
    
	PetscInt					::	nrow1,ncol1,nrow2,ncol2
	PetscInt,allocatable        ::	idx1(:),idx2(:)
	PetscScalar,allocatable     ::	row1(:),row2(:)
	PetscInt					::  ista,iend
	integer						::	i,j
    PetscScalar                 ::  x2,xy,y2,xz,z2,yz

    call MatGetSize(A,nrow1,ncol1,ierr)
    call MatGetOwnershipRange(A,ista,iend,ierr)
  
    !if(ncol /= 3) then
	!	print *, "Error in mat_project_cart2sph: the column of matrix A should be 3" 
	!	stop	
    !endif
    nrow2=nrow1
    ncol2=3
    call mat_create(PX,nrow2,ncol2,ierr)
    call mat_create(PY,nrow2,ncol2,ierr)
    call mat_create(PZ,nrow2,ncol2,ierr)

    allocate(idx1(ncol1),idx2(ncol2),row1(ncol1),row2(ncol2))
    
    do i=ista,iend-1
        do j=1,ncol1
            idx1(j)=j-1
        enddo
        do j=1,ncol2
            idx2(j)=j-1
        enddo
        
        call MatGetRow(A,i,ncol1,idx1,row1,ierr)
       
        x2=row1(1)**2
        y2=row1(2)**2
        z2=row1(3)**2
        xy=row1(1)*row1(2)
        xz=row1(1)*row1(3)
        yz=row1(2)*row1(3)
        call MatRestoreRow(A,i,ncol1,idx1,row1,ierr)
        
        row2(1)= 1-x2
        row2(2)= -xy
        row2(3)= -xz
    	call MatSetValues(PX,1,i,ncol2,idx2,row2,INSERT_VALUES,ierr)
        !print *,"PXrow2=", row2
        
        row2(1)= -xy 
        row2(2)= 1-y2
        row2(3)= -yz
    	call MatSetValues(PY,1,i,ncol2,idx2,row2,INSERT_VALUES,ierr)
        !print *,"PYrow3=", row2
 
        row2(1)= -xz 
        row2(2)= -yz
        row2(3)= 1-z2
    	call MatSetValues(PZ,1,i,ncol2,idx2,row2,INSERT_VALUES,ierr)
        !print *,"PZrow4=", row2

	enddo
    
    call MatAssemblyBegin(PX,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(PX,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PY,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(PY,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyBegin(PZ,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(PZ,MAT_FINAL_ASSEMBLY,ierr)
    deallocate(idx1,idx2,row1,row2)
end subroutine


! -----------------------------------------------------------------------
! Compute the coriolis force of each point in A. The angle is the rotation measured from the equator.
! -----------------------------------------------------------------------
subroutine rbf_coriolis_force(A,angle,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include "mat_math_type.h"
	Mat,			intent(in)	::  A
	PetscReal,      intent(in)  ::  angle	
	Mat,			intent(out)	::  B
	PetscErrorCode,	intent(out)	::	ierr
    
	PetscInt					::	nrow1,ncol1,nrow2,ncol2
	PetscInt,allocatable        ::	idx1(:),idx2(:)
	PetscScalar,allocatable     ::	row1(:),row2(:)
	PetscInt					::  ista,iend
	PetscReal                   ::  omega
    integer						::	i,j

    call MatGetSize(A,nrow1,ncol1,ierr)
    call MatGetOwnershipRange(A,ista,iend,ierr)
  
    nrow2=nrow1
    ncol2=1
    omega=7.292e-5  ! Rotation rate of the earth (1/seconds).
    call mat_create(B,nrow2,ncol2,ierr)

    allocate(idx1(ncol1),idx2(ncol2),row1(ncol1),row2(ncol2))
    
    do i=ista,iend-1
        do j=1,ncol1
            idx1(j)=j-1
        enddo
        do j=1,ncol2
            idx2(j)=j-1
        enddo
        
        call MatGetRow(A,i,ncol1,idx1,row1,ierr)
       
        row2(1)= 2*omega*(row1(1)*sin(angle)+row1(3)*cos(angle))
        call MatRestoreRow(A,i,ncol1,idx1,row1,ierr)
    	
        call MatSetValues(B,1,i,ncol2,idx2,row2,INSERT_VALUES,ierr)
        
	enddo
    
    call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
    deallocate(idx1,idx2,row1,row2)
end subroutine


! -----------------------------------------------------------------------
! Compute the profile of the mountain (multiplied by gravity) 
! -----------------------------------------------------------------------
subroutine rbf_mountain_profile(A,B,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include "mat_math_type.h"
	Mat,			intent(in)	::  A
	Mat,			intent(out)	::  B
	PetscErrorCode,	intent(out)	::	ierr
    
	PetscInt					::	nrow1,ncol1,nrow2,ncol2
	PetscInt,allocatable        ::	idx1(:),idx2(:)
	PetscScalar,allocatable     ::	row1(:),row2(:),row3(:),row4(:)
	PetscInt					::  ista,iend
	PetscReal                   ::  omega
    integer						::	i,j
    PetscReal		            ::  pi,lam_c,thm_c,mR,hm0,r2,g
    PetscScalar     alpha
    
    !Parameters for the mountain:
    alpha=1.0
    pi=4*atan(alpha)
    lam_c = -pi/2;
    thm_c = pi/6;
    mR = pi/9;
    hm0 = 2000;
    g = 9.80616       ! Gravitational constant (m/s^2).
 
    call MatGetSize(A,nrow1,ncol1,ierr)
    call MatGetOwnershipRange(A,ista,iend,ierr)
  
    nrow2=nrow1
    ncol2=1
    call mat_zeros(B,nrow2,ncol2,ierr)

    allocate(idx1(ncol1),idx2(ncol2),row1(ncol1),row2(ncol2))
   
    do i=ista,iend-1
        do j=1,ncol1
            idx1(j)=j-1
        enddo
        do j=1,ncol2
            idx2(j)=j-1
        enddo
        
        call MatGetRow(A,i,ncol1,idx1,row1,ierr)
       
        r2= (row1(1)-lam_c)**2+(row1(2)-thm_c)**2
        call MatRestoreRow(A,i,ncol1,idx1,row1,ierr)
        if(r2 < mR**2) then
            row2(1)=g*hm0*(1-sqrt(r2)/mR)
            call MatSetValues(B,1,i,ncol2,idx2,row2,INSERT_VALUES,ierr)
        endif
	enddo
    
    call MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY,ierr)
	call MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY,ierr)
    deallocate(idx1,idx2,row1,row2)
end subroutine



! -----------------------------------------------------------------------
! Compute the coefficient materix DPX, DPY, and DPZ.
! -----------------------------------------------------------------------
subroutine rbf_matrix_fd_hyper(nodes,ep,fdsize,order,dims,DPx,DPy,DPz,L,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include "mat_math_type.h"
    
    Mat,        intent(in)          ::  nodes
    PetscReal,  intent(in)          ::  ep
    PetscInt,   intent(in)          ::  fdsize,order,dims
    Mat,        intent(out)         ::  DPX,DPY,DPZ,L
    PetscErrorCode, intent(out)     ::  ierr
    integer                         ::  i,j,k

	PetscInt					    ::	N,ncol
	PetscInt					    ::	nrow1,ncol1,nrow2,ncol2
	PetscInt,allocatable            ::	idx1(:),idx2(:)
	PetscScalar,allocatable         ::	row1(:),row2(:)
	PetscInt,allocatable         	::  imat(:)	
	PetscInt,allocatable         	::	rowint1(:),rowint2(:)
	PetscInt					    ::  ista,iend
    PetscScalar                     ::  alpha
   
    Vec     ::  weightsDx, weightsDy, weightsDz, weightsL
    Vec     ::  int_i, int_j 
    Mat     ::  A,NN
    Vec     ::  B
    !IS      ::  idx

    call MatGetSize(nodes,N,ncol,ierr)

    call vec_create(weightsDx,N*fdsize,ierr)
    call vec_duplicate(weightsDx,weightsDy,ierr)
    call vec_duplicate(weightsDx,weightsDz,ierr)
    call vec_duplicate(weightsDx,weightsL,ierr)
    call vec_duplicate(weightsDx,int_i,ierr)
    call vec_duplicate(weightsDx,int_j,ierr)
    
    call mat_ones(A,fdsize+1,fdsize+1,ierr)
    
    alpha=0.0
    call mat_setvalue(A,fdsize,fdsize,alpha,ierr)
    call vec_create(B,fdsize+1,ierr)
    !call mat_view(A,ierr)
    
    call rbf_knnsearch(nodes,fdsize,NN,ierr)

    call MatGetOwnershipRange(NN,ista,iend,ierr)                                       
    allocate(row1(fdsize),rowint1(fdsize),row2(fdsize))
    allocate(imat(fdsize),idx1(fdsize))
            
    do i=ista,iend-1
        call MatGetRow(NN,i,ncol1,idx1,row1,ierr)           
		imat=int(row1)
        !print *,"row1=",rowint1    
        call MatRestoreRow(NN,i,ncol1,idx1,row1,ierr)           
		
		!call VecSetValue(
        
    enddo   
    deallocate(row1,rowint1,row2,imat,idx1) 
 
    !call mat_view(NN,ierr)
    call mat_destroy(NN,ierr)
    call mat_destroy(A,ierr)
    call vec_destroy(B,ierr)
    call vec_destroy(weightsDx,ierr)
    call vec_destroy(weightsDy,ierr)
    call vec_destroy(weightsDz,ierr)
    call vec_destroy(weightsL,ierr)
    call vec_destroy(int_i,ierr)
    call vec_destroy(int_j,ierr)

end subroutine


subroutine rbf_knnsearch(nodes,fdsize,nn,ierr)
	implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscmat.h>
#include "mat_math_type.h"
 
    Mat,			intent(in)	::  nodes	
	PetscInt,	    intent(in)	::  fdsize	
    ! the nearest neighbours
    Mat,			intent(out)	::  nn	
	PetscErrorCode,	intent(out)	::	ierr
    PetscInt                    ::  nrow,ncol	
    character*100   filename
    
    call MatGetSize(nodes,nrow,ncol,ierr)
    
    call mat_create(nn,nrow,fdsize,ierr)
    
    filename="nn.md002.00009.txt"
    call mat_load(filename,nn,ierr)

end subroutine



end module
