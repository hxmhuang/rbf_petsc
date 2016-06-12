subroutine rbfmatrix_fd_hypre(atm,ep,fdsize,order,dim,DPx,DPy,DPz,L,ierr)
    use dm
    use rbf
    use atm_type
    implicit none
	
    type(AtmType),	intent(inout) :: atm
    real(kind=8),	intent(in)    :: ep
    integer,		intent(in)    :: fdsize,order,dim
    type(Matrix),	intent(out)   :: DPx,DPy,DPz,L 
    integer,intent(out)           :: ierr
    type(Matrix)        		  :: A,B,idx,imat
    type(Matrix)        		  :: weightsDx,weightsDy,weightsDz,weightsL
    type(Matrix)        		  :: ind_i,ind_j 
	integer 					  :: N,j
 
    N=atm%pts%x%nrow	
    
    weightsDx=dm_zeros(N,N)		
    weightsDy=dm_zeros(N,N)		
    weightsDz=dm_zeros(N,N)		
    weightsL =dm_zeros(N,N)		

    !ind_i =dm_zeros(N*fdsize,1)		
    !ind_j =dm_zeros(N*fdsize,1)
    
    A=dm_ones(fdsize+1,fdsize+1)
    call dm_setvalue(A,fdsize,fdsize,0,ierr)	
    B=dm_zeros(fdsize+1,1)	
	
	call knnsearch(idx,ierr) 
	call dm_view(idx,ierr)	
   
	do j=1,N
		!imat= dm_getrow(idx,j)
		
	    !dp=atm%pts%x * dm_getsubmatrix(		
		




	enddo
	call dm_destroy(weightsDx,ierr) 
	call dm_destroy(weightsDy,ierr) 
	call dm_destroy(weightsDz,ierr) 
	call dm_destroy(weightsL,ierr) 
	call dm_destroy(A,ierr) 
	call dm_destroy(B,ierr) 
	call dm_destroy(idx,ierr) 
	!call dm_view(A,ierr)	
end subroutine

! -----------------------------------------------------------------------
!Matlab example:
!rbf = @(ep,rd2) exp(-ep^2*rd2)
! -----------------------------------------------------------------------
subroutine rbf(ep,rd2,A,ierr)
    use dm
	implicit none
	real(kind=8),	intent(in)		::  ep		
	type(Matrix),	intent(in)		::	rd2 
	type(Matrix),	intent(out)		::	A 
	integer,		intent(out)		::	ierr
	A=dm_exp((-1.0)*(ep**2)*rd2)
	ierr=0
end subroutine

! -----------------------------------------------------------------------
!Matlab example:
!drbf = @(ep,rd2) -2*ep^2*exp(-ep^2*rd2)
! -----------------------------------------------------------------------
subroutine drbf(ep,rd2,A,ierr)
    use dm
	implicit none
	real(kind=8),	intent(in)		::  ep		
	type(Matrix),	intent(in)		::	rd2 
	type(Matrix),	intent(out)		::	A 
	integer,		intent(out)		::	ierr
	A=(-2.0)*(ep**2)*dm_exp((-1.0)*(ep*ep)*rd2)
	ierr=0
end subroutine


! -----------------------------------------------------------------------
!At present, we load the neighbour from a file
! -----------------------------------------------------------------------
!subroutine knnsearch(nodes,fdsize,nn,ierr)
subroutine knnsearch(nn,ierr)
	use dm
	implicit none
    !type(Matrix),	intent(in)	::  nodes	
	!integer,	    intent(in)	::  fdsize	
    ! the nearest neighbours
    type(Matrix),	intent(out)	::  nn	
	integer,		intent(out)	::	ierr
    character*100   filename
    
    filename="nn.md002.00009.txt"
    call dm_load(filename,nn,ierr)
end subroutine



