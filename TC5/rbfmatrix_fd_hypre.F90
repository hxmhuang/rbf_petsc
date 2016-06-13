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
    type(Matrix)        		  :: ximat,ximat1,ximat2,ximat3
    type(Matrix)        		  :: x,x1,x2,x3,xj
    type(Matrix)        		  :: px,py,pz
    type(Matrix)        		  :: p1,p2,p3
    type(Matrix)        		  :: rd2,rd2v 
    type(Matrix)        		  :: weightsDx,weightsDy,weightsDz,weightsL
    type(Matrix)        		  :: ind_i,ind_j 
	integer 					  :: N,j
	integer 					  :: ista,iend
    
	N=atm%pts%x%nrow	
    ista=atm%pts%x%ista
    iend=atm%pts%x%iend

    weightsDx=dm_zeros(N,N)		
    weightsDy=dm_zeros(N,N)		
    weightsDz=dm_zeros(N,N)		
    weightsL =dm_zeros(N,N)		

    !ind_i =dm_zeros(N*fdsize,1)		
    !ind_j =dm_zeros(N*fdsize,1)
    
    A=dm_ones(fdsize+1,fdsize+1)
    B=dm_zeros(fdsize+1,1)	
	
	call knnsearch(idx,ierr) 
	!call dm_view(idx,ierr)	
  
    x1=atm%pts%x
    x2=atm%pts%y
    x3=atm%pts%z
   	px=atm%pts%p_u
   	py=atm%pts%p_v
   	pz=atm%pts%p_w
 
    x = x1 .hj. x2 .hj. x3
	!do j=0,N-1
	do j=0,0
 		imat	= dm_trans(dm_getrow(idx,j))
        ximat  	= dm_submatrix(x,imat,(0 .to. 2)) 

!    	ximat1 	= dm_getcol(ximat,0) 
!       ximat2 	= dm_getcol(ximat,1) 
!       ximat3 	= dm_getcol(ximat,2) 
!       xj		= dm_getrow(x,j)
!      	xj1		= dm_getcol(xj,0)
!      	xj2		= dm_getcol(xj,1)
!      	xj3		= dm_getcol(xj,2)
!       dp=xj*ximat1 		

 		rd2= 2*(1-dm_xyt(ximat,ximat))
 		rd2= (rd2>0) .em. rd2
 		rd2v= dm_getcol(rd2,0)
			
		A=rbf_guassian(ep,rd2)
		
		A= A .hj. (dm_ones(A%nrow,1))
		A= A .vj. (dm_ones(1,A%ncol))
    	call dm_setvalue(A,A%nrow-1,A%ncol-1,0,ierr)	
		call dm_view(A,ierr)	
!		p1=	dm_submatrix(px, j, dm_seqs(3,1))


	enddo
 	call dm_destroy(weightsDx,ierr) 
 	call dm_destroy(weightsDy,ierr) 
    call dm_destroy(weightsDz,ierr) 
    call dm_destroy(weightsL,ierr) 
    call dm_destroy(A,ierr) 
    call dm_destroy(B,ierr) 
    call dm_destroy(idx,ierr) 
    call dm_destroy(imat,ierr) 
    call dm_destroy(ximat,ierr) 
!	call dm_destroy(ximat1,ierr) 
!	call dm_destroy(ximat2,ierr) 
!	call dm_destroy(ximat3,ierr) 
    call dm_destroy(x,ierr) 
 	call dm_destroy(x1,ierr) 
 	call dm_destroy(x2,ierr) 
 	call dm_destroy(x3,ierr) 
!	call dm_destroy(xj,ierr) 
 	call dm_destroy(px,ierr) 
 	call dm_destroy(py,ierr) 
 	call dm_destroy(pz,ierr) 
    call dm_destroy(rd2,ierr) 
    call dm_destroy(rd2v,ierr) 
!   !call dm_view(A,ierr)	
end subroutine


