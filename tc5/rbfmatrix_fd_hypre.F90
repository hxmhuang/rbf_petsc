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
    type(Matrix)        		  :: A,B,idx
    type(Matrix)        		  :: dp,x,ximat
    type(Matrix)        		  :: rd2,rd2v,drbf_rd2v,ep2r2
    type(Matrix)        		  :: weights
    type(Matrix)        		  :: hypre_result
    type(Matrix)        		  :: m_0
	integer 					  :: N,i,j,myrank
	integer 					  :: ista,iend
   	real(kind=8)				  :: values(3),rmat(fdsize)
	integer						  :: imat(fdsize)
   	real(kind=8)				  :: warray(fdsize)
	
	N=atm%pts%x%nrow	
    ista=atm%pts%x%ista
    iend=atm%pts%x%iend

    DPx=dm_zeros(N,N)
    DPy=dm_zeros(N,N)
    DPz=dm_zeros(N,N)
    L  =dm_zeros(N,N)
    
    A=dm_ones(fdsize+1,fdsize+1,.false.)
    B=dm_zeros(fdsize+1,1,.false.)	
	
	call knnsearch(idx,x,ierr) 
    !x = atm%pts%x .hj. atm%pts%y .hj. atm%pts%z 
	call dm_comm_rank(myrank,ierr)
	
	m_0	  =dm_zeros(1,1,.false.)	
	dp	  =dm_zeros(3,1,.false.)	
    !do j=0,0
    do j=ista,iend-1
		if(myrank==0) print *, ">Rbfmatrix_fd_hypre: current step is",j 
    	call dm_getvalues(idx,(/j/),(/(i,i=0,fdsize-1)/),rmat,ierr)
		imat=int(rmat)

		ximat  	= dm_getsub(x,imat,(/0,1,2/)) 
    	
 		rd2= 2*(1-dm_xyt(ximat,ximat))
   		rd2= (rd2>0) .em. rd2
  		rd2v= dm_getcol(rd2,0)
 		A=rbf_guassian(ep,rd2)
    	drbf_rd2v=drbf_guassian(ep,rd2v)
 		
    	A= A .hj. (dm_ones(A%nrow,1,A%isGlobal))
 		A= A .vj. (dm_ones(1,A%ncol,A%isGlobal))
     	call dm_setvalue(A,A%nrow-1,A%ncol-1,0,ierr)	
    	
    	!compute DPx		
    	!call dm_view(atm%pts%p_u,ierr)
		call dm_getvalues(atm%pts%p_u,(/j/),(/0,1,2/),values,ierr)		
		call dm_setvalues(dp,(/0,1,2/),(/0/),values,ierr)	
    	B= -ximat * dp .em. drbf_rd2v
    	B=B .vj. m_0
    	weights= A .inv. B
    	call dm_getvalues(weights,(/(i,i=0,fdsize-1)/),(/0/),warray,ierr)		
    	call dm_setvalues(DPx,(/j/),imat,warray,ierr)
    		
    	!compute DPy		
    	call dm_getvalues(atm%pts%p_v,(/j/),(/0,1,2/),values,ierr)		
 		call dm_setvalues(dp,(/0,1,2/),(/0/),values,ierr)	
    	B= -ximat * dp .em. drbf_rd2v
    	B=B .vj. m_0
    	weights= A .inv. B
    	call dm_getvalues(weights,(/(i,i=0,fdsize-1)/),(/0/),warray,ierr)		
    	call dm_setvalues(DPy,(/j/),imat,warray,ierr)
 
    	!compute DPz		
    	call dm_getvalues(atm%pts%p_w,(/j/),(/0,1,2/),values,ierr)		
 		call dm_setvalues(dp,(/0,1,2/),(/0/),values,ierr)	
    	B= -ximat * dp .em. drbf_rd2v
    	B=B .vj. m_0
    	weights= A .inv. B
    	call dm_getvalues(weights,(/(i,i=0,fdsize-1)/),(/0/),warray,ierr)		
    	call dm_setvalues(DPz,(/j/),imat,warray,ierr)

    	!There is a bug to put a MAT_XTYPE_IMPLICIT matrix into a function directly.
		!So we have to use the ep2r2 variable 
    	!hypre_result=hypre(ep**2 * rd2v,dim,order)
    	ep2r2=ep**2 * rd2v
     	hypre_result=hypre(ep2r2,dim,order)
    	
    	B=ep**(2*order) * hypre_result .em. (rbf_guassian(ep,rd2v))
    	B=B .vj. m_0 
    	weights= A .inv. B
    	call dm_getvalues(weights,(/(i,i=0,fdsize-1)/),(/0/),warray,ierr)		
    	call dm_setvalues(L,(/j/),imat,warray,ierr)
    enddo
    
 	!call dm_view(DPx,ierr)
 	!call dm_view(DPy,ierr)
 	!call dm_view(DPz,ierr)
 	!call dm_view(L,ierr)
    call dm_destroy(A,ierr) 
    call dm_destroy(B,ierr) 
    call dm_destroy(idx,ierr) 
    call dm_destroy(x,ierr) 
    call dm_destroy(dp,ierr) 
 	call dm_destroy(weights,ierr) 
    call dm_destroy(ximat,ierr) 
    call dm_destroy(rd2,ierr) 
    call dm_destroy(rd2v,ierr) 
    call dm_destroy(drbf_rd2v,ierr) 
    call dm_destroy(ep2r2,ierr) 
    call dm_destroy(hypre_result,ierr) 
    call dm_destroy(m_0,ierr) 
end subroutine


