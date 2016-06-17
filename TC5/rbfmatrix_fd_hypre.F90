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
    type(Matrix)        		  :: x,ximat
    type(Matrix)        		  :: rd2,rd2v,drbf_rd2v,ep2r2
    type(Matrix)        		  :: weights,weightsDx,weightsDy,weightsDz,weightsL
    type(Matrix)        		  :: ind_i,ind_j 
    type(Matrix)        		  :: hypre_result
    type(Matrix)        		  :: m_0,m_0to2,m_0toN 
	integer 					  :: N,j,myrank
	integer 					  :: ista,iend
    
	N=atm%pts%x%nrow	
    ista=atm%pts%x%ista
    iend=atm%pts%x%iend

    !weightsDx=dm_zeros(N,N)
    !weightsDy=dm_zeros(N,N)
    !weightsDz=dm_zeros(N,N)
    !weightsL =dm_zeros(N,N)
	!ind_i =dm_zeros(N*fdsize,1)		
    !ind_j =dm_zeros(N*fdsize,1)
    
    A=dm_ones(fdsize+1,fdsize+1)
    B=dm_zeros(fdsize+1,1)	
	
	call knnsearch(idx,ierr) 
  
    x = atm%pts%x .hj. atm%pts%y .hj. atm%pts%z 
	call dm_comm_rank(myrank,ierr)
	
	ind_i=dm_zeros(0,1)
	ind_j=dm_zeros(0,1)
	weightsDx=dm_zeros(0,1)
	weightsDy=dm_zeros(0,1)
	weightsDz=dm_zeros(0,1)
	weightsL=dm_zeros(0,1)

	m_0to2=dm_m2n(0,2)
	m_0toN=dm_m2n(0,fdsize-1)
	m_0	  =dm_zeros(1,1)	
 	!do j=0,N-1
 	do j=0,30
!   !do j=0,0
 		if(myrank==0) print *,">Rbfmatrix_df_hypre: current", j, " step"
  		imat	= dm_trans(dm_getrow(idx,j))
  		ximat  	= dm_submatrix(x,imat,m_0to2) 
    	
    	ind_i = ind_i .vj. dm_constants(fdsize,1,j)	
    	ind_j = ind_j .vj. imat 

 		rd2= 2*(1-dm_xyt(ximat,ximat))
  		rd2= (rd2>0) .em. rd2
  		rd2v= dm_getcol(rd2,0)
    	
 		A=rbf_guassian(ep,rd2)
    	drbf_rd2v=drbf_guassian(ep,rd2v)
 		
    	A= A .hj. (dm_ones(A%nrow,1))
 		A= A .vj. (dm_ones(1,A%ncol))
     	call dm_setvalue(A,A%nrow-1,A%ncol-1,0,ierr)	
    	B=(-1.0)*ximat * dm_trans(dm_getrow(atm%pts%p_u,j)) .em. drbf_rd2v
    	B=B .vj. m_0
    	weights= A .inv. B
    	weightsDx= weightsDx .vj. (dm_submatrix(weights, m_0toN, m_0)) 
    	
     	B=(-1.0)*ximat * dm_trans(dm_getrow(atm%pts%p_v,j)) .em. drbf_rd2v
     	B=B .vj. m_0 
    	weights= A .inv. B
    	weightsDy= weightsDy .vj. (dm_submatrix(weights, m_0toN, m_0)) 
    		
    	B=(-1.0)*ximat * dm_trans(dm_getrow(atm%pts%p_w,j)) .em. drbf_rd2v
    	B=B .vj. m_0 
    	weights= A .inv. B
    	weightsDz= weightsDz .vj. (dm_submatrix(weights, m_0toN, m_0)) 
    	
    	!There is a bug to put a MAT_XTYPE_IMPLICIT matrix into a function directly.
    	!hypre_result=hypre(ep**2 * rd2v,dim,order)
    	ep2r2=ep**2 * rd2v
     	hypre_result=hypre(ep2r2,dim,order)
    	
    	B=ep**(2*order) * hypre_result .em. (rbf_guassian(ep,rd2v))
    	B=B .vj. m_0 
    	weights= A .inv. B
    	weightsL= weightsL .vj. (dm_submatrix(weights, m_0toN, m_0)) 

    enddo
	
	DPx=dm_sparse(ind_i,ind_j,weightsDx,N,N)
	DPy=dm_sparse(ind_i,ind_j,weightsDy,N,N)
	DPz=dm_sparse(ind_i,ind_j,weightsDz,N,N)
	L  =dm_sparse(ind_i,ind_j,weightsL,N,N)

	!call dm_view(DPx,ierr)
	!call dm_view(DPy,ierr)
	!call dm_view(DPz,ierr)
	!call dm_view(L,ierr)
    call dm_destroy(A,ierr) 
    call dm_destroy(B,ierr) 
    call dm_destroy(idx,ierr) 
    call dm_destroy(x,ierr) 
	call dm_destroy(ind_i,ierr) 
 	call dm_destroy(ind_j,ierr) 
 	call dm_destroy(weights,ierr) 
 	call dm_destroy(weightsDx,ierr) 
 	call dm_destroy(weightsDy,ierr) 
    call dm_destroy(weightsDz,ierr) 
    call dm_destroy(weightsL,ierr) 
    call dm_destroy(imat,ierr) 
    call dm_destroy(ximat,ierr) 
    call dm_destroy(rd2,ierr) 
    call dm_destroy(rd2v,ierr) 
    call dm_destroy(drbf_rd2v,ierr) 
    call dm_destroy(ep2r2,ierr) 
    call dm_destroy(hypre_result,ierr) 
    call dm_destroy(m_0,ierr) 
    call dm_destroy(m_0to2,ierr) 
    call dm_destroy(m_0toN,ierr) 
end subroutine


