subroutine new_rbfmatrix_fd_hypre(atm,ep,fdsize,order,dim,DPx,DPy,DPz,L,ierr)
    use dm
    use rbf
    use atm_type
    implicit none
	
    type(AtmType),	intent(inout) :: atm
    real(kind=8),	intent(in)    :: ep
    integer,		intent(in)    :: fdsize,order,dim
    type(Matrix),	intent(out)   :: DPx,DPy,DPz,L 
    integer,intent(out)           :: ierr
    type(Matrix)        		  :: idx
    type(Matrix)        		  :: x,xt,ximat
    type(Matrix)        		  :: rd2,rbf_rd2,rd2v,drbf_rd2v,ep2r2
    type(Matrix)        		  :: weights,weightsDx,weightsDy,weightsDz,weightsL
    type(Matrix)        		  :: hypre_result
    type(Matrix)        		  :: m_0,m_0to2,m_0toN 
    type(Matrix)        		  :: Mask,M1,M2 
    type(Matrix)        		  :: dp1,dp2,dp3
    type(Matrix)        		  :: K0,K1,K2,K3,K4,K5,K6
    type(Matrix)        		  :: ExtM,ExtV
    type(Matrix)        		  :: RightV1,RightV2,RightV3 
	integer 					  :: N,i,j,myrank
	integer 					  :: ista,iend
	real(kind=8),allocatable	  :: v(:),imat(:)   
	integer,allocatable 		  :: idxm(:),idxn(:) 

    x = atm%pts%x .hj. atm%pts%y .hj. atm%pts%z
	xt= dm_trans(x) 
	N=x%nrow	
    ista=x%ista
    iend=x%iend
    call dm_comm_rank(myrank,ierr)

    call knnsearch(idx,ierr) 

	Mask=dm_zeros(N,N)
	DPx =dm_zeros(N+1,0)
	DPy =dm_zeros(N+1,0)
	DPz =dm_zeros(N+1,0)
	L   =dm_zeros(N+1,0)

	allocate(v(fdsize),idxm(1),idxn(fdsize),imat(fdsize))
	do j=ista,iend-1
		idxm(1)=j
		do i=1,fdsize
			idxn(i)=i-1
		enddo
		call dm_getvalues(idx,1,idxm,fdsize,idxn,imat,ierr)
		v=1.0
		call dm_setvalues(Mask,1,idxm,fdsize,int(imat),v,ierr)		
	enddo
	deallocate(v,idxm,idxn,imat)

	!print *, ">New Mask="
	!call dm_view(Mask, ierr)
	
	dp1= dm_trans( ((-1.0) * atm%pts%p_u * xt) .em.  Mask)
	dp2= dm_trans( ((-1.0) * atm%pts%p_v * xt) .em.  Mask)
	dp3= dm_trans( ((-1.0) * atm%pts%p_w * xt) .em.  Mask)

	rd2= 2*(1-dm_xyt(x,x))
	rd2= (rd2>0) .em. rd2
	rbf_rd2=rbf_guassian(ep,rd2)	
	
	!rd2v = drbf(ep,rd2)
	rd2v= (-2.0)* ep**2 * rbf_rd2
	! Add a small value to keep a non-zero diagnoal, otherwize the solver will report an error. 
	! TODO:  We need to solve this issue later.
	m_0 = 1e-30+dm_zeros(1,1) 
	!m_0 = dm_ones(1,1) 

	m_0toN = dm_m2n(0,N-1) 
	! If the load is not balance in each process, this loop will block!!! (iend-ista)
	! TODO:  We need to solve this issue later.
	
	!do j=ista,iend-1
	do j=ista,ista
		K0=dm_getcol(rd2v,j)
		RightV1= dm_getcol(dp1,j) .em. K0	
		RightV2= dm_getcol(dp2,j) .em. K0	
		RightV3= dm_getcol(dp3,j) .em. K0	
	
		M1=dm_getrow(Mask,j)
		M2=dm_trans(M1)
		
		K1=dm_rep(M2,1,N)
		K2=dm_rep(M1,N,1)
		K3=K1 .em. K2
		call dm_diag_set(K3,1,ierr)	
		K4= K3 .em. rbf_rd2
		ExtM= (K4 .hj. M2) .vj. (M1 .hj. m_0)
	
		ExtV= RightV1 .vj. m_0 
		weights= ExtM .inv. ExtV
		DPx=DPx .hj. weights
		print *, ">DPx="
		call dm_view(DPx, ierr)
	
		ExtV= RightV2 .vj. m_0 
		weights= ExtM .inv. ExtV
		DPy=DPy .hj. weights
		print *, ">DPy="
		call dm_view(DPy, ierr)
	
		ExtV= RightV3 .vj. m_0 
		weights= ExtM .inv. ExtV
		DPz=DPz .hj. weights
		print *, ">DPz="
		call dm_view(DPz, ierr)
			
	enddo
 	
	call dm_destroy(x,ierr)
	call dm_destroy(xt,ierr)
	call dm_destroy(Mask,ierr)
	call dm_destroy(dp1,ierr)
	call dm_destroy(dp2,ierr)
	call dm_destroy(dp3,ierr)
	call dm_destroy(idx,ierr)
    call dm_destroy(rd2,ierr) 
    call dm_destroy(rbf_rd2,ierr) 
    call dm_destroy(rd2v,ierr) 
    call dm_destroy(m_0,ierr) 
    call dm_destroy(m_0toN,ierr) 
    
	call dm_destroy(K0,ierr) 
	call dm_destroy(K1,ierr) 
	call dm_destroy(K2,ierr) 
	call dm_destroy(K3,ierr) 
	call dm_destroy(K4,ierr) 
	call dm_destroy(ExtM,ierr) 
	call dm_destroy(ExtV,ierr) 
	call dm_destroy(weights,ierr) 
    
	call dm_destroy(RightV1,ierr) 
	call dm_destroy(RightV2,ierr) 
	call dm_destroy(RightV3,ierr) 
	
	call dm_destroy(M1,ierr) 
	call dm_destroy(M2,ierr) 
!   !weightsDx=dm_zeros(N,N)
!   !weightsDy=dm_zeros(N,N)
!   !weightsDz=dm_zeros(N,N)
!   !weightsL =dm_zeros(N,N)
!   !ind_i =dm_zeros(N*fdsize,1)		
!   !ind_j =dm_zeros(N*fdsize,1)
!   
!   A=dm_ones(fdsize+1,fdsize+1)
!   B=dm_zeros(fdsize+1,1)	
!   
!   call knnsearch(idx,ierr) 
! 
!   x = atm%pts%x .hj. atm%pts%y .hj. atm%pts%z 
!   call dm_comm_rank(myrank,ierr)
!   
!   ind_i=dm_zeros(0,1)
!   ind_j=dm_zeros(0,1)
!   weightsDx=dm_zeros(0,1)
!   weightsDy=dm_zeros(0,1)
!   weightsDz=dm_zeros(0,1)
!   weightsL=dm_zeros(0,1)

!   m_0to2=dm_m2n(0,2)
!   m_0toN=dm_m2n(0,fdsize-1)
!   m_0	  =dm_zeros(1,1)	
!	do j=0,N-1
!	!do j=0,9
!   !do j=0,0
!		if(myrank==0) print *,">Rbfmatrix_df_hypre: current", j, " step"
! 		imat	= dm_trans(dm_getrow(idx,j))
! 		ximat  	= dm_submatrix(x,imat,m_0to2) 
!   	
!   	ind_i = ind_i .vj. dm_constants(fdsize,1,j)	
!   	ind_j = ind_j .vj. imat 

!		rd2= 2*(1-dm_xyt(ximat,ximat))
! 		rd2= (rd2>0) .em. rd2
! 		rd2v= dm_getcol(rd2,0)
!   	
!		A=rbf_guassian(ep,rd2)
!   	drbf_rd2v=drbf_guassian(ep,rd2v)
!		
!   	A= A .hj. (dm_ones(A%nrow,1))
!		A= A .vj. (dm_ones(1,A%ncol))
!    	call dm_setvalue(A,A%nrow-1,A%ncol-1,0,ierr)	
!   	B=(-1.0)*ximat * dm_trans(dm_getrow(atm%pts%p_u,j)) .em. drbf_rd2v
!   	B=B .vj. m_0
!   	weights= A .inv. B
!   	weightsDx= weightsDx .vj. (dm_submatrix(weights, m_0toN, m_0)) 
!   	
!    	B=(-1.0)*ximat * dm_trans(dm_getrow(atm%pts%p_v,j)) .em. drbf_rd2v
!    	B=B .vj. m_0 
!   	weights= A .inv. B
!   	weightsDy= weightsDy .vj. (dm_submatrix(weights, m_0toN, m_0)) 
!   		
!   	B=(-1.0)*ximat * dm_trans(dm_getrow(atm%pts%p_w,j)) .em. drbf_rd2v
!   	B=B .vj. m_0 
!   	weights= A .inv. B
!   	weightsDz= weightsDz .vj. (dm_submatrix(weights, m_0toN, m_0)) 
!   	
!   	!There is a bug to put a MAT_XTYPE_IMPLICIT matrix into a function directly.
!   	!hypre_result=hypre(ep**2 * rd2v,dim,order)
!   	ep2r2=ep**2 * rd2v
!    	hypre_result=hypre(ep2r2,dim,order)
!   	
!   	B=ep**(2*order) * hypre_result .em. (rbf_guassian(ep,rd2v))
!   	B=B .vj. m_0 
!   	weights= A .inv. B
!   	weightsL= weightsL .vj. (dm_submatrix(weights, m_0toN, m_0)) 
!   	
!   enddo
!   
!   DPx=dm_sparse(ind_i,ind_j,weightsDx,N,N)
!   DPy=dm_sparse(ind_i,ind_j,weightsDy,N,N)
!   DPz=dm_sparse(ind_i,ind_j,weightsDz,N,N)
!   L  =dm_sparse(ind_i,ind_j,weightsL,N,N)
!   !call dm_view(DPx,ierr)
!   !call dm_view(DPy,ierr)
!   !call dm_view(DPz,ierr)
!   !call dm_view(L,ierr)
!   call dm_destroy(idx,ierr) 
!   call dm_destroy(x,ierr) 
!   call dm_destroy(ind_i,ierr) 
!	call dm_destroy(ind_j,ierr) 
!	call dm_destroy(weights,ierr) 
!	call dm_destroy(weightsDx,ierr) 
!	call dm_destroy(weightsDy,ierr) 
!   call dm_destroy(weightsDz,ierr) 
!   call dm_destroy(weightsL,ierr) 
!   call dm_destroy(imat,ierr) 
!   call dm_destroy(ximat,ierr) 
!   call dm_destroy(drbf_rd2v,ierr) 
!   call dm_destroy(ep2r2,ierr) 
!   call dm_destroy(hypre_result,ierr) 
!   call dm_destroy(m_0,ierr) 
!   call dm_destroy(m_0to2,ierr) 
!   call dm_destroy(m_0toN,ierr) 

end subroutine


