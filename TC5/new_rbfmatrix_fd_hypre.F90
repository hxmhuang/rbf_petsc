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
    type(Matrix)        		  :: m_0,m_tiny,m_0to2,m_0toN,m_0toFD 
    type(Matrix)        		  :: Mask,M1,M2 
    type(Matrix)        		  :: dp1,dp2,dp3
    type(Matrix)        		  :: K0,K1,K2,K3,K4,K5,K6
    type(Matrix)        		  :: ExtM,ExtV
    type(Matrix)        		  :: RightV1,RightV2,RightV3,RightV4 
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
	m_0 = dm_zeros(1,1) 
	m_tiny = 1e-20+m_0 
	!m_0 = dm_ones(1,1) 

	m_0toN = dm_m2n(0,N-1) 
	m_0toFD= dm_m2n(0,fdsize-1) 
	! If the load is not balance in each process, this loop will block!!! (iend-ista)
	! TODO:  We need to solve this issue later.
	
	do j=ista,iend-1
	!do j=ista,ista
		K0=dm_getcol(rd2v,j)
		RightV1= dm_getcol(dp1,j) .em. K0	
		RightV2= dm_getcol(dp2,j) .em. K0	
		RightV3= dm_getcol(dp3,j) .em. K0	
	
		M1=dm_getrow(Mask,j)
		M2=dm_trans(M1)
		
		K1=dm_rep(M2,1,N)
		K2=dm_rep(M1,N,1)
		K3=K1 .em. K2
		call dm_setdiag(K3,1,ierr)	
		K4= K3 .em. rbf_rd2
		ExtM= (K4 .hj. M2) .vj. (M1 .hj. m_tiny)
	
		ExtV= RightV1 .vj. m_0 
		weights= ExtM .inv. ExtV
		DPx=DPx .hj. weights
	
		ExtV= RightV2 .vj. m_0 
		weights= ExtM .inv. ExtV
		DPy=DPy .hj. weights
	
		ExtV= RightV3 .vj. m_0 
		weights= ExtM .inv. ExtV
		DPz=DPz .hj. weights
	
		K5= dm_getcol(rd2,j) .em. M2
		ep2r2=ep**2 * K5 
    	hypre_result=hypre(ep2r2,dim,order)
	   	RightV4= (ep**(2*order) * hypre_result) .em. (rbf_guassian(ep,K5)) .em. M2
		ExtV= RightV4 .vj. m_0
		weights= ExtM .inv. ExtV
		L=L .hj. weights
 
	enddo
	print *, DPx%nrow,DPx%ncol
	DPx= dm_trans( dm_getsub(DPx,m_0toN,m_0toFD) )	
!	DPy= dm_trans( dm_getsub(DPy,m_0toN,m_0toFD) )	
!	DPz= dm_trans( dm_getsub(DPz,m_0toN,m_0toFD) )	
!	L  = dm_trans( dm_getsub(L,  m_0toN,m_0toFD) )	
 	
	print *, ">DPx="
	call dm_view(DPx, ierr)
!	print *, ">DPy="
!	call dm_view(DPy, ierr)
!	print *, ">DPz="
!	call dm_view(DPz, ierr)
!	print *, ">L="
!	call dm_view(L, ierr)
	


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
    call dm_destroy(m_tiny,ierr) 
    call dm_destroy(m_0toN,ierr) 
    call dm_destroy(m_0toFD,ierr) 
    
	call dm_destroy(K0,ierr) 
	call dm_destroy(K1,ierr) 
	call dm_destroy(K2,ierr) 
	call dm_destroy(K3,ierr) 
	call dm_destroy(K4,ierr) 
	call dm_destroy(K5,ierr) 
	call dm_destroy(ExtM,ierr) 
	call dm_destroy(ExtV,ierr) 
	call dm_destroy(weights,ierr) 
    
	call dm_destroy(RightV1,ierr) 
	call dm_destroy(RightV2,ierr) 
	call dm_destroy(RightV3,ierr) 
	call dm_destroy(RightV4,ierr) 
	
	call dm_destroy(ep2r2,ierr) 
	call dm_destroy(hypre_result,ierr) 

	call dm_destroy(M1,ierr) 
	call dm_destroy(M2,ierr) 
end subroutine


