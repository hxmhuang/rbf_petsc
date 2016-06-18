subroutine test_case_5_cart_rk4_fd(nfile,ep,fdsize,order,dim,gamma,dt,tend,ierr) 
    use dm
    use rbf
    use atm_type
    implicit none
    
    type(Matrix),intent(in)    	:: nfile
    real(kind=8),intent(in)    	:: ep,tend,gamma
    integer,intent(in)         	:: fdsize,order,dim,dt
    integer,intent(out) 		:: ierr
    type(AtmType)       		:: atm
    type(Matrix)				:: DPx,DPy,DPz,L,H,K 
    type(Matrix)				:: d1,d2,d3,d4 
    type(Matrix)				:: gradghm 
    !real(kind=8)    			:: timeday 
	integer 					:: i,myrank

    call setupT5(nfile,atm,ierr)
	
	call rbfmatrix_fd_hypre(atm,ep,fdsize,order,dim,DPx,DPy,Dpz,L,ierr)	
	L  = gamma * L

	call computeInitialCondition(atm,H,ierr)

	gradghm= (DPx*atm%ghm .hj. DPy*atm%ghm .hj. DPz*atm%ghm)*(1/atm%a)

	call dm_comm_rank(myrank,ierr)
	!do i=1,floor(tend*24*3600/dt)	
	do i=1,0			
 		if(myrank==0) print *,">test_case_5_cart_rk4_rd: the current time is", i*dt/3600.0, "hours"
 		K=H
 		call evalCartRhs_fd(d1,K,DPx,DPy,DPz,L,atm,gradghm,ierr)	
    	K=H+0.5*dt*d1
    	call evalCartRhs_fd(d2,K,DPx,DPy,DPz,L,atm,gradghm,ierr)	
    	K=H+0.5*dt*d2
    	call evalCartRhs_fd(d3,K,DPx,DPy,DPz,L,atm,gradghm,ierr)	
    	K=H+dt*d3
    	call evalCartRhs_fd(d4,K,DPx,DPy,DPz,L,atm,gradghm,ierr)	
    	H=H+1/6.0*dt*(d1+2*d2+2*d3+d4)
    	call dm_destroy(d1,ierr)
    	call dm_destroy(d2,ierr)
    	call dm_destroy(d3,ierr)
    	call dm_destroy(d4,ierr)
	enddo
	
	call dm_destroy(gradghm,ierr)	
	call dm_destroy(DPx,ierr)
    call dm_destroy(DPy,ierr)
    call dm_destroy(DPz,ierr)
    call dm_destroy(L,ierr)
    call dm_destroy(H,ierr)
    call dm_destroy(K,ierr)
    call atm_destroy(atm,ierr)
end subroutine
