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
    type(Matrix)				:: DPx,DPy,DPz,L 

    call setupT5(nfile,atm,ierr)
	
	call rbfmatrix_fd_hypre(atm,ep,fdsize,order,dim,DPx,DPy,Dpz,L,ierr)	
		
    call dm_destroy(DPx,ierr)
    call dm_destroy(DPy,ierr)
    call dm_destroy(DPz,ierr)
    call dm_destroy(L,ierr)
    call atm_destroy(atm,ierr)
end subroutine
