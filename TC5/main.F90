! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! This program implemented the TC5 case based on distributed matrix. 
! -----------------------------------------------------------------------
program main 
    use dm 
	use rbf
    implicit none
    
    integer             :: fd,order,dim,dt
    real(kind=8)		:: ep,tend,gamma
    character(len=100) 	:: filename 
    logical         	:: debug
    type(Matrix)        :: fname
    integer             :: mysize,myrank 
    integer             :: ierr
    debug = .false.

    fd=3 
    tend=15
    order=4
    dim=2
	ep=2.7
    dt=900
    gamma=-2.98e-17
    filename="md002.00009"
    !filename="md059.03600"
    
    call dm_init(ierr)
    
    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
 	
    if(myrank==0) then 
		print *, "==============Input paramenters==========="
		print *, "ep=",ep,"debug=",debug
	endif 
    
    call dm_load(filename,fname,ierr)
	
    call test_case_5_cart_rk4_fd(fname,ep,fd,order,dim,gamma,dt,tend,ierr) 


	call dm_destroy(fname,ierr)	

    call dm_finalize(ierr)
end program


