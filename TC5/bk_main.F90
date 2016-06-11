! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! This program implemented the TC5 case based on distributed matrix. 
! -----------------------------------------------------------------------
program main 
    use dm 
	use rbf
    implicit none
    type(Matrix)    	:: nodes,pts,PX,PY,PZ,F 
    type(Matrix)		:: DM_data,DM_eval,IM,EM
    type(Matrix)		:: u,x,rhx,exact,s,norm 
    integer         	:: myrank, mysize 
    logical         	:: debug 
    integer         	:: ierr	

    real(kind=8)		:: ep,tend,gamma,pi,lam_c,thm_c,mR,hm0
	integer				:: fdsize, order, dims, dt, N	
    real(kind=8)		:: alpha,angle,radius,gh0,u0,g
    character(len=100) 	:: filename 
 
    debug = .false.
    alpha=1.0

    fdsize=31
    tend=15
    order=4
    dims=2
	ep=2.7
    dt=900
    gamma=-2.98e-17
    N=4
    filename="md001.00004"
    
    angle=0         	! Angle of rotation measured from the equator.
    u0=20           	! Speed of rotation in meters/second
    radius=6.37122e6    ! Mean radius of the earth (meters).
    g = 9.80616       	! Gravitational constant (m/s^2).
    gh0 = g*5960      	! Initial condition for the geopotential field (m^2/s^2).

    call dm_init(ierr)
    
    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
 	
    call dm_load(filename,nodes,ierr)
 	call dm_cart2sph(nodes,pts,ierr)	
	call rbf_projection(nodes,PX,PY,PZ,ierr)	   	
	call rbf_coriolis_force(nodes,angle,F,ierr)	   	
	if(myrank==0) then 
		print *, "==============Input paramenters==========="
		print *, "ep=",ep,"debug=",debug
	endif 

    if(debug) then
        if(myrank==0) print *, ">nodes="
        call dm_view(nodes,ierr)
        if(myrank==0) print *, ">pts="
        call dm_view(pts,ierr)
        if(myrank==0) print *, ">PX="
        call dm_view(PX,ierr)
        if(myrank==0) print *, ">PY="
        call dm_view(PY,ierr)
        if(myrank==0) print *, ">PZ="
        call dm_view(PZ,ierr)
        if(myrank==0) print *, ">F="
        call dm_view(F,ierr)
 	endif

	call dm_destroy(nodes,ierr)	
 	call dm_destroy(pts,ierr)	
 	call dm_destroy(PX,ierr)	
 	call dm_destroy(PY,ierr)	
 	call dm_destroy(PZ,ierr)	
 	call dm_destroy(F,ierr)	

    call dm_finalize(ierr)
end program
