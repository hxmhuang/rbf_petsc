program main 
    use dm 
	use rbf
    implicit none
    type(Matrix)    :: dsites,epoints,rhs,ctrs,DM_data,IM,DM_eval,EM,s,exact
    integer         :: myrank, mysize 
    integer         :: m,n 
    integer         :: meval,neval 
    real(kind=8)    :: ep
    real(kind=8)    :: maxerr,rmserr
    logical         :: debug 
    integer         :: ierr
    debug=.false.

    call dm_init(ierr)
    
    call dm_comm_rank(myrank,ierr)
    call dm_comm_size(mysize,ierr)
    
    call dm_option_int('-m',m,ierr)
    call dm_option_int('-n',n,ierr)
    call dm_option_int('-meval',meval,ierr)
    call dm_option_int('-neval',neval,ierr)
    call dm_option_real('-ep',ep,ierr)
    call dm_option_bool('-debug',debug,ierr)
 	
   	call rbf_createpoints(dsites,m,n,ierr)
   	call rbf_createpoints(epoints,meval,neval,ierr)

   	call rbf_testfunction(dsites,rhs,ierr)
  	call rbf_testfunction(epoints,exact,ierr) 
    ctrs=dsites
    call rbf_distancematrix(dsites,ctrs,DM_data,ierr)
    call rbf_distancematrix(epoints,ctrs,DM_eval,ierr)
    IM=rbf_guassian(ep,DM_data)
    EM=rbf_guassian(ep,DM_eval)
   	s=EM*(IM .inv. rhs)
	maxerr=dm_norm_inf(s-exact)
	rmserr=dm_norm_2(s-exact)/neval 
	

	if(myrank==0) then 
		print *, "==============Input paramenters==========="
		print *, "m=",m,",n=",n,"meval=",meval,"neval=",neval,"ep=",ep,"debug=",debug
	endif 

	if(debug) then
		if(myrank==0) print *, ">distes="
        call dm_view(dsites,ierr)
        if(myrank==0) print *, ">epoints="
        call dm_view(epoints,ierr)
		if(myrank==0) print *, ">rhs="
		call dm_view(rhs,ierr)
		if(myrank==0) print *, ">exact="
		call dm_view(exact,ierr)
		if(myrank==0) print *, ">DM_data="
		call dm_view(DM_data,ierr)
		if(myrank==0) print *, ">DM_eval="
		call dm_view(DM_eval,ierr)
		if(myrank==0) print *, ">IM="
		call dm_view(IM,ierr)
		if(myrank==0) print *, ">EM="
		call dm_view(EM,ierr)
		if(myrank==0) print *, ">s="
		call dm_view(s,ierr)
 	endif

	if(myrank==0) then
		print *, "==============Test norm==================="
    	print *, ">RMS 	  error:",rmserr
    	print *, ">Maximum error:",maxerr
	endif

 	call dm_destroy(dsites,ierr)
 	call dm_destroy(epoints,ierr)
    call dm_destroy(rhs,ierr)
    call dm_destroy(ctrs,ierr)
    call dm_destroy(DM_data,ierr)
    call dm_destroy(DM_eval,ierr)
    call dm_destroy(IM,ierr)
    call dm_destroy(EM,ierr)
    call dm_destroy(s,ierr)
    call dm_destroy(exact,ierr)

    call dm_finalize(ierr)
end program
