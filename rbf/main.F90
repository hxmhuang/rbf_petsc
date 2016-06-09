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

    ierr=dm_init()
    
    myrank=dm_comm_rank()
    
    mysize=dm_comm_size()
    
    m=dm_option_int('-m')
    n=dm_option_int('-n')
    meval=dm_option_int('-meval')
    neval=dm_option_int('-neval')
    ep=dm_option_real('-ep')
    debug=dm_option_bool('-debug')

    if(myrank==0) then 
       print *, "==============Input paramenters==========="
        print *, "m=",m,",n=",n,"meval=",meval,"neval=",neval,"ep=",ep,"debug=",debug
     endif 
	
    
 	if(myrank==0) print *, "==============Test rbf_createpoints============"
 	call rbf_createpoints(dsites,m,n,ierr)
 	call rbf_createpoints(epoints,meval,neval,ierr)
    if(debug) then
        if(myrank==0) print *, ">distes="
        ierr=dm_view(dsites)
        if(myrank==0) print *, ">epoints="
        ierr=dm_view(epoints)
 	endif


	if(myrank==0) print *, "==============Test rbf_testfucntion============"
 	call rbf_testfunction(dsites,rhs,ierr)
	call rbf_testfunction(epoints,exact,ierr) 
    if(debug) then
         if(myrank==0) print *, ">rhs="
         ierr=dm_view(rhs)
         if(myrank==0) print *, ">exact="
         ierr=dm_view(exact)
 	endif


	if(myrank==0) print *, "==============Test rbf_distancematrix=========="
    ctrs=dsites
    call rbf_distancematrix(dsites,ctrs,DM_data,ierr)
    call rbf_distancematrix(epoints,ctrs,DM_eval,ierr)
    if(debug) then
         if(myrank==0) print *, ">DM_data="
         ierr=dm_view(DM_data)
         if(myrank==0) print *, ">DM_eval="
         ierr=dm_view(DM_eval)
 	endif

	if(myrank==0) print *, "==============Test rbf_guassian==============="
    call rbf_guassian(ep,DM_data,IM,ierr)
    call rbf_guassian(ep,DM_eval,EM,ierr)
    if(debug) then
         if(myrank==0) print *, ">IM="
         ierr=dm_view(IM)
         if(myrank==0) print *, ">EM="
         ierr=dm_view(EM)
 	endif

	if(myrank==0) print *, "==============Test rbf_solve============="
   	s=EM*(IM .inv. rhs)
    if(debug) then
         if(myrank==0) print *, ">s="
         ierr=dm_view(s)
 	endif

	if(myrank==0) print *, "==============Test norm================="
	maxerr=dm_norm_inf(s-exact)
	rmserr=dm_norm_2(s-exact)/neval 
    if(myrank==0) print *, ">RMS 	  error:",rmserr
    if(myrank==0) print *, ">Maximum error:",maxerr

 	ierr=dm_destroy(dsites)
 	ierr=dm_destroy(epoints)
    ierr=dm_destroy(rhs)
    ierr=dm_destroy(ctrs)
    ierr=dm_destroy(DM_data)
    ierr=dm_destroy(DM_eval)
    ierr=dm_destroy(IM)
    ierr=dm_destroy(EM)
    ierr=dm_destroy(s)
    ierr=dm_destroy(exact)

    ierr=dm_finalize()
end program
