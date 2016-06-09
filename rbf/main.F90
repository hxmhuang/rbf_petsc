program main 
    use dm 
	use rbf
    implicit none
    type(Matrix)    :: A,B,C,D,E,F,G,H 
    integer         :: myrank, mysize 
    integer         :: m,n 
    integer         :: meval,neval 
    real(kind=8)    :: ep,alpha
    logical         :: debug 
    integer         :: ierr
    character(len=50):: filename
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
 	call rbf_createpoints(A,m,n,ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
 	endif
    ierr=dm_destroy(A)


	if(myrank==0) print *, "==============Test rbf_testfucntionD============"
 	call rbf_createpoints(A,m,n,ierr)
 	call rbf_testfunction(A,B,ierr)
    if(debug) then
        if(myrank==0) print *, ">A="
        ierr=dm_view(A)
         if(myrank==0) print *, ">B="
         ierr=dm_view(B)
 	endif
    ierr=dm_destroy(A)
     ierr=dm_destroy(B)





    call PetscFinalize(ierr)
end program
