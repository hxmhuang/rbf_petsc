subroutine evalCartRhs_fd(F,H,DPx,DPy,DPz,L,atm,gradghm,ierr) 
    use dm
    use rbf
    use atm_type
    implicit none
 
    type(Matrix),	intent(in)    	:: H,DPx,DPy,DPz,L,gradghm 
    type(AtmType),	intent(in)    	:: atm
    type(Matrix),	intent(out)	    :: F
	integer,		intent(out)		:: ierr
	
    type(Matrix)					:: x,y,z
    type(Matrix)					:: Tx,Ty,Tz 
    type(Matrix)					:: H0,H1,H2,H3 
    type(Matrix)					:: Tx0,Tx1,Tx2,Tx3
    type(Matrix)					:: Ty0,Ty1,Ty2,Ty3
    type(Matrix)					:: Tz0,Tz1,Tz2,Tz3
    type(Matrix)					:: F0,F1,F2,F3
    type(Matrix)					:: p,q,s
	
	real(kind=8)					:: g,a,gh0 
	x = atm%pts%x  
	y = atm%pts%y  
	z = atm%pts%z 
	f = atm%f 
	g = atm%g 
	a = atm%a  
	gh0 = atm%gh0
	
	Tx=DPx*H*(1/a)
	Ty=DPy*H*(1/a)
	Tz=DPz*H*(1/a)
	
	H0=dm_getcol(H,0)	
	H1=dm_getcol(H,1)	
	H2=dm_getcol(H,2)	
	H3=dm_getcol(H,3)	
	
	Tx0=dm_getcol(Tx,0)	
	Tx1=dm_getcol(Tx,1)	
	Tx2=dm_getcol(Tx,2)	
	Tx3=dm_getcol(Tx,3)	
	
	Ty0=dm_getcol(Ty,0)	
	Ty1=dm_getcol(Ty,1)	
	Ty2=dm_getcol(Ty,2)	
	Ty3=dm_getcol(Ty,3)	
	
	Tz0=dm_getcol(Tz,0)	
	Tz1=dm_getcol(Tz,1)	
	Tz2=dm_getcol(Tz,2)	
	Tz3=dm_getcol(Tz,3)	

	p = (-1.0)*( (H0 .em. Tx0) + (H1 .em. Ty0) + (H2 .em. Tz0) + (f .em. ((y .em. H2) - (z .em. H1))) + Tx3 )

	q = (-1.0)*( (H0 .em. Tx1) + (H1 .em. Ty1) + (H2 .em. Tz1) + (f .em. ((z .em. H0) - (x .em. H2))) + Ty3)  
	
	s = (-1.0)*( (H0 .em. Tx2) + (H1 .em. Ty2) + (H2 .em. Tz2) + (f .em. ((x .em. H1) - (y .em. H0))) + Tz3) 
	!print *, "=======p========"
	!call dm_view(p,ierr)
	!F=dm_zeros(H%nrow,H%ncol)
	F0= (dm_getcol(atm%pts%p_u,0) .em. p) +(dm_getcol(atm%pts%p_u,1) .em. q) + (dm_getcol(atm%pts%p_u,2) .em. s)
		
	F1= (dm_getcol(atm%pts%p_v,0) .em. p) +(dm_getcol(atm%pts%p_v,1) .em. q) + (dm_getcol(atm%pts%p_v,2) .em. s)

	F2= (dm_getcol(atm%pts%p_w,0) .em. p) +(dm_getcol(atm%pts%p_w,1) .em. q) + (dm_getcol(atm%pts%p_w,2) .em. s)
	
	F3= (-1.0)*((H0 .em. (Tx3-dm_getcol(gradghm,0))) &
			   +(H1 .em. (Ty3-dm_getcol(gradghm,1))) &
			   +(H2 .em. (Tz3-dm_getcol(gradghm,2))) &
			   +(( H3 + gh0-atm%ghm) .em. (Tx0+Ty1+Tz2)))
	
	F= F0 .hj. F1 .hj. F2 .hj. F3
	F=F + L*H		
	
	call dm_destroy(x,ierr)
	call dm_destroy(y,ierr)
	call dm_destroy(z,ierr)

    call dm_destroy(Tx0,ierr)
    call dm_destroy(Tx1,ierr)
    call dm_destroy(Tx2,ierr)
    call dm_destroy(Tx3,ierr)

    call dm_destroy(Ty0,ierr)
    call dm_destroy(Ty1,ierr)
    call dm_destroy(Ty2,ierr)
    call dm_destroy(Ty3,ierr)

    call dm_destroy(Tz0,ierr)
    call dm_destroy(Tz1,ierr)
    call dm_destroy(Tz2,ierr)
    call dm_destroy(Tz3,ierr)

    call dm_destroy(H0,ierr)
    call dm_destroy(H1,ierr)
    call dm_destroy(H2,ierr)
    call dm_destroy(H3,ierr)

    call dm_destroy(F0,ierr)
    call dm_destroy(F1,ierr)
    call dm_destroy(F2,ierr)
    call dm_destroy(F3,ierr)
    
    call dm_destroy(p,ierr)
    call dm_destroy(q,ierr)
    call dm_destroy(s,ierr)

	call dm_destroy(Tx,ierr)
	call dm_destroy(Ty,ierr)
	call dm_destroy(Tz,ierr)
	
end subroutine	
