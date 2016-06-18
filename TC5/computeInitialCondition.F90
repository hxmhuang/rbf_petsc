subroutine computeInitialCondition(atm,H,ierr)
	use dm
    use rbf
    use atm_type
    implicit none
    
    type(AtmType),intent(in)   	:: atm
    type(Matrix),intent(out)    :: H 
    integer,intent(out) 		:: ierr
	
	type(Matrix)				:: uc,gh
	type(Matrix)				:: x,y,z 
	real(kind=8)				:: alpha,a,omega,u0	
	real(kind=8)				:: b1,b2	
	alpha=atm%alpha
	a=atm%a
	omega=atm%omega
	u0=atm%u0
	
	x=atm%pts%x
	y=atm%pts%y
	z=atm%pts%z

	b1 = cos(alpha)
	b2 = sin(alpha)

	gh = (-1.0)*(a*omega*u0+(u0**2)/2.0) * dm_squ(z*b1-x*b2)				
	uc = u0*( (-1.0)*b1*y .hj. (x*b1+z*b2) .hj. (-1.0)*b2*y )
	H= uc .hj. gh		
	
	call dm_destroy(x,ierr)
	call dm_destroy(y,ierr)
	call dm_destroy(z,ierr)
	call dm_destroy(uc,ierr)
	call dm_destroy(gh,ierr)

end subroutine
