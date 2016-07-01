subroutine setupT5(nfile,atm,ierr)
    use dm
    use rbf
    use atm_type
    implicit none

    type(Matrix),	intent(in)    :: nfile
    type(AtmType),	intent(inout) :: atm
    integer,		intent(out)   :: ierr
    type(Matrix)        		  :: nodes,work 
    type(Matrix)	    		  ::  x,y,z,x2,xy,y2,xz,z2,yz,r2,id
    real(kind=8)         		  :: pi

    nodes=nfile

 	x=dm_getcol(nodes,0)
 	y=dm_getcol(nodes,1)
 	z=dm_getcol(nodes,2)

    atm%pts%x=x
    atm%pts%y=y
    atm%pts%z=z
    
    call dm_cart2sph(nodes,work,ierr)

    atm%pts%la=dm_getcol(work,0)
    atm%pts%th=dm_getcol(work,1)
    atm%pts%r =dm_getcol(work,2)

 	x2= 1-dm_squ(x)
    y2= 1-dm_squ(y)
    z2= 1-dm_squ(z)
    xy=	(-1.0)*x .em. y
    xz= (-1.0)*x .em. z
    yz= (-1.0)*y .em. z 
    
    atm%pts%p_u= x2 .hj. xy .hj. xz	
    atm%pts%p_v= xy	.hj. y2 .hj. yz	
    atm%pts%p_w= xz .hj. yz .hj. z2	
        
    atm%alpha   = 0.0
    atm%u0      = 20.0
    atm%a       = 6.37122e6
    atm%omega   = 7.292e-5
    atm%g       = 9.80616
    atm%gh0     = atm%g*5960.0
    atm%f       = 2.0*atm%omega*(x*sin(atm%alpha)+z*cos(atm%alpha))

    pi        = 4.0*atan(real(1.0,kind=8))
    atm%lam_c = -pi/2.0
    atm%thm_c = pi/6.0
    atm%mR    = pi/9.0
    atm%hm0   = 2000.0

     
    atm%ghm=dm_zeros(nodes%nrow,1)
    r2=dm_squ(atm%pts%la-atm%lam_c)+dm_squ(atm%pts%th-atm%thm_c)

    id= r2 < (atm%mR**2)
    atm%ghm=(atm%g*atm%hm0*(1-dm_sqrt(r2)*(1.0/atm%mR))) .em. id
    
!   call dm_view(atm%pts%x,ierr)
!   call dm_view(atm%pts%y,ierr)
!   call dm_view(atm%pts%z,ierr)
!   call dm_view(atm%pts%la,ierr)
!   call dm_view(atm%pts%th,ierr)
!   call dm_view(atm%pts%r,ierr)
!   call dm_view(atm%pts%p_u,ierr)
!   call dm_view(atm%pts%p_v,ierr)
!   call dm_view(atm%pts%p_w,ierr)
!   call dm_view(atm%f,ierr)
!   call dm_view(r2,ierr)
!   call dm_view(id,ierr)
!   call dm_view(atm%ghm,ierr)
    
    call dm_destroy(nodes,ierr)
    call dm_destroy(work,ierr)
	call dm_destroy(x,ierr)	
	call dm_destroy(y,ierr)	
	call dm_destroy(z,ierr)	
    call dm_destroy(x2,ierr)
    call dm_destroy(y2,ierr)
    call dm_destroy(z2,ierr)
	call dm_destroy(xy,ierr)
	call dm_destroy(xz,ierr)
	call dm_destroy(yz,ierr)
	call dm_destroy(r2,ierr)	
	call dm_destroy(id,ierr)	

end subroutine
