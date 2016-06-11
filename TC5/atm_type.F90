module atm_type
    use dm
    implicit none
    type Points
        type(Matrix) :: x
        type(Matrix) :: y
        type(Matrix) :: z
        type(Matrix) :: la
        type(Matrix) :: th
        type(Matrix) :: r
        type(Matrix) :: p_u 
        type(Matrix) :: p_v 
        type(Matrix) :: p_w 
   end type Points

   type AtmType
        type(Points) :: pts
        real(kind=8) :: alpha
        real(kind=8) :: u0 
        real(kind=8) :: a 
        real(kind=8) :: omega 
        real(kind=8) :: g 
        real(kind=8) :: gh0 
        real(kind=8) :: lam_c 
        real(kind=8) :: thm_c 
        real(kind=8) :: mR 
        real(kind=8) :: hm0 
        type(Matrix) :: f 
        type(Matrix) :: ghm 
   end type AtmType

contains

subroutine atm_destroy(atm,ierr)
    implicit none
    type(AtmType),intent(inout) :: atm
    integer,intent(out)         :: ierr
    call dm_destroy(atm%pts%x,ierr)
    call dm_destroy(atm%pts%y,ierr)
    call dm_destroy(atm%pts%z,ierr)
    call dm_destroy(atm%pts%la,ierr)
    call dm_destroy(atm%pts%th,ierr)
    call dm_destroy(atm%pts%r,ierr)
    call dm_destroy(atm%pts%p_u,ierr)
    call dm_destroy(atm%pts%p_v,ierr)
    call dm_destroy(atm%pts%p_w,ierr)
    call dm_destroy(atm%f,ierr)
    call dm_destroy(atm%ghm,ierr)
end subroutine

end module

