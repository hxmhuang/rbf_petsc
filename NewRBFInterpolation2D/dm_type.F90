#include <petsc/finclude/petscmatdef.h>
#include "mat_type.h"

module dm_type 
    implicit none
    type Matrix 
        Mat  x         
        integer xtype
    end type Matrix 

!    type MatrixIm 
!        Mat  x         
!    end type MatrixIm 

end module 


