CFLAGS	         =  
FFLAGS	         =-Wno-tabs
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = src/ksp/ksp/examples/tutorials/
MANSEC           = KSP
CLEANFILES       = main *.o *.mod 
NP               = 1
OBJ				 = matrixalgebra.o particle.o main.o 

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main: ${OBJ}  chkopts
	-${FLINKER} -o main ${OBJ}  ${PETSC_KSP_LIB}
#	${RM} *.mod *.o

#----------------------------------------------------------------------------
#runex1:
#	-@${MPIEXEC} -n 1 ./ex1 -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always > ex1_1.tmp 2>&1;	  
#runex1_2:
#	-@${MPIEXEC} -n 1 ./ex1 -pc_type sor -pc_sor_symmetric -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always >\
#	   ex1_2.tmp 2>&1;   
#runex1_3:
#	-@${MPIEXEC} -n 1 ./ex1 -pc_type eisenstat -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always >\
#	   ex1_3.tmp 2>&1;   
#NP = 1
#M  = 4
#N  = 5
#MDOMAINS = 2
#NDOMAINS = 1
#OVERLAP=1
#runex8:
#	-@${MPIEXEC} -n ${NP} ./ex8 -m $M -n $N -user_set_subdomains -Mdomains ${MDOMAINS} -Ndomains ${NDOMAINS} -overlap ${OVERLAP} -print_error ${ARGS}

#clean:

#runex8_1:
#	-@${MPIEXEC} -n 1 ./ex8 -print_error -ksp_view 

run:
	make clean
	make main
	-@${MPIEXEC} -n 3 ./main

#include ${PETSC_DIR}/lib/petsc/conf/test
