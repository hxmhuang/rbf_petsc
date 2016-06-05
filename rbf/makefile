CFLAGS	         =  
FFLAGS	         =-Wno-tabs
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = src/ksp/ksp/examples/tutorials/
MANSEC           = KSP
CLEANFILES       = main*.o *.mod 
NP               = 1
OBJ				 = dm_type.o dm_mat.o dm.o 
OBJMAIN			 = ${OBJ} main.o 

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main: ${OBJMAIN}  chkopts
	-${FLINKER} -o main ${OBJMAIN}  ${PETSC_KSP_LIB}

small:
	make clean
	make main 
	-@${MPIEXEC} -n 4 ./main -m 3 -n 2 -ep 3.1 -debug -ksp_type bcgs -pc_type bjacobi -sub_ksp_type preonly -sub_pc_type sor

middle:
	make clean
	make main 
	-@${MPIEXEC} -n 4 ./main -m 3 -n 2 -ep 3.1 -log_view -ksp_type bcgs -pc_type bjacobi -sub_ksp_type preonly -sub_pc_type sor

big:
	make clean
	make main 
	-@${MPIEXEC} -n 16 ./main -m 300 -n 300 -ep 3.1 -log_view -ksp_type bcgs -pc_type bjacobi -sub_ksp_type preonly -sub_pc_type sor



#include ${PETSC_DIR}/lib/petsc/conf/test
