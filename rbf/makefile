CFLAGS	         =  
FFLAGS	         =-Wno-tabs
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = src/ksp/ksp/examples/tutorials/
MANSEC           = KSP
CLEANFILES       = main*.o *.mod 
NP               = 1
OBJ				 = dm_type.o dm_mat.o dm.o rbf.o 
OBJMAIN			 = ${OBJ} main.o 

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main: ${OBJMAIN}  chkopts
	-${FLINKER} -o main ${OBJMAIN}  ${PETSC_KSP_LIB}

small:
	make clean
	make main 
	-@${MPIEXEC} -n 4 ./main -m 3 -n 3 -meval 4 -neval 4 -ep 11.1 -debug -ksp_type bcgs -pc_type bjacobi -sub_ksp_type preonly -sub_pc_type sor

middle:
	make clean
	make main 
	-@${MPIEXEC} -n 8 ./main -m 14 -n 14 -meval 40 -neval 40 -ep 11.1 -log_view -ksp_type bcgs -pc_type bjacobi -sub_ksp_type preonly -sub_pc_type sor

big:
	make clean
	make main 
	-@${MPIEXEC} -n 16 ./main -m 30 -n 30 -meval 50 -neval 50 -ep 11.1 -log_view -ksp_type bcgs -pc_type bjacobi -sub_ksp_type preonly -sub_pc_type sor

huge:
	make clean
	make main 
	-@${MPIEXEC} -n 16 ./main -m 40 -n 40 -meval 50 -neval 50 -ep 11.1 -log_view -ksp_type bcgs -pc_type bjacobi -sub_ksp_type preonly -sub_pc_type sor


mainsjob:
	qsub job_mains.job 

mainmjob:
	qsub job_mainm.job 

mainbjob:
	qsub job_mainb.job 

mainhjob:
	qsub job_mainh.job 




#include ${PETSC_DIR}/lib/petsc/conf/test
