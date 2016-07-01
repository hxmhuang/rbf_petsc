CFLAGS	         =  
FFLAGS	         =-Wno-tabs
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = src/ksp/ksp/examples/tutorials/
MANSEC           = KSP
CLEANFILES       = main *.o *.mod 
NP               = 1
OBJ				 = dm_type.o dm_mat.o dm.o rbf.o atm_type.o setupT5.o rbfmatrix_fd_hypre.o computeInitialCondition.o evalCartRhs_fd.o test_case_5_cart_rk4_fd.o
OBJMAIN			 = ${OBJ} main.o 

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main: ${OBJMAIN}  chkopts
	-${FLINKER} -o main ${OBJMAIN}  ${PETSC_KSP_LIB}

tiny:
	make clean
	make main 
	-@${MPIEXEC} -n 1 ./main -ep 2.7 -debug -ksp_type gmres -pc_type asm -sub_ksp_type preonly -sub_pc_type ilu 

small:
	make clean
	make main 
	-@${MPIEXEC} -n 3 ./main -ep 2.7 -debug -ksp_type gmres -pc_type asm -sub_ksp_type preonly -sub_pc_type ilu 

middle:
	make clean
	make main 
	-@${MPIEXEC} -n 4 ./main -ep 2.7 -log_view -ksp_type gmres -pc_type asm -sub_ksp_type preonly -sub_pc_type ilu 

big:
	make clean
	make main 
	-@${MPIEXEC} -n 12 ./main -ep 2.7 -log_view -ksp_type gmres -pc_type asm -sub_ksp_type preonly -sub_pc_type ilu 

mainsjob:
	qsub job_mains.job 

mainmjob:
	qsub job_mainm.job 

mainbjob:
	qsub job_mainb.job 

mainhjob:
	qsub job_mainh.job 




#include ${PETSC_DIR}/lib/petsc/conf/test
