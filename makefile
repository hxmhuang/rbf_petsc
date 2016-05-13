CFLAGS	         =  
FFLAGS	         =-Wno-tabs
CPPFLAGS         =
FPPFLAGS         =
LOCDIR           = src/ksp/ksp/examples/tutorials/
MANSEC           = KSP
CLEANFILES       = main test *.o *.mod 
NP               = 1
OBJ				 = matrix.o vector.o rbf.o  
OBJMAIN			 = ${OBJ} main.o 
OBJTEST			 = ${OBJ} test.o 

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main: ${OBJMAIN}  chkopts
	-${FLINKER} -o main ${OBJMAIN}  ${PETSC_KSP_LIB}
#	${RM} *.mod *.o
test: ${OBJTEST}  chkopts
	-${FLINKER} -o test ${OBJTEST}  ${PETSC_KSP_LIB}

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

mains:
	make clean
	make main
	-@${MPIEXEC} -n 4 ./main -m 3 -n 3 -debug -mat_composite_merge 

mainm:
	make clean
	make main
	-@${MPIEXEC} -n 8 ./main -m 30 -n 30 -log_view -mat_composite_merge 

mainb:
	make clean
	make main
	-@${MPIEXEC} -n 16 ./main -m 100 -n 100 -log_view -mat_composite_merge 


small:
	make clean
	make test 
	-@${MPIEXEC} -n 4 ./test -m 3 -n 2 -debug -mat_composite_merge  

middle:
	make clean
	make test 
	-@${MPIEXEC} -n 8 ./test -m 100 -n 100 -log_view -mat_composite_merge 

big:
	make clean
	make test 
	-@${MPIEXEC} -n 16 ./test -m 900 -n 900 -log_view -mat_composite_merge 

smalljob:
	make clean
	make test
	qsub job_small.job 

middlejob:
	make clean
	make test
	qsub job_middle.job 

bigjob:
	make clean
	make test
	qsub job_big.job 

hugejob:
	make clean
	make test
	qsub job_huge.job 

#include ${PETSC_DIR}/lib/petsc/conf/test
