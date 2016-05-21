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
PETSCFLAG        =-ksp_type bcgs -pc_type bjacobi -sub_ksp_type preonly -sub_pc_type sor

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
	-@${MPIEXEC} -n 4 ./main -debug

mainm:
	make clean
	make main
	-@${MPIEXEC} -n 8 ./main -ep 6.1 -m 14 -n 14 -meval 40 -neval 40 -log_view -mat_composite_merge -ksp_type bcgs -pc_type bjacobi -sub_ksp_type preonly -sub_pc_type sor

mainb:
	make clean
	make main
	-@${MPIEXEC} -n 16 ./main -ep 6.1 -m 30 -n 30 -meval 50 -neval 50 -log_view -mat_composite_merge -ksp_type bcgs -pc_type bjacobi -sub_ksp_type preonly -sub_pc_type sor

mainh:
	make clean
	make main
	-@${MPIEXEC} -n 16 ./main -ep 6.1 -m 40 -n 40 -meval 50 -neval 50 -log_view -mat_composite_merge -ksp_type bcgs -pc_type bjacobi -sub_ksp_type preonly -sub_pc_type sor


small:
	make clean
	make test 
	-@${MPIEXEC} -n 4 ./test -m 3 -n 2 -debug -mat_composite_merge -ksp_type bcgs -pc_type bjacobi -sub_ksp_type preonly -sub_pc_type sor

middle:
	make clean
	make test 
	-@${MPIEXEC} -n 8 ./test -m 100 -n 100 -log_view -mat_composite_merge -ksp_type bcgs -pc_type bjacobi -sub_ksp_type preonly -sub_pc_type sor

big:
	make clean
	make test 
	-@${MPIEXEC} -n 16 ./test -m 1000 -n 1000 -log_view -mat_composite_merge -ksp_type bcgs -pc_type bjacobi -sub_ksp_type preonly -sub_pc_type sor

testsjob:
	qsub job_tests.job 

testmjob:
	qsub job_testm.job 

testbjob:
	qsub job_testb.job 

testhjob:
	qsub job_testh.job 

mainsjob:
	qsub job_mains.job 

mainmjob:
	qsub job_mainm.job 

mainbjob:
	qsub job_mainb.job 

mainhjob:
	qsub job_mainh.job 



#include ${PETSC_DIR}/lib/petsc/conf/test
