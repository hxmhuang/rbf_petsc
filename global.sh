#!/bin/csh -f

set ksptype=gmres
set pctype=bjacobi
set subpctype=jacobi
set level=0

#foreach ksptype (gmres bcgs tfqmr)
#	foreach pctype (bjacobi asm)
#		foreach subpctype (jacobi sor ilu)
foreach ksptype (gmres)
	foreach pctype (asm)
		foreach subpctype (ilu)
 			if ($subpctype == ilu) then
 				foreach level (0 1 2)
 					echo "********Beginning new run*********"
 					echo $ksptype $pctype $subpctype $level
 					echo "***********************************"
 					mpirun -n 16 ./main -ep 6.1 -m 30 -n 30 -meval 50 -neval 50 \
 						-log_view  \
 						-pc_type $pctype -ksp_type $ksptype \
 						-sub_ksp_type preonly -sub_pc_type $subpctype \
 						-sub_pc_ilu_levels $level 
# 						-ksp_monitor -sles_view -optionsleft \
# 						 2>&1 | tee log_$ksptype_$pctype_$subpctype_$level.log 
 				end
 			else
				echo "********Beginning new run*********"
				echo $ksptype $pctype $subpctype
				echo "***********************************"
				mpirun -n 16 ./main -ep 6.1 -m 30 -n 30 -meval 50 -neval 50 \
					-log_view  \
					-pc_type $pctype -ksp_type $ksptype \
					-sub_ksp_type preonly -sub_pe_type $subpctype 
					#-ksp_monitor -sles_view -optionsleft \
#					-ksp_monitor -sles_view -optionsleft 
#					2>&1 | tee log_$ksptype_$pctype_$subpctype.log 
			endif
		end
	end
end
