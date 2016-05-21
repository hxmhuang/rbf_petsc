#!/bin/csh -f

set ksptype=gmres
set pctype=asm
set subpctype=ilu
set level=0

#foreach ksptype (gmres bcgs)
#	foreach pctype (bjacobi asm)
#		foreach subpctype (jacobi sor ilu)
foreach ksptype (bcgs)
	foreach pctype (bjacobi)
		foreach subpctype (sor)
 			if ($subpctype == ilu) then
 				foreach level (0)
 				#foreach level (0)
 					echo "********Beginning new run*********"
 					echo $ksptype $pctype $subpctype $level
 					echo "***********************************"
 					mpirun -n 16 ./main -ep 6.1 -m 40 -n 40 -meval 50 -neval 50 \
 						-log_view  \
 						-pc_type $pctype -ksp_type $ksptype \
 						-sub_ksp_type preonly -sub_pc_type $subpctype 
 					#	-sub_pc_ilu_levels $level 
 				    #		-pc_factor_levels $level 
# 						-ksp_monitor -sles_view -optionsleft \
# 						 2>&1 | tee log_$ksptype_$pctype_$subpctype_$level.log 
 				end
 			else
				echo "********Beginning new run*********"
				echo $ksptype $pctype $subpctype
				echo "***********************************"
				mpirun -n 16 ./main -ep 6.1 -m 40 -n 40 -meval 50 -neval 50 \
					-log_view  \
					-pc_type $pctype -ksp_type $ksptype \
					-sub_ksp_type preonly -sub_pc_type $subpctype 
					#-ksp_monitor -sles_view -optionsleft \
#					-ksp_monitor -sles_view -optionsleft 
#					2>&1 | tee log_$ksptype_$pctype_$subpctype.log 
			endif
		end
	end
end
