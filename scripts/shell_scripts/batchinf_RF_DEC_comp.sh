export LSF_DOCKER_VOLUMES="/storage1/fs1/michael.landis/Active:/storage1/fs1/michael.landis/Active"
export JOBDIR="/storage1/fs1/michael.landis/Active/RFBS_RIS"


QMAT_LIST=$(seq 1 1)

 DG=(   "true"   "true"     )
 DL=(   "false"  "true"    )
 SG=(   "true"   "true"     )
 SL=(   "true"   "true"     )
 bB=(   "false"  "false"    )
DEC=(   "false"  "false"    )
 
RF_LIST=$(seq 0 0) 
bRF=(   "true"    "false"  )

GL_LIST=$(seq 0 0) 
bGL=(   "true"    "false"  ) 

DS_LIST=$(seq 0 0) 
bDS=(   "true"    "false"  ) 

MISSING_LIST=$(seq 0 0)   
UP=(    "1.00" "1.00" "1.00" "1.00" "1.00" "1.00" "1.00" "1.00" "1.00" "1.00"   )
SGCOV=( "0.33" "0.33" "0.66" "0.33" "0.66" "0.33" "0.66" "0.33" "0.66" "1.00"   )
SGPER=( "0.25" "0.25" "0.25" "0.50" "0.50" "0.75" "0.75" "1.00" "1.00" "1.00"   )
UT=(   "false" "true" "true" "true" "true" "true" "true" "true" "true" "true"  )

PRIOR_LIST=$(seq 0 2)
PRIOR=("0.5" "1.0" "2.0")

PRIORONLY_LIST=$(seq 1 1)
PRIORONLY=("true" "false")

TREE_LIST=$(seq 0 2)
TREE_SIZE=("50" "150" "500" )

RUN_LIST=$(seq 1 110)


for i in ${RUN_LIST[@]}
do

	for x in ${QMAT_LIST[@]}
	do

	for t in ${TREE_LIST[@]}
	do
	
	for r in ${RF_LIST[@]}
	do
	
	for g in ${GL_LIST[@]}
	do
		
	for d in ${DS_LIST[@]}
	do
	
	for m in ${MISSING_LIST[@]}
	do
	for p in ${PRIOR_LIST[@]}
	do
	for o in ${PRIORONLY_LIST[@]}
	do

	
			NAME="${DG[$x]}_${DL[$x]}_${SG[$x]}_${SL[$x]}_${bB[$x]}_${bRF[$r]}_${bGL[$g]}_${bDS[$d]}_${DEC[$x]}_${UP[$m]}_${SGCOV[$m]}_${SGPER[$m]}_${UT[$m]}_${PRIOR[$p]}_${PRIORONLY[$o]}_${TREE_SIZE[$t]}_$i"
		
			bsub -G compute-michael.landis \
			-cwd /storage1/fs1/michael.landis/Active/RFBS_RIS/ \
			-o /storage1/fs1/michael.landis/Active/RFBS_RIS/stdout/$NAME  \
			-J $NAME \
			-q general \
			-g /m.seanwmchugh/rf_DEC_comp \
			-n 1 -M 5GB -R "rusage [mem=5GB] span[hosts=1]" \
			-a 'docker(sswiston/rb_tp:4)' /bin/bash /storage1/fs1/michael.landis/Active/RFBS_RIS/rf_DEC_comp.sh
			
	done		
	done	
	done
	done 
	done
	done
	done
	done
done







