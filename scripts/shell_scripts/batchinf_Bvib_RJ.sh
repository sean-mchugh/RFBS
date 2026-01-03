export LSF_DOCKER_VOLUMES="/storage1/fs1/michael.landis/Active:/storage1/fs1/michael.landis/Active"
export JOBDIR="/storage1/fs1/michael.landis/Active/Sean/RFBS"



BATCHES=$(seq 0 0)
INCBATCHES=$(seq 0 0)
EXCLBATCHES=$(seq 0 0)
ADMATBATCHES=$(seq 0 0)
RJ_ANABATCHES=$(seq 0 1)
RJ_CLADOBATCHES=$(seq 1 1)


 
    DG=("true" "true" )
    DL=("true" "false" )
    SG=("true"  "true"  )
    SL=("true"  "true"  )
    bB=("false" "false" )
   bRF=("true"  "true"  )
   bGL=("true"  "true"  )
   bDS=("true" "true"   )
   
   
   DEC=("false" "true" )
  
  
  IncF=("true"  "false")
  ExcF=("true"  "false") 
AdjMAT=("true"  "false") 

RJ_ANA=("true" "false")
RJ_CLADO=("true" "false")


RUN_LIST=$(seq 1 5)
EXCIND_LIST=$(seq 1 2)
INCIND_LIST=$(seq 1 2)



for x in ${BATCHES[@]}
do


	for y in ${INCBATCHES[@]}
	do
	
	for z in ${EXCLBATCHES[@]}
	do
	
	for a in ${ADMATBATCHES[@]}
	do
	
		for b in ${RJ_ANABATCHES[@]}
		do
		
		for c in ${RJ_CLADOBATCHES[@]}
		do
		
		
		
			for e in ${EXCIND_LIST[@]}
			do
		
		
			for n in ${INCIND_LIST[@]}
			do
			
			for i in ${RUN_LIST[@]}
			do
		
				
				NAME="${DG[$x]}_${DL[$x]}_${SG[$x]}_${SL[$x]}_${bB[$x]}_${bRF[$x]}_${bGL[$x]}_${bDS[$x]}_${DEC[$x]}_${IncF[$y]}_${ExcF[$z]}_${AdjMAT[$a]}_${RJ_ANA[$b]}_${RJ_CLADO[$c]}_${e}_${n}_${i}"	
				bsub -G compute-michael.landis \
				-cwd /storage1/fs1/michael.landis/Active/Sean/RFBS/ \
				-o /storage1/fs1/michael.landis/Active/Sean/RFBS/outfiles/emp/viburnum/stdout/$NAME  \
				-J $NAME \
				-q general \
				-g /m.seanwmchugh/Bvib \
				-n 1 -M 2GB -R "rusage [mem=2GB] span[hosts=1]" \
				-a 'docker(sswiston/rb_tp:4)' /bin/bash /storage1/fs1/michael.landis/Active/Sean/RFBS/scripts/shell_scripts/inf_Bvib_RJ.sh
				
				
			done
			done
			done
		done
		done
	done
	done
	done
done

