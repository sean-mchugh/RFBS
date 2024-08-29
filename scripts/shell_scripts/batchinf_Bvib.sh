export LSF_DOCKER_VOLUMES="/storage1/fs1/michael.landis/Active:/storage1/fs1/michael.landis/Active"
export JOBDIR="/storage1/fs1/michael.landis/Active/RFBS_RIS"



BATCHES=$(seq 0 0)
INCBATCHES=$(seq 0 0)
EXCLBATCHES=$(seq 1 1)
ADMATBATCHES=$(seq 0 1)

 
    DG=("true" "true" )
    DL=("false" "false" )
    SG=("true"  "true"  )
    SL=("true"  "true"  )
    bB=("false" "false" )
   bRF=("true"  "true"  )
   bGL=("true"  "true"  )
   bDS=("true" "true" )
   
   
   DEC=("false" "true" )
  
  
  IncF=("true"  "false")
  ExcF=("true"  "false") 
AdjMAT=("true"  "false") 

RUN_LIST=$(seq 1 10)
EXCIND_LIST=$(seq 1 1)
INCIND_LIST=$(seq 1 1)

for x in ${BATCHES[@]}
do


for y in ${INCBATCHES[@]}
do

for z in ${EXCLBATCHES[@]}
do

for a in ${ADMATBATCHES[@]}
do


	for e in ${EXCIND_LIST[@]}
	do


	for n in ${INCIND_LIST[@]}
	do
	
	for i in ${RUN_LIST[@]}
	do

		
		NAME="${DG[$x]}_${DL[$x]}_${SG[$x]}_${SL[$x]}_${bB[$x]}_${bRF[$x]}_${bGL[$x]}_${bDS[$x]}_${DEC[$x]}_${IncF[$y]}_${ExcF[$z]}_${AdjMAT[$a]}_${e}_${n}_${i}"	
		bsub -G compute-michael.landis \
		-cwd /storage1/fs1/michael.landis/Active/RFBS_RIS/ \
		-o /storage1/fs1/michael.landis/Active/RFBS_RIS/Bvib_stdout/$NAME  \
		-J $NAME \
		-q general \
		-g /m.seanwmchugh/Bvib \
		-n 1 -M 5GB -R "rusage [mem=5GB] span[hosts=1]" \
		-a 'docker(sswiston/rb_tp:4)' /bin/bash /storage1/fs1/michael.landis/Active/RFBS_RIS/inf_Bvib.sh
		
		
	done
	done
	done
done
done
done
done





BATCHES=$(seq 0 0)
INCBATCHES=$(seq 0 0)
EXCLBATCHES=$(seq 1 1)
ADMATBATCHES=$(seq 0 1)

 
    DG=("true" "true" )
    DL=("false" "false" )
    SG=("true"  "true"  )
    SL=("true"  "true"  )
    bB=("false" "false" )
   bRF=("true"  "true"  )
   bGL=("true"  "true"  )
   bDS=("true" "true" )
   
   
   DEC=("false" "true" )
  
  
  IncF=("true"  "false")
  ExcF=("true"  "false") 
AdjMAT=("true"  "false") 

RUN_LIST=$(seq 1 10)
EXCIND_LIST=$(seq 1 1)
INCIND_LIST=$(seq 2 2)

for x in ${BATCHES[@]}
do


for y in ${INCBATCHES[@]}
do

for z in ${EXCLBATCHES[@]}
do

for a in ${ADMATBATCHES[@]}
do


	for e in ${EXCIND_LIST[@]}
	do


	for n in ${INCIND_LIST[@]}
	do
	
	for i in ${RUN_LIST[@]}
	do

		
		NAME="${DG[$x]}_${DL[$x]}_${SG[$x]}_${SL[$x]}_${bB[$x]}_${bRF[$x]}_${bGL[$x]}_${bDS[$x]}_${DEC[$x]}_${IncF[$y]}_${ExcF[$z]}_${AdjMAT[$a]}_${e}_${n}_${i}"	
		bsub -G compute-michael.landis \
		-cwd /storage1/fs1/michael.landis/Active/RFBS_RIS/ \
		-o /storage1/fs1/michael.landis/Active/RFBS_RIS/Bvib_stdout/$NAME  \
		-J $NAME \
		-q general \
		-g /m.seanwmchugh/Bvib \
		-n 1 -M 5GB -R "rusage [mem=5GB] span[hosts=1]" \
		-a 'docker(sswiston/rb_tp:4)' /bin/bash /storage1/fs1/michael.landis/Active/RFBS_RIS/inf_Bvib.sh
		
		
	done
	done
	done
done
done
done
done




