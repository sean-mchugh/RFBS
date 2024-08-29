PATH=$PATH:/julia-1.9.0/bin

PATH=$PATH:/opt/R/4.2.3/lib/R
#mkdir /julia/something
#echo $UID


#export R_HOME="$R_HOME:/opt/R/4.2.3/lib/R" 
#export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/R/4.2.3/lib"
export HOME=“/julia”
export JULIA_DEPOT_PATH="/storage1/fs1/michael.landis/Active/RFBS_RIS/julia_packages"
export JULIA_CPU_THREADS=1
#export JULIA_DEPOT_PATH="/root/.julia"
#export JULIA_DEPOT_PATH="/julia/.julia"


IFS="_"
read -ra arr <<< $LSB_JOBNAME
IFS=" "


 DG=${arr[0]}
 DL=${arr[1]}
 SG=${arr[2]}
 SL=${arr[3]}
 bB=${arr[4]}
bRF=${arr[5]}
bGL=${arr[6]}
bDS=${arr[7]}
DEC=${arr[8]}

UP=${arr[9]}
SGCOV=${arr[10]}
SGPER=${arr[11]}
UT=${arr[12]}
PRIOR=${arr[13]}
PRIORONLY=${arr[14]}
TREE=${arr[15]}

RUN=${arr[16]}

#julia scripts/RFBSsimestim_DECestim_comp_MCMC.jl ${DG} ${DL} ${SG} ${SL} ${bB} ${bRF} ${bGL} ${bDS} ${DEC} ${UP} ${SGCOV} ${SGPER} ${UT} ${RUN}

julia scripts/RFBSsimestim_comp_MCMC.jl ${DG} ${DL} ${SG} ${SL} ${bB} ${bRF} ${bGL} ${bDS} ${DEC} ${UP} ${SGCOV} ${SGPER} ${UT} ${PRIOR} ${PRIORONLY} ${TREE} ${RUN}
julia scripts/DECestim_comp_MCMC.jl ${DG} ${DL} ${SG} ${SL} ${bB} ${bRF} ${bGL} ${bDS} ${DEC} ${UP} ${SGCOV} ${SGPER} ${UT} ${PRIOR} ${PRIORONLY} ${TREE} ${RUN}
