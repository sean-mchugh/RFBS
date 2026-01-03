PATH=$PATH:/julia-1.9.0/bin

PATH=$PATH:/opt/R/4.2.3/lib/R
#mkdir /julia/something
#echo $UID


#export R_HOME="$R_HOME:/opt/R/4.2.3/lib/R" 
#export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/R/4.2.3/lib"
export HOME=“/julia”
export JULIA_DEPOT_PATH="/storage1/fs1/michael.landis/Active/Sean/RFBS/julia_packages"
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

IncF=${arr[9]}
ExcF=${arr[10]}
AdjMat=${arr[11]}
RJ_ANA=${arr[12]}
RJ_CLADO=${arr[13]}
ExcInd=${arr[14]}  #choose which index to pull for exlcusion treatment whether using leafing germination etc, vector of treatment names are in file
IncInd=${arr[15]}  #choose which index to pull for exlcusion treatment whether using leafing germination etc, vector of treatment names are in file
RUN=${arr[16]}

#julia scripts/RFBSsimestim_DECestim_comp_MCMC.jl ${DG} ${DL} ${SG} ${SL} ${bB} ${bRF} ${bGL} ${bDS} ${DEC} ${UP} ${SGCOV} ${SGPER} ${UT} ${RUN}

#julia scripts/Bvib_parameterized_clado_realfun_test.jl ${DG} ${DL} ${SG} ${SL} ${bB} ${bRF} ${bGL} ${bDS} ${DEC} ${IncF} ${ExcF} ${AdjMat} ${RUN}   

julia scripts/analysis_scripts/Bvib_emp_25_6_29.jl ${DG} ${DL} ${SG} ${SL} ${bB} ${bRF} ${bGL} ${bDS} ${DEC} ${IncF} ${ExcF} ${AdjMat} ${RJ_ANA} ${RJ_CLADO}  ${ExcInd} ${IncInd} ${RUN}   
