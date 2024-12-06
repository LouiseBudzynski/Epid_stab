#!/bin/bash

T=8
LAM=1.0
GAMMA=$1
SIG0=$2
N=$3
MAXITER=$4
FIXSEEDS=$5
NBRUN=$6

#local
#nohup nice -n 19 julia main_epids_1graph_stab.jl ${GAMMA} ${SIG0} ${N} ${MAXITER} ${FIXSEEDS} ${NBRUN} > fileres_epid_1graph_stab_T${T}_LAM${LAM}_GAM${GAMMA}_N${N}_iter${MAXITER}_SIG0${SIG0}_fixseeds${FIXSEEDS}_run${NBRUN} & 

#Sybil
#nohup nice -n 19 julia ../julia-1.8.5/bin/main_epids_1graph_stab.jl ${GAMMA} ${SIG0} ${N} ${MAXITER} ${NBRUN} > fileres_epid_1graph_stab_T${T}_LAM${LAM}_GAM${GAMMA}_N${N}_iter${MAXITER}_SIG0${SIG0}_run${NBRUN} & 

#idra
nohup julia main_epids_1graph_stab.jl ${GAMMA} ${SIG0} ${N} ${MAXITER} ${NBRUN} > fileres_epid_1graph_stab_T${T}_LAM${LAM}_GAM${GAMMA}_N${N}_iter${MAXITER}_SIG0${SIG0}_run${NBRUN} & 

sleep 2

