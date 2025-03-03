#!/bin/bash

T=8
LAM=1.0
GAMMA=$1
SIG0=$2
N=$3
MAXITER=$4
NBRUN=$5

#local0.2 
nohup nice -n 19 julia main_epids_observ.jl ${GAMMA} ${SIG0} ${N} ${MAXITER} ${NBRUN} > fileres_epid_observ_T${T}_LAM${LAM}_GAM${GAMMA}_N${N}_iter${MAXITER}_SIG0${SIG0}_run${NBRUN} & 

#Sybil

#idra

sleep 2

