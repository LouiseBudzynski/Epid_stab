#!/bin/bash

T=8
LAM=1.0
GAMMA=$1
N=$2
MAXITER=$3
NBRUN=$4

nohup julia main_epids_1graph_stab.jl ${GAMMA} ${N} ${MAXITER} ${NBRUN} > fileres_epid_1graph_stab_T${T}_LAM${LAM}_GAM${GAMMA}_N${N}_iter${MAXITER}_run${NBRUN} & 

sleep 2

