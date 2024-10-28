#!/bin/bash

T=8
LAM=1.0
GAMMA=$1
N=$2
MAXITER=$3
NBRUN=$4

nohup nice -n 19 ../julia-1.8.5/bin/julia main_epids_conv.jl ${GAMMA} ${N} ${MAXITER} ${NBRUN} > fileres_epid_conv_T${T}_LAM${LAM}_GAM${GAMMA}_N${N}_iter${MAXITER}_run${NBRUN} & 
#nohup nice -n 19 julia main_epids_conv.jl ${GAMMA} ${N} ${MAXITER} ${NBRUN} > fileres_epid_conv_T${T}_LAM${LAM}_GAM${GAMMA}_N${N}_iter${MAXITER}_run${NBRUN} & 

sleep 2

