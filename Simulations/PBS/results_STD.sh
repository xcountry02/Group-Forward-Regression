#!/bin/bash

#===========================================================================
#      Program:  results_STD.sh
#
#       Author:  K. Michels
#     Language:  Bash
#   To Compile:  ./results_STD.sh nodes p0 n iter csv
#
#-----------------------------------------------------------------------------
#
#  Description:  This bash script mass creates PBS scripts using create_STD.sh.
#
#        Input:  $1 -- nodes -- Number of mpi nodes needed for job
#                $2 -- p0 -- Number of groups
#                $3 -- n -- sample size
#                $4 -- iter -- Number of iterations for each data set
#                $5 -- csv -- 1 -- Output comma separate or not
#                             0 -- Don't output comma separated
#
#       Output:  To get simulation results, have csv=1.  This submits a job
#                onto HPC.
#
#   Known Bugs:  None; all operations work as they should.
#
#===========================================================================

nodes=$1
p0=$2
n=$3
iter=$4
csv=$5

for i in $(seq 3 5) # ngrp
do
    for j in $(seq 0 1) # variance AR1 or IID
    do
        for k in $(seq 1 3) # sigma
        do
            ./create_STD.sh $nodes $p0 $n $i $j $k $iter $csv
        done
    done
done
