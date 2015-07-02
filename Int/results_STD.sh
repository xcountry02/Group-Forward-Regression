#!/bin/bash

#===========================================================================
#      Program:  subset_DS.sh
#
#       Author:  K. Michels
#     Language:  Bash
#   To Compile:  ./results_DS.sh node p0 n iter
#
#-----------------------------------------------------------------------------
#
#  Description:  This bash script works on a subset of all_DS.sh
#
#        Input:  $1 -- nodes -- Number of mpi nodes needed for job
#                $2 -- p0 -- Number of groups
#                $3 -- n -- Sample size
#                $4 -- iter -- Number of iterations for each data set
#                $5 -- csv -- 1 -- Output comma separate or not
#                             0 -- Don't output comma separated
#
#       Output:  To get simulation results, have csv=1, the results files 
#                that are created can be inputted into the simulation data 
#                results.  Enough information is obtained to summarize the 
#                method well.
#
#   Known Bugs:  None; all operations work as they should.
#
#===========================================================================

nodes=$1
p0=$2
n=$3
iter=$4
csv=$5

# I need this for three levels of sigma {1, 2, 3}, two levels of variance {0, 1}
#   and 3 levels of ngrp {3, 4, 5}.  Also I want to report csv, so csv = 1.  
#   This program does it all, it creates the data set, runs mpi analysis on it,
#   reports the analysis, then deletes the data set

for g in $(seq 3 5) # ngrp
do
  for v in $(seq 0 1) # variance AR1 or IID
  do
    for s in "0.33" "1" "3" #SNR 
    do
      for l in $(seq 0 $(($iter - 1)))
      do
        # Make sure to have the directory to the correct print function
        #
        # echo /gsfs1/xdisk/kamichels/gFRmpi/Program/Data/Int/print $p0 $n $g $s $l $v
        /gsfs1/xdisk/kamichels/gFRmpi/Program/Data/Int/print $p0 $n $g $s $l $v 0
        (mpirun -np $nodes ./gFRmpiINT STD Y$l.$p0.$g.$n.${s}_${v} X$l.$p0.$g.$n.${s}_${v} 0 $csv $v $s > "$p0$n$g$s$m$i$v$l.out") >& myerror
        rm X$l.$p0.$g.$n.${s}_${v} Y$l.$p0.$g.$n.${s}_${v}
      done
    done
  done
done
