#!/bin/bash

#===========================================================================
#      Program:  subset_DS.sh
#
#       Author:  K. Michels
#     Language:  Bash
#   To Compile:  ./subset_DS.sh p0 n iter oracle
#
#-----------------------------------------------------------------------------
#
#  Description:  This bash script works on a subset of all_DS.sh
#
#        Input:  $1 -- p0 -- Number of iterations for each data set
#                $2 -- n -- Number of iterations for each data set
#                $3 -- iter -- Number of iterations for each data set
#                $4 -- oracle -- bool -- 0 -- print data sets
#                                        1 -- create summary result cards
#
#       Output:  If oracle==0 the output are data sets, o/w the output is
#                two files that contain all of the necessary
#                information for the simulation.  These files have a very
#                specific format and they are referred to as result cards.
#
#   Known Bugs:  None; all operations work as they should.
#
#===========================================================================

p0=$1
n=$2
iter=$3
ora=$4

for i in $(seq 3 5) #ngrp
do
    for j in $(seq 0 1) #variance AR1 or IID
    do
        for k in $(seq 1 3) #sigma
        do
            for l in $(seq 0 $(($iter - 1)))
            do
                ./print $p0 $n $i $k $l $j 0
            done
        done
    done
done
