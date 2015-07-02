#!/bin/bash

#===========================================================================
#      Program:  all_DS.sh
#
#       Author:  K. Michels
#     Language:  Bash
#   To Compile:  ./all_DS.sh iter oracle
#
#-----------------------------------------------------------------------------
#
#  Description:  This bash script takes in two parameters and it will go
#                through every possible combination. 
#
#        Input:  $1 -- iter -- Number of iterations for each data set
#                $2 -- oracle -- bool -- 0 -- print data sets
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

iter=$1
ora=$2

for p in 1000 10000 100000 #total number of groups
do
    for n in 100 500 #sample size
    do
        for g in $(seq 3 5) #ngrp
        do
            for v in $(seq 0 1) #variance AR1 or IID
            do
                for s in $(seq 1 3) #sigma
                do
                    for l in $(seq 0 $(($iter - 1)))
                    do
                        ./print $p $n $g $s $l $v $ora
                        # echo ./print $p $n $g $s $l $v $ora
                    done
                done
            done
        done
    done
done
