#!/bin/bash


#===========================================================================
#      Program:  gather.sh
#
#       Author:  K. Michels
#     Language:  Bash
#   To Compile:  ./gather.sh p0 n n_grp var sigma mbic ibic iter
#
#-----------------------------------------------------------------------------
#
#  Description:  This bash script simply takes the first 3 lines of each
#                result file (from analysis), and gathers them into one
#                file.
#
#        Input:  $1 -- p0 -- Number of groups
#                $2 -- n -- Sample size
#                $3 -- n_grp -- group size
#                $4 -- var -- var_ind (0 or 1)
#                $5 -- sigma -- Sigma value 
#                $6 -- mbic -- main BIC value
#                $6 -- ibic -- interaction BIC value
#                $7 -- iter -- Number of iteractions for each data set
#
#       Output:  Creates "results.*" file
#
#   Known Bugs:  None; all operations work as they should.
#
#===========================================================================

# ./gather.sh p0 n n_grp var sigma mbic ibic iter

p0=$1
n=$2
ngrp=$3
var=$4
sigma=$5
mbic=$6
ibic=$7
iter=$8

for i in $(seq 0 $(($iter - 1)))
do
    if [ $i -eq 0 ]
    then
        head -3 $p0$n$ngrp$sigma$mbic$ibic$var$i.out > "results.$p0.$n.$ngrp.$sigma.$mbic.$ibic.$var.out"
    else 
        head -3 $p0$n$ngrp$sigma$mbic$ibic$var$i.out >> "results.$p0.$n.$ngrp.$sigma.$mbic.$ibic.$var.out"
    fi
done
