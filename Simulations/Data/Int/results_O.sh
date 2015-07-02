#!/bin/bash

# ./results_O.sh

# I need to do this for three levels of sigma and two levels of variance

for p in "100000"
do
  for n in "600"
  do
    for g in "5" # ngrp
    do
      for v in $(seq 0 1) # variance AR1 or IID
      do
        for s in "3" "1" "0.33" # SNR
        do
          ./print $p $n $g $s 0 $v 1 >& myerror
        done
      done
    done
  done
done
