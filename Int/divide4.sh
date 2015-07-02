#!/bin/bash

# ./gather.sh p0 n n_grp var sigma mbic ibic iter

filelist=ls

for file in `ls`
do

    f=${file:0:1}
    if [ $f -eq "1" ]
    then
        line1=`sed '1q;d' $file`
        line2=`sed '2q;d' $file`
        line3=`sed '3q;d' $file`
        line3=$(echo "scale=9; ${line3}/4" | bc)
        echo $line1 > "new/$file"
        echo $line2 >> "new/$file"
        echo $line3 >> "new/$file"
    fi

done

#p0=$1
#n=$2
#ngrp=$3
#var=$4
#sigma=$5
#mbic=$6
#ibic=$7
#iter=$8
#
#for i in $(seq 0 $(($iter - 1)))
#do
#    if [ $i -eq 0 ]
#    then
#        head -3 $p0$n$ngrp$sigma$mbic$ibic$var$i.out > "results.$p0.$n.$ngrp.$sigma.$mbic.$ibic.$var.out"
#    else 
#        head -3 $p0$n$ngrp$sigma$mbic$ibic$var$i.out >> "results.$p0.$n.$ngrp.$sigma.$mbic.$ibic.$var.out"
#    fi
#done
