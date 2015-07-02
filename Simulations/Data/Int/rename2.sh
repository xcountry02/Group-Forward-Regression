#!/bin/bash

# ./qdel.sh sub1 sub2

filelist=ls

for file in `ls`
do
    echo $file
    c1=${file:27:38}
    c=${c1:0:1}
    f=${file:0:27}
    l=${file:28:38}

    c=$((c + 4))

    mv $file $f$c$l

done

filelist=ls

for file in `ls`
do

    c1=${file:27:38}
    c=${c1:0:1}
    f=${file:0:27}
    l=${file:28:38}

    c=$((c - 3))

    mv $file $f$c$l

done
