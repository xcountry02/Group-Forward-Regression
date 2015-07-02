#!/bin/bash

# ./qdel.sh sub1 sub2

for i in $(seq $1 $2)
do
    qdel $i
done
