#!/bin/bash

#===========================================================================
#      Program:  rename.sh
#
#       Author:  K. Michels
#     Language:  Bash
#   To Compile:  ./rename.sh
#
#-----------------------------------------------------------------------------
#
#  Description:  The job names on HPC can only have a certain length,
#                so I need to rename the file appropriately so it can
#                be used in gather.sh.
#
#        Input:  None 
#
#       Output:  Correct name
#
#   Known Bugs:  None; all operations work as they should.
#
#===========================================================================

filelist=ls
ext=".out"

for file in `ls`
do

    f=${file:0:2}
    len=${#file}
    len=$(($len - 10)) 
    newfile=${file:2:$len}
    if [ "$f" = "15" ]
    then
        first=100000
        mv $file $first$newfile$ext
    elif [ "$f" = "14" ]
    then
        first=10000
        mv $file $first$newfile$ext
#echo $first$newfile$ext
    elif [ "$f" = "13" ]
    then
        first=1000
        mv $file $first$newfile$ext
    fi
done
