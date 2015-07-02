#!/bin/bash

#===========================================================================
#      Program:  create_STD.sh
#
#       Author:  K. Michels
#     Language:  Bash
#   To Compile:  ./create_STD.sh nodes p0 n n_grp var sigma iter csv
#
#-----------------------------------------------------------------------------
#
#  Description:  This bash script creates a PBS script, and submits it.  Make
#                sure to specify which directory data sets are in.
#
#        Input:  $1 -- nodes -- Number of mpi nodes needed for job
#                $2 -- p0 -- Number of groups
#                $3 -- n -- sample size
#                $4 -- n_grp -- group size
#                $5 -- var -- variance index (0 or 1) (AR(1) or IID)
#                $6 -- sigma -- Sample variance (>0)
#                $7 -- iter -- Number of iterations for each data set
#                $8 -- csv -- 1 -- Output comma separate or not
#                             0 -- Don't output comma separated
#
#       Output:  To get simulation results, have csv=1.  This submits a job
#                onto HPC.
#
#   Known Bugs:  None; all operations work as they should.
#
#===========================================================================

# The job output will have multiple results for multiple BIC.  When the 
# analysis is done, need to separate the results into multiple files.
# Also this bash script can be used for both MAIN and INT.  Make sure to
# change the directory and executable depending on that.

nodes=$1
p0=$2
n=$3
ngrp=$4
var=$5
sigma=$6
csv=$8
iter=$7

fp0=${p0:0:1}
p0len=$((${#p0} - 1))

for i in $(seq 0 $(($iter - 1)))
do
    echo "#!/bin/bash" > "STDmpi$(printf "%03d" "$i").sh"
    echo "#PBS -N $fp0$p0len$n$ngrp$sigma$var$i" >> "STDmpi$(printf "%03d" "$i").sh"
    echo "#PBS -W group_list=haozhang" >> "STDmpi$(printf "%03d" "$i").sh"
    echo "#PBS -q standard" >> "STDmpi$(printf "%03d" "$i").sh"
### Use the following for memory intensive mpi ###
    echo "#PBS -l jobtype=cluster_only" >> "STDmpi$(printf "%03d" "$i").sh"
    echo "#PBS -l select=4:ncpus=12:mem=46gb:pcmem=4gb" >> "STDmpi$(printf "%03d" "$i").sh"
    echo "#PBS -l pvmem=46gb" >> "STDmpi$(printf "%03d" "$i").sh"
    echo "#PBS -l walltime=2:00:00" >> "STDmpi$(printf "%03d" "$i").sh"
    echo "#PBS -l cput=96:00:00" >> "STDmpi$(printf "%03d" "$i").sh"
### Use the following fo less memory intensive mpi (faster) ###
#   echo "#PBS -l jobtype=small_mpi" >> "STDmpi$(printf "%03d" "$i").sh"
#   echo "#PBS -l select=4:ncpus=12:mem=31gb" >> "STDmpi$(printf "%03d" "$i").sh"
#   echo "#PBS -l walltime=00:10:00" >> "STDmpi$(printf "%03d" "$i").sh"
#   echo "#PBS -l cput=08:00:00" >> "STDmpi$(printf "%03d" "$i").sh"
    echo "module load openmpi" >> "STDmpi$(printf "%03d" "$i").sh"
    echo "cd /gsfs1/xdisk/kamichels/gFRmpi/Program/Main" >> "STDmpi$(printf "%03d" "$i").sh"
    echo "mpirun -np $nodes ./gFRmpi STD Y$i.$p0.$ngrp.$n.${sigma}_${var} X$i.$p0.$ngrp.$n.${sigma}_${var} 0 $csv $var $sigma" >> "STDmpi$(printf "%03d" "$i").sh"
done

for i in $(seq 0 $(($iter - 1)))
do
    qsub "STDmpi$(printf "%03d" "$i").sh"
done
