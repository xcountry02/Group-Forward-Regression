#!/bin/bash
#PBS -N ATmpi60
#PBS -W group_list=haozhang
#PBS -q standard
#PBS -l jobtype=cluster_only
#PBS -l select=5:ncpus=12:mem=46gb:pcmem=4gb
#PBS -l pvmem=46gb
#PBS -l walltime=6:00:00
#PBS -l cput=360:00:00

module load openmpi
cd /gsfs1/xdisk/kamichels/gFRmpi/final_prog/
mpirun -np 60 ./gFRmpi PM genotypeFT10.ped genotypeFT10.map 0 0 0
