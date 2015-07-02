#!/bin/bash
#PBS -N TRUEdataMPI30
#PBS -W group_list=haozhang
#PBS -q standard
#PBS -l jobtype=cluster_only
#PBS -l select=6:ncpus=12:mem=46gb:pcmem=4gb
#PBS -l pvmem=46gb
#PBS -l walltime=48:00:00
#PBS -l cput=3456:00:00

module load openmpi
cd /gsfs1/xdisk/kamichels/gFRmpi/final_prog/
mpirun -np 30 -bynode ./gFRmpi PM genotype.ped genotype.map 0 0 0
