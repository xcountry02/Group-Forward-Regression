#!/bin/bash
#PBS -N intATmpi60
#PBS -W group_list=haozhang
#PBS -q standard
#PBS -l jobtype=cluster_only
#PBS -l select=5:ncpus=12:mem=46gb:pcmem=4gb
#PBS -l pvmem=46gb
#PBS -l walltime=24:00:00
#PBS -l cput=1440:00:00

module load openmpi
cd /gsfs1/xdisk/kamichels/gFRmpi/final_prog/INT2/REAL
mpirun -np 50 -bynode ./gFRmpiINT PM genotypeFT10.ped genotypeFT10.map 0 0 0 0 1 0
