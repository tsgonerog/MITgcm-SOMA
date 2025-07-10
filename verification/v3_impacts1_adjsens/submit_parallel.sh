#!/bin/bash

#SBATCH -J parV3imp
#SBATCH -o parV3imp.%j.out
#SBATCH -e parV3imp.%j.err
#SBATCH -N 2
#SBATCH -n 4
#SBATCH -t 48:00:00

#SBATCH --mail-user=tanvirshahriar@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

#module load gnu12/12.2.0
#module load openmpi4/4.1.5
#module load phdf5/1.14.0
#module load netcdf-fortran/4.6.0
#module load netcdf/4.9.0
#module load prun

# Prepare run directory and execute model
cd run_tap_parallel
rm *
ln -s ../input_tap/* .
../input_tap/prepare_run
ln -s ../build_tap_parallel/mitgcmuv_tap_adj .
#prun ./mitgcmuv_tap_adj > output_tap_adj.txt 2>&1
mpiexec -n 4 ./mitgcmuv_tap_adj > output_tap_adj.txt 2>&1

