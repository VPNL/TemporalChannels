#!/bin/bash
#
#SBATCH --job-name=c1l
#
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2GB
#SBATCH -p hns,normal

module load matlab/R2017a
matlab -nodisplay < optimize_1ch_lin.m
