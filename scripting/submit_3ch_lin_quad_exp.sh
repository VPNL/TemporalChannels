#!/bin/bash
#
#SBATCH --job-name=c3lqe
#
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2GB

module load matlab/R2017a
matlab -nodisplay < optimize_3ch_lin_quad_exp.m