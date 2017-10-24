#!/bin/bash
#
#SBATCH --job-name=o3ch_exp_quad_exp
#
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4GB

module load matlab/R2017a
matlab -nodisplay < optimize_3ch_exp_quad_exp.m
