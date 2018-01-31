#!/bin/bash
#
#SBATCH --job-name=c3eqco
#
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4GB
#SBATCH -p hns,normal

module load matlab/R2017a
matlab -nodisplay < optimize_3ch_exp_quad_crect_opt.m
