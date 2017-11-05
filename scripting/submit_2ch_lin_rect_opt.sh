#!/bin/bash
#
#SBATCH --job-name=c2lro
#
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2GB

module load matlab/R2017a
matlab -nodisplay < optimize_2ch_lin_rect_opt.m
