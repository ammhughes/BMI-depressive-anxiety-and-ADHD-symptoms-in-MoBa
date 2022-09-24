#!/bin/bash
#
# job name (to identify your job)
#SBATCH --job-name=shell_imp_100
#
#Project:
#SBATCH --account=p471_tsd
#
## wall clock limit (job is terminated if limit exceeded)
#SBATCH --time=50:00:00
#
## upper bound on memory (job is terminated if limit exceeded)
#SBATCH --mem-per-cpu=2G
#
## number of (parallel) tasks (e.g. 8 cores)
#SBATCH --nodes=1 --cpus-per-task=8

#This script just calls the .do file once, but asks for 100 imputations
#To run this script: sbatch scripts/shell_imputation_x100.sh

#cd to where the do file is and other files are 
cd /cluster/p/p471/cluster/projects/mh_bmi/fullsample/imputation

#load stata
module avail Stata
module load Stata/16
module list

#define the program
MYEX="stata-mp do scripts/imputation_colossus_x100.do ${1}"

#Run the program:
$MYEX

