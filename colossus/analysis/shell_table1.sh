#!/bin/bash
#
# job name (to identify your job)
#SBATCH --job-name=shell_table1
#
#Project:
#SBATCH --account=p471
#
## wall clock limit (job is terminated if limit exceeded)
#SBATCH --time=50:00:00
#
## upper bound on memory (job is terminated if limit exceeded)
#SBATCH --mem-per-cpu=2G
#
## number of (parallel) tasks (e.g. 8 cores)
#SBATCH --nodes=1 --cpus-per-task=8

#cd to where the do file is and other files are 
cd /cluster/p/p471/cluster/projects/mh_bmi/fullsample/analysis

#To run this script: sbatch scripts/shell_table1.sh

#load stata
module avail Stata
module load Stata/16
module list

#define the program
MYEX="stata-mp do scripts/table1.do ${1}"
#MYEX="stata-mp do scripts/table1.do"

#Run the program:
$MYEX

