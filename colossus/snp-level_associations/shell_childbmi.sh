#!/bin/bash
#
# job name (to identify your job)
#SBATCH --job-name=shell_childbmi
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

#This script repeatedly calls the same .do file, but subbing in a different snp each time
#To run this script: sbatch /cluster/p/p471/cluster/projects/mh_bmi/fullsample/snp-level_associations/scripts/shell_childbmi.sh

#cd to where the do file is and other files are 
cd /cluster/p/p471/cluster/projects/mh_bmi/fullsample/snp-level_associations

#load stata
module avail Stata
module load Stata/16
module list

#define the program
MYEX="stata-mp do scripts/innermost_childbmi.do ${1}"
#MYEX="stata-mp do scripts/innermost_childbmi.do snp"

echo "is this passing {snp}"
echo "is this passing ${snp}"
echo "is this passing $snp"
echo "is this passing ${1}"

#Run the program:
$MYEX

