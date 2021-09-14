#!/bin/bash
#
# Job name (this is what the cluster will show as running - keep short):
#SBATCH --job-name=childbmi
#
# Project:
#SBATCH --account=p471
#
# Wall clock limit (change based on resources required):
#SBATCH --time=16:00:00 
#
#SBATCH --ntasks=1
#
# Output filename customization
# This specification is Jobname_User_JobID
#SBATCH --output=./reports/%x_%u_%j.out
#SBATCH --error=./adultbmi_error.txt
#
#Number of CPUs per task. 
#SBATCH --cpus-per-task=8
#
# Max memory usage (change based on resources required):
#SBATCH --mem-per-cpu=32G

#to run this script, sbatch 
#/tsd/p471/data/durable/projects/mh_bmi/scripts/colossus/manual_adultbmi_prs.sh

# Name used for output file location and filestems (keep short)
# Filepath for inputs
INM="/cluster/p/p471/cluster/"
# Filepath for outputs
OUTM="/cluster/p/p471/cluster/projects/mh_bmi/"

## NOT using prcise. instead, extract using imported list of snps, pre-clumped in mrbase. this way, we know exctly which snps in the pgs for sensitivity checks.
head  ${OUTM}/depression_snplist.txt
wc -l  ${OUTM}/depression_snplist.txt
#38
head  ${OUTM}/adhd_snplist.txt
wc -l  ${OUTM}/adhd_snplist.txt
#9
head  ${OUTM}/asd_snplist.txt
wc -l  ${OUTM}/asd_snplist.txt
#8
head  ${OUTM}/anxiety_snplist.txt
wc -l  ${OUTM}/anxiety_snplist.txt
#6

module add plink2/2.00a2LM

plink2 --bfile ${INM}/projects/moba_interim_release_post_imp_qc/qcd_genetic_data_Feb22/merge.no_batch.noX.geno.mind.hwe \
	--extract ${OUTM}/depression_snplist.txt \
	--export A \
	--out ${OUTM}/depression_snps
	
plink2 --bfile ${INM}/projects/moba_interim_release_post_imp_qc/qcd_genetic_data_Feb22/merge.no_batch.noX.geno.mind.hwe \
	--extract ${OUTM}/adhd_snplist.txt \
	--export A \
	--out ${OUTM}/adhd_snps
	
plink2 --bfile ${INM}/projects/moba_interim_release_post_imp_qc/qcd_genetic_data_Feb22/merge.no_batch.noX.geno.mind.hwe \
	--extract ${OUTM}/asd_snplist.txt \
	--export A \
	--out ${OUTM}/asd_snps
	
plink2 --bfile ${INM}/projects/moba_interim_release_post_imp_qc/qcd_genetic_data_Feb22/merge.no_batch.noX.geno.mind.hwe \
	--extract ${OUTM}/anxiety_snplist.txt \
	--export A \
	--out ${OUTM}/anxiety_snps

#copy over to tsd