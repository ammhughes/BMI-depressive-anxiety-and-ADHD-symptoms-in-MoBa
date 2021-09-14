#!/bin/bash
#
# Job name (this is what the cluster will show as running - keep short):
#SBATCH --job-name=ea3
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
#SBATCH --error=./childbmi_error.txt
#
#Number of CPUs per task. 
#SBATCH --cpus-per-task=8
#
# Max memory usage (change based on resources required):
#SBATCH --mem-per-cpu=32G

#to run this script, sbatch 
#/tsd/p471/data/durable/projects/mh_bmi/scripts/colossus/manual_ea3_prs.sh

# Name used for output file location and filestems (keep short)
# Filepath for inputs
INM="/cluster/p/p471/cluster/"
# Filepath for outputs
OUTM="/cluster/p/p471/cluster/projects/mh_bmi/"

## NOT using prcise because it won't work, and because we want to know exctly which snps in there for sensitivity checks.

cp /tsd/p471/data/durable/projects/parental_education/data/lee_clumped_snps_ex_moba.txt ${OUTM}/ea3_snps.txt

head ${OUTM}/ea3_snps.txt

grep 'rs' ${OUTM}/ea3_snps.txt > ${OUTM}/temp.txt
head  ${OUTM}/temp.txt
awk '{print $1}' ${OUTM}/temp.txt > ${OUTM}/ea3_snplist.txt
head  ${OUTM}/ea3_snplist.txt
wc -l ${OUTM}/ea3_snplist.txt
#311
rm ${OUTM}/temp.txt

module add plink2/2.00a2LM

plink2 --bfile ${INM}/projects/moba_interim_release_post_imp_qc/qcd_genetic_data_Feb22/merge.no_batch.noX.geno.mind.hwe \
	--extract ${OUTM}/ea3_snplist.txt \
	--export A \
	--out ${OUTM}/ea3_snps

	
#then make a freq file! need MAF to compare aligment of snps in moba vs gwas
#already done this, as below:
#plink2 --bfile ${INM}/projects/moba_interim_release_post_imp_qc/qcd_genetic_data_Feb22/merge.no_batch.noX.geno.mind.hwe --freq --out ${OUTM}/mobq_qc_freq
#head  ${OUTM}/mobq_qc_freq.afreq
	
#copy over to tsd