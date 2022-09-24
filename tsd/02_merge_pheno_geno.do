*MERGE COMPLETE PHENOTYPE DATA TO EACH PHENOTYPE'S GENETIC FILE - including the pgs, individual snps, and genetic covariates (PCs, batch vars etc)

clear all
set maxvar 20000
cd "N:\durable\projects\mh_bmi\fullsample\"
sysdir set PLUS "N:\durable\people\Mandy\Statafiles\ado\plus"

use "data\phenotypes_bmi_analyses.dta", clear
capture drop _merge
count
*114,307

describe M_ID_2306
count if M_ID_2306==""
describe F_ID_2306
count if F_ID_2306==""

*check how many are missing BARN_NR (because not in the birth registry file, and no genetic id for a child from which this could be derived)
count if BARN_NR==.
*475
*these are a subset of the people not present in the birth registry file:
fre present_MBRN
*annoying - can't use them. DROP THEM.
*FIRST EXCLUSION: people not in the birth registry file:
drop if present_MBRN!=1
count
*113,691

**************************************************
*MERGE IN GENETIC DATA.

*adultbmi genetic data:
merge 1:1 PREG_ID_2306 BARN_NR using "scratch\reshaped_linkable_adultbmi_genetic_pgs.dta"
rename _merge _merge_adultbmi
*matched: 76,469  (_merge==3)

count
*142,089
*too many, because includes genetic records not merged because only parent(s) genotyped and hence missing BARN_NR

count if BARN_NR==.
*28,286
count if BARN_NR==. & _merge_adultbmi==2
*28,286
count if BARN_NR!=. & _merge_adultbmi==2
*112

*the big group are trios for which only mother and/or father data available
foreach var in M_ID_2306 F_ID_2306 c_iid m_iid f_iid {
count if `var'!="" & BARN_NR==.
}
*these are partial trios, observations for which there was no BARN_NR derivable from information in the original genetic .cov file because that trio only had one or both parents, but no genetic record for the child
count if c_iid==""
count if c_iid=="" & (m_iid=="" & f_iid=="")
count if c_iid=="" & !(m_iid=="" & f_iid=="")
*could drop these, save them separately, and re-merge with 1:m on just PREG_ID_2306? 
*would be of pretty restricted use, since can't do any analysis with families where no genetic info on the child.

*DROP THEM - CAN SEPARATELY MERGE ON THAT INFORMATION WITHOUT USING BARN_NR LATER IF NEEDED.
*however, DON'T do this by whether c_iid=="", as need to keep kids without any genetic data for complete-case descriptives.
drop if BARN_NR==. & _merge_adultbmi==2
drop if BARN_NR!=. & _merge_adultbmi==2
*28,286 and 112

*how many trios?
tab complete_trio
*42,208, but haven't yet recoded families where the father withdrew consent

*repeat for childbmi:
merge 1:1 PREG_ID_2306 BARN_NR using "scratch\reshaped_linkable_childbmi_genetic_pgs.dta"
rename _merge _merge_childbmi
*exclusions as above:
drop if _merge_childbmi==2

*repeat for outcomes:
foreach outcome in depression adhd asd {
merge 1:1 PREG_ID_2306 BARN_NR using "scratch\reshaped_linkable_`outcome'_genetic_pgs.dta"
rename _merge _merge_`outcome'
*exclusions as above:
drop if _merge_`outcome'==2
}
*repeat for EA4, both versions:

foreach outcome in ea4 ea4excl23andme {
merge 1:1 PREG_ID_2306 BARN_NR using "scratch\reshaped_linkable_`outcome'_genetic_pgs.dta"
rename _merge _merge_`outcome'
*exclusions as above:
drop if _merge_`outcome'==2
}

count
*113,691

**************************************************************************************************************
*new thing: MERGE IN THE FULL GENOME PGS FOR ADULT AND CHILD BMI.
merge 1:1 PREG_ID_2306 BARN_NR using "scratch/reshaped_linkable_full_pgs_adultbmi.dta"
drop if _merge==2
drop _merge

merge 1:1 PREG_ID_2306 BARN_NR using "scratch/reshaped_linkable_full_pgs_childbmi.dta"
drop if _merge==2
drop _merge

*rename the new full gwas pgs vars for consistency with other pgs vars:
capture rename z_c_full_pgs_* c_z_full_pgs_*
capture rename z_m_full_pgs_* m_z_full_pgs_*
capture rename z_f_full_pgs_* f_z_full_pgs_*

**************************************************************************************************************

*NOW DROP ANYONE (WHOLE FAMILY) WHERE THE MOTHER WITHDREW CONSENT:
merge m:1 PREG_ID_2306 using "N:\durable\projects\mh_bmi\fullsample\PDB2306_SV_INFO_v12_mothers_children.dta"
keep if mothers_consent==1
drop _merge

merge m:1 PREG_ID_2306 using "N:\durable\projects\mh_bmi\fullsample\PDB2306_SV_INFO_v12_fathers.dta"
drop _merge

*then remove all the fathers data:
replace complete_trio=. if fathers_consent!=1
fre complete_trio
*42207

*code to missing
replace F_ID_2306="" if fathers_consent!=1
replace f_iid="" if fathers_consent!=1
replace f_fid="" if fathers_consent!=1
*everything else:
foreach var of varlist fathers_bmi_Q1QF fathers_bmi_QF fathers_inc_grp fathers_bmi_FQ2 f_*pgs* FF*  f_rs* {
replace `var'=. if fathers_consent!=1
}

save "data\merged_geno_pheno.dta", replace
