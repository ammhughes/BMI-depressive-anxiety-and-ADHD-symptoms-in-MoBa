

*FOR RUNNING SNP-LEVEL ASSOCIATIONS ON COLOSSUS: MAKE MINI DOSAGE FILES

*CHILDBMI

clear all
set maxvar 30000
sysdir set PLUS "N:/durable/people/Mandy/statafiles/ado/plus"
cd "N:/durable/projects/mh_bmi/fullsample"

*****************************************************************************
*export sorted list of SNPs:
use data/snplevel_info/moba_gwas_snplevelinfo_childbmi.dta, clear
keep SNP
sort SNP

outsheet using scratch/childbmi_sorted_snplist_test.txt, nonames noquote replace
*keep if strpos(SNP,"rs1")>0
*IN NOTEPAD++, WILL NEED TO MANUALLY USE EDIT-->LINE OPERATIONS-->JOIN LINES TO GET ON ONE LINE. CAN THEN COPY AND PASTE INTO RELEVANT SHELL SCRIPT

******************************************************************************
*then open dosage data:
use "scratch/reshaped_linkable_childbmi_genetic_pgs.dta", clear
capture drop _merge
*get the list of snps for child bmi:
merge 1:1 _n using data/snplevel_info/moba_gwas_snplevelinfo_childbmi.dta
drop _merge
list SNP harm_beta_gwas palindromic in 1/322

*save the order of snps in a global for later:
levelsof SNP, local(childbmi_snps)
*change every " to a space (34 to 32)
local childbmi_snplist_temp: subinstr local childbmi_snps "char(34)" "char(32)", all
global childbmi_snplist: subinstr local childbmi_snplist_temp "rs" " rs", all
display $childbmi_snplist

*restrict to full trio:
keep if complete_trio==1

*prep the snp vars for this usage. 
*need to trim off the last part (specification of effect allele) to use easily in a loop storing snp-outcome coefficients
rename *_a *
rename *_t *
rename *_g *
rename *_c *

*then also rename everyone's snp vars to have a person SUFFIX rather than prefix, since a prefix but not a suffix comes with insertion of an uncalled-for space in the loop that stops the loop from running.
rename c_rs* rs*_c
rename m_rs* rs*_m
rename f_rs* rs*_f

*now section off  each snp in turn
*get rid of it with ltrim
foreach snp in $childbmi_snplist {
display "`snp'"
local snp=ltrim("`snp'")
preserve
keep `snp'_c `snp'_m `snp'_f PREG_ID_2306 M_ID_2306 F_ID_2306 c_fid c_iid 
save "scratch/childbmi_singlesnps/dosage/`snp'_dosage.dta", replace
restore
}

***************************************************************************************************************************
*ADULTBMI

clear all
set maxvar 30000
sysdir set PLUS "N:/durable/people/Mandy/statafiles/ado/plus"
cd "N:/durable/projects/mh_bmi/fullsample"

*****************************************************************************
*export sorted list of SNPs:
use data/snplevel_info/moba_gwas_snplevelinfo_adultbmi.dta, clear
keep SNP
sort SNP

outsheet using scratch/adultbmi_sorted_snplist_test.txt, nonames noquote replace
*keep if strpos(SNP,"rs1")>0
*IN NOTEPAD++, WILL NEED TO MANUALLY USE EDIT-->LINE OPERATIONS-->JOIN LINES TO GET ON ONE LINE. CAN THEN COPY AND PASTE INTO RELEVANT SHELL SCRIPT

******************************************************************************
*then open dosage data:
use "scratch/reshaped_linkable_adultbmi_genetic_pgs.dta", clear
capture drop _merge
*get the list of snps for child bmi:
merge 1:1 _n using data/snplevel_info/moba_gwas_snplevelinfo_adultbmi.dta
drop _merge
list SNP harm_beta_gwas palindromic in 1/955

*save the order of snps in a global for later:
levelsof SNP, local(adultbmi_snps)
*change every " to a space (34 to 32)
local adultbmi_snplist_temp: subinstr local adultbmi_snps "char(34)" "char(32)", all
global adultbmi_snplist: subinstr local adultbmi_snplist_temp "rs" " rs", all
display $adultbmi_snplist

*restrict to full trio:
keep if complete_trio==1

*prep the snp vars for this usage. 
*need to trim off the last part (specification of effect allele) to use easily in a loop storing snp-outcome coefficients
rename *_a *
rename *_t *
rename *_g *
rename *_c *

*then also rename everyone's snp vars to have a person SUFFIX rather than prefix, since a prefix but not a suffix comes with insertion of an uncalled-for space in the loop that stops the loop from running.
rename c_rs* rs*_c
rename m_rs* rs*_m
rename f_rs* rs*_f

*now section off  each snp in turn
*get rid of it with ltrim
foreach snp in $adultbmi_snplist {
display "`snp'"
local snp=ltrim("`snp'")
preserve
keep `snp'_c `snp'_m `snp'_f PREG_ID_2306 M_ID_2306 F_ID_2306 c_fid c_iid 
save "scratch/adultbmi_singlesnps/dosage/`snp'_dosage.dta", replace
restore
}