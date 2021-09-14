

*FOR RUNNING SNP-LEVEL ASSOCIATIONS ON COLOSSUS: MAKE MINI DOSAGE FILES

*CHILDBMI

clear all
set maxvar 30000
sysdir set PLUS "N:/data/durable/people/Mandy/statafiles/ado/plus"
cd "N:/data/durable/projects/mh_bmi/"

*****************************************************************************
*export chunked, sorted lists of the ones you HAVEN'T dropped for being palindromic:
use data/moba_gwas_snplevelinfo_childbmi.dta, clear
keep SNP
sort SNP

*outsheet using scratch/childbmi_sorted_snplist_test.txt, nonames noquote replace
*keep if strpos(SNP,"rs1")>0

*for simplicity, make 5 mini-lists, based on first number:
foreach num of numlist 0 2 4 6 8 {
    local nplus1=`num'+1
	display `nplus1'
preserve
keep if strpos(SNP,"rs`num'")>0 | strpos(SNP,"rs`nplus1'")>0
list
outsheet using scratch/childbmi_sorted_snplist_part`num'.txt, nonames noquote replace
restore
}
clear
*IN EACH, IN NOTEPAD++, WILL NEED TO MANUALLY USE LINE OPERATIONS-->JOIN LINES TO GET ON ONE LINE. CAN THEN COPY AND PASTE INTO RELEVANT SHELL SCRIPT

*****************************************************************************


*start with dosage data - trim to single snps in turn:
use "scratch/childbmi_justsnps.dta", clear

*get the list of snps for child bmi:
merge 1:1 _n using data/moba_gwas_snplevelinfo_childbmi.dta
drop _merge
list SNP harm_beta_gwas palindrome in 1/315

*save the order of snps in a global for later:

levelsof SNP, local(childbmi_snps) clean
display "`childbmi_snps'"

/*
levelsof SNP, local(childbmi_snps) clean
global childbmi_snplist="`childbmi_snps'"
display $childbmi_snplist
*error:
*rs10095724 ambiguous abbreviation
*/

levelsof SNP, local(childbmi_snps)
*change every " to a space (34 to 32)
local childbmi_snplist_temp: subinstr local childbmi_snps "char(34)" "char(32)", all
global childbmi_snplist: subinstr local childbmi_snplist_temp "rs" " rs", all
display $childbmi_snplist
*yes

/*this doesn't work, cam't make subinstr involving combos of numbers and spaces
levelsof SNP, local(childbmi_snps)
display "childbmi_snps"
*change every " to a space (34 to 32)
local childbmi_snplist_temp: subinstr local childbmi_snps "char(34)" "char(32)", all
local childbmi_snplist_temp2: subinstr local childbmi_snplist_temp "4rs" "4 rs", all
global childbmi_snplist: subinstr local childbmi_snplist_temp2 "  " ""
display $childbmi_snplist
*no
*/

*it's ok - can get rid of it with ltrim!
*check it works in a loop:
foreach snp in $childbmi_snplist {
	*display "`snp'"
	local snp=ltrim("`snp'")
	display "`snp'"
}

*****************************************************
*first 4:
*rs10095724 rs10116891 rs10133279 rs10182458
*****************************************************

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
keep `snp'_c `snp'_m `snp'_f offspring_id  
save "scratch/childbmi_singlesnps/`snp'_dosage.dta", replace
restore
}

*/

*size of each: 572K
*size of mi wide imputed data: 390M

***************************************************************************************************************************

*ADULTBMI

clear all
set maxvar 30000
sysdir set PLUS "N:/data/durable/people/Mandy/statafiles/ado/plus"
cd "N:/data/durable/projects/mh_bmi/"

*****************************************************************************
*export chunked, sorted lists of the ones you HAVEN'T dropped for being palindromic:
use data/moba_gwas_snplevelinfo_adultbmi.dta, clear
keep SNP
sort SNP

*outsheet using scratch/adultbmi_sorted_snplist_test.txt, nonames noquote replace
*keep if strpos(SNP,"rs1")>0

*for simplicity, make 5 mini-lists, based on first number:
foreach num of numlist 0 2 4 6 8 {
    local nplus1=`num'+1
	display `nplus1'
preserve
keep if strpos(SNP,"rs`num'")>0 | strpos(SNP,"rs`nplus1'")>0
list
outsheet using scratch/adultbmi_sorted_snplist_part`num'.txt, nonames noquote replace
restore
}
*IN EACH, IN NOTEPAD++, WILL NEED TO MANUALLY USE LINE OPERATIONS-->JOIN LINES TO GET ON ONE LINE. CAN THEN COPY AND PASTE INTO RELEVANT SHELL SCRIPT
clear
*****************************************************************************


*start with dosage data - trim to single snps in turn:
use "scratch/adultbmi_justsnps.dta", clear

*get the list of snps for child bmi:
merge 1:1 _n using data/moba_gwas_snplevelinfo_adultbmi.dta
drop _merge
list SNP harm_beta_gwas palindrome in 1/928

*save the order of snps in a global for later:

levelsof SNP, local(adultbmi_snps)
*change every " to a space (34 to 32)
local adultbmi_snplist_temp: subinstr local adultbmi_snps "char(34)" "char(32)", all
global adultbmi_snplist: subinstr local adultbmi_snplist_temp "rs" " rs", all
display $adultbmi_snplist

*this works, but gives leading spaces

*check it works in a loop:
foreach snp in $adultbmi_snplist {
	*display "`snp'"
	local snp=ltrim("`snp'")
	display "`snp'"
}

*****************************************************
*first 4:
*rs10007906 rs10009336 rs1000940 rs10058464
*****************************************************

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

*now section off  each snp in turn:
*now section off  each snp in turn
*get rid of it with ltrim
foreach snp in $adultbmi_snplist {
display "`snp'"
local snp=ltrim("`snp'")
preserve
keep `snp'_c `snp'_m `snp'_f offspring_id  
save "scratch/adultbmi_singlesnps/`snp'_dosage.dta", replace
restore
}
