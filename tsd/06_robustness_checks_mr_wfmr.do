
*ADULTBMI

*ROBUSTNESS CHECKS BASED ON SNP-OUTCOME ASSOCIATIONS CALCULATED ON COLOSSUS

*did this in parallel on colossus because otherwise it would have taken weeks (approx 1hr15 for one SNP with 5 x outcomes and 2 x genetic models)

*transferred the sets of tiny files holding just the SNP-outcome associations back over via WinSCP

clear all
set maxvar 30000
sysdir set PLUS "N:/durable/people/Mandy/statafiles/ado/plus"
cd "N:/durable/projects/mh_bmi/fullsample"

*get the list of snps for child bmi:
use "data/snplevel_info/moba_gwas_snplevelinfo_adultbmi.dta"
list SNP harm_beta_gwas in 1/954

levelsof SNP, local(adultbmi_snps) clean
global adultbmi_snplist="`adultbmi_snps'"
clear


*APPEND SNP-OUTCOME ASSOCIATIONS:

*start with the first in the list:
use "scratch/adultbmi_singlesnps/associations/rs10007906_outcome_assoc.dta", clear

*append the rest:
foreach snp in $adultbmi_snplist {
append using "scratch/adultbmi_singlesnps/associations/`snp'_outcome_assoc.dta"
}
*drop the duplicate of the first
duplicates drop
count
*954 SNPs
save "scratch/adultbmi_singlesnps/associations/adultbmi_all_snps_outcome_assoc.dta", replace

************************************************************************************************************

**MERGE BACK ON THE SNP-LEVEL INFO INCL EXTERNAL BETAS, MATCHING ON SNP

clear all
set maxvar 30000
cd "N:/durable/projects/mh_bmi/fullsample"

use "scratch/adultbmi_singlesnps/associations/adultbmi_all_snps_outcome_assoc.dta", clear

/*not in this version: trim off the un-asked for space whihc has been added: replace SNP = trim(SNP)*/
merge 1:1 SNP using "data/snplevel_info/moba_gwas_snplevelinfo_adultbmi.dta"
drop _merge
list SNP harm_beta_gwas in 1/954

*make vars to save the paramaters from the robustness checks
gen out = ""
gen ivw = .
gen ivw_se = .
gen ivw_p = .
gen egger_slope = .
gen egger_slope_se = .
gen egger_slope_p = .
gen egger_cons = .
gen egger_cons_se = .
gen egger_cons_p = .
gen median = .
gen median_se = .
gen median_p = .
gen modal = .
gen modal_se = .
gen modal_p = .

*export this:
save "scratch/forexport_adultbmi_snp_associations", replace

/*NB: for MR-modal, it says 'moremata required' even though I JUST INSTALLED THIS. It's not seeing it for some reason. 
*UPDATE: FIXED, by following instructions in moremata readme file:
*import, unzip, save in temporary folder, then:
*net from N:\durable\projects\mh_bmi\fullsample\for_mrrobust
*net install moremata, replace
*/

*for classic-MR and WFMR in turn,
*Now the actual checks, using those coefficients.
*MR robust takes the outcome first
foreach m in mr wfmr {
local k=0
foreach out in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr z_inADHD_8yr z_hyADHD_8yr {
local k=`k'+1
qui replace out = "`out'" in `k'
*IVW
mregger `m'_Beta_`out' harm_beta_gwas [aw=1/(se_gwas^2)], ivw
replace ivw = _b[harm_beta_gwas] in `k'
replace ivw_se = _se[harm_beta_gwas] in `k'
replace ivw_p = 2*normal(-abs(ivw/ivw_se)) in `k'
*egger with I^2_GX statistic
mregger `m'_Beta_`out' harm_beta_gwas [aw=1/(`m'_seBeta_`out'^2)], gxse(se_gwas)
replace egger_slope = _b[slope] in `k'
replace egger_slope_se = _se[slope] in `k'
replace egger_slope_p = 2*normal(-abs(egger_slope/egger_slope_se)) in `k'
replace egger_cons = _b[_cons] in `k'
replace egger_cons_se = _se[_cons] in `k'
replace egger_cons_p = 2*normal(-abs(egger_cons/egger_cons_se)) in `k'
*median
mrmedian `m'_Beta_`out' `m'_seBeta_`out'  harm_beta_gwas se_gwas, weighted
replace median = _b[beta] in `k'
replace median_se = _se[beta] in `k'
replace median_p = 2*normal(-abs(median/median_se)) in `k'
*mode
mrmodal `m'_Beta_`out' `m'_seBeta_`out'  harm_beta_gwas se_gwas, weighted
replace modal = _b[beta] in `k'
replace modal_se = _se[beta] in `k'
replace modal_p = 2*normal(-abs(modal/modal_se)) in `k'
}
*Export rounded results as to excel:
preserve
keep out ivw* egger* median* modal*
*actually, drop the SE variables since keeping the p's
capture drop *_se*
drop if out==""
foreach var of varlist ivw-modal_p {
gen r_`var'=round(`var',0.01)
list r_`var' `var'
drop `var'
rename r_`var' `var'
}
export excel using "output/`m'_alloutcomes_imputed_adultbmi_robustness_checks.xlsx", firstrow(variables) replace 
restore
}

*********************************************************************************************************************************************
*********************************************************************************************************************************************

*CHILDBMI

*ROBUSTNESS CHECKS BASED ON SNP-OUTCOME ASSOCIATIONS CALCULATED ON COLOSSUS

*did this for the 314 + 328 SNPs in parallel on colossus because otherwise it would have taken weeks (approx 1hr15 for one SNP with 5 x outcomes and 2 x genetic models)

*transferred the sets of tiny files holding just the SNP-outcome associations back over via WinSCP

clear all
set maxvar 30000
sysdir set PLUS "N:/durable/people/Mandy/statafiles/ado/plus"
cd "N:/durable/projects/mh_bmi/fullsample"

*get the list of snps for child bmi:
use "data/snplevel_info/moba_gwas_snplevelinfo_childbmi.dta"
list SNP harm_beta_gwas palindromic in 1/321

levelsof SNP, local(childbmi_snps) clean
global childbmi_snplist="`childbmi_snps'"
clear

*APPEND SNP-OUTCOME ASSOCIATIONS:

*start with the first in the list:
use "scratch/childbmi_singlesnps/associations/rs10095724_outcome_assoc.dta", clear

*append the rest:
foreach snp in $childbmi_snplist {
append using "scratch/childbmi_singlesnps/associations/`snp'_outcome_assoc.dta"
}
*drop the duplicate of the first
duplicates drop
count
*321 SNPs
save "scratch/childbmi_singlesnps/associations/childbmi_all_snps_outcome_assoc.dta", replace

************************************************************************************************************

**MERGE BACK ON THE SNP-LEVEL INFO INCL EXTERNAL BETAS, MATCHING ON SNP

use "scratch/childbmi_singlesnps/associations/childbmi_all_snps_outcome_assoc.dta", clear

/*not in this version: trim off the un-asked for space whihc has been added: replace SNP = trim(SNP)*/
merge 1:1 SNP using "data/snplevel_info/moba_gwas_snplevelinfo_childbmi.dta"
drop _merge
list SNP harm_beta_gwas palindromic in 1/321


*make vars to save the paramaters from the robustness checks
gen out = ""
gen ivw = .
gen ivw_se = .
gen ivw_p = .
gen egger_slope = .
gen egger_slope_se = .
gen egger_slope_p = .
gen egger_cons = .
gen egger_cons_se = .
gen egger_cons_p = .
gen median = .
gen median_se = .
gen median_p = .
gen modal = .
gen modal_se = .
gen modal_p = .

*export this:
save "scratch/forexport_childbmi_snp_associations", replace

*NB: for MR-modal, it says 'moremata required' even though I JUST INSTALLED THIS. It's not seeing it for some reason. 
*This has happened before, cannot remember how I got around that

sysdir set PLUS "N:/durable/people/Mandy/statafiles/ado/plus"

*for classic-MR and WFMR in turn,
*Now the actual checks, using those coefficients.
*MR robust takes the outcome first
foreach m in mr wfmr {
local k=0
foreach out in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr z_inADHD_8yr z_hyADHD_8yr {
local k=`k'+1
qui replace out = "`out'" in `k'
*IVW
mregger `m'_Beta_`out' harm_beta_gwas [aw=1/(se_gwas^2)], ivw
replace ivw = _b[harm_beta_gwas] in `k'
replace ivw_se = _se[harm_beta_gwas] in `k'
replace ivw_p = 2*normal(-abs(ivw/ivw_se)) in `k'
*egger with I^2_GX statistic
mregger `m'_Beta_`out' harm_beta_gwas [aw=1/(`m'_seBeta_`out'^2)], gxse(se_gwas)
replace egger_slope = _b[slope] in `k'
replace egger_slope_se = _se[slope] in `k'
replace egger_slope_p = 2*normal(-abs(egger_slope/egger_slope_se)) in `k'
replace egger_cons = _b[_cons] in `k'
replace egger_cons_se = _se[_cons] in `k'
replace egger_cons_p = 2*normal(-abs(egger_cons/egger_cons_se)) in `k'
*median
mrmedian `m'_Beta_`out' `m'_seBeta_`out'  harm_beta_gwas se_gwas, weighted
replace median = _b[beta] in `k'
replace median_se = _se[beta] in `k'
replace median_p = 2*normal(-abs(median/median_se)) in `k'
*mode
mrmodal `m'_Beta_`out' `m'_seBeta_`out'  harm_beta_gwas se_gwas, weighted
replace modal = _b[beta] in `k'
replace modal_se = _se[beta] in `k'
replace modal_p = 2*normal(-abs(modal/modal_se)) in `k'
}
*Export rounded results as to excel:
preserve
keep out ivw* egger* median* modal*
*actually, drop the SE variables since keeping the p's
capture drop *_se*
drop if out==""
foreach var of varlist ivw-modal_p {
gen r_`var'=round(`var',0.01)
list r_`var' `var'
drop `var'
rename r_`var' `var'
}
export excel using "output/`m'_alloutcomes_imputed_childbmi_robustness_checks.xlsx", firstrow(variables) replace 
restore
}
