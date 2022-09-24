
*JUST THE INSIDE PART OF THE LOOP - RELEVANT SNP BEING SUBBED IN FROM THE SUPERSHELL

set processors 8

cd "/cluster/p/p471/cluster/projects/mh_bmi/fullsample/snp-level_associations"
sysdir set PLUS "/cluster/p/p471/cluster/people/Mandystatafiles/ado/plus"
set maxvar 20000
sysdir

args snp
local snp="`snp'"
display "`snp'"

*calculate snp-level associations and export the results:
*start with mini dosage file
use "childbmi_singlesnps/dosage/`snp'_dosage.dta", clear
count
*merge onto imputed data, using child's genetic id variables:
merge 1:1 c_iid c_fid using "/cluster/p/p471/cluster/projects/mh_bmi/fullsample/imputation/imp_x100_30Jun2022.dta"
*need to drop those annoying people not in the analytic sample because unknown sex etc - now _merge==1, because :
drop if _merge==1
drop _merge
gen SNP="`snp'"
*for each outcome:
foreach out in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr z_inADHD_8yr z_hyADHD_8yr {
count if `out'!=.
capture drop mr_Beta_`out' mr_seBeta_`out' wfmr_Beta_`out' wfmr_seBeta_`out'
*generate empty variables for SNP-level associations:
*Classic MR
gen mr_Beta_`out'=.
gen mr_seBeta_`out'=.
mi estimate, post: regress `out' `snp'_c KJONN c_yob /* PCs */ c_PC* m_PC* f_PC* ///
/*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1 m_genotyping_center2 f_genotyping_center1  f_genotyping_center2 ///
/*genotyping chip: ommitted (largest) group is 2, 1, 1*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6 ///
, vce (cluster fid)
*store outcome in first line
replace mr_Beta_`out'=_b[`snp'] in 1
replace mr_seBeta_`out'=_se[`snp'] in 1
*WFMR
gen wfmr_Beta_`out'=.
gen wfmr_seBeta_`out'=.
mi estimate, post: regress `out' `snp'_c `snp'_m `snp'_f KJONN c_yob /* PCs */ c_PC* m_PC* f_PC* ///
/*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1 m_genotyping_center2 f_genotyping_center1  f_genotyping_center2 ///
/*genotyping chip: ommitted (largest) group is 2, 1, 1*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6 ///
, vce (cluster fid)
*store outcome in first line
replace wfmr_Beta_`out'=_b[`snp'_c] in 1
replace wfmr_seBeta_`out'=_se[`snp'_c] in 1
}
*keep just the first line with associations in
keep in 1
*and just keep the associations and SEs:
keep SNP *Beta* *seBeta*
save "childbmi_singlesnps/associations/`snp'_outcome_assoc.dta", replace