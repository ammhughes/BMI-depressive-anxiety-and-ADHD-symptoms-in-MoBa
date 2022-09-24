**************************************************************************
*ANALYSIS MODELS WITH IMPUTED DATA: TABLES 1 AND 2. ALL OUTCOMES.
**************************************************************************

*ENTRY POINT: TABLE 1

clear all
cd "N:\durable\projects\mh_bmi\fullsample"
set maxvar 20000
sysdir
sysdir set PLUS "N:\durable\people\Mandy\Statafiles\ado\plus"

use "data/imp_x100_17July2022.dta", clear
count

*for misum, data needs to be in flong format.
***ssc install misum***
mi convert flong

********************************************************************************
*TABLE 1:

capture log close
log using "N:\durable\projects\mh_bmi\fullsample\output\post_imp_descriptives.log", replace

*NB: the variable ADHD_8yr is kind of a full-scale ADHD scale, but it is NOT appropriate for descriptives.
*it was generated like this:
*mi passive: gen ADHD_8yr=z_inADHD_8yr+z_hyADHD_8yr
*and was then standardized again to make z_ADHD_8yr (the one used in analysis).
*but because it's the sum of two things which were standardized, it's on a different scale to everything else.
mi passive: gen descr_ADHD_8yr=inADHD_8yr+hyADHD_8yr

*Means for contin ones
foreach var in c_yob MORS_ALDER FARS_ALDER mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF mothers_bmi_Q1 fathers_bmi_Q1QF moffspring_bmi_Q8y MFQ_8yr SCARED_8yr descr_ADHD_8yr inADHD_8yr hyADHD_8yr {
*change display format
format `var' %12.2fc
misum `var' 
}

foreach var in descr_ADHD_8yr  {
*change display format
format `var' %12.2fc
misum `var' 
}

*then change back to a sensible format?
*mi convert wide
*would need to save it first. easier to reopen as mi wide:
use "data/imp_x100_17July2022.dta", clear

*Categ ones:
*no need to do KJONN PARITET_5 SIVST_2: non missingness so will be the same as complete-case
foreach var in AA1124_v2_flipped f_educ_Q1QF_v2_flipped mumsmokes_Q1 {
mi estimate: proportion `var' 
}

*for those:
foreach var in KJONN PARITET_5 SIVST_2 {
proportion `var' 
}

log close

*********************************************************************************************************************************************************************************************

*ENTRY POINT: DEP, ANXIETY AND ADHD

clear all
cd "N:\durable\projects\mh_bmi\fullsample"
set maxvar 20000
sysdir
sysdir set PLUS "N:\durable\people\Mandy\Statafiles\ado\plus"

use "data/imp_x100_17July2022.dta", clear

pwd
count

*Notes:
*now adjusting for the genetic covariates also in non-genetic models, for comparability
*adjusting for genotyping center and chip but not batch to avoid collider bias issue
*using RESCALED BMI as per 5kg/m2 to have less stupidly small coefficients

*TABLE 2/PLOTS
capture erase "output\tables\imputed_dep_anx_adhd_adultbmi.txt"
capture erase "output\tables\imputed_dep_anx_adhd_adultbmi.xls"
*IV models
eststo clear
foreach var in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr {
*adultbmi
*non-genetic
eststo adultbmi_`var'_pheno: mi estimate, post cmdok: reg `var' rescaled_moffspring_bmi_Q8 rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF KJONN c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1 i.dadsmokes_Q1QF i.PARITET_5 /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*classic MR
eststo adultbmi_`var'_mr: mi estimate, post cmdok:  ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2 m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*WFMR
eststo adultbmi_`var'_wf: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid)
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
}
estout using "output\tables\imputed_dep_anx_adhd_adultbmi.xls", cells ("b(fmt(2)) ci_l(fmt(2)) ci_u(fmt(2)) se(fmt(3)) p(fmt(3))") keep(*bmi*) replace title(BMI to age 8 outcomes) note("comparing non-genetic, classic MR and within-family MR models.")


*for reviewer's response, August 2022: give results of non-genetic models without adjustment for the genetic covariates.
eststo clear
foreach var in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr {
*adultbmi
*non-genetic
eststo adultbmi_`var'_pheno: mi estimate, post cmdok: reg `var' rescaled_moffspring_bmi_Q8 rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF KJONN c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1 i.dadsmokes_Q1QF i.PARITET_5, vce(cluster fid) 
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
}

*for the next bit, log the results:
capture log close
log using "N:\durable\projects\mh_bmi\fullsample\output\nonlinear_tests.log", replace

*in phenotypic model: check for nonlinearity of effects.
/*this doesn't work, mi estimate only uses complete-case data: mi xeq: xtile quintiles_moffspring_bmi_Q8 = rescaled_moffspring_bmi_Q8, nq(5)*/
*alternative:
*need to get cutpoints from the m=0 data to avoid supervarying issues
*these come from the complete_case script:
mi passive: gen quintiles_moffspring_bmi_Q8=0
mi passive: replace quintiles_moffspring_bmi_Q8=1 if rescaled_moffspring_bmi_Q8 >0 & rescaled_moffspring_bmi_Q8<=2.93
mi passive: replace quintiles_moffspring_bmi_Q8=2 if rescaled_moffspring_bmi_Q8 >2.93 & rescaled_moffspring_bmi_Q8<=3.10
mi passive: replace quintiles_moffspring_bmi_Q8=3 if rescaled_moffspring_bmi_Q8 >3.10 & rescaled_moffspring_bmi_Q8<=3.29
mi passive: replace quintiles_moffspring_bmi_Q8=4 if rescaled_moffspring_bmi_Q8 >3.29 & rescaled_moffspring_bmi_Q8<=3.55
mi passive: replace quintiles_moffspring_bmi_Q8=5 if rescaled_moffspring_bmi_Q8 >3.55

eststo clear
foreach var in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr z_inADHD_8yr z_hyADHD_8yr {
mi estimate, post cmdok: reg `var' ib3.quintiles_moffspring_bmi_Q8 KJONN c_yob rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1 i.dadsmokes_Q1QF i.PARITET_5,  vce(cluster fid) 
}
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar

*squared term!
mi passive: gen sq_resc_moffspring_bmi_Q8= rescaled_moffspring_bmi_Q8*rescaled_moffspring_bmi_Q8

eststo clear
foreach var in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr z_inADHD_8yr z_hyADHD_8yr {
mi estimate, post cmdok: reg `var' rescaled_moffspring_bmi_Q8 sq_resc_moffspring_bmi_Q8 KJONN c_yob rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1 i.dadsmokes_Q1QF i.PARITET_5,  vce(cluster fid) 
}
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar

log close

**************************
*ADHD subscales

*IV models: save for export:
capture erase "output\tables\imputed_adhd_subscales_adultbmi.txt"
capture erase "output\tables\imputed_adhd_subscales_adultbmi.xls"
eststo clear
*IV models
eststo clear
foreach var in z_ADHD_8yr z_inADHD_8yr z_hyADHD_8yr {
eststo adultbmi_`var'_pheno: mi estimate, post cmdok: reg `var' rescaled_moffspring_bmi_Q8 rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF KJONN c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1 i.dadsmokes_Q1QF i.PARITET_5 /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*classic MR
eststo adultbmi_`var'_mr: mi estimate, post cmdok:  ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2 m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*WFMR
eststo adultbmi_`var'_wf: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid)
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
}
estout using "output\tables\imputed_adhd_subscales_adultbmi.xls", cells ("b(fmt(2)) ci_l(fmt(2)) ci_u(fmt(2)) se(fmt(3)) p(fmt(3))") keep(*bmi*) replace title(BMI to age 8 outcomes) note("comparing non-genetic, classic MR and within-family MR models.")



*************************************************************************************
*
*ENTRY POINT: HOW MUCH MORE PRECISE WITHOUT CLUSTERING OVER FAMILIES?
/*
clear all
cd "N:\durable\projects\mh_bmi\fullsample"
set maxvar 20000
sysdir
sysdir set PLUS "N:\durable\people\Mandy\Statafiles\ado\plus"
use "data/imp_x100_17July2022.dta"
*/
pwd
count

capture log close
log using "N:\durable\projects\mh_bmi\fullsample\output\nocluster.log", replace


*adultbmi
eststo clear
foreach var in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr  {
*non-genetic
eststo adultbmi_`var'_pheno: mi estimate, post cmdok: reg `var' rescaled_moffspring_bmi_Q8 rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF KJONN c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1 i.dadsmokes_Q1QF i.PARITET_5 /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(robust) 
*classic MR
eststo adultbmi_`var'_mr: mi estimate, post cmdok:  ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2 m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(robust) 
*WFMR
eststo adultbmi_`var'_wf: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(robust)
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
}

log close

*************************************************************************************
*
*ENTRY POINT: ADJUSTMENT FOR PARENTAL EDUCATION. DOES BETWEEN-FAMILY EFFECT DISAPPEAR?

/*
clear all
cd "N:\durable\projects\mh_bmi\fullsample"
set maxvar 20000
sysdir
sysdir set PLUS "N:\durable\people\Mandy\Statafiles\ado\plus"
use "data/imp_x100_17July2022.dta"
*/
pwd
count

*adultbmi
capture erase "output\tables\imputed_all_adultbmi_sep_adj.txt"
capture erase "output\tables\imputed_all_adultbmi_sep_adj.xls"
eststo clear
foreach var in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr z_inADHD_8yr z_hyADHD_8yr {
*non-genetic
eststo adultbmi_`var'_pheno: mi estimate, post cmdok: reg `var' rescaled_moffspring_bmi_Q8 rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF KJONN c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1 i.dadsmokes_Q1QF i.PARITET_5 /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*classic MR
eststo adultbmi_`var'_mr: mi estimate, post cmdok:  ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2 m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*WFMR
eststo adultbmi_`var'_wf: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid)
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
}
estout using "output\tables\imputed_all_adultbmi_sep_adj.xls", cells ("b(fmt(2)) ci_l(fmt(2)) ci_u(fmt(2)) se(fmt(3)) p(fmt(3))") keep(*bmi*) replace title(BMI to age 8 outcomes) note("comparing non-genetic, classic MR and within-family MR models.")



*childbmi, even though redundant as nothing to attenuate:
capture erase "output\tables\imputed_all_childbmi_sep_adj.txt"
capture erase "output\tables\imputed_all_childbmi_sep_adj.xls"
eststo clear
foreach var in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr z_inADHD_8yr z_hyADHD_8yr {
*childbmi
*non-genetic
eststo childbmi_`var'_pheno: mi estimate, post cmdok: reg `var' rescaled_moffspring_bmi_Q8 rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF KJONN c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1 i.dadsmokes_Q1QF i.PARITET_5 /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*classic MR
eststo childbmi_`var'_mr: mi estimate, post cmdok:  ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_childbmi_pgs) KJONN c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2 m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*WFMR
eststo childbmi_`var'_wf: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_childbmi_pgs m_z_childbmi_pgs f_z_childbmi_pgs) KJONN c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid)
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
}
estout using "output\tables\imputed_all_childbmi_sep_adj.xls", cells ("b(fmt(2)) ci_l(fmt(2)) ci_u(fmt(2)) se(fmt(3)) p(fmt(3))") keep(*bmi*) replace title(BMI to age 8 outcomes) note("comparing non-genetic, classic MR and within-family MR models.")


*version with the full gwas BMI PGS (just for adultbmi)

**************************************************************************************
*ENTRY POINT: USE THE FULL GWAS PGS

capture erase "output\tables\imputed_all_adultbmi_full_pgs.txt"
capture erase "output\tables\imputed_all_adultbmi_full_pgs.xls"

eststo clear
foreach var in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr z_inADHD_8yr z_hyADHD_8yr {
*adultbmi
*non-genetic
*classic MR
eststo adultbmi_`var'_mr: mi estimate, post cmdok:  ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_full_pgs_adultbmi) KJONN c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2 m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*WFMR
eststo adultbmi_`var'_wf: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_full_pgs_adultbmi m_z_full_pgs_adultbmi f_z_full_pgs_adultbmi) KJONN c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid)
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
}
estout using "output\tables\imputed_all_adultbmi_full_pgs.xls", cells ("b(fmt(2)) ci_l(fmt(2)) ci_u(fmt(2)) se(fmt(3)) p(fmt(3))") keep(*bmi*) replace title(BMI to age 8 outcomes) note("comparing non-genetic, classic MR and within-family MR models.")


**********************************************************************************************

*ENTRY POINT: log-transformed versions:

/*
clear all
cd "N:\durable\projects\mh_bmi\fullsample"
set maxvar 20000
sysdir
sysdir set PLUS "N:\durable\people\Mandy\Statafiles\ado\plus"
use "data/imp_x100_17July2022.dta"
*/
pwd
count

*PREP FOR LOG-TRANSFORMED MODELS:
*need here to recreate log-transformed ADHD full scale - bit of a fudge though. since not trnasforming root variables, contribution could be different and so not directly comparable to main analyses
mi passive: gen inter_ADHD_8yr=inADHD_8yr + hyADHD_8yr
foreach var in MFQ_8yr SCARED_8yr inter_ADHD_8yr inADHD_8yr hyADHD_8yr {
summ `var', det
hist `var'
*log-transform, adding 1 to everything to deal with the 0s
mi passive: gen ln_`var'=ln(`var'+1)
}
*name of some still too long to store the estimates so rename
foreach varstem in ln_MFQ ln_SCARED ln_inter_ADHD ln_inADHD ln_hyADHD {
mi rename `varstem'_8yr `varstem'
}
mi rename ln_inter_ADHD ln_ADHD

*save for export:
capture erase "output\tables\imputed_dep_anx_adhd_adultbmi_logtransformed.txt"
capture erase "output\tables\imputed_dep_anx_adhd_adultbmi_logtransformed.xls"

eststo clear
foreach var in ln_MFQ ln_SCARED ln_ADHD {
*adultbmi
*non-genetic
eststo adultbmi_`var'_pheno: mi estimate, post cmdok: reg `var' rescaled_moffspring_bmi_Q8 rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF KJONN c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1 i.dadsmokes_Q1QF i.PARITET_5 /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*classic MR
eststo adultbmi_`var'_mr: mi estimate, post cmdok:  ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2 m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*WFMR
eststo adultbmi_`var'_wf: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid)
}
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar

estout using "output\tables\imputed_dep_anx_adhd_adultbmi_logtransformed.xls", cells ("b(fmt(2)) ci_l(fmt(2)) ci_u(fmt(2)) se(fmt(3)) p(fmt(3))") keep(*bmi*) replace title(BMI to age 8 outcomes) note("comparing non-genetic, classic MR and within-family MR models.")

****************************************************************************************

*ENTRY POINT: R2 AND WEAK INSTRUMENT TESTS ACROSS IMPUTATIONS

*MOVED TO A SEPARATE FILE.

***************************************************************************

*ENTRY POINT: FIGURES

/*
clear all
cd "N:\durable\projects\mh_bmi\fullsample"
set maxvar 20000
sysdir
sysdir set PLUS "N:\durable\people\Mandy\Statafiles\ado\plus"
use "data/imp_x100_17July2022.dta"
*/

*Coefplot guide: http://repec.sowi.unibe.ch/stata/coefplot/getting-started.html

*Label them:
label variable z_MFQ_8yr "depressive symptoms"
label variable z_SCARED_8yr "anxiety symptoms"
label variable z_ADHD_8yr "ADHD symptoms"
label variable z_inADHD_8yr "ADHD symptoms (inattention)"
label variable z_hyADHD_8yr "ADHD symptoms (hyperactivity)"
label variable rescaled_moffspring_bmi_Q8y "Child's BMI"
label variable rescaled_mothers_bmi_Q1 "Mother's BMI"
label variable rescaled_fathers_bmi_Q1 "Father's BMI"

*********************
*Repeat for each outcome, pulling graph titles from macros defined here:
global z_MFQ_8yr_title = "Depressive symptoms"
global z_SCARED_8yr_title = "Anxiety symptoms"
global z_ADHD_8yr_title = "ADHD symptoms"

foreach var in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr {
*non-genetic:
mi estimate, post cmdok: reg `var' rescaled_moffspring_bmi_Q8 rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF KJONN c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1 i.dadsmokes_Q1QF i.PARITET_5 /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
estimates store `var'_pheno
*classic MR
mi estimate, post cmdok:  ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2 m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
estimates store `var'_mr
*WFMR
mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid)
estimates store `var'_wfmr
*plot the estimates:
coefplot ///
(`var'_pheno, label (Non-genetic)) ///
(`var'_mr, label (Classic MR)) ///
(`var'_wfmr, label (Within-family MR)), ///
keep (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF) drop(_cons) xline(0) grid(none) msize(small) title("{bf:$`var'_title}", size(medlarge) pos(4) ring(0) color(black)) scheme(s1color) graphregion(margin(l=25)  color(white)) legend(rows(1) size(small)) xsc(r(-0.27 0.67)) xtick(-0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6) xlabel(-0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6) xlabel(,labsize(small)) ylabel(,labsize(small)) legend(region(lwidth(none)))
graph save output/graphs/adultbmi_`var'_coefeplot.gph, replace
}
*combine:
*grc1leg2 output/graphs/adultbmi_z_MFQ_8yr_coefeplot.gph output/graphs/adultbmi_z_SCARED_8yr_coefeplot.gph output/graphs/adultbmi_z_ADHD_8yr_coefeplot.gph, cols(1) xcommon ycommon scheme(s1color) imargin(20 5 0 0) position(6) 
grc1leg2 output/graphs/adultbmi_z_MFQ_8yr_coefeplot.gph output/graphs/adultbmi_z_SCARED_8yr_coefeplot.gph output/graphs/adultbmi_z_ADHD_8yr_coefeplot.gph, cols(1) scheme(s1color) imargin(35 40 3 0) position(6) 
graph save output/graphs/adultbmi_alloutcomes_coefeplot.gph, replace
graph export adultbmi_alloutcomes_coefeplot.tif, replace width(1200)

***************************************************
***FOR CHILD BMI PGS:

*Coefplot guide: http://repec.sowi.unibe.ch/stata/coefplot/getting-started.html

foreach var in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr {
*non-genetic:
mi estimate, post cmdok: reg `var' rescaled_moffspring_bmi_Q8 rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF KJONN c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1 i.dadsmokes_Q1QF i.PARITET_5 /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
estimates store `var'_pheno
*classic MR
mi estimate, post cmdok:  ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_childbmi_pgs) KJONN c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2 m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
estimates store `var'_mr
*WFMR
mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_childbmi_pgs m_z_childbmi_pgs f_z_childbmi_pgs) KJONN c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid)
estimates store `var'_wfmr
*plot the estimates:
coefplot ///
(`var'_pheno, label (Non-genetic)) ///
(`var'_mr, label (Classic MR)) ///
(`var'_wfmr, label (Within-family MR)), ///
keep (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF) drop(_cons) xline(0) grid(none) msize(small) title("{bf:$`var'_title}", size(medlarge) pos(4) ring(0) color(black)) scheme(s1color) graphregion(margin(l=25)  color(white)) legend(rows(1) size(small)) xsc(r(-0.27 0.67)) xtick(-0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6) xlabel(-0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6) xlabel(,labsize(small)) ylabel(,labsize(small)) legend(region(lwidth(none)))
graph save output/graphs/childbmi_`var'_coefeplot.gph, replace
}
*combine:
*grc1leg2 output/graphs/childbmi_z_MFQ_8yr_coefeplot.gph output/graphs/childbmi_z_SCARED_8yr_coefeplot.gph output/graphs/childbmi_z_ADHD_8yr_coefeplot.gph, cols(1) xcommon ycommon scheme(s1color) imargin(0 0 0 0) position(6) 
grc1leg2 output/graphs/childbmi_z_MFQ_8yr_coefeplot.gph output/graphs/childbmi_z_SCARED_8yr_coefeplot.gph output/graphs/childbmi_z_ADHD_8yr_coefeplot.gph, cols(1) scheme(s1color) imargin(35 40 3 0) position(6) 
graph save output/graphs/childbmi_alloutcomes_coefeplot.gph, replace
graph export childbmi_alloutcomes_coefeplot.tif, replace width(1200)

********************************************************************************

*ASSORTATIVE MATING CHECK:

*standardize everything for this
foreach var in rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF  mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF {
mi passive:	egen z_`var'=std(`var')
}

*mother's BMI explaining:
*father's BMI
mi estimate: reg z_rescaled_fathers_bmi_Q1QF z_rescaled_mothers_bmi_Q1 /*genetic covariates: PCs, center and chip but not batch*/ m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*father's depression
mi estimate:  reg z_fhopkins_QF z_rescaled_mothers_bmi_Q1 /*genetic covariates: PCs, center and chip but not batch*/ m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*father's 
mi estimate: reg z_fADHD_QF z_rescaled_mothers_bmi_Q1 /*genetic covariates: PCs, center and chip but not batch*/ m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 

*mother's depression explaining:
*father's BMI
mi estimate: reg z_rescaled_fathers_bmi_Q1QF z_mhopkins_Q1 /*genetic covariates: PCs, center and chip but not batch*/ m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*father's depression
mi estimate: reg z_fhopkins_QF z_mhopkins_Q1 /*genetic covariates: PCs, center and chip but not batch*/ m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*father's ADHD
mi estimate: reg z_fADHD_QF z_mhopkins_Q1 /*genetic covariates: PCs, center and chip but not batch*/ m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 

*mother's ADHD explaining:
*father's BMI:
mi estimate: reg z_rescaled_fathers_bmi_Q1QF z_mADHD_Q6 /*genetic covariates: PCs, center and chip but not batch*/ m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid)
*father's depression
mi estimate: reg fhopkins_QF mADHD_Q6 /*genetic covariates: PCs, center and chip but not batch*/ m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*father's ADHD
mi estimate: reg z_fADHD_QF z_mADHD_Q6 /*genetic covariates: PCs, center and chip but not batch*/ m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
