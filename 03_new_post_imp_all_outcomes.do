
**************************************************************************
*ANALYSIS MODELS WITH IMPUTED DATA: TABLES 1 AND 2. ALL OUTCOMES.
**************************************************************************

*ENTRY POINT: TABLE 1

cd "N:\data\durable\projects\mh_bmi\"
set maxvar 20000
sysdir
sysdir set PLUS "N:\data\durable\people\Mandy\Statafiles\ado\plus"

use "scratch/scalesonly_dep_anxiety_adhd_2July.dta", clear

*UPDATE 12 MAY: RESTRICTED TO PEOPLE IN UPDATED SV_INFO FILE (i.e. exclude the new withdrawls) 
count
*26370

*DO NOT DROP RELATEDS: KING FIDs ARE FOR MORE THAN NUCLEAR FAMILIES, ALREADY GROUP FAMILIES OF SIBLINGS IN PARENTS GEN.
*SO CLUSTERING ON THAT FOR THE SIBLINGS IN CHILD GEN ALREADY HANDLES THIS.

*reg z_MFQ_8yr rescaled_moffspring_bmi_Q8 KJONN, vce(robust) 
*reg z_MFQ_8yr rescaled_moffspring_bmi_Q8 KJONN, vce(cluster fid) 
*mi estimate, cmdok: reg z_MFQ_8yr rescaled_moffspring_bmi_Q8 KJONN, vce(cluster fid) 

********************************************************************************

*TABLE 1:


*NB: the variable ADHD_8yr is kind of a full-scale ADHD scale, but it is NOT appropriate for descriptives.
*it was generated like this:
*mi passive: gen ADHD_8yr=z_inADHD_8yr+z_hyADHD_8yr
*and was then standardized again to make z_ADHD_8yr (the one used in analysis).
*but because it's the sum of two things which were standardized, it's on a different scale to everything else.
mi passive: gen descr_ADHD_8yr=inADHD_8yr+hyADHD_8yr

*for misum, data needs to be in flong format.
***ssc install misum***
mi convert flong

*Means for contin ones
foreach var in MORS_ALDER FARS_ALDER mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF mothers_bmi_Q1 fathers_bmi_Q1QF moffspring_bmi_Q8y MFQ_8yr SCARED_8yr descr_ADHD_8yr inADHD_8yr hyADHD_8yr {
*change display format
format `var' %12.2fc
misum `var' 
}

foreach var in descr_ADHD_8yr  {
*change display format
format `var' %12.2fc
misum `var' 
}

*then change back to a sensible format:
mi convert wide

*Categ ones:
*no need to do KJONN PARITET_5 SIVST_2: non missingness so will be the same as complete-case
foreach var in AA1124_v2_flipped f_educ_Q1QF_v2_flipped mumsmokes_Q1_v2 {
mi estimate: proportion `var' 
}

*for those:
foreach var in KJONN PARITET_5 SIVST_2 {
proportion `var' 
}
*********************************************************************************************************************************************************************************************

*ENTRY POINT: DEP, ANXIETY AND ADHD

cd "N:\data\durable\projects\mh_bmi\"
set maxvar 20000
sysdir
sysdir set PLUS "N:\data\durable\people\Mandy\Statafiles\ado\plus"

use "scratch/scalesonly_dep_anxiety_adhd_2July.dta", clear

*UPDATE 12 MAY: RESTRICTED TO PEOPLE IN UPDATED SV_INFO FILE (i.e. exclude the new withdrawls) 
count
*26370

*DO NOT DROP RELATEDS: KING FIDs ARE FOR MORE THAN NUCLEAR FAMILIES, ALREADY GROUP FAMILIES OF SIBLINGS IN PARENTS GEN.
*SO CLUSTERING ON THAT FOR THE SIBLINGS IN CHILD GEN ALREADY HANDLES THIS.

*TABLE 2/PLOTS

*ADULT AND CHILD BMI PGS:
*r2 in new data: can't get from mi estimate, but look in first imputation:
*mi xeq 1:  regress moffspring_bmi_Q8 c_z_adultbmi_pgs
*mi xeq 1:  regress moffspring_bmi_Q8 c_z_childbmi_pgs

*TEST:
*mi estimate, post cmdok: ivregress 2sls z_MFQ_8yr (moffspring_bmi_Q8y mothers_bmi_Q1 fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*

/*mi estimate: omitted terms vary
    The set of omitted variables or categories is not consistent between m=1 and m=4; this is not allowed.
    To identify varying sets, you can use mi xeq to run the command on individual imputations or you can
	reissue the command with mi estimate, noisily*/
*This happened before. Last year, fixed this by dropping the redundant dummy for the largest batch. This is HARVEST for everyone.
*nb: if this is in wide format, need to specify prefixes since different variables corresponding to diff imputations. went back and did that pre-imputation, and it worked.

*using RESCALED BMI as per 5kh/m2 to have less stupidly small coefficients

*IV models: save for export:
capture erase "output\tables\dep_anxiety_ivreg_estout.txt"
capture erase "output\tables\dep_anxiety_ivreg_estout.xls"
eststo clear
*IV models
eststo clear
foreach var in z_MFQ_8yr z_SCARED_8yr {
*adultbmi
*phenotypic
eststo adultbmi_`var'_pheno: mi estimate, post cmdok: reg `var' rescaled_moffspring_bmi_Q8 KJONN rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1_v2 i.dadsmokes_Q1QF i.PARITET_5,  vce(cluster fid) 
*normal genetic -  iv with pgs
eststo adultbmi_`var'_mr: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_PC* c_batch*, vce(cluster fid)
*with parents genetic -  iv with pgs
eststo adultbmi_`var'_wf: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
*And then childbmi
*phenotypic
eststo childbmi_`var'_pheno: mi estimate, post cmdok: reg `var' rescaled_moffspring_bmi_Q8 KJONN rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1_v2 i.dadsmokes_Q1QF i.PARITET_5, vce(cluster fid)
*normal genetic -  iv with pgs
eststo childbmi_`var'_mr: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_childbmi_pgs) KJONN c_PC* c_batch*, vce(cluster fid)
*with parents genetic -  iv with pgs
eststo childbmi_`var'_wf: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_childbmi_pgs m_z_childbmi_pgs f_z_childbmi_pgs) KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
}
estout using "output\tables\dep_anxiety_ivreg_estout.xls", cells ("b(fmt(2)) ci_l(fmt(2)) ci_u(fmt(2)) se(fmt(3)) p(fmt(3))") keep(*bmi*) replace title(BMI to age 8 outcomes) note("comparing standard and within-family IV models.")


*in phenotypic model: check for nonlinearity of effects.
mi xeq: xtile deciles_moffspring_bmi_Q8 = rescaled_moffspring_bmi_Q8, nq(10)
mi xeq: xtile quintiles_moffspring_bmi_Q8 = rescaled_moffspring_bmi_Q8, nq(5)

eststo clear
foreach var in z_MFQ_8yr z_SCARED_8yr {
mi estimate, post cmdok: reg `var' i.deciles_moffspring_bmi_Q8 KJONN rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1_v2 i.dadsmokes_Q1QF i.PARITET_5,  vce(cluster fid) 
}
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar

eststo clear
foreach var in z_MFQ_8yr z_SCARED_8yr {
mi estimate, post cmdok: reg `var' ib3.quintiles_moffspring_bmi_Q8 KJONN rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1_v2 i.dadsmokes_Q1QF i.PARITET_5,  vce(cluster fid) 
}
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar

**************************
*ADHD and subscales

*IV models: save for export:
capture erase "output\tables\adhd_subscales_ivreg_estout.txt"
capture erase "output\tables\adhd_subscales_ivreg_estout.xls"
eststo clear
*IV models
eststo clear
foreach var in z_ADHD_8yr z_inADHD_8yr z_hyADHD_8yr {
*adultbmi
*phenotypic
eststo adultbmi_`var'_pheno: mi estimate, post cmdok: reg `var' rescaled_moffspring_bmi_Q8 KJONN rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1_v2 i.dadsmokes_Q1QF i.PARITET_5, vce(cluster fid) 
*normal genetic -  iv with pgs
eststo adultbmi_`var'_mr: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_PC* c_batch*, vce(cluster fid)
*with parents genetic -  iv with pgs
eststo adultbmi_`var'_wf: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
*And then childbmi
*phenotypic
eststo childbmi_`var'_pheno: mi estimate, post cmdok: reg `var' rescaled_moffspring_bmi_Q8 KJONN rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1_v2 i.dadsmokes_Q1QF i.PARITET_5, vce(cluster fid) 
*normal genetic -  iv with pgs
eststo childbmi_`var'_mr: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_childbmi_pgs) KJONN c_PC* c_batch*, vce(cluster fid)
*with parents genetic -  iv with pgs
eststo childbmi_`var'_wf: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_childbmi_pgs m_z_childbmi_pgs f_z_childbmi_pgs) KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
}
estout using "output\tables\adhd_subscales_ivreg_estout.xls", cells ("b(fmt(2)) ci_l(fmt(2)) ci_u(fmt(2)) se(fmt(3)) p(fmt(3))") keep(*bmi*) replace title(BMI to age 8 outcomes) note("comparing standard and within-family IV models.")


**********************************************************************************************
*ENTRY POINT: log-transformed versions:

cd "N:\data\durable\projects\mh_bmi\"
set maxvar 20000
sysdir
sysdir set PLUS "N:\data\durable\people\Mandy\Statafiles\ado\plus"

use "scratch/scalesonly_dep_anxiety_adhd_2July.dta", clear


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
capture erase "output\tables\log_transformed_estout.txt"
capture erase "output\tables\log_transformed_estout.xls"
eststo clear

foreach var in ln_MFQ ln_SCARED ln_ADHD ln_inADHD ln_hyADHD {
*adultbmi
*phenotypic
eststo adultbmi_`var'_pheno: mi estimate, post cmdok: reg `var' rescaled_moffspring_bmi_Q8 KJONN rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1_v2 i.dadsmokes_Q1QF i.PARITET_5, vce(cluster fid) 
*normal genetic -  iv with pgs
eststo adultbmi_`var'_mr: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_PC* c_batch*, vce(cluster fid)
*with parents genetic -  iv with pgs
eststo adultbmi_`var'_wf: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
*And then childbmi
*phenotypic
eststo childbmi_`var'_pheno: mi estimate, post cmdok: reg `var' rescaled_moffspring_bmi_Q8 KJONN rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1_v2 i.dadsmokes_Q1QF i.PARITET_5, vce(cluster fid) 
*normal genetic -  iv with pgs
eststo childbmi_`var'_mr: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_childbmi_pgs) KJONN c_PC* c_batch*, vce(cluster fid)
*with parents genetic -  iv with pgs
eststo childbmi_`var'_wf: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_childbmi_pgs m_z_childbmi_pgs f_z_childbmi_pgs) KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
}
estout using "output\tables\log_transformed_estout.xls", cells ("b(fmt(2)) ci_l(fmt(2)) ci_u(fmt(2)) se(fmt(3)) p(fmt(3))") keep(*bmi*) replace title(BMI to age 8 outcomes) note("comparing standard and within-family IV models.")

*************************************************************************************

*ENTRY POINT: alternartive wfmr specification with single instruments:

cd "N:\data\durable\projects\mh_bmi\"
set maxvar 20000
sysdir
sysdir set PLUS "N:\data\durable\people\Mandy\Statafiles\ado\plus"

use "scratch/scalesonly_dep_anxiety_adhd_2July.dta", clear


capture erase "output\tables\triplicate_estimates_estout.txt"
capture erase "output\tables\triplicate_estimates_estout.xls"

eststo clear
foreach var in z_MFQ_8yr z_SCARED_8yr ADHD_8yr inADHD_8yr hyADHD_8yr {
*with parents genetic -  iv with pgs
eststo c_adultbmi_`var'_wf: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y = c_z_adultbmi_pgs) m_z_adultbmi_pgs f_z_adultbmi_pgs KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)
*mother
eststo m_adultbmi_`var'_wf: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_mothers_bmi_Q1 = m_z_adultbmi_pgs) c_z_adultbmi_pgs f_z_adultbmi_pgs KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)
*father
eststo f_adultbmi_`var'_wf: mi estimate, post cmdok: ivregress 2sls `var' (rescaled_fathers_bmi_Q1QF = f_z_adultbmi_pgs) c_z_adultbmi_pgs m_z_adultbmi_pgs KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
}

estout using "output\tables\triplicate_estimates_estout.xls", cells ("b(fmt(2)) ci_l(fmt(2)) ci_u(fmt(2)) se(fmt(3)) p(fmt(3))") keep(*bmi*) replace title(BMI to age 8 outcomes) note("comparing standard and within-family IV models.")



mi estimate, post cmdok: ivregress 2sls z_MFQ_8yr (rescaled_moffspring_bmi_Q8y = c_z_adultbmi_pgs) m_z_adultbmi_pgs f_z_adultbmi_pgs KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)

mi estimate, post cmdok: ivregress 2sls z_MFQ_8yr ///
(rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1  rescaled_fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) ///
KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)


mi estimate, post cmdok: ivregress 2sls z_MFQ_8yr (rescaled_mothers_bmi_Q1 = m_z_adultbmi_pgs) c_z_adultbmi_pgs f_z_adultbmi_pgs KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)

mi estimate, post cmdok: ivregress 2sls z_MFQ_8yr ///
(rescaled_mothers_bmi_Q1 rescaled_moffspring_bmi_Q8y rescaled_fathers_bmi_Q1QF = m_z_adultbmi_pgs c_z_adultbmi_pgs  f_z_adultbmi_pgs) ///
KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)



****************************************************************************************
*ENTRY POINT: R2 AND WEAK INSTRUMENT TESTS ACROSS IMPUTATIONS

cd "N:\data\durable\projects\mh_bmi\"
set maxvar 20000
sysdir
sysdir set PLUS "N:\data\durable\people\Mandy\Statafiles\ado\plus"

use "scratch/scalesonly_dep_anxiety_adhd_2July.dta", clear

*UPDATE 12 MAY: RESTRICT TO PEOPLE IN UPDATED SV_INFO FILE (i.e. exclude the new withdrawls) 
merge m:1 PREG_ID_2306 using "N:\data\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\Statafiles\PDB2306_SV_INFO_v12_UPDATED_April_2021"
*the 51 not matched are the ones that need to be dropped as new exclusions.
drop if _merge!=3
*51 dropped
drop _merge

*correlation of the two scores:
foreach p in c m f {
pwcorr `p'_z_adultbmi_pgs `p'_z_childbmi_pgs
}

*adult BMI PGS:

*how much of BMI does the PGS explain?
gen r2_c_z_adultbmi_pgs=.
gen r2_m_z_adultbmi_pgs=.
gen r2_f_z_adultbmi_pgs=.

forvalues m = 1/100 {
reg _`m'_rescaled_moffspring_bmi_Q8 c_z_adultbmi_pgs, vce(cluster fid)
replace r2_c_z_adultbmi_pgs=e(r2) in `m'
reg _`m'_rescaled_mothers_bmi_Q1 m_z_adultbmi_pgs, vce(cluster fid)
replace r2_m_z_adultbmi_pgs=e(r2) in `m'
reg _`m'_rescaled_fathers_bmi_Q1QF f_z_adultbmi_pgs, vce(cluster fid)
replace r2_f_z_adultbmi_pgs=e(r2) in `m'
}
mean r2_c_z_adultbmi_pgs
mean r2_m_z_adultbmi_pgs
mean r2_f_z_adultbmi_pgs

*after mi estimate, ivreg2 doesn't give you the postestimation parts 
*Manual strategy of ivreg plus estat firststage/estat endogenous doesn't work either. 

*can just do one outcome because the first stage is the same!
foreach var in z_MFQ_8yr {
capture gen `var'_c_adultbmi_cdf=.
capture gen `var'_m_adultbmi_cdf=.
capture gen `var'_f_adultbmi_cdf=.
forvalues m = 1/100 {
	*child
ivregress 2sls _`m'_`var' (_`m'_rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs), vce(cluster fid) 
estat firststage
*r(singleresults) 
mat fstat = r(singleresults)
scalar fs = fstat[1,4] 
display fs
replace `var'_c_adultbmi_cdf=fs in `m'
*mother
ivregress 2sls _`m'_`var' (_`m'_rescaled_mothers_bmi_Q1=m_z_adultbmi_pgs), vce(cluster fid) 
estat firststage
*r(singleresults) 
mat fstat = r(singleresults)
scalar fs = fstat[1,4] 
display fs
replace `var'_m_adultbmi_cdf=fs in `m'
*father
ivregress 2sls _`m'_`var' (_`m'_rescaled_fathers_bmi_Q1=f_z_adultbmi_pgs), vce(cluster fid) 
estat firststage
*r(singleresults) 
mat fstat = r(singleresults)
scalar fs = fstat[1,4] 
display fs
replace `var'_f_adultbmi_cdf=fs in `m'
}
}

foreach var in z_MFQ_8yr  {
mean `var'_c_adultbmi_cdf
mean `var'_m_adultbmi_cdf
mean `var'_f_adultbmi_cdf
}

************************************************
*recalled BMI PGS:
*how much of BMI does the PGS explain?
gen r2_c_z_childbmi_pgs=.
gen r2_m_z_childbmi_pgs=.
gen r2_f_z_childbmi_pgs=.

forvalues m = 1/100 {
reg _`m'_rescaled_moffspring_bmi_Q8 c_z_childbmi_pgs, vce(cluster fid)
replace r2_c_z_childbmi_pgs=e(r2) in `m'
reg _`m'_rescaled_mothers_bmi_Q1 m_z_childbmi_pgs, vce(cluster fid)
replace r2_m_z_childbmi_pgs=e(r2) in `m'
reg _`m'_rescaled_fathers_bmi_Q1QF f_z_childbmi_pgs, vce(cluster fid)
replace r2_f_z_childbmi_pgs=e(r2) in `m'
}
mean r2_c_z_childbmi_pgs
mean r2_m_z_childbmi_pgs
mean r2_f_z_childbmi_pgs

*after mi estimate, ivreg2 doesn't give you the postestimation parts 
*Manual strategy of ivreg plus estat firststage/estat endogenous doesn't work either. 

*use the below to average across imps (in mi wide format!)
*childbmi pgs
*can just do z_ADHD_8yr, not all outcomes, because the first stage is the same!
foreach var in z_MFQ_8yr {
capture gen `var'_c_childbmi_cdf=.
capture gen `var'_m_childbmi_cdf=.
capture gen `var'_f_childbmi_cdf=.
forvalues m = 1/100 {
	*child
ivregress 2sls _`m'_`var' (_`m'_rescaled_moffspring_bmi_Q8=c_z_childbmi_pgs), vce(cluster fid) 
estat firststage
*r(singleresults) 
mat fstat = r(singleresults)
scalar fs = fstat[1,4] 
display fs
replace `var'_c_childbmi_cdf=fs in `m'
*mother
ivregress 2sls _`m'_`var' (_`m'_rescaled_mothers_bmi_Q1=m_z_childbmi_pgs), vce(cluster fid) 
estat firststage
*r(singleresults) 
mat fstat = r(singleresults)
scalar fs = fstat[1,4] 
display fs
replace `var'_m_childbmi_cdf=fs in `m'
*father
ivregress 2sls _`m'_`var' (_`m'_rescaled_fathers_bmi_Q1=f_z_childbmi_pgs), vce(cluster fid) 
estat firststage
*r(singleresults) 
mat fstat = r(singleresults)
scalar fs = fstat[1,4] 
display fs
replace `var'_f_childbmi_cdf=fs in `m'
}
}

foreach var in z_MFQ_8yr  {
mean `var'_c_childbmi_cdf
mean `var'_m_childbmi_cdf
mean `var'_f_childbmi_cdf
}

***************************************************************************

*ENTRY POINT: FIGURES

cd "N:\data\durable\projects\mh_bmi\"
set maxvar 20000
sysdir
sysdir set PLUS "N:\data\durable\people\Mandy\Statafiles\ado\plus"

use "scratch/scalesonly_dep_anxiety_adhd_2July.dta", clear

*UPDATE 12 MAY: RESTRICTED TO PEOPLE IN UPDATED SV_INFO FILE (i.e. exclude the new withdrawls) 
count
*26370

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
*repeated for all with different titles (hence no loop) - separate graphs then combined
mi estimate, post cmdok: reg z_MFQ_8yr rescaled_moffspring_bmi_Q8 KJONN rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1_v2 i.dadsmokes_Q1QF i.PARITET_5, vce(cluster fid)
estimates store z_MFQ_8yr_pheno
mi estimate, post cmdok: ivregress 2sls z_MFQ_8yr (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_PC* c_batch*, vce(cluster fid)
estimates store z_MFQ_8yr_mr
*with parents genetic -  iv with pgs
mi estimate, post cmdok: ivregress 2sls z_MFQ_8yr (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)
estimates store z_MFQ_8yr_wfmr
coefplot ///
(z_MFQ_8yr_pheno, label (Multivariable adjusted)) ///
(z_MFQ_8yr_mr, label (Classic MR)) ///
(z_MFQ_8yr_wfmr, label (Within-families MR)), ///
keep (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF) drop(_cons) xline(0) grid(none) title("{bf:Depressive symptoms}", size(medlarge) pos(11) color(black)) scheme(s1mono) graphregion(margin(l=25)  color(white)) legend(rows(1) size(medsmall)) xsc(r(-0.5 1.0)) xtick(-0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0) xlabel(-0.5 0.0 0.5 1.0) legend(region(lwidth(none)))
graph save output/graphs/adultbmi_z_MFQ_8yr_coefeplot.gph, replace

mi estimate, post cmdok: reg z_SCARED_8yr rescaled_moffspring_bmi_Q8 KJONN rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1_v2 i.dadsmokes_Q1QF i.PARITET_5, vce(cluster fid)
estimates store z_SCARED_8yr_pheno
mi estimate, post cmdok: ivregress 2sls z_SCARED_8yr (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_PC* c_batch*, vce(cluster fid)
estimates store z_SCARED_8yr_mr
*with parents genetic -  iv with pgs
mi estimate, post cmdok: ivregress 2sls z_SCARED_8yr (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)
estimates store z_SCARED_8yr_wfmr
coefplot ///
(z_SCARED_8yr_pheno, label (Multivariable adjusted)) ///
(z_SCARED_8yr_mr, label (Classic MR)) ///
(z_SCARED_8yr_wfmr, label (Within-families MR)), ///
keep (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF) drop(_cons) xline(0) grid(none) title("{bf:Anxiety symptoms}", size(medlarge) pos(11) color(black)) scheme(s1mono) graphregion(margin(l=25)  color(white)) legend(rows(1) size(medsmall)) xsc(r(-0.5 1.0)) xtick(-0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0) xlabel(-0.5 0.0 0.5 1.0) legend(region(lwidth(none)))
graph save output/graphs/adultbmi_z_SCARED_8yr_coefeplot.gph, replace

mi estimate, post cmdok: reg z_ADHD_8yr rescaled_moffspring_bmi_Q8 KJONN rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1_v2 i.dadsmokes_Q1QF i.PARITET_5, vce(cluster fid)
estimates store z_ADHD_8yr_pheno
mi estimate, post cmdok: ivregress 2sls z_ADHD_8yr (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_PC* c_batch*, vce(cluster fid)
estimates store z_ADHD_8yr_mr
*with parents genetic -  iv with pgs
mi estimate, post cmdok: ivregress 2sls z_ADHD_8yr (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)
estimates store z_ADHD_8yr_wfmr
coefplot ///
(z_ADHD_8yr_pheno, label (Multivariable adjusted)) ///
(z_ADHD_8yr_mr, label (Classic MR)) ///
(z_ADHD_8yr_wfmr, label (Within-families MR)), ///
keep (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF) drop(_cons) xline(0) grid(none) title("{bf:ADHD symptoms}", size(medlarge) pos(11) color(black)) scheme(s1mono) graphregion(margin(l=25)  color(white)) legend(rows(1) size(medsmall)) xsc(r(-0.5 1.0)) xtick(-0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0) xlabel(-0.5 0.0 0.5 1.0) legend(region(lwidth(none)))
graph save output/graphs/adultbmi_z_ADHD_8yr_coefeplot.gph, replace

*combine:
grc1leg2 output/graphs/adultbmi_z_MFQ_8yr_coefeplot.gph output/graphs/adultbmi_z_SCARED_8yr_coefeplot.gph output/graphs/adultbmi_z_ADHD_8yr_coefeplot.gph, cols(1) xcommon ycommon scheme(s1mono) imargin(20 5 0 0) position(6) 
graph save output/graphs/adultbmi_alloutcomes_coefeplot.gph, replace
graph export adultbmi_alloutcomes_coefeplot.tif, replace width(1200)

***************************************************
***FOR CHILD BMI PGS:

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

*repeated for all with different titles (hence no loop) - separate graphs then combined
mi estimate, post cmdok: reg z_MFQ_8yr rescaled_moffspring_bmi_Q8 KJONN rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1_v2 i.dadsmokes_Q1QF i.PARITET_5, vce(cluster fid)
estimates store z_MFQ_8yr_pheno
mi estimate, post cmdok: ivregress 2sls z_MFQ_8yr (rescaled_moffspring_bmi_Q8=c_z_childbmi_pgs) KJONN c_PC* c_batch*, vce(cluster fid)
estimates store z_MFQ_8yr_mr
*with parents genetic -  iv with pgs
mi estimate, post cmdok: ivregress 2sls z_MFQ_8yr (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_childbmi_pgs m_z_childbmi_pgs f_z_childbmi_pgs) KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)
estimates store z_MFQ_8yr_wfmr
coefplot ///
(z_MFQ_8yr_pheno, label (Multivariable adjusted)) ///
(z_MFQ_8yr_mr, label (Classic MR)) ///
(z_MFQ_8yr_wfmr, label (Within-families MR)), ///
keep (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF) drop(_cons) xline(0) grid(none) title("{bf:Depressive symptoms}", size(medlarge) pos(11) color(black)) scheme(s1mono) graphregion(margin(l=25)  color(white)) legend(rows(1) size(medsmall)) xsc(r(-0.5 1.0)) xtick(-0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0) xlabel(-0.5 0.0 0.5 1.0) legend(region(lwidth(none)))
graph save output/graphs/childbmi_z_MFQ_8yr_coefeplot.gph, replace

mi estimate, post cmdok: reg z_SCARED_8yr rescaled_moffspring_bmi_Q8 KJONN rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1_v2 i.dadsmokes_Q1QF i.PARITET_5, vce(cluster fid)
estimates store z_SCARED_8yr_pheno
mi estimate, post cmdok: ivregress 2sls z_SCARED_8yr (rescaled_moffspring_bmi_Q8=c_z_childbmi_pgs) KJONN c_PC* c_batch*, vce(cluster fid)
estimates store z_SCARED_8yr_mr
*with parents genetic -  iv with pgs
mi estimate, post cmdok: ivregress 2sls z_SCARED_8yr (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_childbmi_pgs m_z_childbmi_pgs f_z_childbmi_pgs) KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)
estimates store z_SCARED_8yr_wfmr
coefplot ///
(z_SCARED_8yr_pheno, label (Multivariable adjusted)) ///
(z_SCARED_8yr_mr, label (Classic MR)) ///
(z_SCARED_8yr_wfmr, label (Within-families MR)), ///
keep (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF) drop(_cons) xline(0) grid(none) title("{bf:Anxiety symptoms}", size(medlarge) pos(11) color(black)) scheme(s1mono) graphregion(margin(l=25)  color(white)) legend(rows(1) size(medsmall)) xsc(r(-0.5 1.0)) xtick(-0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0) xlabel(-0.5 0.0 0.5 1.0) legend(region(lwidth(none)))
graph save output/graphs/childbmi_z_SCARED_8yr_coefeplot.gph, replace

mi estimate, post cmdok: reg z_ADHD_8yr rescaled_moffspring_bmi_Q8 KJONN rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1_v2 i.dadsmokes_Q1QF i.PARITET_5, vce(cluster fid)
estimates store z_ADHD_8yr_pheno
mi estimate, post cmdok: ivregress 2sls z_ADHD_8yr (rescaled_moffspring_bmi_Q8=c_z_childbmi_pgs) KJONN c_PC* c_batch*, vce(cluster fid)
estimates store z_ADHD_8yr_mr
*with parents genetic -  iv with pgs
mi estimate, post cmdok: ivregress 2sls z_ADHD_8yr (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_childbmi_pgs m_z_childbmi_pgs f_z_childbmi_pgs) KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)
estimates store z_ADHD_8yr_wfmr
coefplot ///
(z_ADHD_8yr_pheno, label (Multivariable adjusted)) ///
(z_ADHD_8yr_mr, label (Classic MR)) ///
(z_ADHD_8yr_wfmr, label (Within-families MR)), ///
keep (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF) drop(_cons) xline(0) grid(none) title("{bf:ADHD symptoms}", size(medlarge) pos(11) color(black)) scheme(s1mono) graphregion(margin(l=25)  color(white)) legend(rows(1) size(medsmall)) xsc(r(-0.5 1.0)) xtick(-0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0) xlabel(-0.5 0.0 0.5 1.0)  legend(region(lwidth(none)))
graph save output/graphs/childbmi_z_ADHD_8yr_coefeplot.gph, replace

*combine:
grc1leg2 output/graphs/childbmi_z_MFQ_8yr_coefeplot.gph output/graphs/childbmi_z_SCARED_8yr_coefeplot.gph output/graphs/childbmi_z_ADHD_8yr_coefeplot.gph, cols(1) xcommon ycommon scheme(s1mono) imargin(20 5 0 0) position(6) 
graph save output/graphs/childbmi_alloutcomes_coefeplot.gph, replace
graph export childbmi_alloutcomes_coefeplot.tif, replace width(1200)
******************************************************************************

/*
*AlL TOGETHER:
grc1leg2 output/graphs/adultbmi_z_MFQ_8yr_coefeplot.gph output/graphs/adultbmi_z_SCARED_8yr_coefeplot.gph output/graphs/adultbmi_z_ADHD_8yr_coefeplot.gph ///
output/graphs/childbmi_z_MFQ_8yr_coefeplot.gph output/graphs/childbmi_z_SCARED_8yr_coefeplot.gph output/graphs/childbmi_z_ADHD_8yr_coefeplot.gph, cols(2) xcommon ycommon scheme(s1mono) imargin(20 5 0 0) position(6) 
graph save output/graphs/adultbmi_childbmi_alloutcomes_coefeplot.gph, replace
graph export adultbmi_childbmi_alloutcomes_coefeplot.tif, replace width(1200)
*/

/*****************************************************************************
*looped for all - all on one graph directly
estimates clear
foreach outcome in z_MFQ_8yr z_SCARED_8yr  z_ADHD_8yr z_inADHD_8yr z_hyADHD_8yr {
reg `outcome' rescaled_moffspring_bmi_Q8 KJONN, vce(cluster fid)
estimates store `outcome'_pheno
ivregress 2sls `outcome' (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_PC* c_batch*, vce(cluster fid)
estimates store `outcome'_mr
*with parents genetic -  iv with pgs
ivregress 2sls `outcome' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*, vce(cluster fid)
estimates store `outcome'_wfmr
}
coefplot ///
(z_MFQ_8yr_pheno, label (Multivariable adjusted) graphregion(margin(l=30)) )  ///
(z_MFQ_8yr_mr, label (Classic MR) graphregion(margin(l=30)) ) ///
(z_MFQ_8yr_wfmr, label (Within-families MR) graphregion(margin(l=30))), bylabel ({bf:Depressive symptoms}) xsc(r(-0.5 0.8))  || /// 
(z_SCARED_8yr_pheno, label (Multivariable adjusted)graphregion(margin(l=30))) ///
(z_SCARED_8yr_mr, label (Classic MR) graphregion(margin(l=30))) ///
(z_SCARED_8yr_wfmr, label (Within-families MR) graphregion(margin(l=30))), bylabel ({bf:Anxiety symptoms}) xsc(r(-0.5 0.8)) || ///
(z_ADHD_8yr_pheno, label (Multivariable adjusted) graphregion(margin(l=30))) ///
(z_ADHD_8yr_mr, label (Classic MR)graphregion(margin(l=30)) ) ///
(z_ADHD_8yr_wfmr, label (Within-families MR) graphregion(margin(l=30))), bylabel ({bf:ADHD symptoms}) xsc(r(-0.5 0.8)) ///
keep (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF) drop(_cons) xline(0) coeflabels(, wrap(20)) scheme(s1mono) grid(none) graphregion(margin(l=30) color(white))  ///
subtitle(, size(medium) margin(vsmall) justification(left) color(black) bcolor(white) bmargin(none) position(11)) ///
byopts(cols(1)) legend(rows(1))
*/

