*COMPLETE-CASE DESCRIPTIVES AND ANALYSIS 

clear all
sysdir set PLUS "N:\durable\people\Mandy\statafiles\ado\plus"
cd "N:\durable\projects\mh_bmi\fullsample"
set maxvar 20000

use "data\merged_geno_pheno.dta", clear


*********************************************************************************************

*IDENTIFY THE FINAL SAMPLE.

*FULL SAMPLE:
count
*114,030

*get rid of people not in the birth registry file:
drop if present_MBRN!=1
count
*113,691

*THE ORIGINAL NUMBER WHO HAVEN'T WITHDRAWN CONSENT.

capture drop analytic_sample
gen analytic_sample=1

*remove people not in the birth registry file:
fre present_MBRN
replace analytic_sample=0 if present_MBRN!=1
*(0 or 339 real changes made, depending on whether you've excluded them)

*DO NOT EXCLUDE BASED ON MISSING SEX INFO - can get this back from genetic files!

*exclude people not present at any of the actual questionnaires
fre present*
gen nQ_present=0
foreach q in Q1 QF Q5 Q6 Q5y Q7y Q8y {
	replace nQ_present=nQ_present+1 if present_`q'==1
}
fre nQ_present

*OK, so there's a chunk of people who never were included in a single questionnaire.
*remove them.
replace analytic_sample=0 if nQ_present==0
*(8,776 real changes made)

fre analytic_sample

*who has genetic information? 
tab complete_trio
*restrict to these:
replace analytic_sample=0 if complete_trio!=1
*(63,966 real changes made)

fre analytic_sample
*40949

*****************************************************************************

*standardize all the outcomes:
foreach var in MFQ_8yr SCARED_8yr ADHD_8yr inADHD_8yr hyADHD_8yr SCQ_8yr S_SCQ_8yr R_SCQ_8yr {
egen z_`var'=std(`var')
}

*rename these otherwise names too long to work in loop storing estimates
rename z_ADHD_8yr z_ADHD_8yr
rename z_inADHD_8yr z_inADHD_8yr
rename z_hyADHD_8yr z_hyADHD_8yr
rename z_S_SCQ_8yr z_S_SCQ_8yr 
rename z_R_SCQ_8yr z_R_SCQ_8yr

*RESCALE BMI as per 5kg/m2 to have less stupidly small coefficients:
gen rescaled_moffspring_bmi_Q8y=moffspring_bmi_Q8/5
gen rescaled_mothers_bmi_Q1=mothers_bmi_Q1/5
gen rescaled_fathers_bmi_Q1QF=fathers_bmi_Q1QF/5

******************************************************************************************
*DESCRIPTIVES FOR ALL PEOPLE IN BIRTH REGISTRY FILE, minus the new withdrawals:

count 
*113,691

foreach var in MORS_ALDER FARS_ALDER mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF mothers_bmi_Q1 fathers_bmi_Q1QF  moffspring_bmi_Q8y MFQ_8yr SCARED_8yr ADHD_8yr inADHD_8yr hyADHD_8yr  {
display "`var'"
format `var' %3.1f
summ `var', format
display ""
}

foreach var in KJONN AA1124_v2_flipped f_educ_Q1QF_v2_flipped PARITET_5 SIVST_2 mumsmokes_Q1  {
fre `var', f(1) nomiss
}

*finally, get min and max birtweight and length IN FULL SAMPLE for truncreg limits in imputation.
foreach var in  mothers_bmi_Q1 mothers_bmi_Q6 fathers_bmi_Q1QF moffspring_bmi_Q5y moffspring_bmi_Q7y moffspring_bmi_Q8y VEKT LENGDE FARS_ALDER {
	summ `var'
}

******************************************************************************************
*DESCRIPTIVES FOR ANALYTIC SAMPLE, UNIMPUTED DATA.

count if analytic_sample==1
*40,949

*VEKT
foreach var in MORS_ALDER FARS_ALDER mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF mothers_bmi_Q1 fathers_bmi_Q1QF  moffspring_bmi_Q8y MFQ_8yr SCARED_8yr ADHD_8yr inADHD_8yr hyADHD_8yr {
display "`var'"
format `var' %3.1f
summ `var' if analytic_sample==1, format
display ""
}

foreach var in KJONN AA1124_v2_flipped f_educ_Q1QF_v2_flipped PARITET_5 SIVST_2 mumsmokes_Q1  {
fre `var' if analytic_sample==1, f(1) nomiss
}

*****************************************************************************
*COMPARISON OF ANALYTIC SAMPLE WITH OTHERS IN THE BIRTH REGISTRY FILE, 
*DESCRIPTIVE TABLE FOR ALL PEOPLE IN BIRTH REGISTRY FILE

*comparison - contin vars
*sd test
foreach var in c_yob VEKT MORS_ALDER FARS_ALDER mothers_bmi_Q1 fathers_bmi_Q1QF moffspring_bmi_Q8y MFQ_8yr SCARED_8yr ADHD_8yr c_yob {
display "`var'"
sdtest `var', by(analytic_sample)
}
*all of them have unequal variance.

*unequal variance t-tests:
foreach var in c_yob VEKT MORS_ALDER FARS_ALDER mothers_bmi_Q1 fathers_bmi_Q1QF moffspring_bmi_Q8y MFQ_8yr SCARED_8yr ADHD_8yr c_yob {
display "`var'"
ttest `var', by(analytic_sample) unequal
}

*comparison - categ vars
foreach var in KJONN PARITET_5 SIVST_2 AA1124_v2_flipped f_educ_Q1QF_v2_flipped mumsmokes_Q1 dadsmokes_Q1QF  {
	tab `var' analytic_sample, chi col
}


*********************************************************************************************
*GET QUINTILE CUTPOINTS FOR LATER ANALYSIS, EVEN IN IMPUTED DATA.
*NEEDS TO BE FROM UNIMPUTED DATA OTHERWISE SUPER-VARYING ISSUES.

pctile quint_cutpoints_bmi_Q8 = rescaled_moffspring_bmi_Q8, nq(5)
tab quint_cutpoints_bmi_Q8

keep if analytic_sample==1

*compare the r2 on the normal pgs and the full one

*adultbmi
regress rescaled_moffspring_bmi_Q8 KJONN c_yob c_z_adultbmi_pgs
regress rescaled_moffspring_bmi_Q8 KJONN c_yob c_z_full_pgs_adultbmi
regress rescaled_moffspring_bmi_Q8 KJONN c_yob c_z_full_pgs_adultbmi_5e08

regress rescaled_mothers_bmi_Q1 m_z_adultbmi_pgs
regress rescaled_mothers_bmi_Q1 m_z_full_pgs_adultbmi
regress rescaled_mothers_bmi_Q1 m_z_full_pgs_adultbmi_5e08

regress rescaled_fathers_bmi_Q1  f_z_adultbmi_pgs
regress rescaled_fathers_bmi_Q1 f_z_full_pgs_adultbmi
regress rescaled_fathers_bmi_Q1 f_z_full_pgs_adultbmi_5e08


*childbmi
regress rescaled_moffspring_bmi_Q8 KJONN c_yob c_z_childbmi_pgs
regress rescaled_moffspring_bmi_Q8 KJONN c_yob c_z_full_pgs_childbmi
regress rescaled_moffspring_bmi_Q8 KJONN c_yob c_z_full_pgs_childbmi_5e08

regress rescaled_mothers_bmi_Q1 m_z_childbmi_pgs
regress rescaled_mothers_bmi_Q1 m_z_full_pgs_childbmi
regress rescaled_mothers_bmi_Q1 m_z_full_pgs_childbmi_5e08

regress rescaled_fathers_bmi_Q1 f_z_childbmi_pgs
regress rescaled_fathers_bmi_Q1 f_z_full_pgs_childbmi
regress rescaled_fathers_bmi_Q1 f_z_full_pgs_childbmi_5e08


****************************************************************************************
/*
*how much missingness?
*bmi
count if !inlist(., rescaled_moffspring_bmi_Q8, rescaled_mothers_bmi_Q1, rescaled_fathers_bmi_Q1QF)
*bmi and all phenotypic covars
count if !inlist(., rescaled_moffspring_bmi_Q8, rescaled_mothers_bmi_Q1, rescaled_fathers_bmi_Q1QF, KJONN, AA1124_v2_flipped, f_educ_Q1QF_v2_flipped, mhopkins_Q1, fhopkins_QF, mADHD_Q6, fADHD_QF, mumsmokes_Q1, dadsmokes_QF, PARITET_5)

*but of course, this doesn't take account of missingness in the outcomes.
foreach var in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr /*z_inADHD_8yr z_hyADHD_8yr*/ {
count if !inlist(., `var', rescaled_moffspring_bmi_Q8, rescaled_mothers_bmi_Q1, rescaled_fathers_bmi_Q1QF, KJONN, AA1124_v2_flipped, f_educ_Q1QF_v2_flipped, mhopkins_Q1, fhopkins_QF, mADHD_Q6, fADHD_QF, mumsmokes_Q1, dadsmokes_QF, PARITET_5)
}
* 5,158
* 5,177
* 5,174
*/

*flag complete-case including all phenotypic covariates
gen cc_sample=1 if !inlist(., rescaled_moffspring_bmi_Q8, rescaled_mothers_bmi_Q1, rescaled_fathers_bmi_Q1QF, KJONN, AA1124_v2_flipped, f_educ_Q1QF_v2_flipped, mhopkins_Q1, fhopkins_QF, mADHD_Q6, fADHD_QF, mumsmokes_Q1, dadsmokes_QF, PARITET_5, c_PC1)

* /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center: omitted (largest group is 3 for everyone)*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip: ommitted (largest) group is 2, 1, 1*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6  ///


*EXPORT RESULTS:
*ADULTBMI
capture erase "output\tables\completecase_all_adultbmi.txt"
capture erase "output\tables\completecase_all_adultbmi.txt"
eststo clear
foreach var in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr z_inADHD_8yr z_hyADHD_8yr {
*ADULTBMI PGS
*non-genetic
eststo adultbmi_`var'_pheno: reg `var' rescaled_moffspring_bmi_Q8 rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF KJONN  c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1 i.dadsmokes_QF i.PARITET_5 /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6 if cc_sample==1, vce(cluster fid) 
*classic MR
eststo adultbmi_`var'_mr:  ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2 m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6  if cc_sample==1, vce(cluster fid) 
*WFMR
eststo adultbmi_`var'_wf: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN  c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6  if cc_sample==1, vce(cluster fid) 
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
*OK - so that weird increase in the ADHD effect isn't because of something that happened in the imputation.
}
estout using "output\tables\completecase_all_adultbmi.xls", cells ("b(fmt(2)) ci_l(fmt(2)) ci_u(fmt(2)) p N") keep(*bmi*) replace title(BMI to age 8 outcomes, using adult BMI PGS) note("comparing non-genetic, classic MR and WFMR models")

*CHILDBMI PGS
capture erase "output\tables\completecase_all_childbmi.txt"
capture erase "output\tables\completecase_all_childbmi.txt"
eststo clear
foreach var in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr z_inADHD_8yr z_hyADHD_8yr {
*non-genetic
eststo childbmi_`var'_pheno: reg `var' rescaled_moffspring_bmi_Q8 rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF KJONN c_yob i.AA1124_v2_flipped i.f_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1 i.dadsmokes_QF i.PARITET_5 /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6 if cc_sample==1, vce(cluster fid) 
*classic MR
eststo childbmi_`var'_mr:  ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_childbmi_pgs) KJONN c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2 m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6  if cc_sample==1, vce(cluster fid) 
*WFMR
eststo childbmi_`var'_wf: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_childbmi_pgs m_z_childbmi_pgs f_z_childbmi_pgs) KJONN c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6  if cc_sample==1, vce(cluster fid) 
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
}
estout using "output\tables\completecase_all_childbmi.xls", cells ("b(fmt(2)) ci_l(fmt(2)) ci_u(fmt(2)) p N") keep(*bmi*) replace title(BMI to age 8 outcomes, using childhood BMI PGS) note("comparing non-genetic, classic MR and WFMR models")

**************************************************

*REVIEWER COMMENTS:

*assortative mating check:

*genotypic assortative mating check, since doesn't require imputed data:
foreach m_pgs in z_adultbmi_pgs z_childbmi_pgs z_depression_pgs z_adhd_pgs  {
foreach f_pgs in z_adultbmi_pgs z_childbmi_pgs z_depression_pgs z_adhd_pgs  {
display ""
display "`m_pgs' `f_pgs'"
regress f_`f_pgs' m_`m_pgs' m_PC* f_PC* m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6
}
}

*regression models:
foreach outcome in z_depression_pgs z_anxiety_pgs z_adhd_pgs z_asd_pgs {
regress f_`outcome' m_z_adultbmi_pgs m_PC* f_PC* m_batch* f_batch*, vce(cluster fid) 
}


*5) Are the actual phenotypes (BMI, depression or ADHD) correlated between the parents? If so, would this not suffice as evidence of cross-trait assortative mating? It is known that the genetic correlation between parents as a result of assortative mating is a function of the correlation in their phenotypes and the heritabilities underlying the two traits (e.g., see Yengo and Visscher 2018). An alternative way to estimate the genetic correlation between parents without using PGS (which is noisy and therefore underpowered) would be to use the phenotypic correlation and heritability estimated using GREML or LDSC. Perhaps this is outside the scope of the paper but I would like to hear the author's thoughts on this.

*needs to be done in imputed data.

*7) Finally, what is the correlation between PGS and genetic PCs/geography in their sample? A correlation might provide evidence to support the point that classic MR effects are inflated due to stratification.

*child's pgs on child's PCs, adjusting for child's genotyping centre and chip:
eststo clear
foreach pgs in z_adultbmi_pgs z_childbmi_pgs  {
eststo `pgs': regress c_`pgs'  /*KJONN c_yob*/ /*genetic covariates: PCs, center and chip but not batch*/ c_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2, vce(cluster fid) 
estimates store `pgs'_PCs
esttab, b(%9.2f) ci(%9.2f) keep(*c_PC*) compress noparen nostar
}
*7 of the 40 show sig associations. 2 would be expected by chance.
*plot these:
*plot the estimates:

global z_adultbmi_pgs_title = "Adult BMI PGS"
global z_childbmi_pgs_title = "Child BMI PGS"

foreach pgs in z_adultbmi_pgs z_childbmi_pgs {
coefplot (`pgs'_PCs), ///
keep (c_PC*) drop(_cons) xline(0) grid(none) title("{bf:$`pgs'_title}", size(med) pos(4) ring(0)) scheme(s1color) graphregion( color(white)) plotregion( color(white)) xtick(-5 -4 -3 -2 -1 0 1 2 3 4 5) xlabel (-5 -4 -3 -2 -1 0 1 2 3 4 5,labsize(small))
graph save `pgs'_PCs.gph, replace
}
graph combine z_adultbmi_pgs_PCs.gph z_childbmi_pgs_PCs.gph, xcommon title("{bf:Unadjusted associations of child's BMI polygenic scores with ancestry PCs}", size(medsmall) color(black) pos(11))  graphregion( color(white)) plotregion( color(white))
graph export output/graphs/association_pgs_PCs.tif, replace width(1200)


*child's pgs on child's PCs and genotyping centre, adjusting for parents' genotype:
eststo clear
foreach pgs in z_adultbmi_pgs z_childbmi_pgs  {
eststo adj_`pgs': regress c_`pgs' m_`pgs' f_`pgs'  /*KJONN c_yob*/  /*genetic covariates: PCs, center and chip but not batch*/ c_PC* c_genotyping_center1 c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 /* m_PC* f_PC* m_genotyping_chip2 m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6*/, vce(cluster fid) 
estimates store adj_`pgs'_PCs
esttab, b(%9.2f) ci(%9.2f) keep(*c_PC*) compress noparen nostar
}

foreach pgs in z_adultbmi_pgs z_childbmi_pgs {
coefplot (adj_`pgs'_PCs), ///
keep (c_PC*) drop(_cons) xline(0) grid(none) title("{bf:$`pgs'_title}", size(med) pos(4) ring(0)) scheme(s1color) graphregion( color(white)) plotregion( color(white)) xtick(-5 -4 -3 -2 -1 0 1 2 3 4 5) xlabel (-5 -4 -3 -2 -1 0 1 2 3 4 5,labsize(small))
graph save adj_`pgs'_PCs.gph, replace
}
graph combine adj_z_adultbmi_pgs_PCs.gph adj_z_childbmi_pgs_PCs.gph, xcommon title("{bf:Adjusted associations of child's BMI polygenic scores with ancestry PCs}", size(medsmall) color(black) pos(11)) graphregion( color(white)) plotregion( color(white))
graph export output/graphs/adj_association_pgs_PCs.tif, replace width(1200)

*******************************************************************
*check for which ea4 pgs to use:

*for kids, using placeholder outcomes:
fre NN239 NN240 NN241

foreach var in NN239 NN240 NN241 {
foreach pgs in ea4 ea4excl23andme {
	regress `var' KJONN c_yob c_z_`pgs'_pgs
	ologit `var' KJONN c_yob c_z_`pgs'_pgs
}
}
*and for parents, using qualifications:
foreach pgs in ea4 ea4excl23andme {
	regress AA1124_v2_flipped m_yob m_z_`pgs'_pgs
	regress f_educ_Q1QF_v2_flipped f_yob f_z_`pgs'_pgs
}
*in all cases, regardless of the specification, R2 is higher with the first score.

