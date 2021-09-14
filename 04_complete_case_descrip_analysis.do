*COMPLETE-CASE DESCRIPTIVES AND ANALYSIS 

sysdir set PLUS "N:\data\durable\people\Mandy\Statafiles\ado\plus"
cd "N:\data\durable\projects\mh_bmi\"
set maxvar 20000

use "data\statafiles\geno_pheno_all.dta", clear

*UPDATE 12 MAY: RESTRICT TO PEOPLE IN UPDATED SV_INFO FILE (i.e. exclude the new withdrawls) 
merge m:1 PREG_ID_2306 using "N:\data\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\Statafiles\PDB2306_SV_INFO_v12_UPDATED_April_2021"
*the 51 not matched are the ones that need to be dropped as new exclusions.
drop if _merge==1
*51 dropped
drop _merge

count
*114,717

*standardize all the outcomes:
foreach var in MFQ_8yr SCARED_8yr RSDBD_ADHD_8yr RSDBD_inADHD_8yr RSDBD_hyADHD_8yr SCQ_8yr SCQ_SCI_8yr SCQ_RRB_8yr {
egen z_`var'=std(`var')
}

*rename these otherwise names too long to work in loop storing estimates
rename z_RSDBD_ADHD_8yr z_ADHD_8yr
rename z_RSDBD_inADHD_8yr z_inADHD_8yr
rename z_RSDBD_hyADHD_8yr z_hyADHD_8yr
rename z_SCQ_SCI_8yr z_S_SCQ_8yr 
rename z_SCQ_RRB_8yr z_R_SCQ_8yr

*RESCALE BMI as per 5kg/m2 to have less stupidly small coefficients:
gen rescaled_moffspring_bmi_Q8y=moffspring_bmi_Q8/5
gen rescaled_mothers_bmi_Q1=mothers_bmi_Q1/5
gen rescaled_fathers_bmi_Q1QF=fathers_bmi_Q1QF/5

******************************************************************************************
*DESCRIPTIVES FOR ALL PEOPLE IN BIRTH REGISTRY FILE, minus the new withdrawals:

count 
*114,717

foreach var in MORS_ALDER FARS_ALDER mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF mothers_bmi_Q1 fathers_bmi_Q1QF  moffspring_bmi_Q8y MFQ_8yr SCARED_8yr RSDBD_ADHD_8yr RSDBD_inADHD_8yr RSDBD_hyADHD_8yr  {
display "`var'"
format `var' %3.1f
summ `var', format
}

foreach var in KJONN AA1124_v2 fathers_educ_Q1QF_v2 PARITET_5 SIVST_2 mumsmokes_Q1  {
fre `var', f(1) nomiss
}

*finally, get min and max birtweight and length IN FULL SAMPLE for truncreg limits in imputation.

*for both, apply the rule of 4sd from the mean:
foreach var in VEKT LENGDE {
summ `var', det
return list
replace `var'=. if `var'<r(mean)-(4*r(sd))
replace `var'=. if `var'>r(mean)+(4*r(sd))
summ `var', det
}

foreach var in  mothers_bmi_Q1 mothers_bmi_Q6 fathers_bmi_Q1QF moffspring_bmi_Q5y moffspring_bmi_Q7y moffspring_bmi_Q8y VEKT LENGDE FARS_ALDER {
	summ `var'
}
/*
(truncreg, ll(8.1) ul(50.2)) mothers_bmi_Q1 ///
(truncreg, ll(13.6) ul(49.6)) mothers_bmi_Q6 ///
(truncreg, ll(13.1) ul(48.3)) fathers_bmi_Q1QF  ///
(truncreg, ll(8.1) ul(28.3)) moffspring_bmi_Q5y ///
(truncreg, ll(6.4) ul(29.3)) moffspring_bmi_Q7y ///
(truncreg, ll(6.6) ul(32.2)) moffspring_bmi_Q8y  ///

(truncreg, ll(1075) ul(5880)) /*birthweight*/ VEKT  ///
(truncreg, ll(39) ul(61)) /*birthlength*/ LENGDE ///
(truncreg, ll(17) ul(60)) /*paternal age at birth*/ FARS_ALDER ///
*/

*********************************************************************************************

*IDENTIFY THE FINAL SAMPLE.

*FULL SAMPLE:
count
*114,717
count if KJONN==.
*831
count if KJONN!=.
*113,886

*exclude people present in the birth registry file but not at any of the actual questionnaires
fre present*
gen nQ_present=0
foreach q in Q1 QF Q3 Q4 Q5 Q6 Q5y Q7y Q8y F2 {
	replace nQ_present=nQ_present+1 if present_`q'==1
}
fre nQ_present
*OK, so there's a chunk of people who never were included in a single questionnaire.
*these are the people that weren't present at the very first one:
tab nQ_present present_Q1, mi

count if KJONN!=. & nQ_present==0
*8,101

*who has genetic information? 
count if KJONN!=. & nQ_present!=0 & c_batch2!=. 
*30,887
*have lost some, relative to the number of children with genetic info: 31,374 (after removing one of the identical twins)
tab fullgeneticinfo
*26,370 trios

*up to here:
capture drop analytic_sample
gen analytic_sample=0
replace analytic_sample=1 if nQ_present!=0 & fullgeneticinfo==1 & KJONN!=.
fre analytic_sample

/*
*SIBLINGS AND OTHER RELATEDS EXCLUDED POST-IMPUTATION.
*can't merge in exclusion lists here in the same way as for actual analysis because of missing mid pid for the non-genotyped people, so made a flag for inclusion in final sample in separate file.

merge 1:1 PREG_ID_2306 BARN_NR using "scratch/in_analytic_sample_flag.dta"

replace analytic_sample=0 if _merge!=3
drop _merge
fre analytic_sample
*24256

*FINAL COUNT FOR ANALYTIC SAMPLE:
count
*24256
*/

******************************************************************************************
*DESCRIPTIVES FOR ANALYTIC SAMPLE, UNIMPUTED DATA.


*VEKT
foreach var in MORS_ALDER FARS_ALDER mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF mothers_bmi_Q1 fathers_bmi_Q1QF  moffspring_bmi_Q8y MFQ_8yr SCARED_8yr RSDBD_ADHD_8yr RSDBD_inADHD_8yr RSDBD_hyADHD_8yr {
display "`var'"
format `var' %3.1f
summ `var' if analytic_sample==1, format
}

foreach var in KJONN AA1124_v2 fathers_educ_Q1QF_v2 PARITET_5 SIVST_2 mumsmokes_Q1  {
fre `var' if analytic_sample==1, f(1) nomiss
}

*****************************************************************************

*COMPARISON OF ANALYTIC SAMPLE WITH OTHERS IN THE BIRTH REGISTRY FILE, 
*DESCRIPTIVE TABLE FOR ALL PEOPLE IN BIRTH REGISTRY FILE

*comparison - contin vars
*sd test
foreach var in VEKT MORS_ALDER FARS_ALDER mothers_bmi_Q1 fathers_bmi_Q1QF moffspring_bmi_Q8y MFQ_8yr SCARED_8yr RSDBD_ADHD_8yr {
display "`var'"
sdtest `var', by(analytic_sample)
}
*unequal variance t-tests:
foreach var in VEKT MORS_ALDER FARS_ALDER mothers_bmi_Q1 fathers_bmi_Q1QF moffspring_bmi_Q8y MFQ_8yr SCARED_8yr RSDBD_ADHD_8yr {
display "`var'"
ttest `var', by(analytic_sample) unequal
}
*has equal variance:
ttest fathers_bmi_Q1QF, by(analytic_sample)

*comparison - categ vars
foreach var in KJONN PARITET_5 SIVST_2 AA1124_v2 fathers_educ_Q1QF_v2 mumsmokes_Q1 dadsmokes_Q1QF  {
	tab `var' analytic_sample, chi col
}


*********************************************************************************************

*ANALYSIS IN COMPLETE-CASE

keep if analytic_sample==1

*prep
foreach var in AA1124_v2 fathers_educ_Q1QF_v2  {
	fre `var'
	gen `var'_flipped=7-`var'
	tab `var' `var'_flipped
}

*how much missingness?
*bmi
count if !inlist(., rescaled_moffspring_bmi_Q8, rescaled_mothers_bmi_Q1, rescaled_fathers_bmi_Q1QF)
*bmi and all phenotypic covars
count if !inlist(., rescaled_moffspring_bmi_Q8, rescaled_mothers_bmi_Q1, rescaled_fathers_bmi_Q1QF, KJONN, AA1124_v2_flipped, fathers_educ_Q1QF_v2_flipped, mhopkins_Q1, fhopkins_QF, mADHD_Q6, fADHD_QF, mumsmokes_Q1, dadsmokes_QF, PARITET_5)

*flag complete-case including all phenotypic covariates
gen cc_sample=1 if !inlist(., rescaled_moffspring_bmi_Q8, rescaled_mothers_bmi_Q1, rescaled_fathers_bmi_Q1QF, KJONN, AA1124_v2_flipped, fathers_educ_Q1QF_v2_flipped, mhopkins_Q1, fhopkins_QF, mADHD_Q6, fADHD_QF, mumsmokes_Q1, dadsmokes_QF, PARITET_5)


*IV models: save for export:
capture erase "output\tables\completecase_allresults.txt"
capture erase "output\tables\completecase_allresults.txt"

eststo clear
foreach var in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr z_inADHD_8yr z_hyADHD_8yr {
*adultbmi
*phenotypic
eststo adultbmi_`var'_pheno: reg `var' rescaled_moffspring_bmi_Q8 KJONN i.AA1124_v2_flipped i.fathers_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1 i.dadsmokes_QF i.PARITET_5 if cc_sample==1, vce(cluster fid) 
*normal genetic -  iv with pgs
eststo adultbmi_`var'_mr:  ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_PC* c_batch* if cc_sample==1, vce(cluster fid) 
*with parents genetic -  iv with pgs
eststo adultbmi_`var'_wf: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch* if cc_sample==1, vce(cluster fid) 
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
*childbmi
*phenotypic
eststo childbmi_`var'_pheno: reg `var' rescaled_moffspring_bmi_Q8 KJONN i.AA1124_v2_flipped i.fathers_educ_Q1QF_v2_flipped mhopkins_Q1 fhopkins_QF mADHD_Q6 fADHD_QF i.mumsmokes_Q1 i.dadsmokes_QF i.PARITET_5 if cc_sample==1, vce(cluster fid) 
*normal genetic -  iv with pgs
eststo childbmi_`var'_mr:  ivregress 2sls `var' (rescaled_moffspring_bmi_Q8=c_z_childbmi_pgs) KJONN c_PC* c_batch* if cc_sample==1, vce(cluster fid) 
*with parents genetic -  iv with pgs
eststo childbmi_`var'_wf: ivregress 2sls `var' (rescaled_moffspring_bmi_Q8y rescaled_mothers_bmi_Q1 rescaled_fathers_bmi_Q1QF = c_z_childbmi_pgs m_z_childbmi_pgs f_z_childbmi_pgs) KJONN c_PC* m_PC* f_PC* c_batch* m_batch* f_batch* if cc_sample==1, vce(cluster fid) 
esttab, b(%9.2f) ci(%9.2f) keep(*bmi*) compress noparen nostar
*eststo clear
}
estout using "output\tables\completecase_allresults.xls", cells ("b(fmt(2)) ci_l(fmt(2)) ci_u(fmt(2)) N") keep(*bmi*) replace title(BMI to age 8 outcomes) note("comparing standard and within-family IV models.")
*Ns:
*MFQ: 3120
*SCARED: 3147
*ADHD: 3016
*inADHD: 3139
*hyADHD: 3148

**************************************************
*also, assortative mating check, since doesn't require imputed data:


foreach m_pgs in z_adultbmi_pgs z_childbmi_pgs z_depression_pgs z_adhd_pgs  {
foreach f_pgs in z_adultbmi_pgs z_childbmi_pgs z_depression_pgs z_adhd_pgs  {
display ""
display "`m_pgs' `f_pgs'"
pwcorr m_`m_pgs' f_`f_pgs'
regress f_`f_pgs' m_`m_pgs' m_PC* f_PC* m_batch* f_batch*
}
}

*correlation of PGSs for BMI and outcomes in parents generation:
foreach outcome in z_adultbmi_pgs z_childbmi_pgs z_depression_pgs z_adhd_pgs  {
display ""
display "m_`outcome'"
pwcorr m_`outcome' f_z_adultbmi_pgs
display ""
display "m_`outcome'"
pwcorr m_`outcome' f_z_childbmi_pgs
display ""
display "m_`outcome'"
pwcorr m_`outcome' f_z_depression_pgs
display ""
display "m_`outcome'"
pwcorr m_`outcome' f_z_adhd_pgs
}
*nada: max 0.01 for everything

*regression models:
foreach outcome in z_depression_pgs z_anxiety_pgs z_adhd_pgs z_asd_pgs {
regress f_`outcome' m_z_adultbmi_pgs m_PC* f_PC* m_batch* f_batch*
}


