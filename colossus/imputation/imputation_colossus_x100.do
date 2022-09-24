
*needs to be called with the shell script, but the supershell isn't necessary

set processors 8

cd /cluster/p/p471/cluster/projects/mh_bmi/fullsample/imputation
sysdir set PLUS "/cluster/p/p471/cluster/people/Mandystatafiles/ado/plus"
set maxvar 20000
sysdir

*load pre-imputation data
use "pre_imputation.dta", clear
count

*ACTUAL IMPUTATION

*in the model below, using information from:
*MBRN Q1 QF Q5 Q6 Q5y Q7y Q8y

capture erase impstats.dta
mi impute chained ///
(ologit) /*SEP: parental education*/ AA1124_v2_flipped f_educ_Q1QF_v2_flipped save_EE584 ///
(truncreg, ll(1) ul(7)) /*parental income, grouped*/ AA1315 fathers_inc_grp  ///
(ologit) /*parental smoking*/ mumsmokes_Q1 dadsmokes_Q1QF  ///
(pmm, knn(10)) /*parental anxiety/dep */ mhopkins_Q1 ///
(pmm, knn(10)) /*parental anxiety/dep */ fhopkins_QF ///
(truncreg, ll(8.1) ul(50.2)) mothers_bmi_Q1 ///
(truncreg, ll(13.6) ul(49.6)) mothers_bmi_Q6 ///
(truncreg, ll(14.9) ul(48.3)) fathers_bmi_Q1QF  ///
(truncreg, ll(8.1) ul(28.3)) moffspring_bmi_Q5y ///
(truncreg, ll(6.4) ul(29.3)) moffspring_bmi_Q7y ///
(truncreg, ll(6.6) ul(32.2)) moffspring_bmi_Q8y  ///
(pmm, knn(10))  /*OUTCOMES: age 8 MFQ score*/ MFQ_8yr ///
(pmm, knn(10)) /*OUTCOMES: age 8 anxiety score*/ SCARED_8yr ///
(pmm, knn(10)) /*OUTCOMES: age 8 anxiety score*/ inADHD_8yr ///
(pmm, knn(10)) /*OUTCOMES: age 8 anxiety score*/ hyADHD_8yr ///
(truncreg, ll(0) ul(24)) /*parental ADHD summ scores, bumped down to start at 0*/ mADHD_Q6 fADHD_QF  ///
(pmm, knn(10)) /*Q5y: ADHD*/ conners_Q5y ///
(truncreg, ll(1075) ul(5880)) /*birthweight*/ VEKT  ///
(truncreg, ll(39) ul(61)) /*birthlength*/ LENGDE ///
(truncreg, ll(17) ul(60)) /*paternal age at birth*/ FARS_ALDER ///
= /*child sex, year of birth, then parity, maternal age, and marital status from birth registry file*/ KJONN c_yob PARITET_5 MORS_ALDER SIVST_2 ///
/*bmi and other polygenic scores */ *z_childbmi_pgs *z_adultbmi_pgs *z_full_pgs_adultbmi *z_full_pgs_childbmi *z_ea4_pgs* *z_depression_pgs *z_adhd_pgs  ///
/* PCs */ c_PC* m_PC* f_PC* ///
/*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1 m_genotyping_center2 f_genotyping_center1  f_genotyping_center2 ///
/*genotyping chip: ommitted (largest) group is 2, 1, 1*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6 /// 
  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6 ///
, add(100) rseed(100) dots savetrace(impstats.dta) augment 

************************************************

*POST-IMPUTATION TRANSFORMATIONS with mi passive:

*RESCALE BMI as per 5kg/m2 to have less stupidly small coefficients:
*not super-varying so can use mi passive:
mi passive: gen rescaled_moffspring_bmi_Q8y=moffspring_bmi_Q8/5
mi passive: gen rescaled_mothers_bmi_Q1=mothers_bmi_Q1/5
mi passive: gen rescaled_fathers_bmi_Q1QF=fathers_bmi_Q1QF/5

*transformations of MFQ and SCARED: do this in terms of mean and SD of complete-case sample, because otherwise they will be super-varying variables requiring working in flong format:
foreach var in MFQ_8yr SCARED_8yr  {
summ `var'
return list 
scalar `var'_mean = r(mean)
scalar `var'_sd = r(sd)
mi passive: gen z_`var'=(`var'-`var'_mean)/`var'_sd
summ z_`var'
*mi xeq 1: summ z_`var'
}

*ADHD main scale
*according to Sean, the correct thing to do is standardize, add, then standardize again. 
*the only way to ensure a) equal contributions of the two halves and b) the scale you want on the final variable.

*standardization: do this in terms of mean and SD of complete-case sample, because otherwise they will be super-varying variables requiring working in flong format:
foreach var in inADHD_8yr hyADHD_8yr {
summ `var'
return list 
scalar `var'_mean = r(mean)
scalar `var'_sd = r(sd)
mi passive: gen z_`var'=(`var'-`var'_mean)/`var'_sd
summ z_`var'
*mi xeq 1: summ z_`var'
}
*then add the standardized ones:
mi passive: gen ADHD_8yr=z_inADHD_8yr+z_hyADHD_8yr

*then standardize that again, again in terms of the mean and sd in m=0:
foreach var in ADHD_8yr {
summ `var'
return list 
scalar `var'_mean = r(mean)
scalar `var'_sd = r(sd)
mi passive: gen z_`var'=(`var'-`var'_mean)/`var'_sd
summ z_`var'
*mi xeq 1: summ z_`var'
}
summ ADHD_8yr

save "imp_x100_17July2022.dta", replace
