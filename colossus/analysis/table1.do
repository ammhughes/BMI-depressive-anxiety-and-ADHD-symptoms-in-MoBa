**************************************************************************
*ANALYSIS MODELS WITH IMPUTED DATA: TABLE 1 
**************************************************************************

*FOR ANALYSIS, NO SUPERSHELL REQUIRED

set processors 8

*********************************************************************************

cd "/cluster/p/p471/cluster/projects/mh_bmi/fullsample/analysis"
set maxvar 20000
sysdir
sysdir set PLUS "/cluster/p/p471/cluster/people/Mandy/statafiles/ado/plus"

use "imp_x100_17July2022.dta", clear

********************************************************************************
*TABLE 1:

*for misum, data needs to be in flong format.
***ssc install misum***
mi convert flong

log using "post_imp_descriptives.log", replace

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

use "imp_x100_17July2022.dta", clear

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