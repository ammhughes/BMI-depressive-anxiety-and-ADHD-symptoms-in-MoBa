****************************************************************************************

*ENTRY POINT: R2 AND WEAK INSTRUMENT TESTS ACROSS IMPUTATIONS

clear all
cd "N:\durable\projects\mh_bmi\fullsample"
set maxvar 20000
sysdir
sysdir set PLUS "N:\durable\people\Mandy\Statafiles\ado\plus"
sysdir set PERSONAL "N:\durable\people\Mandy\statafiles\ado\personal"

use "data/imp_x100_17July2022.dta", clear

**************************************************************************
*MANUALLY LOAD IVREG2
*after putting the mlib file in the new PERSONAL folder, manually run:
do "N:\durable\people\Mandy\statafiles\ado\plus\livreg2.do.txt"
do "N:\durable\people\Mandy\statafiles\ado\plus\ivreg2.ado.txt"
do "N:\durable\people\Mandy\statafiles\ado\plus\ivreg2_p.ado.txt"
**************************************************************************
capture log close
log using "N:\durable\projects\mh_bmi\fullsample\output\r2.log", append


*correlation of the two scores:
foreach p in c m f {
pwcorr `p'_z_adultbmi_pgs `p'_z_childbmi_pgs
}

*for the conditional R2, need to get SWr2 from the output using ", first" after ivreg2.
*also need to rename the exposures to have shorter names - otherwise can't store and display all the results.
mi rename rescaled_moffspring_bmi_Q8 c_bmi_Q8
mi rename rescaled_mothers_bmi_Q1 m_bmi_Q1
mi rename rescaled_fathers_bmi_Q1QF f_bmi_Q1QF

*ADULTBMI
*run this across all imputations, then take the average
*can just do one outcome because the first stage is the same!
foreach var in z_MFQ_8yr {
*vars to store the r2 and F stats in:
capture gen `var'_c_adultbmi_swfs=.
capture gen `var'_m_adultbmi_swfs=.
capture gen `var'_f_adultbmi_swfs=.
capture gen `var'_c_adultbmi_swr2=.
capture gen `var'_m_adultbmi_swr2=.
capture gen `var'_f_adultbmi_swr2=.
forvalues m = 1/100 {
ivreg2 _`m'_`var' (_`m'_c_bmi_Q8 _`m'_m_bmi_Q1 _`m'_f_bmi_Q1QF = c_z_adultbmi_pgs m_z_adultbmi_pgs f_z_adultbmi_pgs) KJONN  c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, first 
*retrieve:
mat list e(first)
mat ivreg2stats=e(first)
*SW R2:
*for child:
scalar c_swr2 = ivreg2stats[14,1] 
display c_swr2
replace `var'_c_adultbmi_swr2=c_swr2 in `m'
*for mother:
scalar m_swr2 = ivreg2stats[14,2] 
display m_swr2
replace `var'_m_adultbmi_swr2=m_swr2 in `m'
*for father:
scalar f_swr2 = ivreg2stats[14,3] 
display f_swr2
replace `var'_f_adultbmi_swr2=f_swr2 in `m'
*SW F stats:
*for child:
scalar c_swfs = ivreg2stats[8,1] 
display c_swfs
replace `var'_c_adultbmi_swf=c_swfs in `m'
*for mother:
scalar m_swfs = ivreg2stats[8,2] 
display m_swfs
replace `var'_m_adultbmi_swf=m_swfs in `m'
*for father:
scalar f_swfs = ivreg2stats[8,3] 
display f_swfs
replace `var'_f_adultbmi_swf=f_swfs in `m'
}
}

*summarize across imputations:
*conditional r2
foreach var in z_MFQ_8yr  {
mean `var'_c_adultbmi_swr2
mean `var'_m_adultbmi_swr2
mean `var'_f_adultbmi_swr2
}
/*
Mean estimation                                         Number of obs = 100
---------------------------------------------------------------------------
                          |       Mean   Std. err.     [95% conf. interval]
--------------------------+------------------------------------------------
z_MFQ_8yr_c_adultbmi_swr2 |   .0172819   .0001309      .0170223    .0175415
--------------------------+------------------------------------------------
z_MFQ_8yr_m_adultbmi_swr2 |   .0317014   .0002519      .0312016    .0322012
--------------------------+------------------------------------------------
z_MFQ_8yr_f_adultbmi_swr2 |   .0301918   .0002565      .0296828    .0307007
---------------------------------------------------------------------------

*/
*conditional F stat
foreach var in z_MFQ_8yr  {
mean `var'_c_adultbmi_swfs
mean `var'_m_adultbmi_swfs
mean `var'_f_adultbmi_swfs
}
/*
Mean estimation                                        Number of obs = 100
---------------------------------------------------------------------------
                          |       Mean   Std. err.     [95% conf. interval]
--------------------------+------------------------------------------------
z_MFQ_8yr_c_adultbmi_swfs |    718.735    5.53818      707.7461     729.724
--------------------------+------------------------------------------------
z_MFQ_8yr_m_adultbmi_swfs |   1338.205   10.98332      1316.412    1359.999
--------------------------+------------------------------------------------
z_MFQ_8yr_f_adultbmi_swfs |   1272.519   11.14962      1250.396    1294.643
---------------------------------------------------------------------------

*/


*CHILDBMI
*for the conditional R2, need to get SWr2 from the output using ", first" after ivreg2.
*run this across all imputations, then take the average
*can just do one outcome because the first stage is the same!
foreach var in z_MFQ_8yr {
*vars to store the r2 and F stats in:
capture gen `var'_c_childbmi_swfs=.
capture gen `var'_m_childbmi_swfs=.
capture gen `var'_f_childbmi_swfs=.
capture gen `var'_c_childbmi_swr2=.
capture gen `var'_m_childbmi_swr2=.
capture gen `var'_f_childbmi_swr2=.
forvalues m = 1/100 {
ivreg2 _`m'_`var' (_`m'_c_bmi_Q8 _`m'_m_bmi_Q1 _`m'_f_bmi_Q1QF = c_z_childbmi_pgs m_z_childbmi_pgs f_z_childbmi_pgs) KJONN  c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, first 
*retrieve:
mat list e(first)
mat ivreg2stats=e(first)
*SW R2:
*for child:
scalar c_swr2 = ivreg2stats[14,1] 
display c_swr2
replace `var'_c_childbmi_swr2=c_swr2 in `m'
*for mother:
scalar m_swr2 = ivreg2stats[14,2] 
display m_swr2
replace `var'_m_childbmi_swr2=m_swr2 in `m'
*for father:
scalar f_swr2 = ivreg2stats[14,3] 
display f_swr2
replace `var'_f_childbmi_swr2=f_swr2 in `m'
*SW F stats:
*for child:
scalar c_swfs = ivreg2stats[8,1] 
display c_swfs
replace `var'_c_childbmi_swf=c_swfs in `m'
*for mother:
scalar m_swfs = ivreg2stats[8,2] 
display m_swfs
replace `var'_m_childbmi_swf=m_swfs in `m'
*for father:
scalar f_swfs = ivreg2stats[8,3] 
display f_swfs
replace `var'_f_childbmi_swf=f_swfs in `m'
}
}

*summarize across imputations:
*conditional r2
foreach var in z_MFQ_8yr  {
mean `var'_c_childbmi_swr2
mean `var'_m_childbmi_swr2
mean `var'_f_childbmi_swr2
}
/*
Mean estimation                                         Number of obs = 100
---------------------------------------------------------------------------
                          |       Mean   Std. err.     [95% conf. interval]
--------------------------+------------------------------------------------
z_MFQ_8yr_c_childbmi_swr2 |   .0220099   .0001292      .0217536    .0222663
--------------------------+------------------------------------------------
z_MFQ_8yr_m_childbmi_swr2 |   .0255556   .0001091      .0253392    .0257721
--------------------------+------------------------------------------------
z_MFQ_8yr_f_childbmi_swr2 |    .022955   .0000975      .0227615    .0231484
---------------------------------------------------------------------------
*/

*conditional F stat
foreach var in z_MFQ_8yr  {
mean `var'_c_childbmi_swfs
mean `var'_m_childbmi_swfs
mean `var'_f_childbmi_swfs
}
/*
Mean estimation                                         Number of obs = 100
---------------------------------------------------------------------------
                          |       Mean   Std. err.     [95% conf. interval]
--------------------------+------------------------------------------------
z_MFQ_8yr_c_childbmi_swfs |    919.772   5.522608      908.8139      930.73
--------------------------+------------------------------------------------
z_MFQ_8yr_m_childbmi_swfs |   1071.798   4.693837      1062.484    1081.111
--------------------------+------------------------------------------------
z_MFQ_8yr_f_childbmi_swfs |   960.1581    4.17484      951.8744    968.4419

*/


********************************************************************************************

*REPEAT FOR THE FULL GWAS VERSIONS - JUST R2.
*hang on  - what's relevant here is the r2 from a model without parental genotype in.

*for the conditional R2, need to get SWr2 from the output using ", first" after ivreg2.
*also need to rename the exposures to have shorter names - otherwise can't store and display all the results.
capture mi rename rescaled_moffspring_bmi_Q8 c_bmi_Q8
capture mi rename rescaled_mothers_bmi_Q1 m_bmi_Q1
capture mi rename rescaled_fathers_bmi_Q1QF f_bmi_Q1QF

*run this across all imputations, then take the average
*adultbmi
foreach var in z_MFQ_8yr {
*vars to store the r2 in:
capture gen `var'_c_full_abmi_swr2=.
capture gen `var'_m_full_abmi_swr2=.
capture gen `var'_f_full_abmi_swr2=.
forvalues m = 1/100 {
*for child:
ivreg2 _`m'_`var' (_`m'_c_bmi_Q8  = c_z_full_pgs_adultbmi) KJONN  c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, first 
*retrieve:
mat list e(first)
mat ivreg2stats=e(first)
*SW R2:
scalar c_swr2 = ivreg2stats[14,1] 
display c_swr2
replace `var'_c_full_abmi_swr2=c_swr2 in `m'
*for mother:
ivreg2 _`m'_`var' (_`m'_m_bmi_Q1  = m_z_full_pgs_adultbmi) KJONN  c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, first 
*retrieve:
mat list e(first)
mat ivreg2stats=e(first)
*SW R2:
scalar m_swr2 = ivreg2stats[14,1] 
display m_swr2
replace `var'_m_full_abmi_swr2=m_swr2 in `m'
*for father:
ivreg2 _`m'_`var' (_`m'_f_bmi_Q1QF  = f_z_full_pgs_adultbmi) KJONN  c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, first 
*retrieve:
mat list e(first)
mat ivreg2stats=e(first)
*SW R2:
scalar f_swr2 = ivreg2stats[14,1] 
display f_swr2
replace `var'_f_full_abmi_swr2=f_swr2 in `m'
}
}
*summarize across imputations:
*conditional r2
foreach var in z_MFQ_8yr  {
mean `var'_c_full_abmi_swr2
mean `var'_m_full_abmi_swr2
mean `var'_f_full_abmi_swr2
}

/*Mean estimation                                          Number of obs = 100
----------------------------------------------------------------------------
                           |       Mean   Std. err.     [95% conf. interval]
---------------------------+------------------------------------------------
z_MFQ_8yr_c_full_abmi_swr2 |   .0430111    .000197      .0426202     .043402
---------------------------+------------------------------------------------
z_MFQ_8yr_m_full_abmi_swr2 |   .0820331   .0000416      .0819505    .0821157
---------------------------+------------------------------------------------
z_MFQ_8yr_f_full_abmi_swr2 |   .0806239   .0000545      .0805156    .0807321
----------------------------------------------------------------------------
*/

*childbmi
foreach var in z_MFQ_8yr {
*vars to store the r2 in:
capture gen `var'_c_full_cbmi_swr2=.
capture gen `var'_m_full_cbmi_swr2=.
capture gen `var'_f_full_cbmi_swr2=.
forvalues m = 1/100 {
*for child:
ivreg2 _`m'_`var' (_`m'_c_bmi_Q8  = c_z_full_pgs_childbmi) KJONN  c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, first 
*retrieve:
mat list e(first)
mat ivreg2stats=e(first)
*SW R2:
scalar c_swr2 = ivreg2stats[14,1] 
display c_swr2
replace `var'_c_full_cbmi_swr2=c_swr2 in `m'
*for mother:
ivreg2 _`m'_`var' (_`m'_m_bmi_Q1  = m_z_full_pgs_childbmi) KJONN  c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, first 
*retrieve:
mat list e(first)
mat ivreg2stats=e(first)
*SW R2:
scalar m_swr2 = ivreg2stats[14,1] 
display m_swr2
replace `var'_m_full_cbmi_swr2=m_swr2 in `m'
*for father:
ivreg2 _`m'_`var' (_`m'_f_bmi_Q1QF  = f_z_full_pgs_childbmi) KJONN  c_yob /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, first 
*retrieve:
mat list e(first)
mat ivreg2stats=e(first)
*SW R2:
scalar f_swr2 = ivreg2stats[14,1] 
display f_swr2
replace `var'_f_full_cbmi_swr2=f_swr2 in `m'
}
}
*summarize across imputations:
*conditional r2
foreach var in z_MFQ_8yr  {
mean `var'_c_full_cbmi_swr2
mean `var'_m_full_cbmi_swr2
mean `var'_f_full_cbmi_swr2
}


/*Mean estimation                                          Number of obs = 100
----------------------------------------------------------------------------
                           |       Mean   Std. err.     [95% conf. interval]
---------------------------+------------------------------------------------
z_MFQ_8yr_c_full_cbmi_swr2 |   .0565768   .0002226      .0561352    .0570184
---------------------------+------------------------------------------------
z_MFQ_8yr_m_full_cbmi_swr2 |   .0333813   .0000308      .0333201    .0334425
---------------------------+------------------------------------------------
z_MFQ_8yr_f_full_cbmi_swr2 |   .0300506   .0000328      .0299854    .0301157
----------------------------------------------------------------------------
*/
*****************************************************************************************

clear all
set maxvar 20000
cd "N:\durable\projects\mh_bmi\fullsample"
sysdir
sysdir set PLUS "N:\durable\people\Mandy\Statafiles\ado\plus"
sysdir set PERSONAL "N:\durable\people\Mandy\statafiles\ado\personal"

use "data/imp_x100_17July2022.dta", clear

*for power calculations: need R2 from a non-genetic model for how much variation in the outcomes is explained by BMI:

*across imputations:
foreach var in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr {
*vars to store r2 in:
capture gen `var'_base_r2=.
capture gen `var'_withbmi_r2=.
capture gen `var'_frombmi_r2=.
*then loop through all imputations:
forvalues m = 1/100 {
*without the relevant variable:
reg _`m'_`var' _`m'_rescaled_mothers_bmi_Q1 _`m'_rescaled_fathers_bmi_Q1QF KJONN c_yob i._`m'_AA1124_v2_flipped i._`m'_f_educ_Q1QF_v2_flipped _`m'_mhopkins_Q1 _`m'_fhopkins_QF _`m'_mADHD_Q6 _`m'_fADHD_QF i._`m'_mumsmokes_Q1 i._`m'_dadsmokes_Q1QF i.PARITET_5 /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*get R2:
ereturn list
replace `var'_base_r2=e(r2) in `m'
*with child's BMI added in:
reg _`m'_`var' _`m'_rescaled_moffspring_bmi_Q8 _`m'_rescaled_mothers_bmi_Q1 _`m'_rescaled_fathers_bmi_Q1QF KJONN c_yob i._`m'_AA1124_v2_flipped i._`m'_f_educ_Q1QF_v2_flipped _`m'_mhopkins_Q1 _`m'_fhopkins_QF _`m'_mADHD_Q6 _`m'_fADHD_QF i._`m'_mumsmokes_Q1 i._`m'_dadsmokes_Q1QF i.PARITET_5 /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*get R2:
ereturn list
replace `var'_withbmi_r2=e(r2) in `m'
*subtract:
replace `var'_frombmi_r2=(`var'_withbmi_r2-`var'_base_r2) in `m'
}
}
list *r2 in 1/101

foreach var in z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr {
summ `var'_base_r2 `var'_withbmi_r2 `var'_frombmi_r2 
}

*just the bmi-attributable bit:
summ *_frombmi_r2 

*****************

*for power calculations, the way we have done them: need r2 for the first stage (child's BMI on child's BMI PGS, both of them) net of PCs etc but NOT the parents PGSs:

*across imputations:
foreach pgs in adultbmi childbmi {
*vars to store r2 in:
capture gen `pgs'_base_r2=.
capture gen `pgs'_withpgs_r2=.
capture gen `pgs'_frompgs_r2=.
*then loop through all imputations:
forvalues m = 1/100 {
*without the pgs:
reg _`m'_rescaled_moffspring_bmi_Q8 KJONN c_yob /*i._`m'_AA1124_v2_flipped i._`m'_f_educ_Q1QF_v2_flipped _`m'_mhopkins_Q1 _`m'_fhopkins_QF _`m'_mADHD_Q6 _`m'_fADHD_QF i._`m'_mumsmokes_Q1 i._`m'_dadsmokes_Q1QF i.PARITET_5*/ /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*get R2:
ereturn list
replace `pgs'_base_r2=e(r2) in `m'
*with relevant pgs added in:
reg _`m'_rescaled_moffspring_bmi_Q8 c_z_`pgs'_pgs KJONN c_yob /*i._`m'_AA1124_v2_flipped i._`m'_f_educ_Q1QF_v2_flipped _`m'_mhopkins_Q1 _`m'_fhopkins_QF _`m'_mADHD_Q6 _`m'_fADHD_QF i._`m'_mumsmokes_Q1 i._`m'_dadsmokes_Q1QF i.PARITET_5*/ /*genetic covariates: PCs, center and chip but not batch*/ c_PC* m_PC* f_PC* /*genotyping center*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1 f_genotyping_center2 /*genotyping chip*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6, vce(cluster fid) 
*get R2:
ereturn list
replace `pgs'_withpgs_r2=e(r2) in `m'
*subtract:
replace `pgs'_frompgs_r2=(`pgs'_withpgs_r2-`pgs'_base_r2) in `m'
}
}
list *r2 in 1/101

foreach pgs in adultbmi childbmi {
summ `pgs'_base_r2 `pgs'_withpgs_r2 `pgs'_frompgs_r2 
}

*just the bmi-attributable bit:
summ *_frompgs_r2 

/*. 
summ *_frompgs_r2 

    Variable |        Obs        Mean    Std. dev.       Min        Max
-------------+---------------------------------------------------------
 _frompgs_r2 |          0
adultbmi_f~2 |        100    .0344899    .0017436   .0301369   .0399378
childbmi_f~2 |        100    .0517473    .0021822   .0457711   .0572893
*/

