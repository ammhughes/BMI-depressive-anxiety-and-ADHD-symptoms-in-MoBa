
*STEP 1: PREPARE PHENOTYPE DATA

clear all
*cd to project folder
cd "N:\durable\projects\mh_bmi\fullsample"

*to use the user-writen files, change the path to the PLUS part:
sysdir
/*
   STATA:  C:\Program Files\Stata16MP\
    BASE:  C:\Program Files\Stata16MP\ado\base\
    SITE:  C:\Program Files\Stata16MP\ado\site\
    PLUS:  C:\Users\p471-ammh\ado\plus\
PERSONAL:  C:\Users\p471-ammh\ado\personal\
OLDPLACE:  c:\ado\
*/
sysdir set PLUS "N:\durable\people\Mandy\Statafiles\ado\plus"
**********************************************************************************************

*PRELIMINARIES:
*UPDATE 21/06/2022: REMOVE THE NEW CONSENT WITHDRAWALS FROM THIS MONTH
*This is done by restricting to people in the updated SV_INFO file.
*needs to be done separately for mothers/children and fathers.
*The correct file was taken from the SPSS folder, using the file called "N:\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\SPSS\PDB2306_SV_INFO_v12.sav"
*this is ALWAYS the name of the most up to date file; when dates are added to the name of an SV_INFO file, that refers to the date when it was taken out of use.
*NB: couldn't import directly into stata, so opened in SPSS and saved as a csv, then imported that:
import delimited "N:\durable\projects\mh_bmi\fullsample\PDB2306_SV_INFO_v12.csv"
rename preg_id_2306 PREG_ID_2306
rename m_id_2306 M_ID_2306
rename f_id_2306 F_ID_2306
save "N:\durable\projects\mh_bmi\fullsample\PDB2306_SV_INFO_v12.dta", replace
*Then split this into one for excluding mothers and children (on PREG_ID_2306) and the one for excluding fathers (on F_ID_2306).
*For mothers/children, apply filter at the end of this file, since less error-prone than applying individually at each merge.
*For fathers, do this twice: filter the first father's questionnaire before merging in, and the second father's questionnaire, before merging in.
use "N:\durable\projects\mh_bmi\fullsample\PDB2306_SV_INFO_v12.dta", clear
preserve
keep PREG_ID_2306 M_ID_2306
count if M_ID_2306==" "
*none missing
count
*112,198
*formerly 112,265
*make a flag variable:
gen mothers_consent=1
save "N:\durable\projects\mh_bmi\fullsample\PDB2306_SV_INFO_v12_mothers_children.dta", replace
restore
keep PREG_ID_2306 F_ID_2306
keep if F_ID_2306!=" "
count
*86,449
*formerly 86,509
*make a flag variable:
gen fathers_consent=1
save "N:\durable\projects\mh_bmi\fullsample\PDB2306_SV_INFO_v12_fathers.dta", replace
*apply these as you go.

**************************************************************
*Start with the BIRTH REGISTRY FILE, which you'll always need to do a pheno-geno linkage.
/*
## Step 1: Read in MoBa Birth Registry file

## This file includes all births in MoBa (in version 12 there are 114,143 children from 112,644 pregnancies)
## key variables: PREG_ID_2306 = main ID for linkage, unique to pregnancy, KJONN=child's sex, BARN_NR=birth order, 
## FAAR=birth year, MORS_ALDER= mothers age
*/
use "N:\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\Statafiles\PDB2306_MBRN_541_v12.dta", clear
keep PREG_ID_2306 BARN_NR KJONN FAAR MORS_ALDER FARS_ALDER PARITET_5 SIVST_2 VEKT LENGDE ZSCORE_BW_GA  
count
*114,143

*FROM BIRTH REGISTRY FILE, to add to imputation:
fre MORS_ALDER PARITET_5
*recode the <17s group to 16 and the >45 group to 46, then can at a pinch use as continous
recode MORS_ALDER 917=16
recode MORS_ALDER 945=46
*PARITET_5 can be used as-is
*paternal age:
fre FARS_ALDER
recode FARS_ALDER 918=17
recode FARS_ALDER 959=60


*birth weight and length: will need for imputation.
*trim at 4sd from the mean:
foreach var in VEKT LENGDE {
summ `var', det
return list
replace `var'=. if `var'<r(mean)-(4*r(sd))
replace `var'=. if `var'>r(mean)+(4*r(sd))
summ `var', det
}

*some people in questionnaires but not in birth registry file.
*make a flag so you can drop them later:
gen present_MBRN=1

**********************************************************************************************
*PREPARE PHENOTYPE DATA FROM QUESTIONNAIRES.

*Height and weight for all members of trios, from various places:

*Q1: by mothers at 17 weeks (no BARN_NR in here as pre-birth!)
merge m:1 PREG_ID_2306 using "N:\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\Statafiles\PDB2306_Q1_v12.dta"
*441 unmatched - not in the birth registry file
*how bizarre
*will be a problem as they don't have BARN_NR:
count if BARN_NR==.
*441
*cannot use them for anything, so drop will need to drop them
*DO THIS AT THE END THOUGH - otherwise will need to repeat at every merge

*presence flag:
gen present_Q1=1 if _merge!=1
fre present_Q1

*clean vars for BMI:
*mother's own height and pre-pregnancy weight
*mother's report of partner's height and weight
foreach var in AA85 AA87 AA88 AA89 {
summ `var', det
return list
replace `var'=. if `var'<r(mean)-(4*r(sd))
replace `var'=. if `var'>r(mean)+(4*r(sd))
summ `var', det
}
*now make BMI:
gen mothers_bmi_Q1=AA85/((AA87/100)^2)
label variable mothers_bmi_Q1 "Q1 mother's pre-preg bmi from self-report h & w"
summ mothers_bmi_Q1, det
*looks acceptable
gen mfathers_bmi_Q1=AA89/((AA88/100)^2)
label variable mfathers_bmi_Q1 "Q1 father's bmi from mother's self-report h & w"
summ mfathers_bmi_Q1, det
*looks acceptable

*Also gather demographic variables for imputation:
fre AA1123 AA1124 AA1125 AA1126 AA1127

*marital status:
fre AA1123 
*remove the >1 box checked people
recode AA1123 0=.
*some tiny groups which might mess with the impuation but let's see. Can merge if necessary

*mother's educ, finished and ongoing, father's " ":
foreach var in AA1124 AA1125 AA1126 AA1127{
recode `var' 0=.
}
*use just completed educ for now? AA1124 AA1126
*actually, a fair few have this missing but filled in the ongoing one:
fre AA1127 if AA1126 ==.
fre AA1125 if AA1124 ==.
*so for those cases, replace the . with the value from the ongoing education variable, minus 1

*mother's education: new version:
clonevar AA1124_v2=AA1124
fre AA1124_v2
*where missing, replace with 'ongoing' value minus 1:
fre AA1125 if AA1124==.
replace AA1124_v2=AA1125-1 if AA1124_v2==.
*nb: now includes 0, for people working towards primary education
*collapse the (now three) tiny categories?
*update: keep 1 and 2 separate (but merge 0 into 1)
recode AA1124_v2 0=1
fre AA1124_v2
*for imputation: flip to have larger baseline groups.
*label for the flipped version:
label define educ_v2_flipped 6"9-year elementary ed" 5"up to 2: FE 1-2 yrs" 4"FE - vocational" 3"FE 3yrs" 2"Higher Ed, <=4yrs" 1"Higher Ed, >4yrs"
gen AA1124_v2_flipped=7-AA1124_v2
label values AA1124_v2_flipped educ_v2_flipped
tab AA1124_v2 AA1124_v2_flipped

*father's education: new version:
clonevar AA1126_v2=AA1126
*where missing, replace with 'ongoing' value minus 1:
fre AA1127 if AA1126==.
replace AA1126_v2=AA1127-1 if AA1126_v2==.
*nb: now includes 0, for people working towards primary education
*collapse the (now three) tiny categories:
*update: keep 1 and 2 separate (but merge 0 into 1)
recode AA1126_v2 0=1
fre AA1126_v2
*for imputation: flip to have larger baseline groups.
gen AA1126_v2_flipped=7-AA1126_v2
label values AA1126_v2_flipped educ_v2_flipped
tab AA1126_v2 AA1126_v2_flipped

*income, financial difficulty
*def add to imputation models:
fre AA1315 AA1316 AA1317
*last one needs fixing
recode AA1317 0=.
*change the dont know in mother's report of partner's income to missing:
recode AA1316  8=.

*there's also stuff about housing type, but leave for now.
fre  AA1318- AA1328 

*also non-Norwegian languages - ethnicity marker?
*Q72-80

*full variables about smoking, alcohol and drugs
fre AA1348-AA1473

*smoking?
fre AA1328-AA1357
*clean them

*make a grouped smoking variable for mums:
capture drop mumsmokes_Q1
gen mumsmokes_Q1=.
label var mumsmokes_Q1 "mother's smoking status at 17 weeks"
label define mumsmokes 0"never" 1"pre-pregnanacy" 2"during, sometimes" 3"during, daily"
*using these two:
fre AA1355 AA1356
label values mumsmokes mumsmokes
replace mumsmokes_Q1=0 if AA1355==1
*then for other groups:
replace mumsmokes=1 if AA1355!=1 & AA1356==1 
replace mumsmokes=2 if AA1356==2 
replace mumsmokes=3 if AA1356==3 
*check it
fre mumsmokes
tab mumsmokes AA1355
tab mumsmokes AA1356
tab AA1355 AA1356 if mumsmokes==.

*mother's smoking: collapse the 3rd and fourth categories?
fre mumsmokes_Q1
gen mumsmokes_Q1_v2=mumsmokes_Q1
recode mumsmokes_Q1_v2 3=2
label define mumsmokes_Q1_v2 0"never" 1"pre-pregnancy" 2"during"
label values mumsmokes_Q1_v2 mumsmokes_Q1_v2
fre mumsmokes_Q1_v2

*and for dads?
tab AA1353 AA1354
gen dadsmokes_Q1=.
label var dadsmokes_Q1 "mother report of dad's smoking status at 17 weeks"
label define dadsmokes 0"not before pregnancy" 1"stopped during pregnanacy" 2"during preganacy"
label values dadsmokes dadsmokes
replace dadsmokes_Q1=0 if AA1353==1
replace dadsmokes_Q1=1 if AA1353==2 & AA1354==1
replace dadsmokes_Q1=2 if AA1354==2
fre dadsmokes
tab AA1353 AA1354 if dadsmokes==0
tab AA1353 AA1354 if dadsmokes==1
tab AA1353 AA1354 if dadsmokes==2
*looks ok.

*parents' psych measures?
*have you had Anorexia/ bulimia/ other eating disorder // depression // anxiety
*pairs: before pregnancy, during pregancy
fre AA806 AA807 AA869 AA870 AA878 AA879
*looks like this was a tick vs no tick, so no way to distinguish no from missing among people present at Q1.
*identify people for whom it's definitely a genuine missing by using the _merge variable: they are unmatched from master, i.e. _merge==1. DO NOT CODE .=2 FOR THESE PEOPLE!
*Most will be no's, so assume that. Label as 2.
label define AA806 1"yes" 2"not reported"
foreach var of varlist AA806 AA807 AA869 AA870 AA878 AA879 {
recode `var' .=2 if _merge!=1
label values `var' AA806
fre `var'
}
*Eating disorder screening questions
fre AA1475 AA1476 AA1477 AA1478 AA1479 AA1480 AA1481 AA1482 AA1483 AA1484 AA1485 AA1486 AA1487 AA1488 
*fix multiple ticks
foreach var of varlist AA1475 AA1476 AA1477 AA1478 AA1479 AA1480 AA1481 AA1482 AA1483 AA1484 AA1485 AA1486 AA1487 AA1488 {
recode `var' 0=.
fre `var'
}
*Satisfaction with life scale
fre AA1527 AA1528 AA1529 AA1530 AA1531 
*fix multiple ticks
foreach var of varlist AA1527 AA1528 AA1529 AA1530 AA1531 {
recode `var' 0=.
fre `var'
}
*Hopkins checklist
fre AA1548 AA1549 AA1550 AA1551 AA1552
*fix multiple ticks
foreach var of varlist AA1548 AA1549 AA1550 AA1551 AA1552 {
recode `var' 0=.
recode `var' 12/56=.
fre `var'
}

*********************************************************************************
*to add to imputation models: 
*social, smoking
fre AA1123 AA1124 AA1126 AA1315 AA1316 AA1317 mumsmokes_Q1 dadsmokes_Q1
*MH
fre AA806 AA807 AA869 AA870 AA878 AA879 ///
AA1475 AA1476 AA1477 AA1478 AA1479 AA1480 AA1481 AA1482 AA1483 AA1484 AA1485 AA1486 AA1487 AA1488 ///
AA1527 AA1528 AA1529 AA1530 AA1531 ///
AA1548 AA1549 AA1550 AA1551 AA1552
*********************************************************************************

*trim:
keep PREG_ID_2306 BARN_NR present_MBRN KJONN FAAR MORS_ALDER FARS_ALDER PARITET_5 SIVST_2 VEKT LENGDE ZSCORE_BW_GA VERSJON_SKJEMA1_TBL1 AA85 AA86 AA87 AA88 AA89 mothers_bmi_Q1 mfathers_bmi_Q1 ///
AA1123 AA1124 AA1125 AA1126 AA1127 AA1315 AA1316 AA1317 AA1318- AA1328 AA1348-AA1473 mumsmokes_Q1 mumsmokes_Q1_v2 dadsmokes_Q1 ///
AA806 AA807 AA869 AA870 AA878 AA879 ///
AA1475 AA1476 AA1477 AA1478 AA1479 AA1480 AA1481 AA1482 AA1483 AA1484 AA1485 AA1486 AA1487 AA1488 ///
AA1527 AA1528 AA1529 AA1530 AA1531 ///
AA1548 AA1549 AA1550 AA1551 AA1552 /// 
AA1124_v2* AA1126_v2* ///
present_Q1 

*save as temp file while you sort out the father's file:
save "temp_phenotypes.dta", replace

count if BARN_NR==.
*441

**********************************************************************************************
*open QF: first fathers questionnaire, 17 week
use "N:\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\Statafiles\PDB2306_QF_v12.dta"
*just keep what you need
keep PREG_ID_2306 FF15 FF16-FF363 FF333 FF334 FF335 FF336 FF214-FF474 ///
FF251 FF252 FF253 FF254 FF255 FF256 FF257 FF258 ///
FF259 FF260 FF261 FF262 FF263 FF264 FF478 FF479 ///
FF266 FF267 FF268  ///
FF146 FF147 FF148 FF386 FF387 FF388 FF389 FF390 FF391 FF392 FF393 FF394 FF395 FF386 FF397 FF398 FF399 FF400 ///
FF480- FF529  ///
FF269 FF270 FF271 FF272 FF273 ///
FF535 FF536 FF537 FF538 FF539 FF540 ///
FF277 FF278 FF279 FF280 FF281 FF282
*presence flag:
gen present_QF=1
*ERASE THIS DATA for any fathers who withdrew consent:
merge 1:1 PREG_ID_2306 using "N:\durable\projects\mh_bmi\fullsample\PDB2306_SV_INFO_v12_fathers.dta"
keep if fathers_consent==1
*264 excluded
*then merge back into growing file:
drop _merge
merge 1:m PREG_ID_2306 using "temp_phenotypes.dta"

*MH: reported health problems: 
*in each case: yes/no for ever had, age started, and age stopped
*15. Sleep problems 
fre FF146 FF147 FF148
*25. ADHD 
fre FF386 FF387 FF388 
*26. Anorexia/bulimia/eating disorders 
fre FF389 FF390 FF391 
*27. Manic depressive illness 
fre FF392 FF393 FF394 
*28. Schizophrenia
fre FF395 FF386 FF397 
*29. Other long-term mental illnesses or health problems 
fre FF398 FF399 FF400

*ever/never:
fre FF146 FF386 FF389 FF392 FF395 FF398 if _merge!=1
*for ever had, was a tick vs no tick, so no way to distinguish no from missing among people present at F1.
*identify people for whom it's definitely a genuine missing by using the _merge variable: they are unmatched from master, i.e. _merge==1. DO NOT CODE .=2 FOR THESE PEOPLE!
*Most will be no's, so assume that. Label as 2.
*label define AA806 1"yes" 2"not reported"
foreach var of varlist FF146 FF386 FF389 FF392 FF395 FF398 {
recode `var' .=2 if _merge!=1
label values `var' AA806
fre `var'
}
*cell counts so small this will be pretty useless for imputation

*FF333 FF334: current height and weight
*FF335 FF336: heaviest and lightest they have ever been
foreach var in FF333 FF334 FF335 FF336 {
summ `var', det
return list
replace `var'=. if `var'<r(mean)-(4*r(sd))
replace `var'=. if `var'>r(mean)+(4*r(sd))
summ `var', det
}
*bmi from current height and weight
gen fathers_bmi_QF=FF334/((FF333/100)^2)
label variable fathers_bmi_QF "QF father's bmi from self-report current h & w"
summ fathers_bmi_QF, det

*****************************************************
*COMBINED BMI report from around this time: from dads or, failing that, mum's report:
gen fathers_bmi_Q1QF=fathers_bmi_QF
replace fathers_bmi_Q1QF=mfathers_bmi_Q1 if fathers_bmi_Q1QF==.
summ fathers_bmi_Q1QF

*for methods section:
summ fathers_bmi_QF
summ mfathers_bmi_Q1
pwcorr fathers_bmi_QF mfathers_bmi_Q1
*0.98
*check for % replaced in pre-imp file (i.e, analytic sample)
*****************************************************

*demographic:
fre FF15 FF16-FF363

*clean:
*marital status:
recode FF15 0=.
*income, grouped:
recode FF341 0=.
*completed education
recode FF16 0=.
*ongoing education:
recode FF17 0=.

*like bmi, smoking, and education, make a variable for father's income beginning with self-report and subbing in mum's report
*check labelling is the same:
fre AA1316 FF341
gen fathers_inc_grp=FF341
replace fathers_inc_grp=AA1316 if fathers_inc_grp==.
label values fathers_inc_grp FF341
fre fathers_inc_grp

*FATHER'S EDUCATION LEVEL: as for smoking and BMI, make a composite variable of dad'd report of mum's where dad'd not available.
*confirm labelling is identical:
fre FF16 FF17 AA1126
clonevar f_educ_Q1QF_v2=FF16
*add the step of incorporating ongoing study:
fre FF17 if FF16==.
replace f_educ_Q1QF_v2=FF17-1 if FF16==.
*recode the new 0s to 1
recode f_educ_Q1QF_v2 0=1
fre f_educ_Q1QF_v2
*flip to have larger baseline group:
gen f_educ_Q1QF_v2_flipped=7-f_educ_Q1QF_v2
*label it
label values f_educ_Q1QF_v2_flipped educ_v2_flipped
tab f_educ_Q1QF_v2 f_educ_Q1QF_v2_flipped
*then sub in mother's report, including from ongoing study, FLIPPED VERSION:
fre AA1126_v2_flipped if f_educ_Q1QF_v2_flipped==.
replace f_educ_Q1QF_v2_flipped=AA1126_v2_flipped if f_educ_Q1QF_v2_flipped==.
fre f_educ_Q1QF_v2_flipped

*there's also stuff on employment status and benefits but could be tricky, so leave for now. Keep in the file

*smoking, alcohol, drugs: Q49 to 61
fre FF214-FF474

*smoking status
tab FF214 FF215
gen dadsmokes_QF=.
label var dadsmokes_QF "father report of dad's smoking status at 17 weeks"
label define dadsmokes_QF 0"not before pregnancy" 1"stopped during pregnanacy" 2"during preganacy"
label values dadsmokes_QF dadsmokes_QF
replace dadsmokes_QF=0 if FF214==1
replace dadsmokes_QF=1 if FF214==2 & FF215==1
replace dadsmokes_QF=2 if FF215==2 | FF215==3
fre dadsmokes_QF
tab FF214 FF215 if dadsmokes_QF==0
tab FF214 FF215 if dadsmokes_QF==1
tab FF214 FF215 if dadsmokes_QF==2
*looks ok.

*Dad'd smoking: as with BMI, make a combined report from dad's own questionnaire or failing that the mum's
gen dadsmokes_Q1QF=dadsmokes_QF
replace dadsmokes_Q1QF=dadsmokes_Q1 if dadsmokes_Q1QF==.
label values dadsmokes_Q1QF dadsmokes_QF
fre dadsmokes_Q1QF

*MH stuff:

*Q66: Hopkins checklist
fre FF251 FF252 FF253 FF254 FF255 FF256 FF257 FF258
*fix multiple ticks
foreach var of varlist FF251 FF252 FF253 FF254 FF255 FF256 FF257 FF258 {
recode `var' 0=.
fre `var'
}
*67-68: Lifetime history of depression:
fre FF259 FF260 FF261 FF262 FF263 FF264 
fre FF478 FF479
*fix multiple ticks
foreach var of varlist FF259 FF260 FF261 FF262 FF263 FF264  {
recode `var' 0=.
fre `var'
}
*For the last 2 items - how many times have you had 3 of these symptoms at once, and for how many weeks,
*the 0s are legit. Leave those.

*69 Rosenberg self-esteem:
fre FF266 FF267 FF268
*fix multiple ticks
foreach var of varlist FF266 FF267 FF268  {
recode `var' 0=.
fre `var'
}
*70 Big 5 Personality
fre FF480- FF529
foreach var of varlist FF480- FF529  {
recode `var' 0=.
fre `var'
}
*71 Satisfaction with life scale
fre FF269 FF270 FF271 FF272 FF273
foreach var of varlist FF269 FF270 FF271 FF272 FF273  {
recode `var' 0=.
fre `var'
}
*72 Adult ADHD
fre FF535 FF536 FF537 FF538 FF539 FF540 
foreach var of varlist FF535 FF536 FF537 FF538 FF539 FF540  {
recode `var' 0=.
fre `var'
}
*78: Differential emotional scale, enjoyment and anger subscales
fre FF277 FF278 FF279 FF280 FF281 FF28
foreach var of varlist FF277 FF278 FF279 FF280 FF281 FF28  {
recode `var' 0=.
fre `var'
}

/*add to imputation: FF15 FF16 FF341 dadsmokes_QF FF251 FF252 FF253 FF254 FF255 FF256 FF257 FF258 ///
FF259 FF260 FF261 FF262 FF263 FF264 FF478 FF479 ///
FF266 FF267 FF268  ///
FF480- FF529  ///
FF269 FF270 FF271 FF272 FF273 ///
FF535 FF536 FF537 FF538 FF539 FF540 ///
FF277 FF278 FF279 FF280 FF281 FF282 */

drop _merge

***********************************************************************************************
*Q3: parental MH vars for imputation

merge m:1 PREG_ID_2306 using "N:\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\Statafiles\PDB2306_Q3_v12.dta", ///
keepus (PREG_ID_2306 CC676 CC677 CC678 CC679 CC680 CC688 CC689 CC690 CC691 CC692 ///
CC1202 CC1203 CC1204 CC1205 CC1206 CC1207 CC1208 CC1209 ///
CC1210 CC1211 CC1212 CC1213 CC1214 CC1215 ///
CC1224 CC1225 CC1226 CC1227 CC1228 ///
CC1229 CC1230 CC1231 CC1232)

*presence flag:
gen present_Q3=1 if _merge!=1
fre present_Q3

*Q52, pt 29: depression at diff stages of preganancy
fre CC676 CC677 CC678 CC679 CC680 
*pt 30: other psychological problem, at different stages
fre CC688 CC689 CC690 CC691 CC692 

*looks like this was a tick vs no tick, so no way to distinguish no from missing among people present at Q1.
*identify people for whom it's definitely a genuine missing by using the _merge variable: they are unmatched from master, i.e. _merge==1. DO NOT CODE .=2 FOR THESE PEOPLE!
foreach var of varlist CC676 CC677 CC678 CC679 CC680 CC688 CC689 CC690 CC691 CC692 {
recode `var' .=2 if _merge!=1
label values `var' AA806
fre `var'
}
*123. Hopkins checklist
fre CC1202 CC1203 CC1204 CC1205 CC1206 CC1207 CC1208 CC1209
*fix multiple ticks
foreach var of varlist CC1202 CC1203 CC1204 CC1205 CC1206 CC1207 CC1208 CC1209 {
recode `var' 0=.
fre `var'
}
*124. Emotion: Enjoyment and Anger: Differential Emotional Scale (DES), Enjoyment and Anger Subscales
fre CC1210 CC1211 CC1212 CC1213 CC1214 CC1215
*fix multiple ticks
foreach var of varlist CC1210 CC1211 CC1212 CC1213 CC1214 CC1215 {
recode `var' 0=.
fre `var'
}
*126 Satisfaction with life scale
fre CC1224 CC1225 CC1226 CC1227 CC1228
*fix multiple ticks
foreach var of varlist CC1224 CC1225 CC1226 CC1227 CC1228 {
recode `var' 0=.
fre `var'
}
*127. The Rosenberg Self-Esteem Scale (RSES)
fre CC1229 CC1230 CC1231 CC1232
*fix multiple ticks
foreach var of varlist CC1229 CC1230 CC1231 CC1232 {
recode `var' 0=.
fre `var'
}
drop _merge
**********************************************************************************************

*Q4 parental MH vars for imputation
merge 1:1 PREG_ID_2306 BARN_NR using "N:\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\Statafiles\PDB2306_Q4_6months_v12.dta", ///
keepus (PREG_ID_2306 BARN_NR DD537 DD538 DD539 ///
DD794 DD795 DD796 DD797 DD798 DD799 ///
DD800 DD801 DD802 DD803 DD804 ///
DD827 DD828 DD829 DD830 DD831 DD832 ///
DD833 DD834 DD835  DD836 ///
DD837 DD838 DD839 DD840 DD841 DD842 DD843 DD844)

*presence flag:
gen present_Q4=1 if _merge!=1
fre present_Q4

*MH probs in last part of pregnancy and after the births
fre DD537 DD538 DD539
*like in Q1, this was a tick vs no tick, so no way to distinguish no from missing.
*Most will be no's, so assume that. Label as 2.
*label define AA806 1"yes" 2"not reported"
foreach var of varlist DD537 DD538 DD539 {
recode `var' .=2 if _merge!=1
label values `var' AA806 
fre `var'
}
*89. Emotion: Enjoyment and Anger: Differential Emotional Scale (DES), Enjoyment and Anger Subscales
fre DD794 DD795 DD796 DD797 DD798 DD799
*fix multiple ticks
foreach var of varlist DD794 DD795 DD796 DD797 DD798 DD799 {
recode `var' 0=.
fre `var'
}
*90. Life Satisfaction: The Satisfaction With Life Scale (SWLS)
fre DD800 DD801 DD802 DD803 DD804
*fix multiple ticks
foreach var of varlist DD800 DD801 DD802 DD803 DD804 {
recode `var' 0=.
fre `var'
}
*92. Postnatal Depression: Edinburgh Postnatal Depression Scale (EPDS)
fre DD827 DD828 DD829 DD830 DD831 DD832
*fix multiple ticks
foreach var of varlist DD827 DD828 DD829 DD830 DD831 DD832 {
recode `var' 0=.
fre `var'
}
*93. Rosenberg Self Esteem Scale: The Rosenberg Self-Esteem Scale (RSES)
fre DD833 DD834 DD835  DD836
*fix multiple ticks
foreach var of varlist DD833 DD834 DD835  DD836 {
recode `var' 0=.
fre `var'
}
*94. Hopkins Checklist
fre DD837 DD838 DD839 DD840 DD841 DD842 DD843 DD844
*fix multiple ticks
foreach var of varlist DD837 DD838 DD839 DD840 DD841 DD842 DD843 DD844 {
recode `var' 0=.
fre `var'
}

drop _merge
**********************************************************************************************

*Q5: Lots of stuff on both parents and kids for imputation

*EE427/EE1005 EE432/997 EE986/406 EE986/406 
 
*EE915 EE953 not in the data - why? 
 
merge 1:1 PREG_ID_2306 BARN_NR using "N:\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\Statafiles\PDB2306_Q5_18months_v12.dta", ///
keepus (PREG_ID_2306 BARN_NR EE583 EE584 EE603 EE604 EE605 EE606 EE607 EE608 EE609 ///
EE886  EE433 EE887  EE888  EE889 EE890 EE891 EE892 EE893 EE894 EE895 EE896 EE960 EE897 ///
EE898 EE884 EE885 ///
EE427 EE1005 EE434 EE429 EE430 EE431 EE998 EE432 EE997 EE433 EE428 EE1006 EE900 EE1000 EE879 EE901 EE882 EE986 EE406 EE1001 EE880 EE881 EE1002 EE899 EE833 EE902 ///
EE435 EE961 EE903 EE904 EE905 EE438 EE439 EE962 EE442 EE446 EE447 EE448 EE963 EE964 EE906 EE440 EE907 EE908 EE909 ///
EE910 EE911 EE912 EE913 EE1007 EE914 ///
EE925 EE926 EE927 EE928 EE929 EE930 EE931 EE932 EE933 EE934 EE935 EE936 EE937 EE938 EE939 EE940 EE941 EE942 EE943 ///
EE628 EE629 EE630 EE631 EE632 EE633 ///
EE634 EE635 EE636 EE637 ///
EE638 EE639 EE640 EE641 EE642 EE643 EE644 EE645 ///
EE671 EE672 EE673 EE674 EE675 EE676 EE677 EE678 EE679 EE680 EE681 EE682 EE683 EE684 EE685 EE686 EE687 EE688 EE689 EE690 EE691 EE692 EE693 EE694 EE695 EE696 ///
EE416 EE417 EE418 EE419 EE420 EE421 EE422 EE423 EE424 EE425 EE426 EE877 EE878)

*presence flag:
gen present_Q5=1 if _merge!=1
fre present_Q5

***Household financial strain: income not asked about
fre EE583 EE584

*simplify EE584, for imputation
recode EE584 4=3
label define EE584_new 1"never" 2"yes, but infrequently" 3"merged: sometimes and often"
label values EE584 EE584_new
fre EE584


*Smoking
*Q89
fre EE603 EE604 EE605 EE606
*Alcohol
*90-91
fre EE607 EE608 EE609

*multiple ticks (not for cigs/day, where 0 means something else)
foreach var of varlist EE583 EE584 EE603 EE605 EE607 EE608 EE609 {
recode `var' 0=.
fre `var'
}

**About the child:

*EAS: Temperament
*fix multiple ticks:
foreach var of varlist EE416 EE417 EE418 EE419 EE420 EE421 EE422 EE423 EE424 EE425 EE426 EE877 EE878 {
	recode `var' 0=.
}

*Q35/36: ESAT
*Selective questions from Early Screening of Autistic Traits Questionnaire (ESAT)
*NB: questions VARY A LOT BETWEEN THE A B AND C VERSIONS...
fre EE886  EE433 EE887 EE888 EE889 EE890 EE891 EE892 EE893 EE894 EE895 EE896 EE960 EE897
*fix multiple ticks
foreach var of varlist EE886 EE433 EE887 EE888 EE889 EE890 EE891 EE892 EE893 EE894 EE895 EE896 EE960 EE897 {
recode `var' 0=.
fre `var'
}
*These three are from other scales – the SCQ-Social communication questionnaire and Communication and Symbolic Behaviour Scales – and seem to have been asked for everyone:
fre EE898 EE884 EE885
*fix multiple ticks
foreach var of varlist EE898 EE884 EE885 {
recode `var' 0=.
fre `var'
}
*35/36. Autistic Traits: M-CHAT
*The Modified Checklist for Autism in Toddlers (M-CHAT)
*NB: questions VARY A LOT BETWEEN THE A B AND C VERSIONS...
*These three are pairs: EE427/EE1005 EE432/997 EE986/406
tab EE427 EE1005, mi
tab EE432 EE997, mi
tab EE986 EE406, mi
*not sure what that means...
gen EE432_EE997=EE997
replace EE432_EE997=EE432 if EE997==.
tab EE432_EE997

fre EE427 EE1005 EE434 EE429 EE430 EE431 EE998 EE432 EE997 EE432_EE997 EE433 EE428 EE1006 EE900 EE1000 EE879 EE901 EE882 EE986 EE406 EE1001 EE880 EE881 EE1002 EE899 EE833 EE902
*fix multiple ticks
foreach var of varlist EE427 EE1005 EE434 EE429 EE430 EE431 EE998 EE432 EE997 EE432_EE997 EE433 EE428 EE1006 EE900 EE1000 EE879 EE901 EE882 EE986 EE406 EE1001 EE880 EE881 EE1002 EE899 EE833 EE902 {
recode `var' 0=.
fre `var'
}
*37. Child Behaviour CheckList (CBCL)
fre EE435 EE961 EE903 EE904 EE905 EE438 EE439 EE962 EE442 EE446 EE447 EE448 EE963 EE964 EE906 EE440 EE907 EE908 EE909
*fix multiple ticks
foreach var of varlist EE435 EE961 EE903 EE904 EE905 EE438 EE439 EE962 EE442 EE446 EE447 EE448 EE963 EE964 EE906 EE440 EE907 EE908 EE909 {
recode `var' 0=.
fre `var'
}
*18 months Child Behaviour CheckList (CBCL): make binary versions
foreach var of varlist EE427 EE1005 EE434 EE429 EE430 EE431 EE998 EE432 EE997 EE432_EE997 EE433 EE428 EE1006 EE900 EE1000 EE879 EE901 EE882 EE986 EE406 EE1001 EE880 EE881 EE1002 EE899 EE833 EE902 {
	gen `var'_bin=`var'
	recode `var'_bin 3=2
}

*40. Mother’s concerns
*EE915 EE953 for other concerns not included. 
fre EE910 EE911 EE912 EE913 EE1007 EE914  
*fix multiple ticks
foreach var of varlist EE910 EE911 EE912 EE913 EE1007 EE914  {
recode `var' 0=.
fre `var'
}

******************************
*About the mum:
*****************************
*69-71. Eating Disorders

fre EE925 EE926 EE927 EE928 EE929 EE930 EE931 EE932 EE933 EE934 EE935 EE936 EE937 EE938 EE939 EE940 EE941 EE942 EE943
*fix multiple ticks
foreach var of varlist EE925 EE926 EE927 EE928 EE929 EE930 EE931 EE932 EE933 EE934 EE935 EE936 EE937 EE938 EE939 EE940 EE941 EE942 EE943  {
recode `var' 0=.
fre `var'
}

*97. Emotion: Enjoyment and Anger: Differential Emotional Scale (DES), Enjoyment and Anger Subscales
fre EE628 EE629 EE630 EE631 EE632 EE633
*fix multiple ticks
foreach var of varlist EE628 EE629 EE630 EE631 EE632 EE633  {
recode `var' 0=.
fre `var'
}

*98. The Rosenberg Self-Esteem Scale: Selective questions from the Rosenberg Self-Esteem Scale (RSES)
fre EE634 EE635 EE636 EE637
foreach var of varlist EE634 EE635 EE636 EE637  {
recode `var' 0=.
fre `var'
}
*99. Depression/Anxiety: Selective items from the (Hopkins) Symptoms Checklist-25 (SCL-25)
fre EE638 EE639 EE640 EE641 EE642 EE643 EE644 EE645
foreach var of varlist EE638 EE639 EE640 EE641 EE642 EE643 EE644 EE645  {
recode `var' 0=.
fre `var'
}
*101-107. World Health Organization’s Quality of Life Instrument
fre EE671 EE672 EE673 EE674 EE675 EE676 EE677 EE678 EE679 EE680 EE681 EE682 EE683 EE684 EE685 EE686 EE687 EE688 EE689 EE690 EE691 EE692 EE693 EE694 EE695 EE696
foreach var of varlist EE671 EE672 EE673 EE674 EE675 EE676 EE677 EE678 EE679 EE680 EE681 EE682 EE683 EE684 EE685 EE686 EE687 EE688 EE689 EE690 EE691 EE692 EE693 EE694 EE695 EE696  {
recode `var' 0=.
fre `var'
}

drop _merge

**********************************************************************************************

*Questionnaire 6, Age 3:
*there is info on the child's height and weight, the age the report relates to and who made the measurement, but is this even relevant?
merge 1:1 PREG_ID_2306 BARN_NR using "N:\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\Statafiles\PDB2306_Q6_3yrs_v12.dta", ///
keepus (GG435 GG501 GG502 GG479-GG488 GG37 GG38 GG39 GG40 GG93 GG94 GG95 GG96 GG578 GG579 GG580 GG581 GG109   GG105 GG106 GG107 GG108 GG101 GG582 GG102 GG583 GG103 GG584 GG104 GG585 GG109 GG110 GG111 GG112 ///
GG222 GG223 GG224 GG225 GG237 GG238 GG239 GG240 GG241 GG242 ///
GG227 GG228 GG229 GG230 ///
GG231 GG232 GG233 GG234 GG235 GG236 ///
GG243 GG244 GG245 GG246 GG247 GG248 GG249 GG250 GG251 GG592 GG252 GG253 GG254 ///
GG255 GG257 GG258 GG259 GG260 GG261 GG262 GG263 GG264 GG265 GG266 GG267 GG268 GG269 GG270 GG271 GG272 GG273 GG274 GG256 GG275 GG276 GG277 GG278 GG279 GG280 GG281 GG282 GG283 GG284 GG285 GG286 GG287 GG288 GG289 GG290 GG291 GG292 GG293 GG294 ///
GG295 GG296 GG297 GG298 ///
GG299 GG300 GG301 GG302 GG303 GG304 GG305 GG306 GG307 GG308 GG309 GG310 GG311 GG312 ///
GG313 GG314 GG315 GG316 GG317 GG318 GG319 GG320 GG321 GG322 GG323 GG324 GG325 GG326 GG327 GG328 GG329 GG330 GG331 GG332 GG333 GG334 GG335 GG336 GG337 GG338 ///
GG339 GG340 GG341 GG342 GG343 GG344 GG345 GG346 GG347 GG348 ///
GG349 GG350 GG351 GG352 GG353 GG354 GG355 GG356 GG357 GG358 GG359 GG360 GG361 GG362 GG363 GG364 GG365 GG366 GG367 GG236 ///
GG380 GG381 GG382 GG594 GG595  ///
GG388 GG389 ///
GG452 GG634 GG635 GG636 GG637 GG453 GG638 GG639 GG640 GG641 GG454 GG642 GG643 GG644 GG645 GG455 GG646 GG647 GG648 GG649 GG456 GG650 GG651 GG652 GG653 GG457 GG654 GG655 GG656 GG657 ///
GG462 GG658 GG659 GG660 GG661 GG662 GG663  ///
GG491 GG492 GG493 GG494  GG495 GG496 GG497 GG498 GG499  GG500 GG501 GG502 ///
GG503 GG504 GG505 GG506 GG507 GG508 ///
GG514 GG515 GG516 GG517 GG518 GG519 GG520 GG521 ///
GG600 GG601 GG602 GG603 GG604 GG605 ///
GG606 GG607 GG608 GG609 GG610 GG611 ///
GG612 GG613 GG614 GG615)

*presence flag:
gen present_Q6=1 if _merge!=1
fre present_Q6

*About the child

*Q3 long-term illness
label define Q6_longtermprob 0"no" 1"yes, now" 2"yes, previously"
*3. Delayed motor development (e.g. sits/walks late) 	
fre GG37 	GG38 	GG39 	GG40
*No
*yes now
*yes previously
*combine:
gen delayed_motor_Q6=.
replace delayed_motor_Q6=0 if GG37==1
replace delayed_motor_Q6=1 if GG38==1
replace delayed_motor_Q6=2 if GG39==1
label values delayed_motor_Q6 Q6_longtermprob 
fre delayed_motor_Q6
tab1 GG37 GG38 GG39 if delayed_motor_Q6==. & _merge!=1
*nothing on those ones

*19. Late or abnormal speech development 	
fre GG93 GG94 GG95 GG96 
gen delayed_speech_Q6=.
replace delayed_speech_Q6=0 if GG93==1
replace delayed_speech_Q6=1 if GG94==1
replace delayed_speech_Q6=2 if GG95==1
label values delayed_speech_Q6 Q6_longtermprob 
fre delayed_speech_Q6
tab1 GG93 GG94 GG95 if delayed_speech_Q6==. & _merge!=1

*Q22 & Q23: ESAT and MoBA-specific ASD screen questions:
*fix multiple ticks
foreach var in GG249 GG250 GG592 GG252 GG253 GG254 {
	recode `var' 0=.
	fre `var'
}
*cell sizes might be a bit too small to be useful

*Trouble relating to others
fre GG578 GG579 GG580 GG581 
gen trouble_relating_Q6=.
replace trouble_relating_Q6=0 if GG578==1
replace trouble_relating_Q6=1 if GG579==1
replace trouble_relating_Q6=2 if GG580==1
label values trouble_relating_Q6 Q6_longtermprob 
fre trouble_relating_Q6
tab1 GG578 GG579 GG580 if trouble_relating_Q6==. & _merge!=1

*hyperactivity
fre GG105 GG106 GG107 GG108 
gen hyperactivity_Q6=.
replace hyperactivity_Q6=0 if GG105==1
replace hyperactivity_Q6=1 if GG106==1
replace hyperactivity_Q6=2 if GG107==1
label values hyperactivity_Q6 Q6_longtermprob 
fre hyperactivity_Q6
tab1 GG105 GG106 GG107 if hyperactivity_Q6==. & _merge!=1

********************************************
*autistic traits: VERSION DIFFERENCES
fre GG101 GG582 GG102 GG583 GG103 GG584 GG104 GG585 
tab1 GG101 GG582
tab1 GG102 GG583
tab1 GG103 GG584 
********************************************

gen autistic_traits_Q6=.
replace autistic_traits_Q6=0 if GG101==1 | GG582==1
replace autistic_traits_Q6=1 if GG102==1 | GG583==1
replace autistic_traits_Q6=2 if GG103==1 | GG584==1
label values autistic_traits_Q6 Q6_longtermprob 
fre autistic_traits_Q6
tab1 GG101 GG582 GG102 GG583 GG103 GG584 GG104 GG585  if autistic_traits_Q6==. & _merge!=1

*other behavioural problems
fre GG109 GG110 GG111 GG112 
gen other_behavioural_Q6=.
replace other_behavioural_Q6=0 if GG109==1
replace other_behavioural_Q6=1 if GG110==1
replace other_behavioural_Q6=2 if GG111==1
label values other_behavioural_Q6 Q6_longtermprob 
fre other_behavioural_Q6
tab1 GG109 GG110 GG111 if other_behavioural_Q6==. & _merge!=1

*17 & 21: ages and stages questionnaire:
fre GG222 GG223 GG224 GG225 GG237 GG238 GG239 GG240 GG241 GG242
*fix multiple ticks
foreach var of varlist GG222 GG223 GG224 GG225 GG237 GG238 GG239 GG240 GG241 GG242  {
recode `var' 0=.
fre `var'
}
*19 non verbal communication checklist
fre GG227 GG228 GG229 GG230
*fix multiple ticks
foreach var of varlist GG227 GG228 GG229 GG230  {
recode `var' 0=.
fre `var'
}
*20: strengths and difficulties questionnaire
fre GG231 GG232 GG233 GG234 GG235 GG236
*fix multiple ticks
foreach var of varlist GG231 GG232 GG233 GG234 GG235 GG236  {
recode `var' 0=.
fre `var'
}
*22. Autistic Traits Part I: Modified Checklist for Autism in Toddlers (M-CHAT)  
fre GG243 GG244 GG245 GG246 GG247 GG248 GG251
*fix multiple ticks
foreach var of varlist GG243 GG244 GG245 GG246 GG247 GG248 GG251  {
recode `var' 0=.
fre `var'
}
*fix multiple ticks
foreach var of varlist GG243 GG244 GG245 GG246 GG247 GG248 GG251  {
recode `var' 0=.
fre `var'
}

************************************************************

*23-25. Social Communication Questionnaire (SCQ) 36 months
foreach var of varlist GG255 GG257 GG258 GG259 GG260 GG261 GG262 GG263 GG264 GG265 GG266 GG267 GG268 GG269 GG270 GG271 GG272 GG273 GG274 GG256 GG275 GG276 GG277 GG278 GG279 GG280 GG281 GG282 GG283 GG284 GG285 GG286 GG287 GG288 GG289 GG290 GG291 GG292 GG293 GG294  {
recode `var' 0=.
fre `var'
}
*same as at Q8yr -items are identical

*reverse code the SCQ variables NN152,NN153,NN154,NN155,NN156,NN157,NN159,NN160,NN161,NN162,NN163,NN164,NN165,NN166,NN167
foreach var in GG258 GG259 GG260 GG261 GG262 GG263 GG265 GG266 GG267 GG268 GG269 GG270 GG271 GG272 GG273 {
fre `var'
gen r`var'=`var'
recode r`var' 2=0
replace r`var'=r`var'+1
}
*inspect
tab GG258 rGG258, nol
*ok.

*missigness counter and summary scores:
*full scale:
gen nitems_SCQ_Q6=0
foreach var of varlist GG255 GG257 rGG258 rGG259 rGG260 rGG261 rGG262 rGG263 GG264 rGG265 rGG266 rGG267 rGG268 rGG269 rGG270 rGG271 rGG272 rGG273 GG274 GG256 GG275 GG276 GG277 GG278 GG279 GG280 GG281 GG282 GG283 GG284 GG285 GG286 GG287 GG288 GG289 GG290 GG291 GG292 GG293 GG294 {
replace nitems_SCQ_Q6=nitems_SCQ_Q6+1 if `var'!=.
}
fre nitems_SCQ_Q6 if present_Q6==1
*social
gen nitems_S_SCQ_Q6=0
foreach var of varlist GG255 GG257 rGG259 GG264 rGG265 GG274 GG256 GG275 GG276 GG277 GG278 GG279 GG280 GG281 GG282 GG283 GG284 GG285 GG286 GG287 GG288 GG290 GG291 GG292 GG293 GG294 {
replace nitems_S_SCQ_Q6=nitems_S_SCQ_Q6+1 if `var'!=.
}
fre nitems_S_SCQ_Q6 if present_Q6==1
*repetitive
gen nitems_R_SCQ_Q6=0
foreach var of varlist rGG258 rGG260 rGG261 rGG262 rGG263 rGG266 rGG267 rGG268 rGG269 rGG270 rGG271 rGG273 {
replace nitems_R_SCQ_Q6=nitems_R_SCQ_Q6+1 if `var'!=.
}
fre nitems_R_SCQ_Q6 if present_Q6==1

*not included in either:
tab1 rGG272 GG289 

*summary scales:
*full scale:
gen cc_SCQ_Q6 =(GG255 + GG257 + rGG258 + rGG259 + rGG260 + rGG261 + rGG262 + rGG263 + GG264 + rGG265 + rGG266 + rGG267 + rGG268 + rGG269 + rGG270 + rGG271 + rGG272 + rGG273 + GG274 + GG256 + GG275 + GG276 + GG277 + GG278 + GG279 + GG280 + GG281 + GG282 + GG283 + GG284 + GG285 + GG286 + GG287 + GG288 + GG289 + GG290 + GG291 + GG292 + GG293 + GG294)
*bump down to start at 0
replace cc_SCQ_Q6=cc_SCQ_Q6-40 
summ cc_SCQ_Q6
*social:
gen cc_S_SCQ_Q6 =(GG255 + GG257 + rGG259 + GG264 + rGG265 + GG274 + GG256 + GG275 + GG276 + GG277 + GG278 + GG279 + GG280 + GG281 + GG282 + GG283 + GG284 + GG285 + GG286 + GG287 + GG288 + GG290 + GG291 + GG292 + GG293 + GG294)
*bump down to start at 0
replace cc_S_SCQ_Q6=cc_S_SCQ_Q6-26
summ cc_S_SCQ_Q6
*repetitive
gen cc_R_SCQ_Q6 =(rGG258 + rGG260 + rGG261 + rGG262 + rGG263 + rGG266 + rGG267 + rGG268 + rGG269 + rGG270 + rGG271 + rGG273)
*bump down to start at 0
replace cc_R_SCQ_Q6=cc_R_SCQ_Q6-12
summ cc_R_SCQ_Q6


*full scale:
egen miss20pc_SCQ_Q6 =rowtotal (GG255 GG257 rGG258 rGG259 rGG260 rGG261 rGG262 rGG263 GG264 rGG265 rGG266 rGG267 rGG268 rGG269 rGG270 rGG271 rGG272 rGG273 GG274 GG256 GG275 GG276 GG277 GG278 GG279 GG280 GG281 GG282 GG283 GG284 GG285 GG286 GG287 GG288 GG289 GG290 GG291 GG292 GG293 GG294)
replace miss20pc_SCQ_Q6=miss20pc_SCQ_Q6*(40/39) if nitems_SCQ_Q6==39
replace miss20pc_SCQ_Q6=miss20pc_SCQ_Q6*(40/38) if nitems_SCQ_Q6==38
replace miss20pc_SCQ_Q6=miss20pc_SCQ_Q6*(40/37) if nitems_SCQ_Q6==37
replace miss20pc_SCQ_Q6=miss20pc_SCQ_Q6*(40/36) if nitems_SCQ_Q6==36
replace miss20pc_SCQ_Q6=miss20pc_SCQ_Q6*(40/35) if nitems_SCQ_Q6==35
replace miss20pc_SCQ_Q6=miss20pc_SCQ_Q6*(40/34) if nitems_SCQ_Q6==34
replace miss20pc_SCQ_Q6=miss20pc_SCQ_Q6*(40/33) if nitems_SCQ_Q6==33
replace miss20pc_SCQ_Q6=miss20pc_SCQ_Q6*(40/32) if nitems_SCQ_Q6==32
replace miss20pc_SCQ_Q6=. if nitems_SCQ_Q6<32
*bump down to start at 0
replace miss20pc_SCQ_Q6=miss20pc_SCQ_Q6-40 
fre miss20pc_SCQ_Q6
tab miss20pc_SCQ_Q6 cc_SCQ_Q6, mi

*social:
egen miss20pc_S_SCQ_Q6 =rowtotal(GG255 GG257 rGG259 GG264 rGG265 GG274 GG256 GG275 GG276 GG277 GG278 GG279 GG280 GG281 GG282 GG283 GG284 GG285 GG286 GG287 GG288 GG290 GG291 GG292 GG293 GG294)
replace miss20pc_S_SCQ_Q6=miss20pc_S_SCQ_Q6*(26/25) if nitems_S_SCQ_Q6==25
replace miss20pc_S_SCQ_Q6=miss20pc_S_SCQ_Q6*(26/24) if nitems_S_SCQ_Q6==24
replace miss20pc_S_SCQ_Q6=miss20pc_S_SCQ_Q6*(26/23) if nitems_S_SCQ_Q6==23
replace miss20pc_S_SCQ_Q6=miss20pc_S_SCQ_Q6*(26/22) if nitems_S_SCQ_Q6==22
replace miss20pc_S_SCQ_Q6=miss20pc_S_SCQ_Q6*(26/21) if nitems_S_SCQ_Q6==21
replace miss20pc_S_SCQ_Q6=. if nitems_S_SCQ_Q6<21
fre miss20pc_S_SCQ_Q6
*bump down to start at 0
replace miss20pc_S_SCQ_Q6=miss20pc_S_SCQ_Q6-26
fre miss20pc_S_SCQ_Q6
tab miss20pc_S_SCQ_Q6 cc_S_SCQ_Q6, mi

*repetitive
egen miss20pc_R_SCQ_Q6 =rowtotal(rGG258 rGG260 rGG261 rGG262 rGG263 rGG266 rGG267 rGG268 rGG269 rGG270 rGG271 rGG273)
replace miss20pc_R_SCQ_Q6=miss20pc_R_SCQ_Q6*(12/11) if nitems_R_SCQ_Q6==11
replace miss20pc_R_SCQ_Q6=miss20pc_R_SCQ_Q6*(12/10) if nitems_R_SCQ_Q6==10
replace miss20pc_R_SCQ_Q6=. if nitems_R_SCQ_Q6<10
fre miss20pc_R_SCQ_Q6
*bump down to start at 0
replace miss20pc_R_SCQ_Q6=miss20pc_R_SCQ_Q6-12
fre miss20pc_R_SCQ_Q6
tab miss20pc_R_SCQ_Q6 cc_R_SCQ_Q6, mi

*rename these to be the default:
foreach varstem in SCQ_Q6 S_SCQ_Q6 R_SCQ_Q6 {
	rename miss20pc_`varstem' `varstem'
}


************************************************************

*26. loss of skills
fre GG295 GG296 GG297 GG298
*fix multiple ticks
foreach var of varlist GG295 GG296 GG297 GG298  {
recode `var' 0=.
fre `var'
}
*27. Temperament
fre GG299 GG300 GG301 GG302 GG303 GG304 GG305 GG306 GG307 GG308 GG309 GG310 GG311 GG312
*fix multiple ticks
foreach var of varlist GG299 GG300 GG301 GG302 GG303 GG304 GG305 GG306 GG307 GG308 GG309 GG310 GG311 GG312 {
recode `var' 0=.
fre `var'
}

*28. Child Behaviour Checklist (CBCL)
fre GG313 GG314 GG315 GG316 GG317 GG318 GG319 GG320 GG321 GG322 GG323 GG324 GG325 GG326 GG327 GG328 GG329 GG330 GG331 GG332 GG333 GG334 GG335 GG336 GG337 GG338
*fix multiple ticks
foreach var of varlist GG313 GG314 GG315 GG316 GG317 GG318 GG319 GG320 GG321 GG322 GG323 GG324 GG325 GG326 GG327 GG328 GG329 GG330 GG331 GG332 GG333 GG334 GG335 GG336 GG337 GG338 GG309 GG310 GG311 GG312 {
recode `var' 0=.
fre `var'
}

*36 months Child Behaviour Checklist (CBCL): make binary versions:
foreach var of varlist GG313-GG338 {
	gen `var'_bin=`var'
	recode `var'_bin 3=2
	tab `var'_bin
}

*29. Part I: Child Behavior and Manner
fre GG339 GG340 GG341 GG342 GG343 GG344 GG345 GG346 GG347 GG348
*fix multiple ticks
foreach var of varlist GG339 GG340 GG341 GG342 GG343 GG344 GG345 GG346 GG347 GG348 {
recode `var' 0=.
fre `var'
}

*29. Part II: The Infant-Toddler Social and Emotional Assessment (ITSEA)
fre GG349 GG350 GG351 GG352 GG353 GG354 GG355 GG356 GG357 GG358 GG359 GG360 GG361 GG362 GG363 GG364 GG365 GG366 GG367 GG236
*fix multiple ticks
foreach var of varlist GG349 GG350 GG351 GG352 GG353 GG354 GG355 GG356 GG357 GG358 GG359 GG360 GG361 GG362 GG363 GG364 GG365 GG366 GG367 GG236 {
recode `var' 0=.
fre `var'
}

*31. Maternal Concerns
*GG596 not there
fre GG380 GG381 GG382 GG594 GG595 
*fix multiple ticks
foreach var of varlist GG380 GG381 GG382 GG594 GG595  {
recode `var' 0=.
fre `var'
}

*34-35. Brushing Teeth
*these variable not there. weird.
*fre GG386 GG387

*36 Is your child ever present in a room where someone smokes?
fre GG388 GG389
recode GG388 0=.
*an 0 in the other one is legit

*******************************

*ABOUT THE MOTHER

*NB: GG448- indicator for whether mother is currently pregnant - NOT IN THE ORIGINAL FILE. DO NOT KNOW WHY.

*52. Life Time History of Major Depression (LTH of MD)
*NB: DIFFERENT QUESTIONS IN DIFFERENT VER, COULD BE DIFFICULT
fre GG452 GG634 GG635 GG636 GG637 GG453 GG638 GG639 GG640 GG641 GG454 GG642 GG643 GG644 GG645 GG455 GG646 GG647 GG648 GG649 GG456 GG650 GG651 GG652 GG653 GG457 GG654 GG655 GG656 GG657
**fix multiple ticks
foreach var of varlist GG452 GG634 GG635 GG636 GG637 GG453 GG638 GG639 GG640 GG641 GG454 GG642 GG643 GG644 GG645 GG455 GG646 GG647 GG648 GG649 GG456 GG650 GG651 GG652 GG653 GG457 GG654 GG655 GG656 GG657  {
recode `var' 0=.
fre `var'
}
*53-58. Health and Health Problems
*Are you pregnant now: NOT HERE! WHY???
*fre GG448 
*these aren't here GG459 GG462 GG460 GG461 GG463
*Have you had any long-term illness or health problems that have occurred during the last 3 years? NB: different items between versions a b c 
fre GG462 GG658 GG659 GG660 GG661 GG662 GG663 
**fix multiple ticks
foreach var of varlist GG462 GG658 GG659 GG660 GG661 GG662 GG663   {
recode `var' 0=.
fre `var'
}

*59-62. Smoking and alcohol
*smoking coded below

*64-67 Eating disorders
fre GG491 GG492 GG493 GG494  GG495 GG496 GG497 GG498 GG499  GG500 GG501 GG502
**fix multiple ticks, except for the last 2
foreach var of varlist GG491 GG492 GG493 GG494  GG495 GG496 GG497 GG498 GG499  GG500   {
recode `var' 0=.
*other multiple categories
recode `var' 4/6=.
fre `var'
}

*68. Adult ADHD
fre GG503 GG504 GG505 GG506 GG507 GG508
**fix multiple ticks
foreach var of varlist GG503 GG504 GG505 GG506 GG507 GG508   {
recode `var' 0=.
fre `var'
}
*70. Hopkins checklist
fre GG514 GG515 GG516 GG517 GG518 GG519 GG520 GG521
**fix multiple ticks
foreach var of varlist GG514 GG515 GG516 GG517 GG518 GG519 GG520 GG521  {
recode `var' 0=.
fre `var'
}
*72. Differential Emotional Scale (DES), Enjoyment and Anger Subscales
fre GG600 GG601 GG602 GG603 GG604 GG605
**fix multiple ticks
foreach var of varlist GG600 GG601 GG602 GG603 GG604 GG605 {
recode `var' 0=.
fre `var'
}
*73. Satisfaction with Life Scale (SWLS)
fre GG606 GG607 GG608 GG609 GG610 GG611
*fix multiple ticks
foreach var of varlist GG606 GG607 GG608 GG609 GG610 GG611 {
recode `var' 0=.
fre `var'
}
*74. Rosenberg Self Esteem Scale (RSES)
*fix multiple ticks
fre GG612 GG613 GG614 GG615
foreach var of varlist GG612 GG613 GG614 GG615 {
recode `var' 0=.
fre `var'
}

*BMI
foreach var in GG501 GG502 {
summ `var', det
return list
replace `var'=. if `var'<r(mean)-(4*r(sd))
replace `var'=. if `var'>r(mean)+(4*r(sd))
summ `var', det
}
*now make BMI:
gen mothers_bmi_Q6=GG501/((GG502/100)^2)
label variable mothers_bmi_Q6 "Q1 mother's bmi a child age 3 from self-report h & w"
summ mothers_bmi_Q6, det

*mother's martial status:
fre GG435
recode GG435 0=.

*also for imputation: q 58 through 62
fre GG479-GG488
*make a new smoking variable
gen mumsmokes_Q6=GG479
recode mumsmokes_Q6 4=.
*for the 20 who double-ticked sometimes and daily, assume sometimes?
*yes, but then need to go back and make consistent
recode mumsmokes_Q6 6=3
*bump down by 1
replace mumsmokes_Q6=mumsmokes_Q6-1
label define mumsmokes_Q6 0"no" 1"sometimes" 2"daily"
label values mumsmokes_Q6 mumsmokes_Q6
fre mumsmokes_Q6

*add to imputation: GG435 mumsmokes_Q6
*like at other Qs, stuff on employment, but leave for now

drop _merge

*****************************************************************************************************************

*Age 5:
*child height and weight, mother's weight and if pregnant

*IMPORTANT: THERE'S LOADS OF STUFF HERE ON CHILD DEVELOPMENT, FROM BASIC STANDALONE Qs LIKE 
*11. Does your child receive, or has received any extra resources in the kindergarten?
*TO MUCH MORE IN-DEPTH STUFF. CLEARLY RELEVANT TO LATER NEURODEVELOPMENTAL OUTCOMES, SO COME BACK AND GO THROUGH ALL THIS TO ADD TO IMP MODELS.

merge 1:1 PREG_ID_2306 BARN_NR using "N:\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\Statafiles\PDB2306_Q5yrs_v12.dta", ///
keepus (PREG_ID_2306 BARN_NR ///
LL12 LL13 LL338 LL339 LL508 AGE_SENT_MTHS_Q5AAR AGE_MTHS_Q5AAR AGE_RETURN_MTHS_Q5AAR ///
LL508 LL340- LL515 ///
LL174 LL175 LL176 LL177 LL178 LL179 LL180 ///
LL264 LL265 LL266 LL267 LL268 LL269 LL270 LL271 LL272 LL273 LL274 LL275 ///
LL190 LL191 LL192 LL193 LL194 LL195 LL196 LL197 LL198 LL199 LL200 LL201 LL202 LL203 LL204 LL205 LL206 LL207 LL208 LL209 LL210 LL211 LL212 ///
LL213 LL214 LL215 LL216 LL217 LL218 LL219 LL220 LL221 LL480 LL481 LL482 LL483 ///
LL222 LL223 LL224 LL225 LL226 LL227 LL228 LL229 LL230 LL231 ///
LL232 LL233 LL234 LL235 LL236 LL237 LL238 LL239 LL240 LL241 LL242 LL243 LL244 LL245 LL246 LL247 LL248 LL249 LL250 LL251 ///
LL252 LL253 LL254 LL255 LL256 LL257 LL258 LL259 LL260 LL261 LL262 LL263 ///
LL276 LL277 LL278 LL279 LL280 LL281 LL282 LL283 LL284 LL285 LL286 LL287 ///
LL288 LL289 LL290 LL292 LL293 LL294 LL295 LL296 LL297 LL298 LL299 LL300 ///
LL301 LL302 LL303 LL304 LL305 LL306 LL307 LL308 LL309 LL310 LL311 LL312 LL313 LL314 LL315 LL316 LL317 LL318 LL319 LL320 LL321 LL322 LL323 LL324 LL325 LL504 LL505 ///
LL326 LL327 LL328 LL329 LL330 LL331 LL332 LL333 LL334 LL335 LL336 LL337 ///
LL361 LL362 LL363 ///
LL363 LL364 LL365 LL366 LL367 LL368 LL369 LL370 ///
LL382 LL383 LL384 LL385 LL386 LL387)

*presence flag:
gen present_Q5y=1 if _merge!=1
fre present_Q5y 

*28 & 36. Ages and Stages Questionnaires (ASQ)
fre LL174 LL175 LL176 LL177 LL178 LL179 LL180 
fre LL264 LL265 LL266 LL267 LL268 LL269 LL270 LL271 LL272 LL273 LL274 LL275
**fix multiple ticks
foreach var of varlist LL174 LL175 LL176 LL177 LL178 LL179 LL180 ///
LL264 LL265 LL266 LL267 LL268 LL269 LL270 LL271 LL272 LL273 LL274 LL275 {
recode `var' 0=.
fre `var'
}
**for imputation, MAKE BINARY VERSIONS FOR THE 3-CATEG ONES - the non-no categories too small to use
label define ASQ_binarized 1"yes" 2"sometimes or not yet"
foreach var of varlist LL174 LL175 LL176 LL177 LL178 LL179 LL180  {
	gen `var'_binary=`var'
	recode `var'_binary 3=2
	label variable `var'_binary "ever had this condition, no/yes"
	label values `var'_binary ASQ_binarized
	fre `var'_binary
}

*31. Checklist of 20 Statements about Language-Related Difficulties (Språk20)
fre LL190 LL191 LL192 LL193 LL194 LL195 LL196 LL197 LL198 LL199 LL200 LL201 LL202 LL203 LL204 LL205 LL206 LL207 LL208 LL209 LL210 LL211 LL212
**fix multiple ticks
foreach var of varlist LL190 LL191 LL192 LL193 LL194 LL195 LL196 LL197 LL198 LL199 LL200 LL201 LL202 LL203 LL204 LL205 LL206 LL207 LL208 LL209 LL210 LL211 LL212 {
recode `var' 0=.
fre `var'
}

*32. Children’s Communication Checklist-2 Coherence Sub-scale (CCC-2 Coherence)
fre LL213 LL214 LL215 LL216 LL217 LL218 LL219 LL220 LL221 LL480 LL481 LL482 LL483
**fix multiple ticks
foreach var of varlist LL213 LL214 LL215 LL216 LL217 LL218 LL219 LL220 LL221 LL480 LL481 LL482 LL483 {
recode `var' 0=.
fre `var'
}
*some of these certainly measure the child’s environment rather than ability, but relevant for imputation:
*33. Preschool Activities: Narrative and Communicative SkiLLs
fre  LL222 LL223

*33. Preschool Activities: Experiences with Letter and Sound Knowledge
fre  LL224 LL225

*Preschool Activities: Literacy SkiLLs
fre LL226 LL227 LL228 LL229 LL230

*33. Preschool Activities: Home Reading
fre LL231 

**fix multiple ticks
foreach var of varlist LL222 LL223 LL224 LL225 LL226 LL227 LL228 LL229 LL230 LL231 {
recode `var' 0=.
fre `var'
}

*35a. Childhood Asperger Syndrome Test (CAST)
fre LL232 LL233 LL234 LL235 LL236 LL237 LL238 LL239 LL240 LL241 LL242 LL243 LL244 LL245 LL246 LL247 LL248 LL249 LL250 LL251
**fix multiple ticks
foreach var of varlist LL232 LL233 LL234 LL235 LL236 LL237 LL238 LL239 LL240 LL241 LL242 LL243 LL244 LL245 LL246 LL247 LL248 LL249 LL250 LL251 {
recode `var' 0=.
fre `var'
}

*35. Conners Parent Rating Scale-Revised, Short Form (CPRS-R (S))
fre LL252 LL253 LL254 LL255 LL256 LL257 LL258 LL259 LL260 LL261 LL262 LL263
**fix multiple ticks
foreach var of varlist LL252 LL253 LL254 LL255 LL256 LL257 LL258 LL259 LL260 LL261 LL262 LL263 {
recode `var' 0=.
fre `var'
}
*missigness counter and summary scores:
*missingness
gen nitems_conners_Q5y=0
foreach var of varlist LL252 LL253 LL254 LL255 LL256 LL257 LL258 LL259 LL260 LL261 LL262 LL263 {
replace nitems_conners_Q5y=nitems_conners_Q5y+1 if `var'!=.
}
fre nitems_conners_Q5y if present_Q5y==1
*complete case:
gen cc_conners_Q5y=(LL252 + LL253 + LL254 + LL255 + LL256 + LL257 + LL258 + LL259 + LL260 + LL261 + LL262 + LL263)
*bump down cc scale to start at 0 
replace cc_conners_Q5y=cc_conners_Q5y-12
fre cc_conners_Q5y
*miss20pc version:
*12 items, so make a version allowing up to 2 missing:
egen miss20pc_conners_Q5y=rowtotal(LL252 LL253 LL254 LL255 LL256 LL257 LL258 LL259 LL260 LL261 LL262 LL263)
replace miss20pc_conners_Q5y=miss20pc_conners_Q5y*(12/11) if nitems_conners_Q5y==11
replace miss20pc_conners_Q5y=miss20pc_conners_Q5y*(12/10) if nitems_conners_Q5y==10
replace miss20pc_conners_Q5y=. if nitems_conners_Q5y<10
*bump down cc scale to start at 0 
replace miss20pc_conners_Q5y=miss20pc_conners_Q5y-12
*check these:
tab miss20pc_conners_Q5y cc_conners_Q5y, mi

*rename this to be the default:
foreach varstem in conners_Q5y {
	rename miss20pc_`varstem' `varstem'
}

*37. Selective items from the Emotionality, Activity and Shyness Temperament Questionnaire (EAS)
fre LL276 LL277 LL278 LL279 LL280 LL281 LL282 LL283 LL284 LL285 LL286 LL287
**fix multiple ticks
foreach var of varlist LL276 LL277 LL278 LL279 LL280 LL281 LL282 LL283 LL284 LL285 LL286 LL287 {
recode `var' 0=.
fre `var'
}
*38. Speech and Language Assessment Scale (SLAS)
fre LL288 LL289 LL290 LL292 LL293 LL294 LL295 LL296 LL297 LL298 LL299 LL300
**fix multiple ticks
foreach var of varlist LL288 LL289 LL290 LL292 LL293 LL294 LL295 LL296 LL297 LL298 LL299 LL300 {
recode `var' 0=.
fre `var'
}
*39. Child Behaviour Checklist (CBCL)
fre LL301 LL302 LL303 LL304 LL305 LL306 LL307 LL308 LL309 LL310 LL311 LL312 LL313 LL314 LL315 LL316 LL317 LL318 LL319 LL320 LL321 LL322 LL323 LL324 LL325 LL504 LL505
**fix multiple ticks
foreach var of varlist LL301 LL302 LL303 LL304 LL305 LL306 LL307 LL308 LL309 LL310 LL311 LL312 LL313 LL314 LL315 LL316 LL317 LL318 LL319 LL320 LL321 LL322 LL323 LL324 LL325 LL504 LL505 {
recode `var' 0=.
fre `var'
}

*Q5yr Child Behaviour Checklist (CBCL): make binary versions
foreach var of varlist LL301 LL302 LL303 LL304 LL305 LL306 LL307 LL308 LL309 LL310 LL311 LL312 LL313 LL314 LL315 LL316 LL317 LL318 LL319 LL320 LL321 LL322 LL323 LL324 LL325 LL504 LL505 {
	gen `var'_bin=`var'
	recode `var'_bin 3=2
	tab `var'_bin
}
*Don't use the last two as version-specific - only there for a subset

*42. Maternal Concerns
fre LL326 LL327 LL328 LL329 LL330 LL331 LL332 LL333 LL334 LL335 LL336 LL337
**fix multiple ticks
*the last is weeks when started so 0 is legit
foreach var of varlist LL326 LL327 LL328 LL329 LL330 LL331 LL332 LL333 LL334 LL335 LL336 {
recode `var' 0=.
fre `var'
}

**********************************
*About the mother

*48-50.  Questions about the mother’s health and health problems
*Most of this section looks fiddly as hell - text entries and differences between versions
*50 Have you ever had any problems with your physical or mental health that has prevented you in your work or social activities with family or friends?  LL361- LL363, last one for degree of MH impairment
fre LL361 LL362 LL363

*51 (Hopkins) Symptoms Checklist-25 (SCL-25)
fre LL363 LL364 LL365 LL366 LL367 LL368 LL369 LL370

*54. Satisfaction with Life Scale (SWLS)
fre LL382 LL383 LL384 LL385 LL386 LL387

**fix multiple ticks
foreach var of varlist LL361 LL362 LL363 ///
LL363 LL364 LL365 LL366 LL367 LL368 LL369 LL370 ///
LL382 LL383 LL384 LL385 LL386 LL387 {
recode `var' 0=.
fre `var'
}

*height and weight
foreach var in LL12 LL13 LL338 LL339  {
summ `var', det
return list
replace `var'=. if `var'<r(mean)-(4*r(sd))
replace `var'=. if `var'>r(mean)+(4*r(sd))
summ `var', det
}
*now make BMI:
*#CAREFUL: only in one of the versions were mothers asked for current height.
*#For the others, use an earlier report: GG502 and AA87.
gen mothers_bmi_Q5y=LL339/((LL338/100)^2)
replace mothers_bmi_Q5y=LL339/((GG502/100)^2) if LL338==. 
replace mothers_bmi_Q5y=LL339/((AA87/100)^2) if LL338==. & GG502==. 
label variable mothers_bmi_Q5y "Q5y mother's bmi at child age 5 from self-report h & w"
summ mothers_bmi_Q5y, det
*not bad subbing, but some strange values for weight in kg at Q_5yr, especially when in context with the same people's weight from Q1:
list mothers_bmi_Q5y LL339 LL338 mothers_bmi_Q1 AA85 AA87 GG502  if mothers_bmi_Q5y<15 & mothers_bmi_Q5y>1
list mothers_bmi_Q5y LL339 LL338 mothers_bmi_Q1 AA85 AA87 GG502  if mothers_bmi_Q5y>22 & mothers_bmi_Q5y<22.5
*don't use this one.

gen moffspring_bmi_Q5y=LL13/((LL12/100)^2)
label variable moffspring_bmi_Q5y "Q5y offspring bmi from mother's self-report h & w"
summ moffspring_bmi_Q5y, det
*looks acceptable

capture drop mumsmokes_Q5y
gen mumsmokes_Q5y=LL509
recode mumsmokes_Q5y 4/7=.
*bump down
replace mumsmokes_Q5y=mumsmokes_Q5y-1
label var mumsmokes_Q5y "mother's smoking status a Q5y"
label define mumsmokes_Q5y 0"never" 1"sometimes" 2"daily"
label values mumsmokes_Q5y mumsmokes_Q5y
tab mumsmokes_Q5y LL509

*and for dads?
capture drop dadsmokes_Q5y
gen dadsmokes_Q5y=LL511
recode dadsmokes_Q5y 4/7=.
*bump down
replace dadsmokes_Q5y=dadsmokes_Q5y-1
label var dadsmokes_Q5y "father's smoking status a Q5y"
label define dadsmokes_Q5y 0"never" 1"sometimes" 2"daily"
label values dadsmokes_Q5y dadsmokes_Q5y
tab dadsmokes_Q5y LL511

*add to imputation: mumsmokes_Q5y  dadsmokes_Q5y plus all the behavioural/developmental stuff!!!

drop _merge

**************************************************************************************************************

*Age 7:
merge 1:1 PREG_ID_2306 BARN_NR using "N:\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\Statafiles\PDB2306_Q7yrs_v12.dta", ///
keepus (JJ408 JJ324 JJ325 JJ329 JJ37 JJ38 JJ39 JJ40 JJ41 JJ42 JJ43 JJ44 JJ45 JJ46 AGE_SENT_MTHS_Q7 AGE_MTHS_Q7 AGE_RETURN_MTHS_Q7 ///
JJ300- JJ305 ///
/*vars about child's specific health problems*/ JJ431 JJ435 JJ436 JJ437 JJ438 JJ441 JJ443 JJ445 JJ446 JJ447 JJ448 JJ449 JJ618 JJ619 JJ620 JJ621 JJ622)

gen present_Q7y=1 if _merge!=1
fre present_Q7y

*Q25: child's health conditions
fre JJ431 JJ435 JJ436 JJ437 JJ441 
*NB: for chronic fatigue literally 1 person has a yes, so don't bother:
fre JJ438

*tick vs no tick.
*identify people for whom it's definitely a genuine missing by using the _merge variable: they are unmatched from master, i.e. _merge==1. DO NOT CODE .=2 FOR THESE PEOPLE!
*Most will be no's, so assume that. Label as 2.
*label define AA806 1"yes" 2"not reported"
foreach var of varlist JJ431 JJ435 JJ436 JJ437 JJ441  {
recode `var' .=2 if _merge!=1
label values `var' AA806
fre `var'
}

*merge info in JJ436 and JJ437
gen JJ436_JJ437=JJ436 
replace JJ436_JJ437=JJ437 if JJ437==1
label values JJ436_JJ437 AA806
label variable JJ436_JJ437 "combined report of autistic features OR aspergers"
fre JJ436_JJ437 

**********************************
*Q26
*ones from q26 will be tricky - different items and structure entirely in version A so will need to find equivalent ones and match up
fre JJ443 JJ445 JJ446 JJ447 JJ448 JJ449 JJ618 JJ619 JJ620 JJ621 JJ622
**********************************

*height and weight
*CAREFUL! 
*THERE ARE TWO HEIGHT VARS THAT NEED TO HAVE INFO COMBINED, OTHERWISE YOU GET A F-TONNE OF MISSINGNESS
*for some reason, there's a version of height in metres and a version in cm, but people only have one of either
list  JJ324 JJ408 JJ325  in 1/50
*sub in the cm reports to the m variable:
summ JJ324
replace JJ324=JJ408/100 if JJ324==.
summ JJ324
*clearly some dodgy values, but that's why we trim:
foreach var in JJ324 JJ325  {
summ `var', det
return list
replace `var'=. if `var'<r(mean)-(4*r(sd))
replace `var'=. if `var'>r(mean)+(4*r(sd))
summ `var', det
}

*now make BMI: NB: using height in metres not cm!
*looks acceptable
gen moffspring_bmi_Q7y=JJ325/((JJ324)^2)
label variable moffspring_bmi_Q7y "Q5y offspring bmi from mother's self-report h & w"
summ moffspring_bmi_Q7y, det
*looks acceptable

*smoking
capture drop mumsmokes_Q7y
gen mumsmokes_Q7y=JJ300
recode mumsmokes_Q7y 4/7=.
*bump down
replace mumsmokes_Q7y=mumsmokes_Q7y-1
label var mumsmokes_Q7y "mother's smoking status at Q7y"
label define mumsmokes_Q7y 0"never" 1"sometimes" 2"daily"
label values mumsmokes_Q7y mumsmokes_Q7y
tab mumsmokes_Q7y JJ300

capture drop dadsmokes_Q7y
gen dadsmokes_Q7y=JJ303
recode dadsmokes_Q7y 4/7=.
*bump down
replace dadsmokes_Q7y=dadsmokes_Q7y-1
label var dadsmokes_Q7y "father's smoking status at Q7y"
label define dadsmokes_Q7y 0"never" 1"sometimes" 2"daily"
label values dadsmokes_Q7y dadsmokes_Q7y
tab dadsmokes_Q7y JJ303

*NB no alcohol or parental MH stuff this year.

**********************************************************************************************
*father's height and weight from second father questionnaire
*merge 1:1 PREG_ID_2306 using PDB2306_Far2_V12.dta
*different id var here - is it linkable? looks like it might be the same but with aprefix to strip off, eg 2 corresponds to F000002
*use PDB2306_Far2_V12.dta
*list F_ID_2306 in 1/50


*merge 1:1 PREG_ID_2306 using PDB2306_Far2_V12.dta
*different id var here - is it linkable? looks like it might be the same but with aprefix to strip off, eg 2 corresponds to F000002
*use PDB2306_Far2_V12.dta
*list F_ID_2306 in 1/50

drop _merge

**********************************************************************************************
*Age 8 questionnaire for child phenotypes

merge 1:1 PREG_ID_2306 BARN_NR using "N:\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\Statafiles\PDB2306_Q8yrs_v12.dta"

gen present_Q8y=1 if _merge!=1
fre present_Q8y

*Q9 Has your child ever had any of the following health problems? 
*no / yes currently / yes in past + seen as specialist (y/n)
*for all of these, first three items are tick vs no tick so missing includes nos.
*code as you did for equivalent questions at Q6, and use the same label
*label define Q6_longtermprob 0"no" 1"yes, now" 2"yes, previously"

*1. Delayed psychomotor development 	
fre NN36 	NN37 	NN38 	NN39 
*No
*yes now
*yes previously
*combine:
gen delayed_motor_Q8yr=.
replace delayed_motor_Q8yr=0 if NN36==1
replace delayed_motor_Q8yr=1 if NN37==1
replace delayed_motor_Q8yr=2 if NN38==1
label values delayed_motor_Q8y Q6_longtermprob 
fre delayed_motor_Q8yr
tab1 NN36 NN37 NN38 if delayed_motor_Q8yr==. & _merge!=1
*nothing on those ones

*2. Delayed or abnormal language development 	
fre NN40 	NN41 	NN42 	NN43
gen delayed_language_Q8y=.
replace delayed_language_Q8y=0 if NN40==1
replace delayed_language_Q8y=1 if NN41==1
replace delayed_language_Q8y=2 if NN42==1
label values delayed_language_Q8y Q6_longtermprob 
fre delayed_language_Q8y
tab1 NN40 NN41 NN42 if delayed_language_Q8y==. & _merge!=1
 
*3. Hyperactivity 	
fre NN44 	NN45 	NN46 	NN47 
gen hyperactivity_Q8y=.
replace hyperactivity_Q8y=0 if NN44==1
replace hyperactivity_Q8y=1 if NN45==1
replace hyperactivity_Q8y=2 if NN46==1
label values hyperactivity_Q8y Q6_longtermprob 
fre hyperactivity_Q8y
tab1 NN44 NN45 NN46 if hyperactivity_Q8y==. & _merge!=1

*4. Concentration or attention difficulties 
fre	NN48 	NN49 	NN50 	NN51 
gen concentration_Q8y=.
replace concentration_Q8y=0 if NN48==1
replace concentration_Q8y=1 if NN49==1
replace concentration_Q8y=2 if NN50==1
label values concentration_Q8y Q6_longtermprob 
fre concentration_Q8y
tab1 NN48 NN49 NN50 if concentration_Q8y==. & _merge!=1

*5. Autistic traits /autism/Asperger’s Syndrome 
fre	NN52 NN53 NN54 NN55 
gen autistic_traits_Q8y=.
replace autistic_traits_Q8y=0 if NN52==1
replace autistic_traits_Q8y=1 if NN53==1
replace autistic_traits_Q8y=2 if NN54==1
label values autistic_traits_Q8y Q6_longtermprob 
fre autistic_traits_Q8y
tab1 NN52 NN53 NN54 if autistic_traits_Q8y==. & _merge!=1

*6. Behavioural problems (difficult and unruly) 
fre NN56 NN57 NN58 NN59 
gen behavioural_problems_Q8y=.
replace behavioural_problems_Q8y=0 if NN56==1
replace behavioural_problems_Q8y=1 if NN57==1
replace behavioural_problems_Q8y=2 if NN58==1
label values behavioural_problems_Q8y Q6_longtermprob 
fre behavioural_problems_Q8y
tab1 NN56 NN57 NN58 if behavioural_problems_Q8y==. & _merge!=1


*7. Emotional difficulties (sad or anxious) 	
fre NN60 NN61 NN62 NN63 
gen emotional_diffic_Q8y=.
replace emotional_diffic_Q8y=0 if NN60==1
replace emotional_diffic_Q8y=1 if NN61==1
replace emotional_diffic_Q8y=2 if NN62==1
label values emotional_diffic_Q8y Q6_longtermprob 
fre emotional_diffic_Q8y
tab1 NN60 NN61 NN62 if emotional_diffic_Q8y==. & _merge!=1

*8. Other 	
fre NN64 NN65 NN66 NN67 
gen other_condition_Q8y=.
replace other_condition_Q8y=0 if NN64==1
replace other_condition_Q8y=1 if NN65==1
replace other_condition_Q8y=2 if NN66==1
label values other_condition_Q8y Q6_longtermprob 
fre other_condition_Q8y
tab1 NN60 NN61 NN62 if other_condition_Q8y==. & _merge!=1

*10. Short Mood and Feelings Questionnaire (SMFQ) 
*(already have this prepared)

*11. Short Norwegian Hierarchical Personality Inventory for Children (NHiPIC-30)
*Note: here also there are changes between versions a/b/c
*versions B & C:
fre NN81 NN82 NN83 NN84 NN85 NN86 NN87 NN88 NN89 NN90 NN91 NN92 NN93 NN94 NN95 NN96 NN97 NN98 NN99 NN100 NN101 NN102 NN103 NN104 NN105 NN106 NN107 NN108 NN109 NN110
*In version A, the five items below differ from those in versions B 
*1.Become easily panic
fre NN368
*2. Will get to the bottom of things
fre NN369
*8. Have energy to spare
fre NN370
*10. Seeking contact with new classmates
fre NN371
*27. Feel at ease with him/herself
fre NN372

*fix multiple ticks:
foreach var of varlist NN81 NN82 NN83 NN84 NN85 NN86 NN87 NN88 NN89 NN90 NN91 NN92 NN93 NN94 NN95 NN96 NN97 NN98 NN99 NN100 NN101 NN102 NN103 NN104 NN105 NN106 NN107 NN108 NN109 NN110 NN368 NN369 NN370 NN371 NN372 {
recode `var' 0=.
fre `var'
}

*12-13. Parent/Teacher Rating Scale for Disruptive Behaviour Disorders (RS-DBD)
fre NN111 NN112 NN113 NN114 NN115 NN116 NN117 NN118 NN119 NN120 NN121 NN122 NN123 NN124 NN125 NN126 NN127 NN128 NN129 NN130 NN131 NN132 NN133 NN134 NN135 NN136 NN137 NN138 NN139 NN140 NN141 NN142 NN143 NN144
*fix multiple ticks:
foreach var of varlist NN111 NN112 NN113 NN114 NN115 NN116 NN117 NN118 NN119 NN120 NN121 NN122 NN123 NN124 NN125 NN126 NN127 NN128 NN129 NN130 NN131 NN132 NN133 NN134 NN135 NN136 NN137 NN138 NN139 NN140 NN141 NN142 NN143 NN144 {
recode `var' 0=.
fre `var'
}

*14. Screen for Child Anxiety Related Disorders (SCARED) 
*(already have this prepared)

*15-17. Social Communication Questionnaire (SCQ) 
*(already have this prepared)

*20. Children’s Communication Checklist-2 (CCC-2)
fre NN211-NN226
*fix multiple ticks:
foreach var of varlist NN211-NN226 {
recode `var' 0=.
fre `var'
}
*21. Checklist of 20 Statements about Language-Related Difficulties (Språk20)
*Note: VERSION DIFFERENCE: ITEM 8, NN374, ONLY IN VERSION C
fre NN227 NN228 NN229 NN230 NN231 NN232 NN233 NN374
*fix multiple ticks:
foreach var of varlist NN227 NN228 NN229 NN230 NN231 NN232 NN233 NN374 {
recode `var' 0=.
fre `var'
}

******************************************************************
*PLACEHOLDER EDUCATION OUTCOME:

*Q25 School assessment
*Reading skills in 1st grade 
*Reading skills in 2nd grade
*Arithmetic in 2nd grade
fre NN239 NN240 NN241

*fix multiple ticks
*also the 'Don’t know / have not talked to the teacher about it' needs to be coded to missing
foreach var of varlist NN239 NN240 NN241 {
recode `var' 0=.
*and the don't know category
recode `var' 4=.
fre `var'
}

**auxillary for imputation:
*26 Is an administrative decision made about your child being eligible for special education? 
*in various subjects:
fre NN242 NN244 NN246 NN248
*fix multiple ticks:
foreach var of varlist NN242 NN244 NN246 NN248 {
recode `var' 0=.
fre `var'
}
*how many hours: 0s are legit!
fre NN243 NN245 NN247 NN249 

********************************************************
*28. Reading and Writing Skills
*Potentially useful but differs quite a lot between versions

*MAJOR VERSION DIFFERENCES.
*28 Enter a cross indicating what your child masters:
*In version C, responses were 3-choice categorical:  Yes/ 2- Partially /3-Not yet
fre NN380 NN381 NN382 NN383 NN384 
*fix multiple ticks
foreach var of varlist NN380 NN381 NN382 NN383 NN384 {
recode `var' 0=.
fre `var'
}
*In version A & B, the questions were slightly different, but moreover binary response,
*but with a don't know option:
fre  NN254 NN255 NN256 NN257
*fix multiple ticks
foreach var of varlist NN254 NN255 NN256 NN257 {
recode `var' 0=.
*also the don't knows, here=3:
recode `var' 3=.
fre `var'
}
*Probably leave these out of imputation.

********************************************************

*Q29: Pronunciation, ability to tell a story,
*ability to communicate his/her own needs in a way understandable to adults and friends?
fre NN258 NN259 NN260 NN261 
*fix multiple ticks
foreach var of varlist NN258 NN259 NN260 NN261 {
recode `var' 0=.
*also the 2+3, 3+4 and 4+5 categories
recode `var' 23/45=.
fre `var'
}

*30-31: Home Reading and Self-reading
*VERSION DIFFERENCES: NN263/NN385 NN264/NN386
fre NN262 NN263 NN385 NN264 NN386 NN265
*fix multiple ticks
foreach var of varlist NN262 NN263 NN385 NN264 NN386 NN265 {
recode `var' 0=.
*also the 2+3, 3+4 etc categories
recode `var' 12	/234=.
fre `var'
}

*34. Difficulties, Impairment and Impact: from Strengths and Difficulties Questionnaire (SDQ)
fre NN388 NN389 NN390 NN391 NN392 NN393 NN394 NN395 NN396
*fix multiple ticks
foreach var of varlist NN388 NN389 NN390 NN391 NN392 NN393 NN394 NN395 NN396 {
recode `var' 0=.
fre `var'
}
**************************
*About the mother
*48 Do you have/ have you had any of the following disorders/illnesses? 
*1-No, never/2-Not now, but in the past/3-Yes, now 	
*Have you been treated for the problem/illness? 
*1-No/2-Yes 
*1. ADHD 
fre NN304 	NN305 
*2. Reading and writing difficulties 	
fre NN306 	NN307 
*3. Anorexia 	
fre NN308 	NN309 
*4. Bulimia 	
fre NN310 	NN311 

*Do you have or have you had any other serious illness or health problem? 	
fre NN312 
*If yes, what was the name of the illness (es)? 	NN313 (txt.) 
*not here

*fix multiple ticks
foreach var of varlist NN304 NN305 NN306 NN307 NN308 NN309 NN310 	NN311 NN312 {
recode `var' 0=.
fre `var'
}
*for purposes of imputation, just use the first of each pair:
fre NN304 NN306 NN308 NN310

*49. Social Phobia
fre NN314 NN315 NN316
*fix multiple ticks
foreach var of varlist NN314 NN315 NN316 {
recode `var' 0=.
fre `var'
}

*50. Satisfaction with Life Scale (SWLS)
fre NN317 NN318 NN319 NN320 NN321
*fix multiple ticks
foreach var of varlist NN317 NN318 NN319 NN320 NN321 {
recode `var' 0=.
fre `var'
}

*51. The Autonomic Nervous System Questionnaire (ANS)
fre NN322 NN323 NN324
*fix multiple ticks
foreach var of varlist NN322 NN323 NN324 {
recode `var' 0=.
fre `var'
}

*52. Depression/Anxiety: Selective items from the (Hopkins) Symptoms Checklist-25 (SCL-25)
fre NN325 NN326 NN327 NN328 NN329 NN330 NN331 NN332
*fix multiple ticks
foreach var of varlist NN325 NN326 NN327 NN328 NN329 NN330 NN331 NN332 {
recode `var' 0=.
fre `var'
}

**********************
*Height and weight
*child height cm and weight kg
tab1 NN24 NN25
*mother height cm and weight kg
tab1 NN283 NN284
*educ:
*across 2nd 3rd 4th grade
tab NN12

*clean vars for child and mother's BMI:
*mother's report of child's height and weight
foreach var in NN24 NN25 NN283 NN284 {
summ `var', det
return list
replace `var'=. if `var'<r(mean)-(4*r(sd))
replace `var'=. if `var'>r(mean)+(4*r(sd))
summ `var', det
}
*now make BMI:
gen moffspring_bmi_Q8y=NN25/((NN24/100)^2)
label variable moffspring_bmi_Q8y "Q8 offspring bmi from mother-report h & w"
summ moffspring_bmi_Q8y, det
*looks acceptable
gen mothers_bmi_Q8y=NN284/((NN283/100)^2)
label variable mothers_bmi_Q8y "Q8 mother's bmi from self-report h & w"
summ mothers_bmi_Q8y, det
*relevant var for whether pregnant:
tab NN302

*Mothers smoking and drinking: Q56 to 60
*smoking
capture drop mumsmokes_Q8y
gen mumsmokes_Q8y=NN344
recode mumsmokes_Q8y 4/7=.
*bump down
replace mumsmokes_Q8y=mumsmokes_Q8y-1
label var mumsmokes_Q8y "mother's smoking status at Q8y"
label define mumsmokes_Q8y 0"never" 1"sometimes" 2"daily"
label values mumsmokes_Q8y mumsmokes_Q8y
tab mumsmokes_Q8y NN344

*mother report of father's smoking:
fre NN347 NN348 NN349
*smoking
capture drop dadsmokes_Q8y
gen dadsmokes_Q8y=NN347
recode dadsmokes_Q8y 4/7=.
*bump down
replace dadsmokes_Q8y=dadsmokes_Q8y-1
label var dadsmokes_Q8y "mothe report of father's smoking status at Q8y"
label define dadsmokes_Q8y 0"never" 1"sometimes" 2"daily"
label values dadsmokes_Q8y dadsmokes_Q8y
tab dadsmokes_Q8y NN347

*use zanthro to get a more valid measure of bmi for children!
list AGE_SENT_MTHS_Q8AAR AGE_MTHS_Q8AAR AGE_RETURN_MTHS_Q8AAR in 1/10
*eh?
tab AGE_MTHS_Q8AAR
list AGE_SENT_MTHS_Q8AAR AGE_MTHS_Q8AAR AGE_RETURN_MTHS_Q8AAR if AGE_MTHS_Q8AAR<70
*This one looks most useful:
tab AGE_RETURN_MTHS_Q8AAR
*use AGE_RETURN_MTHS_Q8AAR and KJONN with zanthro
gen AGE_RETURN_YRS_Q8AAR=AGE_RETURN_MTHS_Q8AAR/12
summ AGE_RETURN_YRS_Q8AAR

gen KJONN_orig=KJONN
label variable KJONN_orig "0 Not specified, 1 Male 2 Female 3 Uncertain"
recode KJONN 3=. 0=.
fre KJONN

capture drop  zanthro_bmi_Q8y
egen zanthro_bmi_Q8y=zanthro(moffspring_bmi_Q8y, ba, WHO), xvar(AGE_RETURN_MTHS_Q8AAR) gender(KJONN) gencode (male=1, female=2) ageunit(month) nocutoff
summ zanthro_bmi_Q8y, det

count if zanthro_bmi_Q8y<-3
*304
count if zanthro_bmi_Q8y<-4
*81
list zanthro_bmi_Q8y NN24 NN25 moffspring_bmi_Q8y AGE_RETURN_MTHS_Q8AA KJONN if zanthro_bmi_Q8y<-4

**********************************************************************************************************************************

*Throughout, remove unusable obs where >1 category selected
/*
######### q8yr 
names(q8yr)
q8yr[q8yr=='More than 1 check box filled in'] <-NA
q8yr <- q8yr %>% mutate_if(is.factor,as.numeric)
*/
*Translate: for each item, recode 0s to missings

**********************************************************************************************************************************

*8. Has the child experienced development/hyperactivity/concentration, autistic/behavioural/emotional difficulties; currently/in the past/if yes, referred to a specialist?

fre NN36 NN37 NN38 NN39 NN40 NN41 NN42 NN43 NN44 NN45 NN46 NN47 NN48 NN49 NN50 NN51 NN53 NN54 NN55 NN56 NN57 NN58 NN59 NN60 NN61 NN62 NN63 NN64 NN65 NN66 NN67

*leave for now.

**********************************************************************************************************************************

*10. Short Mood and Feelings Questionnaire (SMFQ)
/*
*Description of original scale: Short Mood and Feelings Questionnaire (SMFQ)
The Mood and Feelings Quesionnaire (MFQ Angold & Costello, 1987) is a 32-item questionnaire based on DSM-III-R criteria for depression. 
The MFQ consists of a series of descriptive phrases regarding how the subject has been feeling or acting recently. 
A 13-item short form was developed, based on the discriminating ability between the depressed and non-depressed (Angold, et al., 1995). 
Both parent and child-report forms are available. The parent version is used in the MoBa 8-year questionnaire.
*/
foreach var of varlist NN68-NN80 {
tab `var' if `var'==0
recode `var' 0=.
}
*Nothing needs reverse-coding

**********************************************************************************************************************************

 *11. Short Norwegian Hierarchical Personality Inventory for Children (NHiPIC-30)
 *NN81-NN110
 tab NN81
 tab NN81, nol
 /*
 Description of original scale: The Hierarchical Personality Inventory for Children (HiPIC)
The HiPIC (Mervielde & De Fruyt, 1999, Mervielde & De Fruyt, 2002) is a questionnaire measuring the Big Five personality factors in children and adolescents. 
By means of 144 items, the HiPIC assesses five broad personality traits: Extraversion, Benevolence, Neuroticism, Conscientiousness, and Imagination. 
Each HiPIC item refers to a specific overt behaviour and is formulated in the third-person singular. 
Items are rated on a five-point Likert scale ranging from ‘not typical’ (1) to very typical (5). 
This section used the 30-item short form, also referred to as NHiPIC-30 (Vollrath, Hampson and Torgersen, submitted 2013). 
It contains five domain scales with 6 items each: 
Extraversion (items 8, 10, 18, 19, 23, 25), 
Benevolence (items 4, 5, 11, 16, 21, 28), 
Conscientiousness (items 3, 7, 9, 13, 15, 26), 
Neuroticism (items 1, 6, 14, 17, 22, 27), 
Imagination (2, 12, 20, 24, 29, 30).
 */
 *Recode to missing >1 response for all items
foreach var of varlist NN81-NN110  {
tab `var' if `var'==0
recode `var' 0=.
}
/*
 #reverse code the HPIC-30 variables NN82, NN83,NN84,NN90,NN92,NN95,NN99,NN100,NN103,NN104,NN105,NN106,NN107,NN108,NN109,NN110
q8yr[,c("rNN82", "rNN83","rNN84","rNN90","rNN92","rNN95","rNN99","rNN100", "rNN103", "rNN104", "rNN105", "rNN106", "rNN107", "rNN108", "rNN109", "rNN110")]<- 
  abs(q8yr[,c("NN82","NN83", "NN84","NN90","NN92", "NN95", "NN99", "NN100", "NN103","NN104", "NN105", "NN106", "NN107", "NN108", "NN109", "NN110" )]-7)
*/

*Translate
foreach var in NN82 NN83 NN84 NN90 NN92 NN95 NN99 NN100 NN103 NN104 NN105 NN106 NN107 NN108 NN109 NN110 {
gen r`var'=abs(`var'-6)
 }
 *inspect
tab NN82 rNN82

**********************************************************************************************************************************

*12-13. Parent/Teacher Rating Scale for Disruptive Behaviour Disorders (RS-DBD)
/*
Description of original scale: 
Parent/Teacher Rating Scale for Disruptive Behaviour Disorders (RS-DBD) Parent/Teacher Rating Scale for Disruptive Behavior Disorders (RS-DBD; Silva et al., 2005) 
consists of 41 DSM-IV items; with 18 items related to ADHD, 8 items related to Oppositional Defiant (OD), and 15 items to Conduct Disorder (CD). 
The 18 items (items 1-18 of section 13) related to ADHD, the 8 items related to OD (items 19-26 of section 13), and 8 items to CD were selected into use in this section.
Each item is rated on a four-point scale (1 = never/rarely, 2 = sometimes, 3 = often, 4 = very often).
*/
 *Recode to missing >1 response for all items
foreach var of varlist NN111-NN144 {
tab `var' if `var'==0
recode `var' 0=.
tab `var', nol
}
*Nothing needs reverse-coding

**********************************************************************************************************************************

*14.Screen for Child Anxiety Related Disorders (SCARED)
/*
Description of original scale: Screen for Child Anxiety Related Disorders (SCARED)
The SCARED (Birmaher et al., 1997) is a multidimensional questionnaire that purports to measure DSM-defined anxiety symptom. 
It contains 41 items which can be allocated to five separate anxiety subscales. 
Four of these subscales represent anxiety disorders that correspond with DSM categories, namely panic disorder, generalized anxiety disorder, social phobia, 
and separation anxiety. The fifth subscale is school phobia. The SCARED comes in two versions: a parent version and a child version. 
The 5-item short version, as used in the MoBa, was developed in Birmaher et al. (1999). 
Mothers rate how true the statements describe their children using a 3-point scale (i.e. 1= Not true, 2=Sometimes true, 3=True).
*/

 *Recode to missing >1 response for all items
foreach var of varlist NN145-NN149  {
tab `var' if `var'==0
recode `var' 0=.
}
*Nothing needs reverse coding.

**********************************************************************************************************************************
*15-17. Social Communication Questionnaire (SCQ)
/*
*Description of original instrument: Social Communication Questionnaire (SCQ)
The SCQ (Ritter, et al., 2003) is a parental-report Autism screening tool developed to serve as a practical piece of early childhood developmental screenings 
which parallels the Autism Diagnostic Interview-Revised (ADI-R; Lord, et al., 1994). 
It is a 40-question screening form designed for children with an age of 4.0 years (and a mental age of 2.0) which takes less than 10 minutes to complete and score. 
The items are administered in a yes/no response format.
*/

*NB: very uneven missingness between the items. 
*most missing <0.5% but NN172, NN181 missing >2%, NN183 missing 3.22%
*has the effect that social subscale has a lot more missingness than repetitive traits subscale
foreach var of varlist NN150-NN189  {
*display "`var'"
*count if `var'==.
*count if `var'==0
fre `var'
}

 *Recode to missing >1 response for all items
foreach var of varlist NN150-NN189  {
tab `var' if `var'==0
recode `var' 0=.
}
/*R code:
#reverse code the SCQ variables NN152,NN153,NN154,NN155,NN156,NN157,NN159,NN160,NN161,NN162,NN163,NN164,NN165,NN166,NN167
q8yr[,c("rNN152", "rNN153","rNN154","rNN155","rNN156","rNN157","rNN159","rNN160", "rNN161","rNN162", "rNN163", "rNN164", "rNN165","rNN166", "rNN167")]<- 
  abs(q8yr[,c("NN152","NN153", "NN154","NN155","NN156", "NN157", "NN159", "NN160", "NN161","NN162", "NN163", "NN164", "NN165", "NN166", "NN167")]-4)
 */

*Translate (differently, easier way to preserve correct weighting between reverse-coded and other items)
foreach var in NN152 NN153 NN154 NN155 NN156 NN157 NN159 NN160 NN161 NN162 NN163 NN164 NN165 NN166 NN167 {
fre `var'
gen r`var'=`var'
recode r`var' 2=0
replace r`var'=r`var'+1
*check
tab `var' r`var'
}
*inspect
tab NN152 rNN152

*************************************************************

*Question 20: Children's communication checklist
fre NN211 NN212 NN213 NN214 NN215 NN216 NN217 NN218 NN219 NN220 NN221 NN222 NN223 NN224 NN225 NN226

foreach var of varlist NN211-NN226  {
tab `var' if `var'==0
recode `var' 0=.
}

*************************************************************

*Question 21: The checklist of 20 Statements about Language-Related Difficulties (Språk 20)
fre NN227 NN228 NN229 NN230 NN232 NN233 NN374

foreach var of varlist NN227 NN228 NN229 NN230 NN232 NN233 NN374 {
tab `var' if `var'==0
recode `var' 0=.
}

*#############################################################################################################################
*#############################################################################################################################
*# Step 4: create sum scores for all total scales longitudinally 
*#############################################################################################################################

**##8yrs

*alldata$RSDBD_ADHD_8yr <-(alldata$NN119 + alldata$NN120 + alldata$NN121 + alldata$NN122 + alldata$NN123 + alldata$NN124 
*                         + alldata$NN125 + alldata$NN126 + alldata$NN127 + alldata$NN128 + alldata$NN129 + alldata$NN130
*                          + alldata$NN131 + alldata$NN132 + alldata$NN133 + alldata$NN134 + alldata$NN135 + alldata$NN136)
*describe(alldata$RSDBD_ADHD_8yr)

*Translate
gen ADHD_8yr=(NN119 + NN120 + NN121 + NN122 + NN123 + NN124 + NN125 + NN126 + NN127 + NN128 + NN129 + NN130 + NN131 + NN132 + NN133 + NN134 + NN135 + NN136)
*make the scales start at 0 so you can use nbreg or poisson
replace ADHD_8yr =ADHD_8yr-18
sum ADHD_8yr 

*alldata$RSDBD_inADHD_8yr <-(alldata$NN119 + alldata$NN120 + alldata$NN121 + alldata$NN122 + alldata$NN123 + alldata$NN124 + alldata$NN125 + alldata$NN126 + alldata$NN127)
*describe(alldata$RSDBD_inADHD_8yr)					
*Translate				
gen inADHD_8yr =(NN119 + NN120 + NN121 + NN122 + NN123 + NN124 + NN125 + NN126 + NN127)
*make the scales start at 0 so you can use nbreg or poisson
replace inADHD_8yr=inADHD_8yr-9
summ inADHD_8yr
						

*alldata$RSDBD_hyADHD_8yr <-(alldata$NN128 + alldata$NN129 + alldata$NN130 + alldata$NN131 + alldata$NN132 + alldata$NN133 + alldata$NN134 + alldata$NN135 + alldata$NN136)
*describe(alldata$RSDBD_hyADHD_8yr)
*Translate
gen hyADHD_8yr =(NN128 + NN129 + NN130 + NN131 + NN132 + NN133 + NN134 + NN135 + NN136)
*make the scales start at 0 so you can use nbreg or poisson
replace hyADHD_8yr=hyADHD_8yr-9
summ hyADHD_8yr

/*not looking at OD or CD in this analysis
*alldata$RSDBD_OD_8yr <-(alldata$NN137 + alldata$NN138 + alldata$NN139 + alldata$NN140 + alldata$NN141 + alldata$NN142 + alldata$NN143 + alldata$NN144)
*describe(alldata$RSDBD_OD_8yr)
*Translate
gen OD_8yr =(NN137 + NN138 + NN139 + NN140 + NN141 + NN142 + NN143 + NN144)
*make the scales start at 0 so you can use nbreg or poisson
replace OD_8yr=OD_8yr-8
summ OD_8yr

*alldata$RSDBD_CD_8yr <-(alldata$NN111 + alldata$NN112 + alldata$NN113 + alldata$NN114 + alldata$NN115 + alldata$NN116 + alldata$NN117 + alldata$NN118)
*describe(alldata$RSDBD_CD_8yr)
*Translate
gen CD_8yr =(NN111 + NN112 + NN113 + NN114 + NN115 + NN116 + NN117 + NN118)
summ CD_8yr
*make the scales start at 0 so you can use nbreg or poisson
replace CD_8yr=CD_8yr-8
summ CD_8yr
*/

*alldata$SCARED_8yr <-(alldata$NN145 + alldata$NN146 + alldata$NN147 + alldata$NN148 + alldata$NN149)
*describe(alldata$SCARED_8yr)
*Translate
gen SCARED_8yr =(NN145 + NN146 + NN147 + NN148 + NN149)
summ SCARED_8yr
*make the scales start at 0 so you can use nbreg or poisson
replace SCARED_8yr=SCARED_8yr-5
sum SCARED_8yr

*alldata$MFQ_8yr <-(alldata$NN68 + alldata$NN69 + alldata$NN70 + alldata$NN71 + alldata$NN72 + alldata$NN73 + alldata$NN74 + alldata$NN75 + alldata$NN76 + alldata$NN77
*                   + alldata$NN78 + alldata$NN79 + alldata$NN80)
*describe(alldata$MFQ_8yr)
*Translate
gen MFQ_8yr =(NN68 + NN69 + NN70 + NN71 + NN72 + NN73 + NN74 + NN75 + NN76 + NN77 + NN78 + NN79 + NN80)
*make the scales start at 0 so you can use nbreg or poisson
replace MFQ_8yr=MFQ_8yr-13
sum MFQ_8yr

*alldata$SCQ_8yr <-(alldata$NN150 + alldata$NN151 + alldata$rNN152 + alldata$rNN153 + alldata$rNN154 + alldata$rNN155 + alldata$rNN156 + alldata$rNN157 + alldata$NN158 + alldata$rNN159 +
*                    alldata$rNN160 + alldata$rNN161 + alldata$rNN162 + alldata$rNN163 + alldata$rNN164 + alldata$rNN165 + alldata$rNN166 + alldata$rNN167 + alldata$NN168 +
*                     alldata$NN169 + alldata$NN170 + alldata$NN171 + alldata$NN172 + alldata$NN173 + alldata$NN174 + alldata$NN175 + alldata$NN176 + alldata$NN177 + alldata$NN178 +
*                     alldata$NN179 + alldata$NN180 + alldata$NN181 + alldata$NN182 + alldata$NN183 + alldata$NN184 + alldata$NN185 + alldata$NN186 + alldata$NN187)
*describe(alldata$SCQ_8yr)

*Translate
gen SCQ_8yr =(NN150 + NN151 + rNN152 + rNN153 + rNN154 + rNN155 + rNN156 + rNN157 + NN158 + rNN159 + ///
                     rNN160 + rNN161 + rNN162 + rNN163 + rNN164 + rNN165 + rNN166 + rNN167 + NN168 + ///
                     NN169 + NN170 + NN171 + NN172 + NN173 + NN174 + NN175 + NN176 + NN177 + NN178 + ///
                     NN179 + NN180 + NN181 + NN182 + NN183 + NN184 + NN185 + NN186 + NN187 + NN188 + NN189)
*make the scales start at 0 so you can use nbreg or poisson
replace SCQ_8yr=SCQ_8yr-40
sum SCQ_8yr

*SEPARATE SYMPTOM CLUSTERS:
*split SCQ into social and repetitive:
*from the spreadsheet, there are three 'groupings': language, behaviour, social development
*also, there are three 'syndrome scales':  communication / Restricted, repetitive, and stereotyped patterns of behavior /Reciprosal Social Interaction
*and a few items don't map to any of these
*thirdly, DSM-oriented scales: SCI and RRB, with a few items mapping to neither
*SPLIT ON THESE LINES:
*SCI
gen S_SCQ_8yr =(NN150 + NN151 + rNN153 + NN158 + rNN159 + ///
                     NN168 + ///
                     NN169 + NN170 + NN171 + NN172 + NN173 + NN174 + NN175 + NN176 + NN177 + NN178 + ///
                     NN179 + NN180 + NN181 + NN182 + NN183 + NN185 + NN186 + NN187 + NN188 + NN189)
*make the scales start at 0 so you can use nbreg or poisson
replace S_SCQ_8yr=S_SCQ_8yr-26
sum S_SCQ_8yr

*RRB
gen R_SCQ_8yr = (rNN152 + rNN154 + rNN155 + rNN156 + rNN157 + rNN160 + rNN161 + rNN162 + rNN163 + rNN164 + rNN165 + rNN167)
*make the scales start at 0 so you can use nbreg or poisson
replace R_SCQ_8yr=R_SCQ_8yr-12
sum R_SCQ_8yr

********************************************************************************************************

*calculate % item-level missingness for each scale:

gen nitems_MFQ_8yr=0
foreach var of varlist NN68 NN69 NN70 NN71 NN72 NN73 NN74 NN75 NN76 NN77 NN78 NN79 NN80 {
replace nitems_MFQ_8yr=nitems_MFQ_8yr+1 if `var'!=.
}
fre nitems_MFQ_8yr if present_Q8y==1

*adhd 
fre NN119 NN120 NN121 NN122 NN123 NN124 NN125 NN126 NN127 NN128 NN129 NN130 NN131 NN132 NN133 NN134 NN135 NN136 
gen nitems_ADHD_8yr=0
foreach var of varlist NN119 NN120 NN121 NN122 NN123 NN124 NN125 NN126 NN127 NN128 NN129 NN130 NN131 NN132 NN133 NN134 NN135 NN136 {
replace nitems_ADHD_8yr=nitems_ADHD_8yr+1 if `var'!=.
}
fre nitems_ADHD_8yr if present_Q8y==1

*inattention subscale
gen nitems_inADHD_8yr=0
foreach var of varlist NN119 NN120 NN121 NN122 NN123 NN124 NN125 NN126 NN127 {
replace nitems_inADHD_8yr=nitems_inADHD_8yr+1 if `var'!=.
}
fre nitems_inADHD_8yr if present_Q8y==1

*hyperactivity subscale
gen nitems_hyADHD_8yr=0
foreach var of varlist NN128 NN129 NN130 NN131 NN132 NN133 NN134 NN135 NN136 {
replace nitems_hyADHD_8yr=nitems_hyADHD_8yr+1 if `var'!=.
}
fre nitems_hyADHD_8yr if present_Q8y==1

*oppositional defiant
gen nitems_OD_8yr=0
foreach var of varlist NN137 NN138 NN139 NN140 NN141 NN142 NN143 NN144 {
replace nitems_OD_8yr=nitems_OD_8yr+1 if `var'!=.
}
fre nitems_OD_8yr if present_Q8y==1

*conduct disorder
gen nitems_CD_8yr=0
foreach var of varlist NN111 NN112 NN113 NN114 NN115 NN116 NN117 NN118 {
replace nitems_CD_8yr=nitems_CD_8yr+1 if `var'!=.
}
fre nitems_CD_8yr if present_Q8y==1

*asd
gen nitems_SCQ_8yr=0
foreach var of varlist NN150 NN151 rNN152 rNN153 rNN154 rNN155 rNN156 rNN157 NN158 rNN159 rNN160 rNN161 rNN162 rNN163 rNN164 rNN165 rNN166 rNN167 NN168 NN169 NN170 NN171 NN172 NN173 NN174 NN175 NN176 NN177 NN178 NN179 NN180 NN181 NN182 NN183 NN184 NN185 NN186 NN187 NN188 NN189 {
replace nitems_SCQ_8yr=nitems_SCQ_8yr+1 if `var'!=.
}
fre nitems_SCQ_8yr if present_Q8y==1

*social subscale
gen nitems_S_SCQ_8yr=0
foreach var of varlist NN150 NN151 rNN153 NN158 rNN159 NN168 NN169 NN170 NN171 NN172 NN173 NN174 NN175 NN176 NN177 NN178 NN179 NN180 NN181 NN182 NN183 NN185 NN186 NN187 NN188 NN189 {
replace nitems_S_SCQ_8yr=nitems_S_SCQ_8yr+1 if `var'!=.
}
fre nitems_S_SCQ_8yr if present_Q8y==1

*repetitive subscale 
gen nitems_R_SCQ_8yr=0
foreach var of varlist rNN152 rNN154 rNN155 rNN156 rNN157 rNN160 rNN161 rNN162 rNN163 rNN164 rNN165 rNN167 {
replace nitems_R_SCQ_8yr=nitems_R_SCQ_8yr+1 if `var'!=.
}
fre nitems_R_SCQ_8yr if present_Q8y==1

*anxiety
gen nitems_SCARED_8yr=0
foreach var of varlist NN145 NN146 NN147 NN148 NN149 {
replace nitems_SCARED_8yr=nitems_SCARED_8yr+1 if `var'!=.
}
fre nitems_SCARED_8yr if present_Q8y==1

******************************************************************
*Next, make versions allowing up to 20% missing:
*1 item for MFQ, 2 for ADHD, 1 for each of the ADHD subscales, 7 for full SCQ, 4 for SCQ-social, 2 for SCQ-repetititive, 1 for scared
display 13/5
display 18/5
display 9/5
display 40/5
display 26/5
display 12/5
display 5/5

*MFQ
*make a new variable, equal to the complete-case version for complete cases:
egen miss20pc_MFQ_8yr=rowtotal (NN68 NN69 NN70 NN71 NN72 NN73 NN74 NN75 NN76 NN77 NN78 NN79 NN80)
replace miss20pc_MFQ_8yr=miss20pc_MFQ_8yr*(13/12) if nitems_MFQ_8yr==12
replace miss20pc_MFQ_8yr=. if nitems_MFQ_8yr<12 | nitems_MFQ_8yr==.
summ miss20pc_MFQ_8yr
*actually, bump this all down - since lowest value for any item was 1 not 0, mimimum possible value is not 0 which is weird
replace miss20pc_MFQ_8yr=miss20pc_MFQ_8yr-13
summ miss20pc_MFQ_8yr
fre miss20pc_MFQ

*ADHD, and subscales
*make a new variable, equal to the complete-case version for complete cases:
egen miss20pc_ADHD_8yr=rowtotal (NN119 NN120 NN121 NN122 NN123 NN124 NN125 NN126 NN127 NN128 NN129 NN130 NN131 NN132 NN133 NN134 NN135 NN136)
replace miss20pc_ADHD_8yr=miss20pc_ADHD_8yr*(18/17) if nitems_ADHD_8yr==17
replace miss20pc_ADHD_8yr=miss20pc_ADHD_8yr*(18/16) if nitems_ADHD_8yr==16
replace miss20pc_ADHD_8yr=miss20pc_ADHD_8yr*(18/15) if nitems_ADHD_8yr==15
replace miss20pc_ADHD_8yr=. if nitems_ADHD_8yr<15 | nitems_ADHD_8yr==.
summ miss20pc_ADHD_8yr
*actually, bump this all down - since lowest value for any item was 1 not 0, mimimum possible value is not 0 which is weird
replace miss20pc_ADHD_8yr=miss20pc_ADHD_8yr-18
fre miss20pc_ADHD_8yr

*inattention subscale
egen miss20pc_inADHD_8yr=rowtotal (NN119 NN120 NN121 NN122 NN123 NN124 NN125 NN126 NN127)
replace miss20pc_inADHD_8yr=miss20pc_inADHD_8yr*(9/8) if nitems_inADHD_8yr==8
replace miss20pc_inADHD_8yr=. if nitems_inADHD_8yr<8 | nitems_inADHD_8yr==.
summ miss20pc_inADHD_8yr
*actually, bump this all down - since lowest value for any item was 1 not 0, mimimum possible value is not 0 which is weird
replace miss20pc_inADHD_8yr=miss20pc_inADHD_8yr-9
fre miss20pc_inADHD_8yr

*hyperactivity subscale						
egen miss20pc_hyADHD_8yr=rowtotal (NN128 NN129 NN130 NN131 NN132 NN133 NN134 NN135 NN136)
replace miss20pc_hyADHD_8yr=miss20pc_hyADHD_8yr*(9/8) if nitems_hyADHD_8yr==8
replace miss20pc_hyADHD_8yr=. if nitems_hyADHD_8yr<8 | nitems_hyADHD_8yr==.
summ miss20pc_hyADHD_8yr
*actually, bump this all down - since lowest value for any item was 1 not 0, mimimum possible value is not 0 which is weird
replace miss20pc_hyADHD_8yr=miss20pc_hyADHD_8yr-9
fre miss20pc_hyADHD_8yr

*SCQ, and subscales
*make a new variable, equal to the complete-case version for complete cases:
egen miss20pc_SCQ_8yr=rowtotal (NN150 NN151 rNN152 rNN153 rNN154 rNN155 rNN156 rNN157 NN158 rNN159 rNN160 rNN161 rNN162 rNN163 rNN164 rNN165 rNN166 rNN167 NN168 NN169 NN170 NN171 NN172 NN173 NN174 NN175 NN176 NN177 NN178 NN179 NN180 NN181 NN182 NN183 NN184 NN185 NN186 NN187 NN188 NN189)
replace miss20pc_SCQ_8yr=miss20pc_SCQ_8yr*(40/39) if nitems_SCQ_8yr==39 
replace miss20pc_SCQ_8yr=miss20pc_SCQ_8yr*(40/38) if nitems_SCQ_8yr==38 
replace miss20pc_SCQ_8yr=miss20pc_SCQ_8yr*(40/37) if nitems_SCQ_8yr==37 
replace miss20pc_SCQ_8yr=miss20pc_SCQ_8yr*(40/36) if nitems_SCQ_8yr==36 
replace miss20pc_SCQ_8yr=miss20pc_SCQ_8yr*(40/35) if nitems_SCQ_8yr==35 
replace miss20pc_SCQ_8yr=miss20pc_SCQ_8yr*(40/34) if nitems_SCQ_8yr==34 
replace miss20pc_SCQ_8yr=miss20pc_SCQ_8yr*(40/33) if nitems_SCQ_8yr==33 
replace miss20pc_SCQ_8yr=miss20pc_SCQ_8yr*(40/32) if nitems_SCQ_8yr==32 
replace miss20pc_SCQ_8yr=. if nitems_SCQ_8yr<32 | nitems_SCQ_8yr==.
summ miss20pc_SCQ_8yr
*actually, bump this all down - since lowest value for any item was 1 not 0, mimimum possible value is not 0 which is weird
replace miss20pc_SCQ_8yr=miss20pc_SCQ_8yr-40
fre miss20pc_SCQ_8yr

*SCI
egen miss20pc_S_SCQ_8yr=rowtotal (NN150 NN151 rNN153 NN158 rNN159 NN168 NN169 NN170 NN171 NN172 NN173 NN174 NN175 NN176 NN177 NN178 NN179 NN180 NN181 NN182 NN183 NN185 NN186 NN187 NN188 NN189)
replace miss20pc_S_SCQ_8yr=miss20pc_S_SCQ_8yr*(26/25) if nitems_S_SCQ_8yr==25
replace miss20pc_S_SCQ_8yr=miss20pc_S_SCQ_8yr*(26/24) if nitems_S_SCQ_8yr==24
replace miss20pc_S_SCQ_8yr=miss20pc_S_SCQ_8yr*(26/23) if nitems_S_SCQ_8yr==23
replace miss20pc_S_SCQ_8yr=miss20pc_S_SCQ_8yr*(26/22) if nitems_S_SCQ_8yr==22
replace miss20pc_S_SCQ_8yr=miss20pc_S_SCQ_8yr*(26/21) if nitems_S_SCQ_8yr==21
replace miss20pc_S_SCQ_8yr=. if nitems_S_SCQ_8yr<21 | nitems_S_SCQ_8yr==.
summ miss20pc_S_SCQ_8yr
*actually, bump this all down - since lowest value for any item was 1 not 0, mimimum possible value is not 0 which is weird
replace miss20pc_S_SCQ_8yr=miss20pc_S_SCQ_8yr-26
fre miss20pc_S_SCQ_8yr

*RRB 
egen miss20pc_R_SCQ_8yr=rowtotal (rNN152 rNN154 rNN155 rNN156 rNN157 rNN160 rNN161 rNN162 rNN163 rNN164 rNN165 rNN167)
replace miss20pc_R_SCQ_8yr=miss20pc_R_SCQ_8yr*(12/11) if nitems_R_SCQ_8yr==11
replace miss20pc_R_SCQ_8yr=miss20pc_R_SCQ_8yr*(12/10) if nitems_R_SCQ_8yr==10
replace miss20pc_R_SCQ_8yr=. if nitems_R_SCQ_8yr<10 | nitems_R_SCQ_8yr==.
summ miss20pc_R_SCQ_8yr
*actually, bump this all down - since lowest value for any item was 1 not 0, mimimum possible value is not 0 which is weird
replace miss20pc_R_SCQ_8yr=miss20pc_R_SCQ_8yr-12
fre miss20pc_R_SCQ_8yr

*anxiety
egen miss20pc_SCARED_8yr=rowtotal (NN145 NN146 NN147 NN148 NN149)
replace miss20pc_SCARED_8yr=miss20pc_SCARED_8yr*(5/4) if nitems_SCARED_8yr==4
replace miss20pc_SCARED_8yr=. if nitems_SCARED_8yr<4 | nitems_SCARED_8yr==.
summ miss20pc_SCARED_8yr
*actually, bump this all down - since lowest value for any item was 1 not 0, mimimum possible value is not 0 which is weird
replace miss20pc_SCARED_8yr=miss20pc_SCARED_8yr-5
fre miss20pc_SCARED_8yr


*check it:
foreach varstem in MFQ_8yr ADHD_8yr inADHD_8yr hyADHD_8yr SCQ_8yr S_SCQ_8yr R_SCQ_8yr {
	tab miss20pc_`varstem' `varstem', mi
}

***************************************************************************************
*in fact, rename so that these are the default and the complete-case versions are labelled as such.
foreach varstem in MFQ_8yr SCARED_8yr ADHD_8yr inADHD_8yr hyADHD_8yr SCQ_8yr S_SCQ_8yr R_SCQ_8yr {
	rename `varstem' cc_`varstem'
	rename miss20pc_`varstem' `varstem'
}

***************************************************************************************

*******************************************************************************************************
*Now cleaned, make standardized versions of everything:

foreach var in MFQ_8yr ADHD_8yr inADHD_8yr hyADHD_8yr /*OD_8yr CD_8yr*/ SCARED_8yr SCQ_8yr S_SCQ_8yr R_SCQ_8yr {
egen z`var'=std(`var')
summ z`var'
}

*Now combine info on mother's BMI from across waves

count if mothers_bmi_Q8y==.& moffspring_bmi_Q8y!=.
list mothers_bmi_Q1 mothers_bmi_Q6 mothers_bmi_Q5y LL508 mothers_bmi_Q8y NN302 if (mothers_bmi_Q8y==.  | NN302==2) & moffspring_bmi_Q8y!=.
tab LL508 if mothers_bmi_Q5y!=.
*?

gen mothers_bmi_last=mothers_bmi_Q8y
replace mothers_bmi_last=mothers_bmi_Q5y if mothers_bmi_last==. 
replace mothers_bmi_last=mothers_bmi_Q6 if mothers_bmi_last==. 
replace mothers_bmi_last=mothers_bmi_Q1 if mothers_bmi_last==. 
summ mothers_bmi_last if moffspring_bmi_Q8y!=., det

******************************************************************************************
*Educational outcomes and auxillary:
*before record linkage, using Q25 items NN239 NN240 NN241 as outcomes
*Q25: Three questions about school results on national exams
fre NN239 NN240 NN241

*also:
*Q22 How is your child enjoying school?
fre NN234
*school enjoyment: collpase and flip for bigger baseline categ:
gen NN234_new=NN234
label define NN234_new 1"very well" 2"well" 3"neither/poor/very poor"
label values NN234_new NN234_new
recode NN234_new 1/3=1 4=2 5=3
*more than one box ticked to missing:
recode NN234_new 0=.
fre NN234_new
replace NN234_new=4-NN234_new
fre NN234_new

*Q26. Special Education. Is an administrative decision made about your child being eligible for special education?
*In Norwegian language?
fre NN242 NN243
*In arithmetic?
fre NN244 NN245
*In other subjects?
fre NN246 NN247
*Does your child receive any other educational support?
fre NN248 NN249
*Does your child get extra help (e.g. an assistant) at school because of a disability or a developmental problem?
fre NN250 
*These look like grouped versions from NN243 NN245 NN247 NN249
fre NN375 NN376 NN377 NN378 NN379
*Although in the questionnaire, looks like you can only respond with categories. Weird. A versions thing?
*Just use the yes/no for now

*Remove the multiple ticks:
foreach var of varlist NN234 NN239 NN240 NN241 NN242 NN243 NN244 NN245 NN246 NN247 NN248 NN249 NN250 NN375 NN376 NN377 NN378 NN379{
tab `var' if `var'==0
recode `var' 0=.
fre `var'
}

**************************************************************************************************
*For the main outcome measures, NN239 NN240 NN241, 4='don't know/have not talked to teacher about it'
*worrying in terms of parental involvement, but I guess we code these to missing?
foreach var in NN239 NN240 NN241 {
recode `var' 4=.
}
**************************************************************************************************
*Homework? Is this measuring ability, SEP, something else? Does it matter?
*Q27 Three questions about homework, help with homework at home and at school
fre NN251 NN252 NN253

*Remove the multiple ticks, including pairs of adjacent categories:
foreach var of varlist NN251 NN252 NN253 {
tab `var' if `var'==0
recode `var' 0=. 
recode `var' 12/56=.
fre `var'
}

*28 Reading and writing skills
*NB ITEMS DIFFER BETWEEN VERSIONS 
*In version C
fre NN380 NN381 NN382 NN383 NN384 
*In version A & B
fre NN254 NN255 NN256 NN257

*Remove the multiple ticks, including pairs of adjacent categories:
foreach var of varlist NN380 NN381 NN382 NN383 NN384 {
tab `var' if `var'==0
recode `var' 0=. 
fre `var'
}

*Remove the multiple ticks, including pairs of adjacent categories:
*also the 'don't know' categ to missing
foreach var of varlist NN254 NN255 NN256 NN257 {
tab `var' if `var'==0
recode `var' 0=. 
recode `var' 3=. 
fre `var'
}

*29. The Child’s Pronunciation. 2 questions about understandability of the child’s speech; 2 questions about the child’s narrative skills
fre NN258 NN259 NN260 NN261
foreach var of varlist NN258 NN259 NN260 NN261 {
tab `var' if `var'==0
recode `var' 0=. 
recode `var' 23/45=. 
fre `var'
}
*30: Home reading and self-reading
fre NN262 NN263 NN385 NN264 NN386 NN265

*Remove the multiple ticks:
foreach var of varlist NN262 NN263 NN385 NN264 NN386 NN265 {
tab `var' if `var'==0
recode `var' 0=.
recode `var' 12/234=.
fre `var'
}
*For NN265, 5='don't know'
recode NN265 5=.

*self-reading: just flip. later, can collapse if needed
gen NN265_new=NN265
label define NN265_new 1"Books with chapters (almost only text)" 2"Simple stories, pictures & text on every page" 3"Picture book (few words)" 4"Does not like to read him-/herself"
label values NN265_new NN265_new
replace NN265_new=5-NN265_new
fre NN265_new

*34. Difficulties, Impairment and Impact
*Individual questions related to impact/impairments based on The Strengths and Difficulties Questionnaires (SDQ)
fre NN388 NN389 NN390 NN391 NN392 NN393 NN394 NN395 NN396
*multiple ticks:

foreach var of varlist NN388 NN389 NN390 NN391 NN392 NN393 NN394 NN395 NN396 {
tab `var' if `var'==0
recode `var' 0=.
fre `var'
}

****************************************************************************************
**ABOUT THE MUM
****************************************************************************************

*asked again about education - worth having
fre NN271
*no prep needed

*Q42-45: Mum's eating disorder 
*NB: some version differences in Q43, items NN293 NN367
fre NN285-NN301
fre NN293 NN367
*fix multiple ticks
foreach var of varlist NN285-NN301 NN367 {
recode `var' 0=.
fre `var'	
}

*Q 48:Mum's health problems (ADHD, reading/writing difficulties, anorexia, bulimia, other)
*ever had
fre NN304 NN306 NN308 NN310 NN312
*treated?
fre NN305 NN307 NN309 NN311 NN313

*for all but the last, fix multiple ticks: 0s to .
foreach var of varlist NN304-NN312 {
recode `var' 0=.
fre `var'	
}
*for imputation, MAKE BINARY VERSIONS FOR NEVER/EVER HAD - the middle category is too small to use
label define neverever 1"never" 2"yes, now or in past"
foreach var of varlist NN304 NN306 NN308 NN310 NN312 {
	gen `var'_binary=`var'
	recode `var'_binary 3=2
	label variable `var'_binary "ever had this condition, no/yes"
	label values `var'_binary neverever
	fre `var'_binary
}

*Q49: Social Phobia
*worth having in as relevant to anxiety.  GO BACK AND ADD FROM FATHER'S QUESTIONNAIRE ETC
fre NN314 NN315 NN316
*fix multiple ticks
foreach var of varlist NN314 NN315 NN316 {
recode `var' 0=.
fre `var'	
}

*Q50: Satisfaction with life
fre NN317-NN321
*fix multiple ticks
foreach var of varlist NN317-NN321 {
recode `var' 0=.
fre `var'	
}

*51. The Autonomic Nervous System Questionnaire (ANS)
fre NN322 NN323 NN324 
*fix multiple ticks
foreach var of varlist NN322 NN323 NN324  {
recode `var' 0=.
fre `var'	
}

*Hopkins checklist:
fre NN325-NN332
*fix multiple ticks
foreach var of varlist NN325-NN332  {
recode `var' 0=.
fre `var'	
}
********************************************************************************************
drop _merge
********************************************************************************************
count
*114,747

*temp save point
save "temp_phenotypes.dta", replace

****************************************************
*Last thing: fathers second questionnaire:
use "N:\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\Statafiles\PDB2306_Far2_V12.dta", clear
keep F_ID_2306 G__5 G__6 G__7_1 G_21 G_66 ///
/*health conditions*/ G__230_1 G__230_2 G__231_1 G__231_2 G__232_1 G__232_2 G__233_1 G__233_2 G__234_1 G__234_2 ///
/*satisfaction with life scale*/  G_51_1 G_51_2 G_51_2 G_51_3 G_51_4 G_51_5 ///
/*Hopkins 12-item version*/ G_52_1 G_52_2 G_52_3 G_52_4 G_52_5 G_52_6 G_52_7 G_52_8 G_52_9 G_5210 G_5211 G_5212 ///
/*suicide ideation and attempts*/ G_54 G_55 ///
/*psychosis symptoms*/ G_56_1_1 G_56_1_2 G_56_2_1 G_56_2_2 G_56_3_1 G_56_4_1 G_56_4_2 G_56_5_1 G_56_5_2 G_56_6_1 G_56_6_2 G_56_7_1 G_56_7_2 G_56_8_1 G_56_8_2 G_56_9_1 G_56_9_2
gen present_QF2=1 
*ERASE THIS DATA for any fathers who withdrew consent.
*need to merge on F_ID_2306, since PDB2306_Far2_V12.dta doesn't actually contain PREG_ID_2306:
merge 1:m F_ID_2306 using "N:\durable\projects\mh_bmi\fullsample\PDB2306_SV_INFO_v12_fathers.dta"
keep if fathers_consent==1
*63 excluded
*then merge back into growing file:
drop _merge
merge 1:m PREG_ID_2306 using "temp_phenotypes.dta"

*sort out:
*MH stuff:
/*specific relevant health conds and age when */
fre G__230_1 G__230_2 G__231_1 G__231_2 G__232_1 G__232_2 G__233_1 G__233_2 G__234_1 G__234_2
*fix multiple ticks in the yes/no items:
foreach var of varlist G__230_1 G__231_1 G__232_1 G__233_1 G__234_1 {
recode `var' 0=.
fre `var'	
}
/*satisfaction with life scale*/  
fre G_51_1 G_51_2 G_51_2 G_51_3 G_51_4 G_51_5 
*fix multiple ticks:
foreach var of varlist G_51_1 G_51_2 G_51_2 G_51_3 G_51_4 G_51_5  {
recode `var' 0=.
fre `var'	
}
/*Hopkins*/ 
*CAREFUL! IN THE SECOND FATEHR'S QUESTIONNAIRE, IT'S A 12-ITEM VERSION RATHER THAN THE USUAL 8
fre G_52_1 G_52_2 G_52_3 G_52_4 G_52_5 G_52_6 G_52_7 G_52_8 G_52_9 G_5210 G_5211 G_5212
*fix multiple ticks:
foreach var of varlist G_52_1 G_52_2 G_52_3 G_52_4 G_52_5 G_52_6 G_52_7 G_52_8 G_52_9 G_5210 G_5211 G_5212  {
recode `var' 0=.
fre `var'	
}
/*suicide ideation and attempts*/ 
fre G_54 G_55 
*multiple ticks:
recode G_54 0=.
recode G_55 0=.
fre G_54 G_55 

/*psychosis symptoms*/ 
fre G_56_1_1 G_56_1_2 G_56_2_1 G_56_2_2 G_56_3_1 G_56_4_1 G_56_4_2 G_56_5_1 G_56_5_2 G_56_6_1 G_56_6_2 G_56_7_1 G_56_7_2 G_56_8_1 G_56_8_2 G_56_9_1 G_56_9_2
*fix multiple ticks:
foreach var of varlist G_56_1_1 G_56_1_2 G_56_2_1 G_56_2_2 G_56_3_1 G_56_4_1 G_56_4_2 G_56_5_1 G_56_5_2 G_56_6_1 G_56_6_2 G_56_7_1 G_56_7_2 G_56_8_1 G_56_8_2 G_56_9_1 G_56_9_2  {
recode `var' 0=.
fre `var'	
}

*clean auxillary vars: smoking status and income bands
gen dadsmokes_QF2=G_21
recode dadsmokes_QF2 0=. 12/34=.
*bump down
replace dadsmokes_QF2=dadsmokes_QF2-1
label var dadsmokes_QF2 "father's smoking status at F2"
label define dadsmokes_QF2 0"never" 1"ex" 2"current, socially" 3"current, daily"
label values dadsmokes_QF2 dadsmokes_QF2
tab dadsmokes_QF2 G_21
*income
recode G_66 0=.

*finish Dad's BMI:
*clean vars:
foreach var in G__5 G__6 {
summ `var', det
return list
replace `var'=. if `var'<r(mean)-(4*r(sd))
replace `var'=. if `var'>r(mean)+(4*r(sd))
summ `var', det
}
*now make BMI:
gen fathers_bmi_FQ2=G__6/((G__5/100)^2)
label variable fathers_bmi_FQ2 "Father's second Q: bmi from own self-report h & w"
summ fathers_bmi_FQ2, det

drop _merge

********************************************************************************
**PREP FOR PARENTAL HOPKINS AND ADHD SUMMARY SCORES.

*As with kids measures at 8yr, make a complete-case version of the summary var and also an average item response, to use in the imputation.

*ADD UP PARENTAL HOPKINS CHECKLISTS
*to impute summary scores. To use the information from people who answered some items but not all, calculate per-item average and include in the imputation

*mother Q1
*missingness
gen nitems_mhopkins_Q1=0
foreach var of varlist AA1548 AA1549 AA1550 AA1551 AA1552 {
replace nitems_mhopkins_Q1=nitems_mhopkins_Q1+1 if `var'!=.
}
fre nitems_mhopkins_Q1 if present_Q1==1
*complete case:
gen cc_mhopkins_Q1=(AA1548 + AA1549 + AA1550 + AA1551 + AA1552)
*bump down to start at 0:
replace cc_mhopkins_Q1=cc_mhopkins_Q1-5
*allowing 20% missingness
egen miss20pc_mhopkins_Q1=rowtotal(AA1548 AA1549 AA1550 AA1551 AA1552)
replace miss20pc_mhopkins_Q1=miss20pc_mhopkins_Q1*(5/4) if nitems_mhopkins_Q1==4
replace miss20pc_mhopkins_Q1=. if nitems_mhopkins_Q1<4
*bump down
replace miss20pc_mhopkins_Q1=miss20pc_mhopkins_Q1-5
fre miss20pc_mhopkins_Q1
*check
tab miss20pc_mhopkins_Q1 cc_mhopkins_Q1 , mi
*rename to be default:
rename miss20pc_mhopkins_Q1 mhopkins_Q1


*QF
*missingness
gen nitems_fhopkins_QF=0
foreach var of varlist FF251 FF252 FF253 FF254 FF255 FF256 FF257 FF258 {
replace nitems_fhopkins_QF=nitems_fhopkins_QF+1 if `var'!=.
}
fre nitems_fhopkins_QF if present_QF==1 
*complete case:
gen cc_fhopkins_QF=(FF251 + FF252 + FF253 + FF254 + FF255 + FF256 + FF257 + FF258)
*bump down to start at 0:
replace cc_fhopkins_QF=cc_fhopkins_QF-8
*allowing 20% missingenss:
egen miss20pc_fhopkins_QF=rowtotal(FF251 FF252 FF253 FF254 FF255 FF256 FF257 FF258)
replace miss20pc_fhopkins_QF=miss20pc_fhopkins_QF*(8/7) if nitems_fhopkins_QF==7
replace miss20pc_fhopkins_QF=. if nitems_fhopkins_QF<7
*bump down
replace miss20pc_fhopkins_QF=miss20pc_fhopkins_QF-8
fre miss20pc_fhopkins_QF
*check
tab miss20pc_fhopkins_QF cc_fhopkins_QF , mi
*list cc_fhopkins_QF  nitems_fhopkins_QF miss20pc_fhopkins_QF FF251 FF252 FF253 FF254 FF255 FF256 FF257 FF258 if miss20pc_fhopkins_QF<0
*rename to be default:
rename miss20pc_fhopkins_QF fhopkins_QF

*Q6 HOPKINS
*missingess
gen nitems_mhopkins_Q6=0
foreach var of varlist GG514 GG515 GG516 GG517 GG518 GG519 GG520 GG521 {
replace nitems_mhopkins_Q6=nitems_mhopkins_Q6+1 if `var'!=.
}
fre nitems_mhopkins_Q6 if present_Q6==1
*complete case:
gen cc_mhopkins_Q6=(GG514 + GG515 + GG516 + GG517 + GG518 + GG519 + GG520 + GG521)
*bump down to start at 0:
replace cc_mhopkins_Q6=cc_mhopkins_Q6-8
*allowing 20% missingness
egen miss20pc_mhopkins_Q6=rowtotal(GG514 GG515 GG516 GG517 GG518 GG519 GG520 GG521)
replace miss20pc_mhopkins_Q6=miss20pc_mhopkins_Q6*(8/7) if nitems_mhopkins_Q6==7
replace miss20pc_mhopkins_Q6=. if nitems_mhopkins_Q6<7
*bump down
replace miss20pc_mhopkins_Q6=miss20pc_mhopkins_Q6-8
fre miss20pc_mhopkins_Q6
*check
tab miss20pc_mhopkins_Q6 cc_mhopkins_Q6, mi
*rename to be default:
rename miss20pc_mhopkins_Q6 mhopkins_Q6

*Q5yr HOPKINS
*missingess
gen nitems_mhopkins_Q5y=0
foreach var of varlist LL363 LL364 LL365 LL366 LL367 LL368 LL369 LL370 {
replace nitems_mhopkins_Q5y=nitems_mhopkins_Q5y+1 if `var'!=.
}
fre nitems_mhopkins_Q5y if present_Q5y==1
*complete case:
gen cc_mhopkins_Q5y=(LL363 + LL364 + LL365 + LL366 + LL367 + LL368 + LL369 + LL370)
*bump down to start at 0:
replace cc_mhopkins_Q5y=cc_mhopkins_Q5y-8
*allowing 20% missingness
egen miss20pc_mhopkins_Q5y=rowtotal(LL363 LL364 LL365 LL366 LL367 LL368 LL369 LL370)
replace miss20pc_mhopkins_Q5y=miss20pc_mhopkins_Q5y*(8/7) if nitems_mhopkins_Q5y==7
replace miss20pc_mhopkins_Q5y=. if nitems_mhopkins_Q5y<7
*bump down
replace miss20pc_mhopkins_Q5y=miss20pc_mhopkins_Q5y-8
fre miss20pc_mhopkins_Q5y
*check
tab miss20pc_mhopkins_Q5y cc_mhopkins_Q5y , mi
*rename to be default:
rename miss20pc_mhopkins_Q5y mhopkins_Q5y


*Q8yr HOPKINS
*missingess
gen nitems_mhopkins_Q8y=0
foreach var of varlist NN325 NN326 NN327 NN328 NN329 NN330 NN331 NN332  {
replace nitems_mhopkins_Q8y=nitems_mhopkins_Q8y+1 if `var'!=.
}
fre nitems_mhopkins_Q8y if present_Q8y==1
*complete case:
gen cc_mhopkins_Q8y=(NN325 + NN326 + NN327 + NN328 + NN329 + NN330 + NN331 + NN332)
*bump down to start at 0:
replace cc_mhopkins_Q8y=cc_mhopkins_Q8y-8
*allowing 20% missingness
egen miss20pc_mhopkins_Q8y=rowtotal(NN325 NN326 NN327 NN328 NN329 NN330 NN331 NN332)
replace miss20pc_mhopkins_Q8y=miss20pc_mhopkins_Q8y*(8/7) if nitems_mhopkins_Q8y==7
replace miss20pc_mhopkins_Q8y=. if nitems_mhopkins_Q8y<7
*bump down
replace miss20pc_mhopkins_Q8y=miss20pc_mhopkins_Q8y-8
fre miss20pc_mhopkins_Q8y
*check
tab miss20pc_mhopkins_Q8y cc_mhopkins_Q8y , mi
*rename to be default:
rename miss20pc_mhopkins_Q8y mhopkins_Q8y

*F2 HOPKINS
*missingess
gen nitems_fhopkins_QF2=0
foreach var of varlist G_52_1 G_52_2 G_52_3 G_52_4 G_52_5 G_52_6 G_52_7 G_52_8 G_52_9 G_5210 G_5211 G_5212 {
replace nitems_fhopkins_QF2=nitems_fhopkins_QF2+1 if `var'!=.
}
fre nitems_fhopkins_QF2 if present_QF2==1
*complete case:
gen cc_fhopkins_QF2=(G_52_1 + G_52_2 + G_52_3 + G_52_4 + G_52_5 + G_52_6 + G_52_7 + G_52_8 + G_52_9 + G_5210 + G_5211 + G_5212)
*bump down to start at 0:
replace cc_fhopkins_QF2=cc_fhopkins_QF2-12
*allowing 20% missingness
egen miss20pc_fhopkins_QF2=rowtotal(G_52_1 G_52_2 G_52_3 G_52_4 G_52_5 G_52_6 G_52_7 G_52_8 G_52_9 G_5210 G_5211 G_5212)
replace miss20pc_fhopkins_QF2=miss20pc_fhopkins_QF2*(12/11) if nitems_fhopkins_QF2==11
replace miss20pc_fhopkins_QF2=miss20pc_fhopkins_QF2*(12/10) if nitems_fhopkins_QF2==10
replace miss20pc_fhopkins_QF2=. if nitems_fhopkins_QF2<10
*bump down
replace miss20pc_fhopkins_QF2=miss20pc_fhopkins_QF2-12
fre miss20pc_fhopkins_QF2
*check
tab miss20pc_fhopkins_QF2 cc_fhopkins_QF2 , mi
*rename to be default:
rename miss20pc_fhopkins_QF2 fhopkins_QF2

*****************
*parental adhd symptoms
*maternal ADHD from Q6
fre GG503 GG504 GG505 GG506 GG507 GG508 
gen nitems_mADHD_Q6=0
foreach var of varlist GG503 GG504 GG505 GG506 GG507 GG508  {
replace nitems_mADHD_Q6=nitems_mADHD_Q6+1 if `var'!=.
}
fre nitems_mADHD_Q6 if present_Q6==1
*complete case
gen cc_mADHD_Q6=(GG503 + GG504 + GG505 + GG506 + GG507 + GG508)
*for descriptives and imputation, bump down to start at 0
replace cc_mADHD_Q6=cc_mADHD_Q6-6
fre cc_mADHD_Q6
*allowing 20% missingness
egen miss20pc_mADHD_Q6=rowtotal(GG503 GG504 GG505 GG506 GG507 GG508)
replace miss20pc_mADHD_Q6=miss20pc_mADHD_Q6*(6/5) if nitems_mADHD_Q6==5
replace miss20pc_mADHD_Q6=. if nitems_mADHD_Q6<5
*bump down
replace miss20pc_mADHD_Q6=miss20pc_mADHD_Q6-6
fre miss20pc_mADHD_Q6
*check
tab miss20pc_mADHD_Q6 cc_mADHD_Q6 , mi
*rename to be default:
rename miss20pc_mADHD_Q6 mADHD_Q6

*paternal ADHD from QF
fre FF535 FF536 FF537 FF538 FF539 FF540  
gen nitems_fADHD_QF=0
foreach var of varlist FF535 FF536 FF537 FF538 FF539 FF540  {
replace nitems_fADHD_QF=nitems_fADHD_QF+1 if `var'!=.
}
fre nitems_fADHD_QF if present_Q6==1
*complete case
gen cc_fADHD_QF=(FF535 + FF536 + FF537 + FF538 + FF539 + FF540)
*for descriptives and imputation, bump down to start at 0
replace cc_fADHD_QF=cc_fADHD_QF-6
fre cc_fADHD_QF
*allowing 20% missingness
egen miss20pc_fADHD_QF=rowtotal(FF535 FF536 FF537 FF538 FF539 FF540)
replace miss20pc_fADHD_QF=miss20pc_fADHD_QF*(6/5) if nitems_fADHD_QF==5
replace miss20pc_fADHD_QF=. if nitems_fADHD_QF<5
*bump down
replace miss20pc_fADHD_QF=miss20pc_fADHD_QF-6
fre miss20pc_fADHD_QF
*check
tab miss20pc_fADHD_QF cc_fADHD_QF , mi
*rename to be default:
rename miss20pc_fADHD_QF fADHD_QF

*******************************************************************************************************
*NOW DROP ANYONE (WHOLE FAMILY) WHERE THE MOTHER WITHDREW CONSENT:
merge m:1 PREG_ID_2306 using "N:\durable\projects\mh_bmi\fullsample\PDB2306_SV_INFO_v12_mothers_children.dta"
keep if mothers_consent==1

#WHO ARE THOSE OTHER 20?
#HAVEN'T DROPPED AYNONE YET.
*flag them:
gen only_in_SV_INFO=1 if _merge==2
drop _merge

count
*114,307


*SAVE IT:
save "data/phenotypes_bmi_analyses.dta", replace
*erase the temp file:
erase temp_phenotypes.dta


summ *hopkins*
summ *ADHD*
summ *MFQ*
summ *SCARED*
