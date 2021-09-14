
*IMPUTATION

*PRELIMINARIES:
*MERGE COMPLETE PHENOTYPE DATA TO ALL THE PGSs, including genetic covariates (PCs, batch vars etc)

cd "N:\data\durable\projects\mh_bmi\"
set maxvar 20000

use "data\statafiles\phenotypes_bmi_analyses.dta", clear
capture drop _merge

*merge all exposures in turn,
*keeping only children with genetic info 
foreach exp in childbmi adultbmi ea3 {
merge 1:1 PREG_ID_2306 BARN using scratch/reshaped_linkable_`exp'_genetic_pgs.dta
*, keep(match) 
drop _merge
} 

*delineate the usable trios:
gen fullgeneticinfo=0
replace fullgeneticinfo=1 if (c_batch!="" & m_batch!="" & f_batch!="")
tab fullgeneticinfo
*26,644 trios


*merge in PGS for outcomes to check for assortment:
capture drop _merge
foreach exp in depression anxiety adhd asd {
merge 1:1 PREG_ID_2306 BARN using scratch/reshaped_linkable_`exp'_genetic_pgs.dta
*, keep(match) 
drop _merge
} 

*sort out batch varaiables. MAKE DUMMIES OUT OF THESE AGAIN because otherwise can't put in models with mi regress etc, where 'factor-variables not allowed':
tab c_batch, gen(c_batch)
tab m_batch, gen(m_batch)
tab f_batch, gen(f_batch)
*DROP THE REDUNDANT BASELINE CATEG DUMMY, since this messes with imputation models later:
fre c_batch m_batch f_batch
*harvest for everyone, categ 1:
fre c_batch1 m_batch1 f_batch1
drop c_batch1 m_batch1 f_batch1
*drop the categorical ones to avoid confusion
drop c_batch m_batch f_batch

save "data\statafiles\geno_pheno_all.dta", replace

**************************************************************************************************************************
**************************************************************************************************************************

*FINAL PREP TO SUPPORT IMPUTATION 

clear all
cd "N:\data\durable\projects\mh_bmi\"
set maxvar 20000
sysdir
sysdir set PLUS "N:\data\durable\people\Mandy\Statafiles\ado\plus"

use "data\statafiles\geno_pheno_all.dta", clear
}

*****************************************************************************
*TRIM OFF PEOPLE WHO CONTRIBUTE NOTHING AND STOP IT FROM RUNNING

*UPDATE 12 MAY: RESTRICT TO PEOPLE IN UPDATED SV_INFO FILE (i.e. exclude the new withdrawls) 
merge m:1 PREG_ID_2306 using "N:\data\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\Statafiles\PDB2306_SV_INFO_v12_UPDATED_April_2021"
*the 51 not matched are the ones that need to be dropped as new exclusions.
drop if _merge==1
*51 dropped
drop _merge

count
*114,717

fre present*
*is there anyone present in the birth registry file but not at any of the actual questionnaires?
gen nQ_present=0
foreach q in Q1 QF Q3 Q4 Q5 Q6 Q5y Q7y Q8y F2 {
	replace nQ_present=nQ_present+1 if present_`q'==1
}
fre nQ_present
*OK, so there's a chunk of people who never were included in a single questionnaire.
*these are the people that weren't present at the very first one:
tab nQ_present present_Q1, mi

*these people basically contribute no information, and are most likely holding up imputation. 
*DROP THEM.

*for flowchart: how many of these among people with nonmissing sex?
count
*114,717
count if KJONN==.
*831
drop if KJONN==.
count
*113,886
count if KJONN!=. & nQ_present==0
*8,101

drop if nQ_present==0
count
*105,785

*who has genetic information? 
count if c_batch2!=. 
*30,887
*have lost some, relative to the number of children with genetic info: 31,374 (after removing one of the identical twins)
tab fullgeneticinfo
*26,370 trios

keep if fullgenetic==1
*26370 remaining

******************************************************************
*other things in birth registry file:
count if BARN_NR==.

*IMPUTATION DOES NOT WORK IF YOU KEEP EVERYONE
*SO RESTRICT TO FULL TRIOS
keep if fullgeneticinfo==1

****************************************************************************************************

*FINAL PREP OF PHENOTYPES FOR IMPUTATION:

*impute with untransformed variables since you will need them for table 1. then standardize in mi passive

*birth weight and length for imputation: trim at 4sd from the mean
foreach var in VEKT LENGDE {
summ `var', det
return list
replace `var'=. if `var'<r(mean)-(4*r(sd))
replace `var'=. if `var'>r(mean)+(4*r(sd))
summ `var', det
}

*CHECK FORMAT OF INDIVIDUAL ITEMS FOR AGE 8 OUTCOMES:

*RSDBD_ADHD_8yr
fre NN119 NN120 NN121 NN122 NN123 NN124 NN125 NN126 NN127 NN128 NN129 NN130 NN131 NN132 NN133 NN134 NN135 NN136
*(ologit) 

*MFQ_8yr
fre NN68 NN69 NN70 NN71 NN72 NN73 NN74 NN75 NN76 NN77 NN78 NN79 NN80
*(ologit) 

*SCARED_8yr
fre NN145 NN146 NN147 NN148 NN149
*(ologit)

*SCQ_8yr
fre NN150 NN151 rNN152 rNN153 rNN154 rNN155 rNN156 rNN157 NN158 rNN159 rNN160 rNN161 rNN162 rNN163 rNN164 rNN165 rNN166 rNN167 NN168 NN169 NN170 NN171 NN172 NN173 NN174 NN175 NN176 NN177 NN178 NN179 NN180 NN181 NN182 NN183 NN184 NN185 NN186 NN187 NN188 NN189
*(logit) 
				 
*S_SCQ_8yr
fre NN150 NN151 rNN153 NN158 rNN159 NN168 NN169 NN170 NN171 NN172 NN173 NN174 NN175 NN176 NN177 NN178 NN179 NN180 NN181 NN182 NN183 NN185 NN186 NN187 NN188 NN189
*as above

*R_SCQ_8yr
fre rNN152 rNN154 rNN155 rNN156 rNN157 rNN160 rNN161 rNN162 rNN163 rNN164 rNN165 rNN167
*as above

*SCQ_Q6
fre GG255 GG257 rGG258 rGG259 rGG260 rGG261 rGG262 rGG263 GG264 rGG265 rGG266 rGG267 rGG268 rGG269  rGG270 rGG271 rGG272 rGG273 GG274 GG256 GG275 GG276 GG277 GG278 GG279 GG280 GG281 GG282  GG283  GG284 GG285 GG286 GG287 GG288 GG289 GG290 GG291 GG292 GG293 GG294

*S_SCQ_Q6
fre GG255 GG257 rGG259 GG264 rGG265 GG274 GG256 GG275 GG276 GG277 GG278 GG279 GG280 GG281 GG282 GG283 GG284 GG285 GG286 GG287 GG288 GG290 GG291 GG292 GG293 GG294

*R_SCQ_Q6
fre rGG258 rGG260 rGG261 rGG262 rGG263 rGG266 rGG267 rGG268 rGG269 rGG270 rGG271 rGG273


 **********************************************
 *39. Child Behaviour Checklist (CBCL)
 
*Q6 - 18 months:

*make binary versions
foreach var of varlist EE427 EE1005 EE434 EE429 EE430 EE431 EE998 EE432 EE997 EE432_EE997 EE433 EE428 EE1006 EE900 EE1000 EE879 EE901 EE882 EE986 EE406 EE1001 EE880 EE881 EE1002 EE899 EE833 EE902 {
	gen `var'_bin=`var'
	recode `var'_bin 3=2
}

*tab1 EE433 EE898 EE884 EE885 EE434 EE429 EE431 EE432_EE997 EE428

*tab1 EE435 EE961 EE903 EE904 EE905 EE438 EE439 EE962 EE442 EE446 EE447 EE448 EE963 EE964 EE906 EE440 EE907 EE908 EE909

*Q5 36 months: make binary versions:
foreach var of varlist GG313-GG338 {
	gen `var'_bin=`var'
	recode `var'_bin 3=2
	tab `var'_bin
}

*Q5yr: make binary versions
foreach var of varlist LL301 LL302 LL303 LL304 LL305 LL306 LL307 LL308 LL309 LL310 LL311 LL312 LL313 LL314 LL315 LL316 LL317 LL318 LL319 LL320 LL321 LL322 LL323 LL324 LL325 LL504 LL505 {
	gen `var'_bin=`var'
	recode `var'_bin 3=2
	tab `var'_bin
}
*Don't use the last two as version-specific - only there for a subset


**********************
*EAS: Temperament
*fix multiple ticks:
foreach var of varlist EE416 EE417 EE418 EE419 EE420 EE421 EE422 EE423 EE424 EE425 EE426 EE877 EE878 {
	recode `var' 0=.
}
fre EE416 EE417 EE418 EE419 EE420 EE421 EE422 EE423 EE424 EE425 EE426 EE877 EE878

*make the three group domains: 
*shyness
*activity
*emotionality
*?

*or...

*these with truncreg
fre EE416 EE418 EE421 EE422 EE423 EE425 EE426 EE878
*no way normal:
fre EE417 EE419 EE420 EE424 EE877 
recode EE417 4/5=4
recode EE419 4/5=4
recode EE420 4/5=4
recode EE424 4/5=4
recode EE877 1/2=2

fre GG299 GG300 GG301 GG302 GG303 GG304 GG305 GG306 GG307 GG308 GG309 GG310 GG311 GG312

*with truncreg: 
fre GG299 GG300 GG301 GG302 GG304 GG305 GG306 GG308 GG309 GG310 GG311 GG312
*no way normal:
fre GG303 GG307 
recode GG303 4/5=4
recode GG307 4/5=4

*Age 5 
fre LL276 LL277 LL278 LL279 LL280 LL281 LL282 LL283 LL284 LL285 LL286 LL287
*these with truncreg
fre LL276 LL277 LL279 LL280 LL281 LL283 LL284 LL285 LL287
*these: no way vaguely normal, so simplify and use ologit:
fre LL278 LL282 LL286 
recode LL278 1=2
recode LL282 5=4
recode LL286 4/5=3

*rename these otherwise names too long to work in loop storing estimates
rename SCQ_SCI_8yr S_SCQ_8yr 
rename SCQ_RRB_8yr R_SCQ_8yr
rename SCQ_SCI_Q6 S_SCQ_Q6
rename SCQ_RRB_Q6 R_SCQ_Q6


*********************************************

*REGISTRATION AND TRIMMING FOR SPEED

*mi set
mi set wide

*registration: auxillary vars include SEP, smoking, lots of parental health stuff, and child developmental questions from age 5yr and 8yr
fre LL252 LL253 LL254 LL255 LL256 LL257 LL258 LL259 LL260 LL261 LL262 LL263

mi register imputed ///
FARS_ALDER VEKT LENGDE AA1123 FF15 AA1124 fathers_educ_Q1QF AA1315 AA1316 AA1317 FF341 EE583 EE584 G_66 AA806 AA807 AA869 AA870 AA878 AA879 FF146 FF386 FF389 FF392 FF395 FF398 NN322 NN323 NN304_binary NN306_binary NN308_binary NN310_binary G__230_1 G__231_1 G__232_1 G__233_1 G__234_1 G_54 G_55 mumsmokes_Q1 dadsmokes_Q1 dadsmokes_QF mumsmokes_Q6 mumsmokes_Q5y mumsmokes_Q7y dadsmokes_Q7y mumsmokes_Q8y dadsmokes_Q8y dadsmokes_FQ2 dadsmokes_Q1QF EE427 EE1005 EE434 EE429 EE430 EE431 EE998 EE432 EE997 EE433 EE428 EE1006 EE900 EE1000 EE879 EE901 EE882 EE986 EE406 EE1001 EE880 EE881 EE1002 EE899 EE833 EE902 JJ431 JJ435 JJ436_JJ437 JJ441 LL252 LL253 LL254 LL255 LL256 LL257 LL258 LL259 LL260 LL261 LL262 LL263 LL174_binary LL175_binary LL176_binary LL177_binary LL178_binary LL179_binary LL180_binary LL264 LL265 LL266 LL267 LL268 LL269 LL270 LL271 LL272 LL273 LL274 LL275 LL232 LL233 LL234 LL235 LL236 LL237 LL238 LL239 LL240 LL241 LL242 LL243 LL244 LL245 LL246 LL247 LL248 LL249 LL250 LL251 LL213 LL214 LL215 LL216 LL217 LL218 LL222 LL223 LL224 LL225 LL226 LL227 LL228 LL229 LL230 LL231 MFQ_8yr SCARED_8yr RSDBD_ADHD_8yr RSDBD_inADHD_8yr RSDBD_hyADHD_8yr  RSDBD_CD_8yr RSDBD_OD_8yr SCQ_8yr S_SCQ_8yr R_SCQ_8yr NN239 NN240 NN241 NN234 NN265 AA1548 AA1549 AA1550 AA1551 AA1552 FF251 FF252 FF253 FF254 FF255 FF256 FF257 FF258 GG514 GG515 GG516 GG517 GG518 GG519 GG520 GG521 LL363 LL364 LL365 LL366 LL367 LL368 LL369 LL370 NN325 NN326 NN327 NN328 NN329 NN330 NN331 NN332 G_52_1 G_52_2 G_52_3 G_52_4 G_52_5 G_52_6 G_52_7 G_52_8 G_52_9 G_5210 G_5211 G_5212 FF535 FF536 FF537 FF538 FF539 FF540  GG503 GG504 GG505 GG506 GG507 GG508 delayed_motor_Q8yr delayed_language_Q8y hyperactivity_Q8y concentration_Q8y behavioural_problems_Q8y emotional_diffic_Q8y other_condition_Q8y mothers_bmi_Q1 mothers_bmi_Q6 mothers_bmi_Q5y mothers_bmi_Q8y mfathers_bmi_Q1 fathers_bmi_QF fathers_bmi_Q1QF fathers_bmi_FQ2 moffspring_bmi_Q5y  moffspring_bmi_Q7y moffspring_bmi_Q8y mhopkins_Q1 fhopkins_QF mhopkins_Q6 mhopkins_Q5y mhopkins_Q8y fhopkins_F2 mADHD_Q6 fADHD_QF conners_Q5y miss1* ///
GG313-GG338 GG313_bin-GG338_bin LL276 LL277 LL278 LL279 LL280 LL281 LL282 LL283 LL284 LL285 LL286 LL287 LL301-LL505 LL301_bin-LL505_bin ///
GG231-GG236 GG299-GG312 SCQ_Q6 S_SCQ_Q6 R_SCQ_Q6 ///
EE898 EE884 EE885 EE432_EE997 ///
EE435 EE961 EE903 EE904 EE905 EE438 EE439 EE962 EE442 EE446 EE447 EE448 EE963 EE964 EE906 EE440 EE907 EE908 EE909 EE427_bin-EE902_bin ///
EE416 EE417 EE418 EE419 EE420 EE421 EE422 EE423 EE424 EE425 EE426 EE877 EE878 ///
/*RSDBD_ADHD_8yr*/ NN119 NN120 NN121 NN122 NN123 NN124 NN125 NN126 NN127 NN128 NN129 NN130 NN131 NN132 NN133 NN134 NN135 NN136 ///
/*MFQ_8yr*/ NN68 NN69 NN70 NN71 NN72 NN73 NN74 NN75 NN76 NN77 NN78 NN79 NN80 ///
/*SCARED_8yr*/ NN145 NN146 NN147 NN148 NN149 ///
/*SCQ_8yr*/ NN150 NN151 rNN152 rNN153 rNN154 rNN155 rNN156 rNN157 NN158 rNN159 rNN160 rNN161 rNN162 rNN163 rNN164 rNN165 rNN166 rNN167 NN168 NN169 NN170 NN171 NN172 NN173 NN174 NN175 NN176 NN177 NN178 NN179 NN180 NN181 NN182 NN183 NN184 NN185 NN186 NN187 NN188 NN189 ///
/*SCQ_Q6*/ GG255 GG257 rGG258 rGG259 rGG260 rGG261 rGG262 rGG263 GG264 rGG265 rGG266 rGG267 rGG268 rGG269  rGG270 rGG271 rGG272 rGG273 GG274 GG256 GG275 GG276 GG277 GG278 GG279 GG280 GG281 GG282  GG283  GG284 GG285 GG286 GG287 GG288 GG289 GG290 GG291 GG292 GG293 GG294 ///
/*simplified SEP variables */ AA1124_v2 fathers_educ_Q1QF_v2

*regular, i.e. not to be imputed, including age at questionnaire return:
mi register regular MORS_ALDER PARITET_5 SIVST_2 KJONN fid offspring_id* mid pid AGE_RETURN_MTHS_Q8AAR 

/*and the pgs for everyone*/ mi register regular *pgs* 

*also PCs and batch variables: if only running imputation on trios with full genetic data, no missing values in these to impute. 
*register the remaining dummies and the PCs
mi register regular  c_PC* m_PC* f_PC* c_batch* m_batch* f_batch*


*for efficiency, KEEP ONLY THOSE NEEDED, DROP THE REST.
*if in flong format, include: _mi_id
keep _mi_m  _mi_miss PREG_ID_2306 BARN_NR KJONN *present* fullgeneticinfo fid mid pid offspring_id*  ///
MORS_ALDER FARS_ALDER PARITET_5 SIVST_2 VEKT LENGDE AA1123 FF15 AA1124 fathers_educ_Q1QF AA1315 AA1316 AA1317 FF341 EE583 EE584 G_66 AA806 AA807 AA869 AA870 AA878 AA879 FF146 FF386 FF389 FF392 FF395 FF398 NN322 NN323 NN304_binary NN306_binary NN308_binary NN310_binary G__230_1 G__231_1 G__232_1 G__233_1 G__234_1 G_54 G_55 mumsmokes_Q1 dadsmokes_Q1 dadsmokes_QF mumsmokes_Q6 mumsmokes_Q5y mumsmokes_Q7y dadsmokes_Q7y mumsmokes_Q8y dadsmokes_Q8y dadsmokes_FQ2 dadsmokes_Q1QF EE427 EE1005 EE434 EE429 EE430 EE431 EE998 EE432 EE997 EE433 EE428 EE1006 EE900 EE1000 EE879 EE901 EE882 EE986 EE406 EE1001 EE880 EE881 EE1002 EE899 EE833 EE902 JJ431 JJ435 JJ436_JJ437 JJ441 LL174_binary LL175_binary LL176_binary LL177_binary LL178_binary LL179_binary LL180_binary LL264 LL265 LL266 LL267 LL268 LL269 LL270 LL271 LL272 LL273 LL274 LL275 LL232 LL233 LL234 LL235 LL236 LL237 LL238 LL239 LL240 LL241 LL242 LL243 LL244 LL245 LL246 LL247 LL248 LL249 LL250 LL251 LL213 LL214 LL215 LL216 LL217 LL218 LL222 LL223 LL224 LL225 LL226 LL227 LL228 LL229 LL230 LL231 LL252 LL253 LL254 LL255 LL256 LL257 LL258 LL259 LL260 LL261 LL262 LL263 MFQ_8yr SCARED_8yr RSDBD_ADHD_8yr RSDBD_inADHD_8yr RSDBD_hyADHD_8yr  RSDBD_CD_8yr RSDBD_OD_8yr SCQ_8yr S_SCQ_8yr R_SCQ_8yr NN239 NN240 NN241 NN234 NN265 AA1548 AA1549 AA1550 AA1551 AA1552 FF251 FF252 FF253 FF254 FF255 FF256 FF257 FF258 GG514 GG515 GG516 GG517 GG518 GG519 GG520 GG521 LL363 LL364 LL365 LL366 LL367 LL368 LL369 LL370 NN325 NN326 NN327 NN328 NN329 NN330 NN331 NN332 G_52_1 G_52_2 G_52_3 G_52_4 G_52_5 G_52_6 G_52_7 G_52_8 G_52_9 G_5210 G_5211 G_5212 FF535 FF536 FF537 FF538 FF539 FF540  GG503 GG504 GG505 GG506 GG507 GG508 delayed_motor_Q8yr delayed_language_Q8y hyperactivity_Q8y concentration_Q8y behavioural_problems_Q8y emotional_diffic_Q8y other_condition_Q8y mothers_bmi_Q1 mothers_bmi_Q6 mothers_bmi_Q5y mothers_bmi_Q8y mfathers_bmi_Q1 fathers_bmi_QF fathers_bmi_Q1QF fathers_bmi_FQ2 moffspring_bmi_Q5y  moffspring_bmi_Q7y  moffspring_bmi_Q8y mhopkins_Q1 fhopkins_QF mhopkins_Q6 mhopkins_Q5y mhopkins_Q8y fhopkins_F2 mADHD_Q6 fADHD_QF conners_Q5y miss1* ///
GG313-GG338 GG313_bin-GG338_bin LL276 LL277 LL278 LL279 LL280 LL281 LL282 LL283 LL284 LL285 LL286 LL287 LL301-LL505 LL301_bin-LL505_bin SCQ_Q6 S_SCQ_Q6 R_SCQ_Q6 ///
EE898 EE884 EE885 EE432_EE997 ///
EE435 EE961 EE903 EE904 EE905 EE438 EE439 EE962 EE442 EE446 EE447 EE448 EE963 EE964 EE906 EE440 EE907 EE908 EE909 EE427_bin-EE902_bin ///
EE416 EE417 EE418 EE419 EE420 EE421 EE422 EE423 EE424 EE425 EE426 EE877 EE878 ///
GG299 GG300 GG301 GG302 GG303 GG304 GG305 GG306 GG307 GG308 GG309 GG310 GG311 GG312 ///
/*RSDBD_ADHD_8yr*/ NN119 NN120 NN121 NN122 NN123 NN124 NN125 NN126 NN127 NN128 NN129 NN130 NN131 NN132 NN133 NN134 NN135 NN136 ///
/*MFQ_8yr*/ NN68 NN69 NN70 NN71 NN72 NN73 NN74 NN75 NN76 NN77 NN78 NN79 NN80 ///
/*SCARED_8yr*/ NN145 NN146 NN147 NN148 NN149 ///
/*SCQ_8yr*/ NN150 NN151 rNN152 rNN153 rNN154 rNN155 rNN156 rNN157 NN158 rNN159 rNN160 rNN161 rNN162 rNN163 rNN164 rNN165 rNN166 rNN167 NN168 NN169 NN170 NN171 NN172 NN173 NN174 NN175 NN176 NN177 NN178 NN179 NN180 NN181 NN182 NN183 NN184 NN185 NN186 NN187 NN188 NN189 ///
/*SCQ_Q6*/ GG255 GG257 rGG258 rGG259 rGG260 rGG261 rGG262 rGG263 GG264 rGG265 rGG266 rGG267 rGG268 rGG269  rGG270 rGG271 rGG272 rGG273 GG274 GG256 GG275 GG276 GG277 GG278 GG279 GG280 GG281 GG282  GG283  GG284 GG285 GG286 GG287 GG288 GG289 GG290 GG291 GG292 GG293 GG294 ///
/*simplified SEP variables */ AA1124_v2 fathers_educ_Q1QF_v2 ///
/*and the pgs for everyone*/ *pgs* ///
/*plus genetic exclusion flags, and covars for the genetic analysis!*/ AGE_RETURN_MTHS_Q8AAR *PC* *batch*  

*binary variable causing problems:
/* fre AA806 AA807 AA869 AA870 AA878 AA879 FF146 FF386 FF389 FF392 FF395 FF398 NN322 NN323 NN304_binary NN306_binary NN308_binary NN310_binary G__230_1 G__231_1 G__232_1 G__233_1 G__234_1 G_54 G_55 JJ431 JJ435 JJ436_JJ437 JJ441 
fre LL226 LL227 LL228 LL229 LL230 LL174_binary LL175_binary LL176_binary LL177_binary LL178_binary LL179_binary LL180_binary LL264 LL265 LL266 LL267 LL268 LL269 LL270 LL271 LL272 LL273 LL274 LL275 LL232 LL233 LL234 LL235 LL236 LL237 LL238 LL239 LL240 LL241 LL242 LL243 LL244 LL245 LL246 LL247 LL248 LL249 LL250 LL251 */

*in STATA, binary variables in imputation need to be coded 0/1 not 1/2.
foreach var of varlist AA806 AA807 AA869 AA870 AA878 AA879 FF146 FF386 FF389 FF392 FF395 FF398 NN322 NN323 NN304_binary NN306_binary NN308_binary NN310_binary G__230_1 G__231_1 G__232_1 G__233_1 G__234_1 G_54 G_55 JJ431 JJ435 JJ436_JJ437 JJ441 LL226 LL227 LL228 LL229 LL230 LL174_binary LL175_binary LL176_binary LL177_binary LL178_binary LL179_binary LL180_binary LL264 LL265 LL266 LL267 LL268 LL269 LL270 LL271 LL272 LL273 LL274 LL275 LL232 LL233 LL234 LL235 LL236 LL237 LL238 LL239 LL240 LL241 LL242 LL243 LL244 LL245 LL246 LL247 LL248 LL249 LL250 LL251 SIVST_2 EE898 EE884 EE885 EE434 EE429 EE431 EE432_EE997 EE433 EE428 EE427_bin-EE902_bin GG313_bin-GG338_bin LL301_bin-LL505_bin ///
NN150 NN151 rNN152 rNN153 rNN154 rNN155 rNN156 rNN157 NN158 rNN159 rNN160 rNN161 rNN162 rNN163 rNN164 rNN165 rNN166 rNN167 NN168 NN169 NN170 NN171 NN172 NN173 NN174 NN175 NN176 NN177 NN178 NN179 NN180 NN181 NN182 NN183 NN184 NN185 NN186 NN187 NN188 NN189 ///
GG255 GG257 rGG258 rGG259 rGG260 rGG261 rGG262 rGG263 GG264 rGG265 rGG266 rGG267 rGG268 rGG269  rGG270 rGG271 rGG272 rGG273 GG274 GG256 GG275 GG276 GG277 GG278 GG279 GG280 GG281 GG282 GG283 GG284 GG285 GG286 GG287 GG288 GG289 GG290 GG291 GG292 GG293 GG294 {
replace `var'=`var'-1
fre `var'
}

**GO BACK AND CHANGE LABELLING?
*#BETTER - POST-IMPUTATION, ADD 1 TO ALL OF THESE IN MI PASSIVE.

*****************************************************
*simplify further categories of SEP vars:

*financial strain items:
fre AA1317 EE583 EE584 
*leave EE583  since it's so similar to EE584 but less helpful
*simplify EE584
recode EE584 4=3
label define EE584_new 1"never" 2"yes, but infrequently" 3"merged: sometimes and often"
label values EE584 EE584_new
fre EE584

*income variables:
fre AA1315 AA1316 FF341 G_66
*change the dont know in mother's report of partner's income to missing:
recode AA1316  8=.
*hist AA1315 AA1316 FF341 G_66
*At a pinch, can run these as truncreg, limits of 1 and 7. Not happy about it.

*like bmi, smoking, and education, make a variable for father's income beginning with self-report and subbing in mum's report
*check labelling is the same:
fre AA1316 FF341
gen fathers_inc_grp=FF341
replace fathers_inc_grp=AA1316 if fathers_inc_grp==.
label values fathers_inc_grp FF341
fre fathers_inc_grp
mi register imputed fathers_inc_grp


*for parental educ, FLIP THEM to have larger baseline groups.
*post-imputation, for Table 1, can flip back by subtracting from 7 again. Or just leave it!
foreach var in AA1124_v2 fathers_educ_Q1QF_v2  {
	fre `var'
	gen `var'_flipped=7-`var'
	tab `var' `var'_flipped
	mi register imputed `var'_flipped
}

*mother's smoking: collapse the 3rd and fourth categories?
fre mumsmokes_Q1
gen mumsmokes_Q1_v2=mumsmokes_Q1
recode mumsmokes_Q1_v2 3=2
label define mumsmokes_Q1_v2 0"never" 1"pre-pregnancy" 2"during"
label values mumsmokes_Q1_v2 mumsmokes_Q1_v2
fre mumsmokes_Q1_v2
mi register imputed mumsmokes_Q1_v2

*school enjoyment: collpase and flip for bigger baseline categ:
gen NN234_new=NN234
label define NN234_new 1"very well" 2"well" 3"neither/poor/very poor"
label values NN234_new NN234_new
recode NN234_new 1/3=1 4=2 5=3
fre NN234_new
replace NN234_new=4-NN234_new
fre NN234_new
mi register imputed NN234_new

*sefl-reading: just flip. later, can collpase if needed
gen NN265_new=NN265
label define NN265_new 1"Books with chapters (almost only text)" 2"Simple stories, pictures & text on every page" 3"Picture book (few words)" 4"Does not like to read him-/herself"
label values NN265_new NN265_new
replace NN265_new=5-NN265_new
fre NN265_new
mi register imputed NN265_new

save scratch/pre_imputation.dta, replace

BREAK

**********************************************************************************************************************************

*ENTRY POINT: 

*STRATEGY 2: IMPUTE 2 x ADHD SUBSCALES, THEN ADD THE TWO POST-IMPUTATION FOR THE FULL SCALE.

*IMPUTATION MODEL FOR DEPRESSION, ANXIETY AND ADHD (SINGLE-ITEMS) 

sysdir set PLUS "N:\data\durable\people\Mandy\Statafiles\ado\plus"
cd "N:\data\durable\projects\mh_bmi\"
set maxvar 20000

use scratch/pre_imputation.dta, clear
count
*26730

*also check:
lookfor fid
summ mothers_bmi_Q1

**************************************************

*New ADHD prep:
*using these totals now: rename and register:
foreach varstem in ADHD_8yr inADHD_8yr hyADHD_8yr {
mi unregister RSDBD_`varstem'
rename RSDBD_`varstem' `varstem'
mi register imputed `varstem'
summ `varstem'
}
*min and max: 0 for all, 54, 27, 27
*so made correctly to start from 0

**UPDATE: ONLY PUT THE SUBSCALES, NOT THE GFULL ADHD SCALE, IN THE MODEL
**RECREATE THE 
mi unregister ADHD_8yr 
drop ADHD_8yr

**************************************************

*************FOR SPEED: prune stuff we definitely don't need(move upstream later):
*rename this one first:
mi unregister EE584
rename EE584 save_EE584
mi register imputed save_EE584 

drop *miss1*
*drop *VEKT* *LENGDE*
drop *AA1123 *AA1552 /**AA1315*/ *AA1316 *AA1317 *AA806 *AA807 *AA869 *AA870 *AA878 *AA879 *AA1548 *AA1549 *AA1550 *AA1551
drop G_* FF* EE* GG* rGG* LL* JJ* NN* rNN* 
drop *mumsmokes_Q6 *mumsmokes_Q5y *mumsmokes_Q7y *dadsmokes_Q7y *mumsmokes_Q8y *dadsmokes_Q8y *dadsmokes_FQ2 
drop *RSDBD_CD_8yr *RSDBD_OD_8yr 
drop *SCQ* 
drop /* *conners_Q5y */ *delayed_motor_Q8yr *delayed_language_Q8y *hyperactivity_Q8y *concentration_Q8y *behavioural_problems_Q8y *emotional_diffic_Q8y *other_condition_Q8y 
drop *mhopkins_Q6 *mhopkins_Q5y *mhopkins_Q8y *fhopkins_F2  /**fathers_inc_grp*/
*drop *moffspring_bmi_Q5y *moffspring_bmi_Q7y *mothers_bmi_Q5y *mothers_bmi_Q6 *mothers_bmi_Q8y  *fathers_bmi_FQ2 *mfathers_bmi_Q1 *fathers_bmi_QF
*use the _v2_flipped versions so don't need the following:
drop *AA1124 *fathers_educ_Q1QF *AA1124_v2 *fathers_educ_Q1QF_v2

*rename that as too long apparently
mi unregister fathers_educ_Q1QF_v2_flipped
rename fathers_educ_Q1QF_v2_flipped f_educ_Q1QF_v2_flipped
mi register imputed f_educ_Q1QF_v2_flipped

*178 vars left

************************************************************************************************

*ACTUAL IMPUTATION

capture erase impstats.dta
mi impute chained ///
(ologit) /*SEP: parental education*/ AA1124_v2_flipped f_educ_Q1QF_v2_flipped save_EE584 ///
(truncreg, ll(1) ul(7)) /*parental income, grouped*/ AA1315 fathers_inc_grp  ///
(ologit) /*parental smoking*/ mumsmokes_Q1_v2 dadsmokes_Q1QF  ///
(pmm, knn(10)) /*parental anxiety/dep */ mhopkins_Q1 ///
(pmm, knn(10)) /*parental anxiety/dep */ fhopkins_QF ///
(truncreg, ll(8.1) ul(50.2)) mothers_bmi_Q1 ///
(truncreg, ll(13.6) ul(49.6)) mothers_bmi_Q6 ///
(truncreg, ll(13.1) ul(48.3)) fathers_bmi_Q1QF  ///
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
= /*child sex, then parity, maternal age, and marital status from birth registry file*/ KJONN PARITET_5 MORS_ALDER SIVST_2 ///
/*bmi and other polygenic scores*/  *z_childbmi_pgs *z_adultbmi_pgs *z_ea3_pgs *z_depression_pgs *z_adhd_pgs ///
/* PCs and batch dummies */ c_PC* m_PC* f_PC* c_batch* m_batch* f_batch* ///
, add(100) rseed(100) dots savetrace(impstats.dta) augment 

*<10 minutes for one!
*100 set going noon Friday

save "scratch/scalesonly_dep_anxiety_adhd_2July.dta", replace

count
*26370

BREAK

********************************************************************

*ENTRY FOR TRANSFORMATIONS:

cd "N:\data\durable\projects\mh_bmi\"
set maxvar 20000
sysdir
sysdir set PLUS "N:\data\durable\people\Mandy\Statafiles\ado\plus"

use "scratch/scalesonly_dep_anxiety_adhd_2July.dta", clear

*TRANSFORMATIONS:

*RESCALE BMI as per 5kg/m2 to have less stupidly small coefficients:
*not super-varying so can use mi passive:
mi passive: gen rescaled_moffspring_bmi_Q8y=moffspring_bmi_Q8/5
mi passive: gen rescaled_mothers_bmi_Q1=mothers_bmi_Q1/5
mi passive: gen rescaled_fathers_bmi_Q1QF=fathers_bmi_Q1QF/5

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

*test model:

mi estimate, post cmdok: ivregress 2sls z_MFQ_8yr (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_PC* c_batch*, vce(cluster fid)
mi estimate, post cmdok: ivregress 2sls z_ADHD_8yr (rescaled_moffspring_bmi_Q8=c_z_adultbmi_pgs) KJONN c_PC* c_batch*, vce(cluster fid)

save "scratch/scalesonly_dep_anxiety_adhd_2July.dta", replace

BREAK
