
*PRE-IMPUTATION

clear all
cd "N:\durable\projects\mh_bmi\fullsample"
set maxvar 20000
sysdir
sysdir set PLUS "N:\durable\people\Mandy\Statafiles\ado\plus"
use "data\merged_geno_pheno.dta", clear

*****************************************************************************

*FINAL PREP TO SUPPORT IMPUTATION 

*TRIM OFF PEOPLE WHO CONTRIBUTE NOTHING AND STOP IT FROM RUNNING
*restrict again to those present in birth registry file:
drop if present_MBRN!=1
*113,691
fre present*
*is there anyone present in the birth registry file but not at any of the actual questionnaires used in imputation?
gen nQ_present=0
foreach q in Q1 QF Q5 Q6 Q5y Q7y Q8y {
	replace nQ_present=nQ_present+1 if present_`q'==1
}
fre nQ_present
*OK, so there's a chunk of people who never were included in a single relevant questionnaire.
*these are the people that weren't present at the very first one:
tab nQ_present present_Q1, mi
*these people basically contribute no information, and are most likely holding up imputation. 
*DROP THEM, but at a later stage as less vital than the genetic info.

*for flowchart: 
*PRESENT IN BIRTH REGISTRY FILE, consent not withdrawn:
count
*113,691
*with full genetic information:
keep if complete_trio==1
count
*42,207
*with nonmissing sex? *not an issue in the sample of full trios
*even if it were, could use c_sex from genetic file
tab KJONN c_sex if complete==1, mi
*present at >=1 questionnaire:
keep if nQ_present>0 
count
*40,949

******************************************************************
*IMPUTATION DOES NOT WORK IF YOU KEEP EVERYONE
*SO RESTRICTING TO FULL TRIOS
****************************************************************************************************

*FINAL PREP OF PHENOTYPES:

*impute with untransformed variables since you will need them for table 1. then standardize in mi passive

*CHECK FORMAT OF INDIVIDUAL ITEMS FOR AGE 8 OUTCOMES:

*ADHD_8yr
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

*AUXILIARY VARIABLES:

*financial strain items:
fre AA1317 EE583 EE584 
*leave EE583 since it's so similar to EE584 but less helpful

*income variables:
fre AA1315 AA1316 FF341 G_66
*hist AA1315 AA1316 FF341 G_66

*39. Child Behaviour Checklist (CBCL)
*18 months, 36 months, Q5yr: use binary versions (made in phenotype prep file)

*EAS: Temperament
fre EE416 EE417 EE418 EE419 EE420 EE421 EE422 EE423 EE424 EE425 EE426 EE877 EE878
*these with truncreg?
fre EE416 EE418 EE421 EE422 EE423 EE425 EE426 EE878
*no way normal:
fre EE417 EE419 EE420 EE424 EE877 
recode EE417 4/5=4
recode EE419 4/5=4
recode EE420 4/5=4
recode EE424 4/5=4
recode EE877 1/2=2

fre GG299 GG300 GG301 GG302 GG303 GG304 GG305 GG306 GG307 GG308 GG309 GG310 GG311 GG312
*with truncreg?
fre GG299 GG300 GG301 GG302 GG304 GG305 GG306 GG308 GG309 GG310 GG311 GG312
*no way normal:
fre GG303 GG307 
recode GG303 4/5=4
recode GG307 4/5=4

fre LL276 LL277 LL278 LL279 LL280 LL281 LL282 LL283 LL284 LL285 LL286 LL287
*these with truncreg?
fre LL276 LL277 LL279 LL280 LL281 LL283 LL284 LL285 LL287
*these: no way vaguely normal, so simplify and use ologit:
fre LL278 LL282 LL286 
recode LL278 1=2
recode LL282 5=4
recode LL286 4/5=3

****************************************************************************************

*check mutually-adjusted batch group and genotyping chip effects

*INFO ON THE SELECTION CRITERIA OF ALL THE BATCHES IS HERE:
*https://github.com/folkehelseinstituttet/mobagen/wiki/Projects-that-have-contributed-to-MoBa-Genetics#selection-criterias-for-samples-in-harvest

*z_MFQ_8yr z_SCARED_8yr z_ADHD_8yr
****************************************************************************
*NB: for children (not mums and dads), none were genotyped using:
*genotyping_chip==Illumina_HumanOmniExpress-24v1.0
*so this dummy has 0 people in the nomissing group, and is omitted from imputation and analysis models. 
*depending on how the dummy has been coded, this is either the 4th of 5th dummy
*largest (baseline) group to omit for children is:
*genotyping_chip==Illumina_GlobalScreeningArray_MD_v.3.0

*so here, chip dummies to include for children are: c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 
****************************************************************************
*try center and chip together:
foreach var in c_z_adultbmi_pgs c_z_childbmi_pgs  {
regress `var' c_genotyping_center2 c_genotyping_center3 c_genotyping_chip1 c_genotyping_chip2 c_genotyping_chip3 c_genotyping_chip5 c_PC*
}

foreach var in c_z_adultbmi_pgs c_z_childbmi_pgs MFQ_8yr SCARED_8yr ADHD_8yr  {
regress `var' c_genotyping_chip1-c_genotyping_chip6 c_PC*
}
*nothing with chip for BMI, but definitely chip effects for the outcomes. is this just because the chips proxy batch?

foreach var in c_z_adultbmi_pgs c_z_childbmi_pgs MFQ_8yr SCARED_8yr ADHD_8yr  {
regress `var' c_batch* c_PC*
}

**********************************************************************************************

*REGISTRATION AND TRIMMING FOR SPEED

*mi set
mi set wide

*registration: auxillary vars include SEP, smoking, lots of parental health stuff, and child developmental questions from age 5yr and 8yr
fre LL252 LL253 LL254 LL255 LL256 LL257 LL258 LL259 LL260 LL261 LL262 LL263

mi register imputed ///
FARS_ALDER VEKT LENGDE AA1123 FF15 AA1315 AA1316 AA1317 FF341 EE583 EE584 G_66 AA806 AA807 AA869 AA870 AA878 AA879 FF146 FF386 FF389 FF392 FF395 FF398 NN322 NN323 NN304_binary NN306_binary NN308_binary NN310_binary G__230_1 G__231_1 G__232_1 G__233_1 G__234_1 G_54 G_55 mumsmokes_Q1 dadsmokes_Q1 dadsmokes_QF mumsmokes_Q6 mumsmokes_Q5y mumsmokes_Q7y dadsmokes_Q7y mumsmokes_Q8y dadsmokes_Q8y dadsmokes_QF2 dadsmokes_Q1QF EE427 EE1005 EE434 EE429 EE430 EE431 EE998 EE432 EE997 EE433 EE428 EE1006 EE900 EE1000 EE879 EE901 EE882 EE986 EE406 EE1001 EE880 EE881 EE1002 EE899 EE833 EE902 JJ431 JJ435 JJ436_JJ437 JJ441 LL252 LL253 LL254 LL255 LL256 LL257 LL258 LL259 LL260 LL261 LL262 LL263 LL174_binary LL175_binary LL176_binary LL177_binary LL178_binary LL179_binary LL180_binary LL264 LL265 LL266 LL267 LL268 LL269 LL270 LL271 LL272 LL273 LL274 LL275 LL232 LL233 LL234 LL235 LL236 LL237 LL238 LL239 LL240 LL241 LL242 LL243 LL244 LL245 LL246 LL247 LL248 LL249 LL250 LL251 LL213 LL214 LL215 LL216 LL217 LL218 LL222 LL223 LL224 LL225 LL226 LL227 LL228 LL229 LL230 LL231 MFQ_8yr SCARED_8yr ADHD_8yr inADHD_8yr hyADHD_8yr /* CD_8yr OD_8yr */ SCQ_8yr S_SCQ_8yr R_SCQ_8yr NN239 NN240 NN241 NN234_new NN265_new AA1548 AA1549 AA1550 AA1551 AA1552 FF251 FF252 FF253 FF254 FF255 FF256 FF257 FF258 GG514 GG515 GG516 GG517 GG518 GG519 GG520 GG521 LL363 LL364 LL365 LL366 LL367 LL368 LL369 LL370 NN325 NN326 NN327 NN328 NN329 NN330 NN331 NN332 G_52_1 G_52_2 G_52_3 G_52_4 G_52_5 G_52_6 G_52_7 G_52_8 G_52_9 G_5210 G_5211 G_5212 FF535 FF536 FF537 FF538 FF539 FF540  GG503 GG504 GG505 GG506 GG507 GG508 delayed_motor_Q8yr delayed_language_Q8y hyperactivity_Q8y concentration_Q8y behavioural_problems_Q8y emotional_diffic_Q8y other_condition_Q8y mothers_bmi_Q1 mothers_bmi_Q6 mothers_bmi_Q5y mothers_bmi_Q8y mfathers_bmi_Q1 fathers_bmi_QF fathers_bmi_Q1QF fathers_bmi_FQ2 moffspring_bmi_Q5y  moffspring_bmi_Q7y moffspring_bmi_Q8y mhopkins_Q1 fhopkins_QF mhopkins_Q6 mhopkins_Q5y mhopkins_Q8y fhopkins_QF2 mADHD_Q6 fADHD_QF conners_Q5y  ///
GG313-GG338 GG313_bin-GG338_bin LL276 LL277 LL278 LL279 LL280 LL281 LL282 LL283 LL284 LL285 LL286 LL287 LL301-LL505 LL301_bin-LL505_bin ///
GG231-GG236 GG299-GG312 SCQ_Q6 S_SCQ_Q6 R_SCQ_Q6 ///
EE898 EE884 EE885 EE432_EE997 ///
EE435 EE961 EE903 EE904 EE905 EE438 EE439 EE962 EE442 EE446 EE447 EE448 EE963 EE964 EE906 EE440 EE907 EE908 EE909 EE427_bin-EE902_bin ///
EE416 EE417 EE418 EE419 EE420 EE421 EE422 EE423 EE424 EE425 EE426 EE877 EE878 ///
/*ADHD_8yr*/ NN119 NN120 NN121 NN122 NN123 NN124 NN125 NN126 NN127 NN128 NN129 NN130 NN131 NN132 NN133 NN134 NN135 NN136 ///
/*MFQ_8yr*/ NN68 NN69 NN70 NN71 NN72 NN73 NN74 NN75 NN76 NN77 NN78 NN79 NN80 ///
/*SCARED_8yr*/ NN145 NN146 NN147 NN148 NN149 ///
/*SCQ_8yr*/ NN150 NN151 rNN152 rNN153 rNN154 rNN155 rNN156 rNN157 NN158 rNN159 rNN160 rNN161 rNN162 rNN163 rNN164 rNN165 rNN166 rNN167 NN168 NN169 NN170 NN171 NN172 NN173 NN174 NN175 NN176 NN177 NN178 NN179 NN180 NN181 NN182 NN183 NN184 NN185 NN186 NN187 NN188 NN189 ///
/*SCQ_Q6*/ GG255 GG257 rGG258 rGG259 rGG260 rGG261 rGG262 rGG263 GG264 rGG265 rGG266 rGG267 rGG268 rGG269  rGG270 rGG271 rGG272 rGG273 GG274 GG256 GG275 GG276 GG277 GG278 GG279 GG280 GG281 GG282  GG283  GG284 GG285 GG286 GG287 GG288 GG289 GG290 GG291 GG292 GG293 GG294 ///
/*simplified SEP variables */ AA1124_v2_flipped f_educ_Q1QF_v2_flipped fathers_inc_grp

*regular, i.e. not to be imputed, including age at questionnaire return:
mi register regular MORS_ALDER PARITET_5 SIVST_2 KJONN fid AGE_RETURN_MTHS_Q8AAR c_yob m_yob f_yob

/*and the pgs for everyone */ mi register regular *pgs* 
*covars for the genetic analysis: if only running imputation on trios with full genetic data, no missing values in these to impute. 
mi register regular ///
/*genetic id vars and PCs*/ fid c_iid c_fid m_iid m_fid f_iid f_fid *PC* ///
/*genotyping center (largest group=3 for everyone)*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1  f_genotyping_center2 ///
/*genotyping chip: omitted (largest) group is 2, 1, 1*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6
*genotyping batch, omitting largest group for each person and the two groups who for children has 0 people in?
*leave for now 


*for efficiency, KEEP ONLY THOSE NEEDED, DROP THE REST.
*if in flong format, include: _mi_id
keep _mi_m  _mi_miss PREG_ID_2306 BARN_NR KJONN *present* c_yob m_yob f_yob ///
MORS_ALDER FARS_ALDER PARITET_5 SIVST_2 VEKT LENGDE AA1123 FF15 AA1315 AA1316 AA1317 FF341 EE583 EE584 G_66 AA806 AA807 AA869 AA870 AA878 AA879 FF146 FF386 FF389 FF392 FF395 FF398 NN322 NN323 NN304_binary NN306_binary NN308_binary NN310_binary G__230_1 G__231_1 G__232_1 G__233_1 G__234_1 G_54 G_55 mumsmokes_Q1 dadsmokes_Q1 dadsmokes_QF mumsmokes_Q6 mumsmokes_Q5y mumsmokes_Q7y dadsmokes_Q7y mumsmokes_Q8y dadsmokes_Q8y dadsmokes_QF2 dadsmokes_Q1QF EE427 EE1005 EE434 EE429 EE430 EE431 EE998 EE432 EE997 EE433 EE428 EE1006 EE900 EE1000 EE879 EE901 EE882 EE986 EE406 EE1001 EE880 EE881 EE1002 EE899 EE833 EE902 JJ431 JJ435 JJ436_JJ437 JJ441 LL174_binary LL175_binary LL176_binary LL177_binary LL178_binary LL179_binary LL180_binary LL264 LL265 LL266 LL267 LL268 LL269 LL270 LL271 LL272 LL273 LL274 LL275 LL232 LL233 LL234 LL235 LL236 LL237 LL238 LL239 LL240 LL241 LL242 LL243 LL244 LL245 LL246 LL247 LL248 LL249 LL250 LL251 LL213 LL214 LL215 LL216 LL217 LL218 LL222 LL223 LL224 LL225 LL226 LL227 LL228 LL229 LL230 LL231 LL252 LL253 LL254 LL255 LL256 LL257 LL258 LL259 LL260 LL261 LL262 LL263 MFQ_8yr SCARED_8yr ADHD_8yr inADHD_8yr hyADHD_8yr  /*CD_8yr OD_8yr*/ SCQ_8yr S_SCQ_8yr R_SCQ_8yr NN239 NN240 NN241 NN234_new NN265_new AA1548 AA1549 AA1550 AA1551 AA1552 FF251 FF252 FF253 FF254 FF255 FF256 FF257 FF258 GG514 GG515 GG516 GG517 GG518 GG519 GG520 GG521 LL363 LL364 LL365 LL366 LL367 LL368 LL369 LL370 NN325 NN326 NN327 NN328 NN329 NN330 NN331 NN332 G_52_1 G_52_2 G_52_3 G_52_4 G_52_5 G_52_6 G_52_7 G_52_8 G_52_9 G_5210 G_5211 G_5212 FF535 FF536 FF537 FF538 FF539 FF540  GG503 GG504 GG505 GG506 GG507 GG508 delayed_motor_Q8yr delayed_language_Q8y hyperactivity_Q8y concentration_Q8y behavioural_problems_Q8y emotional_diffic_Q8y other_condition_Q8y mothers_bmi_Q1 mothers_bmi_Q6 mothers_bmi_Q5y mothers_bmi_Q8y mfathers_bmi_Q1 fathers_bmi_QF fathers_bmi_Q1QF fathers_bmi_FQ2 moffspring_bmi_Q5y  moffspring_bmi_Q7y  moffspring_bmi_Q8y mhopkins_Q1 fhopkins_QF mhopkins_Q6 mhopkins_Q5y mhopkins_Q8y fhopkins_QF2 mADHD_Q6 fADHD_QF conners_Q5y ///
GG313-GG338 GG313_bin-GG338_bin LL276 LL277 LL278 LL279 LL280 LL281 LL282 LL283 LL284 LL285 LL286 LL287 LL301-LL505 LL301_bin-LL505_bin SCQ_Q6 S_SCQ_Q6 R_SCQ_Q6 ///
EE898 EE884 EE885 EE432_EE997 ///
EE435 EE961 EE903 EE904 EE905 EE438 EE439 EE962 EE442 EE446 EE447 EE448 EE963 EE964 EE906 EE440 EE907 EE908 EE909 EE427_bin-EE902_bin ///
EE416 EE417 EE418 EE419 EE420 EE421 EE422 EE423 EE424 EE425 EE426 EE877 EE878 ///
GG299 GG300 GG301 GG302 GG303 GG304 GG305 GG306 GG307 GG308 GG309 GG310 GG311 GG312 ///
/*ADHD_8yr*/ NN119 NN120 NN121 NN122 NN123 NN124 NN125 NN126 NN127 NN128 NN129 NN130 NN131 NN132 NN133 NN134 NN135 NN136 ///
/*MFQ_8yr*/ NN68 NN69 NN70 NN71 NN72 NN73 NN74 NN75 NN76 NN77 NN78 NN79 NN80 ///
/*SCARED_8yr*/ NN145 NN146 NN147 NN148 NN149 ///
/*SCQ_8yr*/ NN150 NN151 rNN152 rNN153 rNN154 rNN155 rNN156 rNN157 NN158 rNN159 rNN160 rNN161 rNN162 rNN163 rNN164 rNN165 rNN166 rNN167 NN168 NN169 NN170 NN171 NN172 NN173 NN174 NN175 NN176 NN177 NN178 NN179 NN180 NN181 NN182 NN183 NN184 NN185 NN186 NN187 NN188 NN189 ///
/*SCQ_Q6*/ GG255 GG257 rGG258 rGG259 rGG260 rGG261 rGG262 rGG263 GG264 rGG265 rGG266 rGG267 rGG268 rGG269  rGG270 rGG271 rGG272 rGG273 GG274 GG256 GG275 GG276 GG277 GG278 GG279 GG280 GG281 GG282  GG283  GG284 GG285 GG286 GG287 GG288 GG289 GG290 GG291 GG292 GG293 GG294 ///
/*simplified SEP variables */ AA1124_v2_flipped f_educ_Q1QF_v2_flipped fathers_inc_grp ///
/*and the pgs for everyone */ *pgs* ///
/*admin vars*/ AGE_RETURN_MTHS_Q8AAR ///
/*plus covars for the genetic analysis!*/ ///
/*genetic id vars and PCs*/ fid c_iid c_fid m_iid m_fid f_iid f_fid *PC* ///
/*genotyping center (largest group=3 for everyone)*/ c_genotyping_center1  c_genotyping_center2 m_genotyping_center1  m_genotyping_center2 f_genotyping_center1  f_genotyping_center2 ///
/*genotyping chip: omitted (largest) group is 2, 1, 1*/ c_genotyping_chip1 c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip6 m_genotyping_chip2  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip2 f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6
*genotyping batch, omitting largest group for each person and the two groups who for children has 0 people in?
*leave for now 

*binary variable causing problems: in STATA, binary variables in imputation need to be coded 0/1 not 1/2.
*fix all these:
foreach var of varlist AA806 AA807 AA869 AA870 AA878 AA879 FF146 FF386 FF389 FF392 FF395 FF398 NN322 NN323 NN304_binary NN306_binary NN308_binary NN310_binary G__230_1 G__231_1 G__232_1 G__233_1 G__234_1 G_54 G_55 JJ431 JJ435 JJ436_JJ437 JJ441 LL226 LL227 LL228 LL229 LL230 LL174_binary LL175_binary LL176_binary LL177_binary LL178_binary LL179_binary LL180_binary LL264 LL265 LL266 LL267 LL268 LL269 LL270 LL271 LL272 LL273 LL274 LL275 LL232 LL233 LL234 LL235 LL236 LL237 LL238 LL239 LL240 LL241 LL242 LL243 LL244 LL245 LL246 LL247 LL248 LL249 LL250 LL251 SIVST_2 EE898 EE884 EE885 EE434 EE429 EE431 EE432_EE997 EE433 EE428 EE427_bin-EE902_bin GG313_bin-GG338_bin LL301_bin-LL505_bin ///
NN150 NN151 rNN152 rNN153 rNN154 rNN155 rNN156 rNN157 NN158 rNN159 rNN160 rNN161 rNN162 rNN163 rNN164 rNN165 rNN166 rNN167 NN168 NN169 NN170 NN171 NN172 NN173 NN174 NN175 NN176 NN177 NN178 NN179 NN180 NN181 NN182 NN183 NN184 NN185 NN186 NN187 NN188 NN189 ///
GG255 GG257 rGG258 rGG259 rGG260 rGG261 rGG262 rGG263 GG264 rGG265 rGG266 rGG267 rGG268 rGG269  rGG270 rGG271 rGG272 rGG273 GG274 GG256 GG275 GG276 GG277 GG278 GG279 GG280 GG281 GG282 GG283 GG284 GG285 GG286 GG287 GG288 GG289 GG290 GG291 GG292 GG293 GG294 {
replace `var'=`var'-1
fre `var'
}

**********************************************************************************************************************************

*STRATEGY 2: IMPUTE 2 x ADHD SUBSCALES, THEN ADD THE TWO POST-IMPUTATION FOR THE FULL SCALE.

mi unregister ADHD_8yr 
drop ADHD_8yr

*************FOR SPEED: prune stuff we definitely don't need(move upstream later):
*rename this one first:
mi unregister EE584
rename EE584 save_EE584
mi register imputed save_EE584 

drop *AA1123 *AA1552 /**AA1315*/ *AA1316 *AA1317 *AA806 *AA807 *AA869 *AA870 *AA878 *AA879 *AA1548 *AA1549 *AA1550 *AA1551
drop G_* FF* EE* GG* rGG* LL* JJ* NN* rNN* 
drop *mumsmokes_Q6 *mumsmokes_Q5y *mumsmokes_Q7y *dadsmokes_Q7y *mumsmokes_Q8y *dadsmokes_Q8y *dadsmokes_QF2 
capture drop *CD_8yr *OD_8yr 
drop *SCQ* 
drop /* *conners_Q5y */ *delayed_motor_Q8yr *delayed_language_Q8y *hyperactivity_Q8y *concentration_Q8y *behavioural_problems_Q8y *emotional_diffic_Q8y *other_condition_Q8y 
drop *mhopkins_Q6 *mhopkins_Q5y *mhopkins_Q8y *fhopkins_QF2  
*drop /*complete-case versions of outcomes, since using versions allowing a small amount of item-level missingness*/ cc_*
*drop the test versions of scores from precise:
drop *_5e08

*190 vars left
*40,949 obs if restricting to people with some relevant questionnaire information

save "scratch/pre_imputation.dta", replace

***AT THIS POINT, TRANSFER THE SAVED FILE TO COLOSSUS AND RUN THE IMPUTATION MODEL THERE.