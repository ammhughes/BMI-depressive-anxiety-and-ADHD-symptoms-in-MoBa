/*# Covariate file #
The FID and IID columns match the FID (column 1) and IID (column 2) of the fam file.
The ID_2306 column contains the F_ID_2306 for fathers, the M_ID_2306 for mothers, and the PREG_ID_BARN_NR for children created using the following command in R: paste(child$PREG_ID, child$BARN_NR, sep="_").
For the majority of individuals the IID is the SENTRIXID, however, for the majority of individuals who have been genotyped multiple times, "_D1", "_D2", "_T1", "_T2", and "_T3" have been added to the end of the IID.
The first 20 principle components are included in the file.
Please note there are significant genotyping batch effects. Therefore, we highly recommend including these as covarites in all analyses. Any column with _num at the end of the column name contains a numeric version of the column with the same name (e.g., imputation_batch and imputation_batch_num).
*/
*DUPLICATES CHECK:
/*
sysdir set PLUS "N:\durable\people\Mandy\Statafiles\ado\plus"
cd "N:\durable\projects\mh_bmi\fullsample"

use "MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.dta", clear

*how many duplicates?
duplicates report sentrixid
duplicates report iid
*none?
*however, there is evidence of duplicates:
count if strpos(iid, "_D1")!=0
count if strpos(iid, "_D2")!=0
count if strpos(iid, "_T1")!=0
count if strpos(iid, "_T2")!=0
count if strpos(iid, "_T3")!=0
*you idiot - of course there wouldn't be duplicate values of sentrixid, these are ids for genotyped samples, so each of a person's duplicates would have a different code

*how to de-dupe?
*strip off the unique part of the iid:
split iid, p("_D1" "_D2" "_T1" "_T2" "_T3")
list iid iid1 if strpos(iid, "_D1")!=0
list iid iid1  if strpos(iid, "_D2")!=0
list iid iid1  if strpos(iid, "_T1")!=0
list iid iid1  if strpos(iid, "_T2")!=0
list iid iid1  if strpos(iid, "_T3")!=0
duplicates report iid1

*still no duplicates?
*try this another way - how many contain the stem of a random person who had an iid with _T1?
count if strpos(iid, "201048710077_R06C02")!=0
*1.
*hang on - has this already been de-duped? it must have been
drop iid1
*need to rename sentrixid iid so you can merge it with the snp-level data...
*no! this has been de-duped, so don't need to worry about that, but merging on (renamed) sentrixid to the .raw files will mean anyone who was one of a set of duplicates will be unmatched. so keep the iid as it is, and just use that.
*/

*RESHAPING AND REFORMATTING

sysdir set PLUS "N:\durable\people\Mandy\Statafiles\ado\plus"
cd "N:\durable\projects\mh_bmi\fullsample"

use "MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.dta", clear

*pre-reshape: make a simplified batch variable by collapsing the smallest groups
fre genotyping_batch

*current numerical version of batch variable is unhelpful because unlablled - recreate a labelled version:
encode genotyping_batch, gen(genotyping_batch_num_labelled)
drop genotyping_batch_num
rename genotyping_batch_num_labelled genotyping_batch_num
*biggest group is 19.norment_nov2020_1066, followed by: 6.norment_aug2020_1029, and 7.norment_aug2020_996 

*SORT OUT DUMMY VARIABLES.
*identify the redundant one (largest, baseline category) in each set
*genotyping centre: make dummies:
tab genotyping_center, gen(genotyping_center)
fre genotyping_center*
*biggest is group 3, deCODE_Genetics_Rekjavik_Iceland

*genotyping chip:
tab genotyping_chip, gen(genotyping_chip)
fre genotyping_chip*
*biggest is 2=Illumina_GlobalScreeningArray_MD_v.3.0

*plate id:
fre plate_id 
*cannot include that in models. 

*genotyping batch: way too many categories to include in models, also makes no sense as adjusting for the outcome (at least for ADHD)
*largest=norment_aug2020_1029
tab genotyping_batch, gen(genotyping_batch)
fre genotyping_batch*
*nb: post-rehape, there only 0s and missing for groups 11 and 12 of original categorical variable, norment_jan2015 and norment_jun2015.
*hence there are only 21 real dummies for children, but 23 for mums and dads.

*imputation release
*for children:
fre imputation_release1_hce imputation_release1_omni imputation_release1_gsa imputation_release2 imputation_release3 imputation_release4

*imputation batch: looks like dummies have already been made, but given the name stem 'imputation release'
foreach var of varlist imputation_release1_hce-imputation_release4 {
tab  imputation_batch `var'
}

*adjust for center and chip together, but not batch?
tab genotyping_center genotyping_chip

*rename the PCs.
rename pc* PC*

save "scratch/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov_updated.dta", replace

********************************************************************************

*CONTINUE TO MAKE RESHAPED VERSIONS:

/*shorten long variables whose names will extend in the reshape?
rename genotyping* geno*
rename imputation* imp*
*/
*get back PREG_ID from the ID_2306 column, where this is equivalent to PREG_ID_BARN_NR?
fre role
*hang on - need to do this separately by all the roles, then merge to the non-genetic data using F_ID_2306, M_ID_2306 and (PREG_ID BARN_NR) for the three groups in turn.
*children:
preserve
list id_2306 if role=="Child"
*split off into two at the underscore:
keep if role=="Child"
split id_2306, p("_")
*variables created as string: 
*id_23061  id_23062
rename id_23061 PREG_ID_2306
*make this numeric so it works in a merge:
destring PREG_ID_2306, gen(PREG_ID_2306_num)
*check it
list PREG_ID_2306 PREG_ID_2306_num in 1/10
drop PREG_ID_2306
rename PREG_ID_2306_num PREG_ID_2306
rename id_23062 BARN_NR
*make this numeric so it works in a merge:
destring BARN_NR, gen(BARN_NR_num)
*check it
list BARN_NR BARN_NR_num in 1/10
drop BARN_NR
rename BARN_NR_num BARN_NR
*drop that
drop id_2306
*rename everything else to have an identifying prefix:
rename sentrixid c_sentrixid
rename fid c_fid
rename iid c_iid
rename yob c_yob
rename sex c_sex
rename genotyping* c_genotyping* 
rename imputation* c_imputation* 
rename plate* c_plate*
rename PC* c_PC*
rename batch* c_batch*
drop role
save "scratch/child_genetic_cov.dta", replace
restore

*fathers:
preserve
keep if role=="Father"
rename id_2306 F_ID_2306
*rename everything else to have an identifying prefix:
rename sentrixid f_sentrixid
rename fid f_fid
rename iid f_iid
rename yob f_yob
rename sex f_sex
rename genotyping* f_genotyping* 
rename imputation* f_imputation* 
rename plate* f_plate*
rename PC* f_PC*
rename batch* f_batch*
drop role
save "scratch/father_genetic_cov.dta", replace
restore

*mothers:
preserve
keep if role=="Mother"
rename id_2306 M_ID_2306
*rename everything else to have an identifying prefix:
rename sentrixid m_sentrixid
rename fid m_fid
rename iid m_iid
rename yob m_yob
rename sex m_sex
rename genotyping* m_genotyping* 
rename imputation* m_imputation* 
rename plate* m_plate*
rename PC* m_PC*
rename batch* m_batch*
drop role
save "scratch/mother_genetic_cov.dta", replace
restore

*then merge these all together using the SV_INFO file.
use "N:\durable\projects\mh_bmi\fullsample\PDB2306_SV_INFO_v12.dta", clear
count
merge 1:m PREG_ID_2306 using "scratch/child_genetic_cov.dta"
rename _merge _merge1
*N increased here because of multiple births, presumably (1:m merge)
count
merge m:1 M_ID_2306 using "scratch/mother_genetic_cov.dta"
rename _merge _merge2
merge m:1 F_ID_2306 using "scratch/father_genetic_cov.dta"
rename _merge _merge3
*for each role, a group of PREG_ID_2306 values not merged from master (i.e., non-genetic SV_INFO file)
*nobody is in the genetic covariate file but not in the SV_INFO file
*of the unmerged, there are the least for mothers, more for children, most for fathers. as expected.
count
*how many PREG_ID_2306 values not represented in the genetic data even partally, i.e. under any role?
count if _merge1==1 & _merge2==1 & _merge3==1
*8,141
*presumably families who were never genotyped even in part.

*check that the fid is consistent within a family (where missing for nobody):
count if  c_fid==m_fid & !(c_fid=="" | m_fid =="")
count if  c_fid!=m_fid & !(c_fid=="" | m_fid =="")
count if  c_fid==f_fid & !(c_fid=="" | f_fid =="")
count if  c_fid!=f_fid & !(c_fid=="" | f_fid =="")
*for good measure:
count if  m_fid==f_fid & !(m_fid=="" | f_fid =="")
count if  m_fid!=f_fid & !(m_fid=="" | f_fid =="")
*good.

*make a flag for complete trios:
gen complete_trio=1 if c_iid!="" & m_iid!="" & f_iid!=""
fre complete_trio
*N=42228

*having checked that, make a single variable for fid, initially equal to that of the child: 
gen fid=c_fid 
count if fid==""
*sub in for where child isn't genotyped:
replace fid=m_fid if fid==""
replace fid=f_fid if fid==""
count if fid==""
*AS EXPECTED - THESE ARE THE PREGNANCIES FOR WHICH NONE OF THE THREE PEOPLE WERE IN THE GENOTYPE FILE.

*drop these people from the genetic covariate file:
drop if _merge1==1 & _merge2==1 & _merge3==1
count
*N=104,816

save "data/genetic/reshaped_genetic_cov.dta", replace
