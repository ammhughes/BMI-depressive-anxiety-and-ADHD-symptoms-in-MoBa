*EA3 - SCORE MADE IN STATA

*snps have been clumped in mrbase from the set common to the gwas and the moba_qc_bim file.
*961 clumped with parameters pval.exposure<5e-08,clump_r2 = 0.01
*snps extracted March 2021 from new qc'd dataset with txt in this script: /tsd/p471/data/durable/projects/mh_bmi/scripts/colossus/manual_ea3_prs.sh
*(had to paste it into the command line because colossus not accepting the job in the queue)

*************************************************************

*set wd to correct project folder:
*sysdir set PLUS "N:\data\durable\people\Mandy\Statafiles\ado\plus"
cd "N:\data\durable\projects\mh_bmi\"

*get EA3 snps from the parental education project folder:

//Import GWAS results GWS clumped excluding MoBa and 23andMe
import delimited "N:\data\durable\projects\parental_education\data\lee_clumped_snps_ex_moba.txt", delimiter(space) bindquote(nobind) varnames(1) asfloat clear 
gen n=_n
save "data\lee_clumped_snps_ex_moba",replace

/*
//Import GWAS results GWS clumped including MoBa and 23andMe
import delimited "N:\data\durable\projects\parental_education\data\lee_clumped_snps_incl_moba.txt", delimiter(space) bindquote(nobind) varnames(1) asfloat clear 
gen n=_n
save "data\lee_clumped_snps_incl_moba",replace

//Import GWS non-independent GWAS results including MoBa and 23andMe
import delimited "N:\data\durable\projects\parental_education\data\lee_clumped_snps_incl_moba_non_indep.txt", delimiter(space) bindquote(nobind) varnames(1) asfloat clear 
gen n=_n
save "data\lee_clumped_snps_incl_moba_non_indep",replace
*/

*************************************************************
*sysdir set PLUS "N:\data\durable\people\Mandy\Statafiles\ado\plus"
cd "N:\data\durable\projects\mh_bmi\"

*FIRST STEP: COMPARE ALIGNMENT IN MOBA DATA VS GWAS AND HARMONIZE.

*UPDATE: need to compare non-effect alleles as well, and allele frequencies, because palindromic snps.
*this info isn't in the dosage (.raw) file, so you need the moba .afreq file (v similar to bim but with a allele frequency).
*made this in plink 22/03/21 and added it to FINAL QC VERSION folder off colossus

/*
import delimited  "data/mobq_qc_freq.afreq", clear
rename id SNP
rename ref ref_allele_moba
rename alt other_allele_moba
rename alt_freqs other_allele_freq_moba
save "data/mobq_qc_freq_afreq.dta", replace
*/

*then, you need the rsid, A1, A2 and beta etc from the gwas - version excluding MoBa and 23andme.
*import and prep:
use "data\lee_clumped_snps_ex_moba.dta", clear
rename snp SNP
rename effect_alleleexposure effect_allele_gwas
rename other_alleleexposure other_allele_gwas
rename eafexposure effect_allele_freq_gwas
rename betaexposure beta_gwas
rename seexposure se_gwas
rename pvalexposure p_gwas
drop mr_keepexposure pval_originexposure idexposure 
save "data/formatted_clumped_moba_ea3_snps_allinfo", replace

*merge and keep common ones:
merge 1:1 SNP using "data/mobq_qc_freq_afreq.dta"
keep if _merge==3
drop _merge
*311 left

*HARMONIZE. 
*Table 1 of this is a good ref for the steps, although doesn't give a rule of thumb for when eaf is too close to 50% to make a guess: https://academic.oup.com/ije/article/45/6/1717/3072174

*make allele freq vars for the other allele in the gwas, and the ref allele in moba, by subtraction (assuming none are triallelic):
*check distribution of eaf from the gwas:
summ effect_allele_freq_gwas
*this runs from 0.02 to 0.98. fine. 
gen other_allele_freq_gwas=1-effect_allele_freq_gwas
label variable other_allele_freq_gwas "other allele frequency in gwas"

*check distribution of other allele in moba:
summ other_allele_freq_moba
*hmmm. ok, it looks like the 'other' allele is always the minor allele, since this runs 0.02-0.49. So the compliment, ref_allele_freq_moba, will always be >0.50
gen ref_allele_freq_moba=1-other_allele_freq_moba
label variable ref_allele_freq_moba "reference allele frequency in moba"
summ ref_allele_freq_moba
***************IN MOBA DATA, AT LEAST IN THE FREQ FILE, REF ALLELE IS THE MAJOR ALLELE.******************

*define same or reverse strand alignment
capture drop strand_alignment
gen strand_alignment=.
capture label define strand_alignment 1"same" 2"reverse"
label values strand_alignment strand_alignment
replace strand_alignment=1 if effect_allele_gwas==ref_allele_moba & other_allele_gwas==other_allele_moba
replace strand_alignment=2 if other_allele_gwas==ref_allele_moba & effect_allele_gwas==other_allele_moba
tab strand_alignment

*however, some of these could have had alignment misclassified because they are palindromes.
*possible palindromes:
*where the pair of possible alleles are A/T or C/G, and 
*where the allele frequencies are different enough that they could could plausibly refer to different alleles, i.e.:
*for ones we think are "same" aligned, where strand_alignment==1 & (abs(effect_allele_freq_gwas-ref_allele_freq_moba))>0.05, or 
*for ones we think are "reverse" aligned, where strand_alignment==2 & (abs(effect_allele_freq_gwas-other_allele_freq_moba))>0.05 
*is this too conservative? not enough?
list effect_allele_gwas other_allele_gwas ref_allele_moba other_allele_moba effect_allele_freq_gwas ref_allele_freq_moba strand_alignment if ///
((effect_allele_gwas=="A" & other_allele_gwas=="T") | (effect_allele_gwas=="T" & other_allele_gwas=="A") | (effect_allele_gwas=="C" & other_allele_gwas=="G") | (effect_allele_gwas=="G" & other_allele_gwas=="C")) ///
& strand_alignment==1 & (abs(effect_allele_freq_gwas-ref_allele_freq_moba))>0.05
*none.
list effect_allele_gwas other_allele_gwas ref_allele_moba other_allele_moba effect_allele_freq_gwas ref_allele_freq_moba strand_alignment if ///
((effect_allele_gwas=="A" & other_allele_gwas=="T") | (effect_allele_gwas=="T" & other_allele_gwas=="A") | (effect_allele_gwas=="C" & other_allele_gwas=="G") | (effect_allele_gwas=="G" & other_allele_gwas=="C")) ///
& strand_alignment==2 & (abs(effect_allele_freq_gwas-other_allele_freq_moba))>0.05
*none.
*also, ones where the frequencies are close enough to 50% each that it could be wrong even whilst passing the above filter:
list effect_allele_gwas other_allele_gwas ref_allele_moba other_allele_moba effect_allele_freq_gwas ref_allele_freq_moba strand_alignment if ///
((effect_allele_gwas=="A" & other_allele_gwas=="T") | (effect_allele_gwas=="T" & other_allele_gwas=="A") | (effect_allele_gwas=="C" & other_allele_gwas=="G") | (effect_allele_gwas=="G" & other_allele_gwas=="C")) ///
& (effect_allele_freq_gwas>0.40 & effect_allele_freq_gwas<0.60)
*12.

capture drop palindrome
gen palindrome=0
replace palindrome=1 if ///
((effect_allele_gwas=="A" & other_allele_gwas=="T") | (effect_allele_gwas=="T" & other_allele_gwas=="A") | (effect_allele_gwas=="C" & other_allele_gwas=="G") | (effect_allele_gwas=="G" & other_allele_gwas=="C")) ///
& ((strand_alignment==1 & abs(effect_allele_freq_gwas-ref_allele_freq_moba)>0.05) ///
| (strand_alignment==2 & abs(effect_allele_freq_gwas-other_allele_freq_moba)>0.05) ///
| (effect_allele_freq_gwas>0.40 & effect_allele_freq_gwas<0.60))

*check this by plotting allele frequencies:
twoway scatter effect_allele_freq_gwas ref_allele_freq_moba if strand_alignment==1
*here they're all identical
twoway scatter effect_allele_freq_gwas ref_allele_freq_moba if strand_alignment==2
*here they all SUM to 1

*make a harmonized version of the beta:
gen harm_beta_gwas=beta_gwas
replace harm_beta_gwas=-(harm_beta_gwas) if strand_alignment==2
list harm_beta_gwas beta_gwas strand_alignment in 1/50
*remove possible palindromes:
*replace harm_beta_gwas=. if palindrome==1
*DROP THEM: otherwise, missing values of beta mess with the pgs calculation
drop if palindrome==1

*make a harmonized version of the allele frequency corresponding to the harmonized beta:
gen harm_effect_allele_freq_gwas=effect_allele_freq_gwas 
replace harm_effect_allele_freq_gwas=1-effect_allele_freq_gwas if strand_alignment==2

*correlate allele frequencies:
*for ones originally same-aligned:
twoway scatter harm_effect_allele_freq_gwas ref_allele_freq_moba if strand_alignment==1
twoway scatter harm_effect_allele_freq_gwas other_allele_freq_moba if strand_alignment==1
*for ones originally reverse-aligned:
twoway scatter harm_effect_allele_freq_gwas ref_allele_freq_moba if strand_alignment==2
twoway scatter harm_effect_allele_freq_gwas other_allele_freq_moba if strand_alignment==2
*seems ok!

sort SNP
save data/moba_gwas_snplevelinfo_ea3.dta, replace

********************************************

/*NEXT CHECK:
*Check that the ref allele in the raw file (one included in the actual snp variable names) is the same as the ref allele in the .freq file (one to which freq refers, i.e., MAJOR ALLELE, since all frequencies >0.50):

import delimited "N:\data\durable\projects\mh_bmi\genetic\ea3_snps.raw", varnames(nonames) clear

*then just strip off the colnames:
keep in 1
*drop the stuff which isn't about snps (FID IID PAT MAT SEX PHENOTYPE):
drop v1-v6
*transpose:
sxpose, clear force
*sysdir

*make two vars for rsid and named allele with split 
rename _var1 snp
split snp, p("_")
rename snp1 SNP
rename snp2 named_allele
*drop original
drop snp
save scratch/rawtemp.dta, replace

*get freq file:
import delimited "N:/data/durable/projects/moba_interim_release_post_imp_qc/FINAL QC VERSION/mobq_qc_freq.afreq", clear
*(use from final qc folder, because keeping a copy in the mh_bmi folder takes up all your disk space)
rename id SNP
rename ref named_allele
rename alt other_allele
rename alt_freqs other_allele_freq_moba
save scratch/freqtemp.dta, replace

*merge on SNP and named_allele to formatted raw file. if ref allele is the same, all obs from .raw should match
merge 1:1 SNP named_allele using scratch/rawtemp.dta
*yep. all good.

clear
erase scratch/freqtemp.dta
erase scratch/rawtemp.dta
*/

********************************************

*IMPORT IID SNP-LEVEL INFO
import delimited "genetic\ea3_snps.raw", clear

drop phenotype
*rename these to make more sense:
rename pat pid
rename mat mid

*some of the rsid variables are string vars, some are byte vars.
*the string ones are ones which have some obs as 'NA', SNP-level missingness which varies between people. 
*need to fix that and have them all as byte vars
foreach var of varlist rs* {
   capture confirm string variable `var'
   if !_rc {
      replace `var'="" if `var'=="NA"
	  destring `var', gen(`var'_num)
	  drop `var'
	  rename `var'_num `var'
   }
}

save "scratch/ea3_genetic.dta", replace

*****************

*NOW MAKE THE PGS.

use scratch/ea3_genetic.dta, clear

*'merge' in snplevel info, now JUST THE ONES IN MOBA AS WELL AS GWAS, even though not really a merge
*in the other file, THE DATA NEEDS TO BE SORTED BY SNP SO THAT THIS COMES OUT IN THE RIGHT ORDER.
merge 1:1 _n using data/moba_gwas_snplevelinfo_ea3.dta

*save the order of snps in a global for later:
levelsof SNP, local(ea3_snps)
display "ea3_snps"
local ea3_snplist_temp: subinstr local ea3_snps "char(34)" "char(32)", all
local ea3_snplist_temp2: subinstr local ea3_snplist_temp "rs" "  rs", all
global ea3_snplist: subinstr local ea3_snplist_temp2 "  " ""
display $ea3_snplist

*check it works in a loop:
foreach snp in $ea3_snplist {
	display "`snp'"
}

/*NB: there is a little bit of SNP-level missingness for some people, enough that you lose about 600 trios
you don't notice this when making scores with PRSise because it does something to fill these in as specified here: 
https://www.prsice.info/step_by_step/#prs-calculation, specifically:
"The --missing option is used to determine how PRSice handle the missingness.
When not specified, the Minor Allele Frequency (MAF) in the target sample will be used as the genotype as the sample with missing genotype.
If --missing SET_ZERO is set, the SNP for the missing samples will be excluded. Alternatively, if --missing CENTER is set, all PRS calculated will be minused by the MAF of the SNP (therefore, missing samples will have PRS of 0).
*because in the scripts we don't specify anything with the missing option, it does the default. It's not clear what that is - that must be a typo or something because 'the Minor Allele Frequency (MAF) in the target sample will be used as the genotype as the sample with missing genotype' doesn't mean anything. The later comment 'NOTE: Missingness imputation is usually based on the target samples. If you would like to impute the missingness using the reference sample, you can use --use-ref-maf parameter to specify all MAF to be calculated using the reference samples' suggests it uses the MAF in the target sample to impute these values, which makes sense. But for now, here, obviously can't do that,
so do the equivalent of --missing SET_ZER0.
*/

*UPDATE 12 MAY: NEW VERSION OF PGS adding as if for mean genotype for snp-level missigness.
capture drop ea3_pgs*
gen ea3_pgs=0 
*if batch_num!=.
*nb: keeping as . for people who weren't genotyped (FIX THIS LATER)
local k=0
foreach var in $ea3_snplist {
display "`var'"
fre `var'
local k=`k'+1
local relevantbeta=harm_beta_gwas in `k'
display `relevantbeta'
replace ea3_pgs=ea3_pgs+(`relevantbeta'*0) if `var'==0
replace ea3_pgs=ea3_pgs+(`relevantbeta'*1) if `var'==1
replace ea3_pgs=ea3_pgs+(`relevantbeta'*2) if `var'==2
*for people for whom the snp is missing, ADD BETA*MEAN VALUE
sum `var'
return list
replace ea3_pgs=ea3_pgs+(`relevantbeta'*r(mean)) if `var'==.
}

*standardize and rename it for consistency:
egen z_ea3_pgs=std(ea3_pgs)
summ z_ea3_pgs

drop _merge
save "scratch/ea3_genetic_pgs.dta", replace

********************************************************************

*31/03/21 add in the new PCs:
/*
import delimited "N:\data\durable\projects\moba_interim_release_post_imp_qc\FINAL QC VERSION\final.1kg-pca_updated.eigenvec", delimiter(space) clear 
rename v1 fid
rename v2 iid
local k=0
foreach var of varlist v3-v22 {
	local k=`k'+1
	rename `var' PC`k'
}
*duplicate iid values:
duplicates report iid
duplicates tag iid, gen(dup)
list fid iid dup if dup!=0
*these look like the people from the reference panel, no? drop them:
drop if dup!=0
drop dup
count
*93,211
save "data/final.1kg-pca_updated.eigenvec.dta", replace
clear
*/
use "scratch/ea3_genetic_pgs.dta", clear
capture drop _merge
merge 1:1 iid using "data/final.1kg-pca_updated.eigenvec.dta"
drop _merge
save "scratch/ea3_genetic_pgs.dta", replace

*******************************************************************************************

*PART 2:
*PREP AND RESHAPE to merge into phenotype data

cd "N:\data\durable\projects\mh_bmi\"

*start with the full linker file:
*keep everything, but rename:
*get pregnancy id, role etc
use "N:\data\durable\projects\moba_interim_release_post_imp_qc\scratch\cleaned_linkage_file",clear
rename SENTRIXID iid
rename Role role
rename BATCH batch
save "data\2019_10_01_MoBaGenetics_all_2306.dta", replace

*merge to the genetic data including the pgs:
merge 1:1 iid using "scratch/ea3_genetic_pgs.dta"
keep if _merge==3
drop _merge
save "scratch/linkable_ea3_genetic_pgs.dta", replace

*NB: if you don't either a) remove relateds, or b) specify BARN_NR in the i variable you get a duplicates problem due to twins
*count if role==""
duplicates report PREG_ID_2306 role
duplicates tag PREG_ID_2306 role, gen(dup)
*examine:
list PREG_ID_2306 BARN_NR fid iid pid mid sex role if dup==1
*so those are the genuine twins! fine.
*nb: also some other halves of twins:
list PREG_ID_2306 BARN_NR fid iid pid mid sex role if BARN_NR==2
*drop one of each of the problem pairs:
duplicates drop PREG_ID_2306 role, force
drop dup
*count
*93,208, of which:
tab role
/*
       role |      Freq.     Percent        Cum.
------------+-----------------------------------
      Child |     31,374       33.66       33.66
     Father |     30,239       32.44       66.10
     Mother |     31,595       33.90      100.00
------------+-----------------------------------
      Total |     93,208      100.00
*/

*Now reshape it:
/**NB: something goes seriously wrong whe using the reshape command. the right people aren't linked within families, because the parental and child scores don't even correlate
*first drop the snp-level gwas info:
*drop SNP-snp_present
*reshape wide iid BARN_NR batch PDB SampleType Retrieval_ID fid pid mid sex rs* *pgs*, i(PREG_ID_2306) j(role) string
*/

*do this Neil's way to see if it avoids the fuckup:

*first drop the snp-level gwas info:
drop SNP-harm_effect_allele_freq_gwas
*then these bits:
drop SampleType PDB SampleType Retrieval_ID
//Need create three files, one for offspring, one for mother and one for fathers
preserve
keep if role=="Child"
*do these renames later to make the mendelian error check easier
*renpfix rs c_rs
*rename ea3_pgs c_ea3_pgs
*rename z_ea3_pgs c_z_ea3_pgs
rename batch c_batch
renpfix PC c_PC
save "scratch/offspring",replace
restore
preserve
keep if role=="Father"
renpfix rs f_rs
rename ea3_pgs f_ea3_pgs
rename z_ea3_pgs f_z_ea3_pgs
rename batch f_batch
renpfix PC f_PC
save "scratch/father",replace
restore
preserve
keep if role=="Mother"
renpfix rs m_rs
rename ea3_pgs m_ea3_pgs
rename z_ea3_pgs m_z_ea3_pgs
rename batch m_batch
renpfix PC m_PC
save "scratch/mother",replace
restore

//Next we need to merge each of these files togeather
use  "scratch/offspring",clear
*31,374
duplicates report PREG_ID_2306
*none
desc using  "scratch/father",
*30,239 
desc using  "scratch/mother",
*31,595
//Just keep ids, sex, and genetic data
keep PREG_ID_2306 BARN_NR fid iid pid mid sex *rs* *ea3* *batch* *PC*
*drop rsi*

//Match in fathers genotypes
rename iid offspring_id
*fix the 0s in pid!!!! otherwise you get duplicates issue
replace pid="" if pid=="0"
rename pid iid
joinby iid using  "scratch/father",unmatched(master) _merge(_merge2)
tab _merge2
drop _merge*
drop role 
drop pid
rename iid pid
//Match in mothers genotypes
*fix the 0s in pid!!!! otherwise you get duplicates issue
replace mid="" if mid=="0"
rename mid iid
joinby iid using  "scratch/mother",unmatched(master) _merge(_merge2)
tab _merge2
drop role
drop mid
drop _merge*
rename iid mid

//Check for Mendelian errors
ds rs*
foreach i in `r(varlist)'{
	reg `i' f_`i' m_`i',ro
}
*and with the scores:
pwcorr z_ea3_pgs m_z_ea3_pgs  f_z_ea3_pgs 

*and the PCs:
forvalues i= 1(1)20 {
    regress c_PC`i' m_PC`i' f_PC`i'
}

*then rename the child's genetic vars with equivalent prefix
renpfix rs c_rs
rename ea3_pgs c_ea3_pgs
rename z_ea3_pgs c_z_ea3_pgs

*UPDATE 12 MAY: RESTRICT TO PEOPLE IN UPDATED SV_INFO FILE (i.e. exclude the new consent withdrawls) 
merge m:1 PREG_ID_2306 using "N:\data\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\Statafiles\PDB2306_SV_INFO_v12_UPDATED_April_2021"
*the 11 not matched are the ones that need to be dropped as new exclusions.
drop if _merge!=3
drop _merge
count
*31,363

save "scratch/reshaped_linkable_ea3_genetic_pgs.dta", replace


*****************************************************************************************************
*suspect there are space issues, so erase intermediate files:
capture erase "scratch/ea3_genetic_pgs.dta"
*not that one, it's useful later
*capture erase "scratch/linkable_ea3_genetic_pgs.dta"
******************************************************************************************************
*****************************************************************************************************

*MERGE TO PHENOTYPE DATA:

cd "N:\data\durable\projects\mh_bmi\"
use "scratch/reshaped_linkable_ea3_genetic_pgs.dta", replace

*OK! Now merge this in the with the rest of the pre-imputation phenotype data and check relationships in the complete-case.

merge 1:1 PREG_ID_2306 BARN using "data\statafiles\phenotypes_bmi_analyses.dta"
*no point in keeping ones without genetic info:
keep if _merge==3

*delineate the usable trios:
gen fullgeneticinfo=0
replace fullgeneticinfo=1 if (c_batch!="" & m_batch!="" & f_batch!="")
tab fullgeneticinfo
*26,644 trios

save "data\statafiles\geno_pheno_ea3.dta", replace
******************************************************************************************************

*COMPLETE-CASE ANALYSIS ENTRY POINT:

sysdir set PLUS "N:\data\durable\people\Mandy\statafiles\ado\plus"
cd "N:\data\durable\projects\mh_bmi\"

use "data\statafiles\geno_pheno_ea3.dta", clear

*COMPLETE-CASE ANALYSIS:

*how much of the education variables does the pgs predict?

*Recap of cleaned placeholder educational outcomes
*before record linkage, using Q25 items NN239 NN240 NN241: Three questions about school results on national exams
*nb: I coded 4='don't know/have not talked to teacher about it' to missing
fre NN239 NN240 NN241

*also:
*Q22 How is your child enjoying school?
fre NN234

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

*28 Reading and writing skills
*NB ITEMS DIFFER BETWEEN VERSIONS 
*In version C
fre NN380 NN381 NN382 NN383 NN384 
*In version A & B
fre NN254 NN255 NN256 NN257

*30: Home reading and self-reading
fre NN262 NN263 NN385 NN264 NN386 NN265

************************************************************************************

*None of these are normal of even continuous, so can't get an R2. But look at average PGS across categories:
*3 x exam results, school enjoyment, self-reading
foreach var in NN239 NN240 NN241 NN234 NN265 {
    tab `var', summ (c_z_ea3_pgs)
}

*OK, seems to go in the right direction at least.

fre NN239 NN240 NN241 NN234 NN265
*for imputation, also change the last two:
*school enjoyment: collpase and flip for bigger baseline categ:
gen NN234_new=NN234
label define NN234_new 1"very well" 2"well" 3"neither/poor/very poor"
label values NN234_new NN234_new
recode NN234_new 1/3=1 4=2 5=3
fre NN234_new
replace NN234_new=5-NN234_new
fre NN234_new


