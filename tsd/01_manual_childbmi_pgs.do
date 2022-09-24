*CHILD BMI - SCORE MADE IN STATA

*using the set of snps clumped (externally, in mrbase) from Tom's full GWAS results, subset to those in the final moba .bim file.

*************************************************************

*HARMONIZATION OF SNP-LEVEL INFO
*need to compare also non-effect alleles and allele frequencies
*this info isn't in the dosage (.raw) file, so you need the moba .afreq file (v similar to bim but with a allele frequency).
*made this in plink 22/03/21 and added it to FINAL QC VERSION folder off colossus

sysdir set PLUS "N:\durable\people\Mandy\Statafiles\ado\plus"
cd "N:\durable\projects\mh_bmi\fullsample"

/*
import delimited "data/genetic/mobq_qc_freq.afreq"
rename id SNP
rename ref ref_allele_moba
rename alt other_allele_moba
rename alt_freqs other_allele_freq_moba
save "data/genetic/mobq_qc_freq_afreq.dta", replace
*/

*then, you need the rsid, A1, A2 and beta etc from the correct table of the gwas:
*import and prep:
import delimited "data/snplevel_info/clumped_moba_childbmi_snps_allinfo.txt", case(preserve) clear 
*import delimited "genetic/childbmi_UKB_deduped_table_S9_childbmisnps.txt", case(preserve) clear 
rename effect_allele effect_allele_gwas
rename other_allele other_allele_gwas
rename eaf effect_allele_freq_gwas
rename betaexposure beta_gwas
rename seexposure se_gwas
rename pvalexposure p_gwas
*drop the confusing things which you're not using:
drop CHISQ_LINREG P_LINREG CHISQ_BOLT_LMM_INF P_BOLT_LMM_INF CHISQ_BOLT_LMM idexposure
save "data/snplevel_info/formatted_childbmisnps", replace

*merge and keep common ones:
merge 1:1 SNP using "data/genetic/mobq_qc_freq_afreq.dta"
keep if _merge==3
drop _merge
*321 

*HARMONIZE. 
*Table 1 of this is a good ref for the steps, although doesn't give a rule of thumb for when eaf is too close to 50% to make a guess: https://academic.oup.com/ije/article/45/6/1717/3072174

*make allele freq vars for the other allele in the gwas, and the ref allele in moba, by subtraction (assuming none are triallelic):
*check distribution of eaf from the gwas:
summ effect_allele_freq_gwas
*this runs from 0.11 to 0.98. fine. 
gen other_allele_freq_gwas=1-effect_allele_freq_gwas
label variable other_allele_freq_gwas "other allele frequency in gwas"

*check distribution of other allele in moba:
summ other_allele_freq_moba
*hmmm. ok, it looks like the 'other' allele is always the minor allele, since this runs 0.02-0.49. So the compliment, ref_allele_freq_moba, will always be >0.50
gen ref_allele_freq_moba=1-other_allele_freq_moba
label variable ref_allele_freq_moba "reference allele frequency in moba"
summ ref_allele_freq_moba
***************IN MOBA DATA, AT LEAST IN THE FREQ FILE, REF ALLELE IS THE MAJOR ALLELE.******************

*define alignment of gwas effect allele w/ respect to moba reference allele:
capture drop eaf_alignment
gen eaf_alignment=.
capture label define eaf_alignment 1"harmonized" 2"alternative"
label values eaf_alignment eaf_alignment
replace eaf_alignment=1 if effect_allele_gwas==ref_allele_moba & other_allele_gwas==other_allele_moba
replace eaf_alignment=2 if other_allele_gwas==ref_allele_moba & effect_allele_gwas==other_allele_moba
fre eaf_alignment
*both of these categories mean that the GWAS and MoBa genotyped in the same direction (both on forward strand, or both on reverse strand) 
*are there any SNPs which don't meet these conditions, indicating that for any SNPs, the GWAS used the forward strand and MoBa the reverse or vice-versa? 
*no, because the total of nonmissing observations is equal to all the SNPs.

*check by plotting allele frequencies:
twoway scatter effect_allele_freq_gwas ref_allele_freq_moba if eaf_alignment==1
*here they're all identical
twoway scatter effect_allele_freq_gwas ref_allele_freq_moba if eaf_alignment==2
*here they all SUM to 1


*check if the alignment of gwas effect allele w/ respect to moba reference allele differs between palindromic and not:
gen palindromic=1 if ///
((effect_allele_gwas=="A" & other_allele_gwas=="T") | (effect_allele_gwas=="T" & other_allele_gwas=="A") | (effect_allele_gwas=="C" & other_allele_gwas=="G") | (effect_allele_gwas=="G" & other_allele_gwas=="C")) 
recode palindromic .=0
tab eaf_alignment palindromic, row
*as expected

**************************************************************************
*UPDATE after speaking with Sean: 
*not necessary to chuck out any SNPs because of A/T C/G combinations, regardless of eaf.
*whether genotyping was on forward or reverse strand is decided at the GWAS level, so the check above (for alignment of gwas effect allele w/ respect to moba reference allele, which gave no missing values for any SNPs) is sufficient.
*DO NOT remove palindromes
**************************************************************************

*make a harmonized version of the beta:
gen harm_beta_gwas=beta_gwas
replace harm_beta_gwas=-(harm_beta_gwas) if eaf_alignment==2
list harm_beta_gwas beta_gwas eaf_alignment in 1/50

*make a harmonized version of the allele frequency corresponding to the harmonized beta:
gen harm_effect_allele_freq_gwas=effect_allele_freq_gwas 
replace harm_effect_allele_freq_gwas=1-effect_allele_freq_gwas if eaf_alignment==2

*correlate allele frequencies:
*for ones originally same-aligned:
twoway scatter harm_effect_allele_freq_gwas ref_allele_freq_moba if eaf_alignment==1
twoway scatter harm_effect_allele_freq_gwas other_allele_freq_moba if eaf_alignment==1
*for ones originally reverse-aligned:
twoway scatter harm_effect_allele_freq_gwas ref_allele_freq_moba if eaf_alignment==2
twoway scatter harm_effect_allele_freq_gwas other_allele_freq_moba if eaf_alignment==2
*seems ok!

sort SNP
save data/snplevel_info/moba_gwas_snplevelinfo_childbmi.dta, replace

********************************************

/*NEXT CHECK:
*Check that the ref allele in the raw file (one included in the actual snp variable names) is the same as the ref allele in the .freq file (one to which freq refers, i.e., MAJOR ALLELE, since all frequencies >0.50):

import delimited "data\genetic\childbmi_snps.raw", varnames(nonames) clear

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
import delimited "data/genetic/mobq_qc_freq.afreq", clear
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

****************************************************************************************
*IMPORT IID SNP-LEVEL INFO
import delimited "data/genetic/childbmi_snps.raw", clear

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

save "scratch/childbmi_genetic.dta", replace

*****************
*NOW MAKE THE PGS.

use "scratch/childbmi_genetic.dta", clear

*'merge' in snplevel info, now JUST THE ONES IN MOBA AS WELL AS GWAS, even though not really a merge
*in the other file, THE DATA NEEDS TO BE SORTED BY SNP SO THAT THIS COMES OUT IN THE RIGHT ORDER.
merge 1:1 _n using "data/snplevel_info/moba_gwas_snplevelinfo_childbmi.dta"

list SNP harm_beta_gwas in 1/322
*yep

*save the order of snps in a global for later:
levelsof SNP, local(childbmi_snps)
display "childbmi_snps"
local childbmi_snplist_temp: subinstr local childbmi_snps "char(34)" "char(32)", all
local childbmi_snplist_temp2: subinstr local childbmi_snplist_temp "rs" "  rs", all
global childbmi_snplist: subinstr local childbmi_snplist_temp2 "  " ""
display $childbmi_snplist

*check it works in a loop:
foreach snp in $childbmi_snplist {
	display "`snp'"
}
/*NB: there is a little bit of SNP-level missingness for some people
this is also the case for the adult bmi snps - but PRSise does something to fill these in as specified here: 
https://www.prsice.info/step_by_step/#prs-calculation, specifically:
"The --missing option is used to determine how PRSice handle the missingness.
When not specified, the Minor Allele Frequency (MAF) in the target sample will be used as the genotype as the sample with missing genotype.
If --missing SET_ZERO is set, the SNP for the missing samples will be excluded. Alternatively, if --missing CENTER is set, all PRS calculated will be minused by the MAF of the SNP (therefore, missing samples will have PRS of 0).
*because in the scripts we don't specify anything with the missing option, it does the default. It's not clear what that is - that must be a typo or something because 'the Minor Allele Frequency (MAF) in the target sample will be used as the genotype as the sample with missing genotype' doesn't mean anything. The later comment 'NOTE: Missingness imputation is usually based on the target samples. If you would like to impute the missingness using the reference sample, you can use --use-ref-maf parameter to specify all MAF to be calculated using the reference samples' suggests it uses the MAF in the target sample to impute these values, which makes sense. But for now, here, obviously can't do that,
so do the equivalent of --missing SET_ZER0.
*/

*to quantify snp-level missingness:
capture gen snp_missingness=.
capture gen name_snp_missingness=""
local k=0
foreach var in $childbmi_snplist {
local k=`k'+1
display `var'
fre `var'
count if `var'==.
local cases = r(N)
replace snp_missingness=`cases' in `k'
replace name_snp_missingness="`var'" in `k'
}
summ snp_missingness, det
*as proportion:
gen prop_snp_missingness=snp_missingness/207409
summ prop_snp_missingness, det


*NEW VERSION OF PGS adding as if for mean genotype for snp-level missigness.
capture drop childbmi_pgs*
gen childbmi_pgs=0 
local k=0
foreach var in $childbmi_snplist {
display "`var'"
fre `var'
local k=`k'+1
local relevantbeta=harm_beta_gwas in `k'
display `relevantbeta'
replace childbmi_pgs=childbmi_pgs+(`relevantbeta'*0) if `var'==0
replace childbmi_pgs=childbmi_pgs+(`relevantbeta'*1) if `var'==1
replace childbmi_pgs=childbmi_pgs+(`relevantbeta'*2) if `var'==2
*for people for whom the snp is missing, ADD BETA*MEAN VALUE
sum `var'
return list
replace childbmi_pgs=childbmi_pgs+(`relevantbeta'*r(mean)) if `var'==.
}

*standardize and rename it for consistency:
egen z_childbmi_pgs=std(childbmi_pgs)
summ z_childbmi_pgs

drop _merge
save "scratch/childbmi_genetic_pgs.dta", replace

************************************************************************************************************
*ADD GENETIC COVARIATES:
*In new release, PCs are part of the covariate file now, so don't need an individual step to add them.

use "scratch/childbmi_genetic_pgs.dta", clear
*save the order of snps in a global for later (mendelian error checks):
levelsof SNP, local(childbmi_snps)
display "childbmi_snps"
local childbmi_snplist_temp: subinstr local childbmi_snps "char(34)" "char(32)", all
local childbmi_snplist_temp2: subinstr local childbmi_snplist_temp "rs" "  rs", all
global childbmi_snplist: subinstr local childbmi_snplist_temp2 "  " ""
display $childbmi_snplist

*check it works in a loop:
foreach snp in $childbmi_snplist {
	display "`snp'"
}
*can now drop snp-level information
drop SNP- prop_snp_missingness

*then stitch on covariate file, prepared version. Can either:
*1.add the long-format covariate file and reshape all the genetic information (pgs, snps and covariates together) before adding to phenotype data, or 
*2.save the pgsand snps separately for each person, then stitch on for each person in to the wide-format covariate file. the first is quicker, overall.
capture drop _merge
merge 1:1 iid using "data/genetic/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov_updated.dta"
drop _merge
*now proceed to reshape for everyone, as you did for the generic covariate file but now including the pgs and snp data:

count
tab role
/*
       Role |      Freq.     Percent        Cum.
------------+-----------------------------------
      Child |     76,530       36.90       36.90
     Father |     53,315       25.71       62.60
     Mother |     77,564       37.40      100.00
------------+-----------------------------------
      Total |    207,409      100.00
*/

*make three smaller files, with the ipd phenotype-specific genetic data for children, mothers and fathers respectively, in which all variables have an identifying prefix.

**********************
*Then separate this out into the three people by role, then merge to the non-genetic data using F_ID_2306, M_ID_2306 and (PREG_ID BARN_NR) for the three groups in turn.
fre role
*children:
preserve
*drop snp-level information:
capture drop CHR- prop_snp_missingness
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
*including individual snp info:
rename rs* c_rs*
*and the score:
rename *pgs* c_*pgs*
drop role
save "scratch/c_childbmi_genetic_pgs.dta", replace
restore

*mothers:
preserve
*drop snp-level information:
capture drop CHR- prop_snp_missingness
list id_2306 if role=="Mother"
*split off into two at the underscore:
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
*including individual snp info:
rename rs* m_rs*
*and the score:
rename *pgs* m_*pgs*
drop role
save "scratch/m_childbmi_genetic_pgs.dta", replace
restore

preserve
*drop snp-level information:
capture drop CHR- prop_snp_missingness
list id_2306 if role=="Father"
*split off into two at the underscore:
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
*including individual snp info:
rename rs* f_rs*
*and the score:
rename *pgs* f_*pgs*
drop role
save "scratch/f_childbmi_genetic_pgs.dta", replace
restore

*now attach all three to the SV_INFO file,
*then merge these all together using the SV_INFO file.
use "N:\durable\projects\mh_bmi\fullsample\PDB2306_SV_INFO_v12.dta", clear
count
*112,265
merge 1:m PREG_ID_2306 using "scratch/c_childbmi_genetic_pgs.dta"
rename _merge _merge1
*N increased here because of multiple births, presumably (1:m merge)
count
*112,996
merge m:1 M_ID_2306 using "scratch/m_childbmi_genetic_pgs.dta"
rename _merge _merge2
merge m:1 F_ID_2306 using "scratch/f_childbmi_genetic_pgs.dta"
rename _merge _merge3
*for each role, a group of PREG_ID_2306 values not merged from master (i.e., non-genetic SV_INFO file)
*nobody is in the genetic covariate file but not in the SV_INFO file
*of the unmerged, there are the least for mothers, more for children, most for fathers. as expected.
count
**112,996
*how many PREG_ID_2306 values not represented in the genetic data even partally, i.e. under any role?
count if _merge1==1 & _merge2==1 & _merge3==1
*8,150
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
*N=42272

*having checked that, make a single variable for fid, initially equal to that of the child: 
gen fid=c_fid 
count if fid==""
*36,436
*sub in for where child isn't genotyped:
replace fid=m_fid if fid==""
replace fid=f_fid if fid==""
count if fid==""
*8,150
*AS EXPECTED - THESE ARE THE PREGNANCIES FOR WHICH NONE OF THE THREE PEOPLE WERE IN THE GENOTYPE FILE.

*drop these people from the genetic covariate file:
drop if _merge1==1 & _merge2==1 & _merge3==1
count
*N=104,816

*TO SUM UP: ON RHS OF ANALYSIS MODELS (AND HENCE, IMPUTATION MODELS), YOU NEED:
*NOT genotyping center, plate, or imputation release
*DO include genotyping_chip and genotyping_batch_simple
*chip:
fre *_genotyping_chip
*biggest is genotyping_chip==Illumina_GlobalScreeningArray_MD_v.3.0
fre *_genotyping_chip2
*omit baseline group 2 for everyone
*include:
fre c_genotyping_chip1  c_genotyping_chip3 c_genotyping_chip4 c_genotyping_chip5 c_genotyping_chip6 m_genotyping_chip1  m_genotyping_chip3 m_genotyping_chip4 m_genotyping_chip5 m_genotyping_chip6 f_genotyping_chip1  f_genotyping_chip3 f_genotyping_chip4 f_genotyping_chip5 f_genotyping_chip6

*genotyping batch:
fre *_genotyping_batch_num
*children: 
*omit largest = c_batch_norment_aug2020_1029 (also omit: c_batch_norment_jan2015 c_batch_norment_jun2015, because no observations in these groups)
*mothers: 
*omit largest = largest, norment_aug2020_996
*omit largest: 18 norment_may2016

*SAVE IT:
*drop _merge
save "scratch/reshaped_linkable_childbmi_genetic_pgs.dta", replace

************************************************
*MENDELIAN ERRORS CHECK:

use "scratch/reshaped_linkable_childbmi_genetic_pgs.dta", clear

//Check for Mendelian errors
*still have the ordered snplist saved as a global macro:
foreach snp in $childbmi_snplist {
display "`snp'"
}
*infuriatingly, the leading spaces for all but the first make it is impossible to use this in lots of ways, including for mendelian error checks.
*a solution - inelegant but effective - is to recreate a local macro from the contents of the global via a variable, and use that:
gen childbmi_snp=""
local k=0
foreach snp in $childbmi_snplist {
local k=`k'+1
display "`snp'"
replace childbmi_snp="`snp'" in `k'
}

*then make a local macro, and loop through with that.
*save the coefficients in a single new variable to go with the variable "snp", and check if any look like Mendelian errors:
gen m_childbmi_snp_coeff=.
gen f_childbmi_snp_coeff=.
levelsof childbmi_snp, local(levels) clean
local k=0
foreach snp of local levels {
local k=`k'+1
reg c_`snp' f_`snp' m_`snp'
replace m_childbmi_snp_coeff=_b[m_`snp'] in `k'
replace f_childbmi_snp_coeff=_b[f_`snp'] in `k'
}
*then summ these:
summ m_childbmi_snp_coeff, det
summ f_childbmi_snp_coeff, det
*looks ok

*ditto for the scores:
pwcorr c_z_childbmi_pgs m_z_childbmi_pgs  f_z_childbmi_pgs 

*and the PCs:
forvalues i= 1(1)20 {
    regress c_PC`i' m_PC`i' f_PC`i'
}

save "scratch/reshaped_linkable_childbmi_genetic_pgs.dta", replace

