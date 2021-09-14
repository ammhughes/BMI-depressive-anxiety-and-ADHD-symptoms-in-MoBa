 *OUTCOMES: DEPRESSION, ANXIETY, ADHD AND ASD - SCORES MADE IN STATA

*snps have been clumped in mrbase from the set common to the gwas and the moba_qc_bim file.
*snps extracted from new qc'd dataset with txt in this script: /tsd/p471/data/durable/projects/mh_bmi/scripts/colossus/outcomes_snp_extraction.sh
*(had to paste it into the command line because colossus not accepting the job in the queue)

*************************************************************
*sysdir set PLUS "N:\data\durable\people\Mandy\Statafiles\ado\plus"
cd "N:\data\durable\projects\mh_bmi\"

*FIRST STEP: COMPARE ALIGNMENT IN MOBA DATA VS GWAS AND HARMONIZE.

*UPDATE: need to compare non-effect alleles as well, and allele frequencies, because palindromic snps.
*this info isn't in the dosage (.raw) file, so you need the moba .afreq file (v similar to bim but with a allele frequency).
*made this in plink 22/03/21 and added it to FINAL QC VERSION folder off colossus

*sysdir set PLUS "N:\data\durable\people\Mandy\Statafiles\ado\plus"
cd "N:\data\durable\projects\mh_bmi\"

/*
import delimited "N:/data/durable/projects/moba_interim_release_post_imp_qc/FINAL QC VERSION/mobq_qc_freq.afreq"
*(use from final qc folder, because keeping a copy in the mh_bmi folder takes up all your disk space)
rename id SNP
rename ref ref_allele_moba
rename alt other_allele_moba
rename alt_freqs other_allele_freq_moba
save "data/mobq_qc_freq_afreq.dta", replace
*/

*then, you need the rsid, A1, A2 and beta etc from the correct table of the gwas.

*IMPORT and initial prep for each:

*depression:
import delimited "genetic/clumped_moba_depression_snps_allinfo.txt", case(preserve) clear 
rename effect_alleleexposure effect_allele_gwas
rename other_alleleexposure other_allele_gwas
*rename eafexposure effect_allele_freq_gwas
*the eafexposure variable is fake - just NA values, because the gwas didn't release it
drop eafexposure
rename betaexposure beta_gwas
rename seexposure se_gwas
rename pvalexposure p_gwas
drop idexposure samplesizeexposure ncaseexposure ncontrolexposure mr_keepexposure pval_originexposure unitsexposure_dat data_sourceexposure exposure unitsexposure
save "data/formatted_clumped_moba_depression_snps_allinfo", replace

*anxiety:
*then, you need the rsid, A1, A2 and beta etc from the correct table of the gwas:
*import and prep:
import delimited "genetic/formatted_anxiety_snps.txt", case(preserve) clear 
rename effect_allele effect_allele_gwas
rename other_allele other_allele_gwas
rename eaf effect_allele_freq_gwas
rename betaexposure beta_gwas
*no SE given! that will be a problem down the line for robustness checks - need to write to authors?
*rename seexposure se_gwas
rename pvalexposure p_gwas
save "data/formatted_clumped_moba_anxiety_snps_allinfo", replace

*adhd:
*then, you need the rsid, A1, A2 and beta etc from the correct table of the gwas:
*import and prep:
import delimited "genetic/clumped_moba_adhd_snps_allinfo.txt", case(preserve) clear 
rename effect_allele effect_allele_gwas
rename other_allele other_allele_gwas
*no eaf given!
*rename eaf effect_allele_freq_gwas
rename betaexposure beta_gwas
rename seexposure se_gwas
rename pvalexposure p_gwas
drop idexposure
save "data/formatted_clumped_moba_adhd_snps_allinfo", replace

*asd:
*then, you need the rsid, A1, A2 and beta etc from the correct table of the gwas:
*import and prep:
import delimited "genetic/clumped_moba_asd_snps_allinfo.txt", case(preserve) clear 
rename effect_allele effect_allele_gwas
rename other_allele other_allele_gwas
*no eaf given!
*rename eaf effect_allele_freq_gwas
rename betaexposure beta_gwas
rename seexposure se_gwas
rename pvalexposure p_gwas
drop idexposure
save "data/formatted_clumped_moba_asd_snps_allinfo", replace

**************************************************************************************************************

*HARMONIZATION: LOOP OVER OUTCOMES FOR THIS BIT.
 
*DEPRESSION, ADHD AND ASD GWAS DID NOT GIVE EAF. 
*for anxiety, moot point because not A/T or C/G pairs among the 6 SNPs
*hence, harmonization process for outcomes is modified to not make reference to allele frequencies

foreach outcome in depression anxiety adhd asd {
use "data/formatted_clumped_moba_`outcome'_snps_allinfo", clear
*merge and keep common ones:
merge 1:1 SNP using "data/mobq_qc_freq_afreq.dta"
keep if _merge==3
drop _merge

*define same or reverse strand alignment
capture drop strand_alignment
gen strand_alignment=.
capture label define strand_alignment 1"same" 2"reverse"
label values strand_alignment strand_alignment
replace strand_alignment=1 if effect_allele_gwas==ref_allele_moba & other_allele_gwas==other_allele_moba
replace strand_alignment=2 if other_allele_gwas==ref_allele_moba & effect_allele_gwas==other_allele_moba
tab strand_alignment

*however, some of these could have had alignment misclassified because they are palindromes.
*possible palindromes: where the pair of possible alleles are A/T or C/G. for depression, adhd and asd cannot check using eaf as don't have this info.
list SNP effect_allele_gwas other_allele_gwas ref_allele_moba other_allele_moba strand_alignment if ///
((effect_allele_gwas=="A" & other_allele_gwas=="T") | (effect_allele_gwas=="T" & other_allele_gwas=="A") | (effect_allele_gwas=="C" & other_allele_gwas=="G") | (effect_allele_gwas=="G" & other_allele_gwas=="C")) 

capture drop palindrome
gen palindrome=0
replace palindrome=1 if ((effect_allele_gwas=="A" & other_allele_gwas=="T") | (effect_allele_gwas=="T" & other_allele_gwas=="A") | (effect_allele_gwas=="C" & other_allele_gwas=="G") | (effect_allele_gwas=="G" & other_allele_gwas=="C"))

*make a harmonized version:
gen harm_beta_gwas=beta_gwas
replace harm_beta_gwas=-(harm_beta_gwas) if strand_alignment==2
list harm_beta_gwas beta_gwas strand_alignment 
*remove possible palindromes:
*replace harm_beta_gwas=. if palindrome==1
*DROP THEM: otherwise, missing values of beta mess with the pgs calculation
drop if palindrome==1

sort SNP
save data/moba_gwas_snplevelinfo_`outcome'.dta, replace
}

*possible palindromic SNPs dropped:
*2 for depression
*	SNP	effect~s	other_~s	ref_al~a	o~e_moba	strand~t							
*19.	rs34215985	C	G	G	C	reverse	
*27.	rs62099069	A	T	T	A	reverse	

*1 for adhd
*SNP	effect~s	other_~s	ref_al~a	o~e_moba	strand~t				
*2.	rs11591402	A	T	T	A	reverse	

*1 for asd
*	SNP	effect~s	other_~s	ref_al~a	o~e_moba	strand~t							
*1.	rs10099100	C	G	G	C	reverse	

********************************************************************************************************************************************************************************

/*NEXT (ALMOST CERTAINALY UNECCESSARY) CHECK:
*Check that the ref allele in the raw file (one included in the actual snp variable names) is the same as the ref allele in the .freq file (one to which freq refers, i.e., MAJOR ALLELE, since all frequencies >0.50):

foreach outcome in depression anxiety adhd asd {
import delimited "genetic/`outcome'_snps.raw", varnames(nonames) clear
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
}
*38, 6, 9, 8 SNPs matched so all good.*/

********************************************

*IMPORT IID SNP-LEVEL INFO

foreach outcome in depression anxiety adhd asd {
import delimited "genetic/`outcome'_snps.raw", clear

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

save "scratch/`outcome'_genetic.dta", replace
}

******************************************************************************************************

*NOW MAKE THE PGS.

*'merge' in snplevel info, now JUST THE ONES IN MOBA AS WELL AS GWAS, even though not really a merge
*in the other file, THE DATA NEEDS TO BE SORTED BY SNP SO THAT THIS COMES OUT IN THE RIGHT ORDER.
foreach outcome in depression anxiety adhd asd {
    
use "scratch/`outcome'_genetic.dta"
merge 1:1 _n using data/moba_gwas_snplevelinfo_`outcome'.dta

*save the order of snps in a global for later:
levelsof SNP, local(`outcome'_snps)
display "`outcome'_snps"
local `outcome'_snplist_temp: subinstr local `outcome'_snps "char(34)" "char(32)", all
local `outcome'_snplist_temp2: subinstr local `outcome'_snplist_temp "rs" "  rs", all
global `outcome'_snplist: subinstr local `outcome'_snplist_temp2 "  " ""
display $`outcome'_snplist

*check it works in a loop:
foreach snp in $`outcome'_snplist {
	display "`snp'"
}
/*NB: there is a little bit of SNP-level missingness for some people, enough that you lose about 600 trios
you don't notice this when making scores with PRSise because it does something to fill these in as specified here: 
https://www.prsice.info/step_by_step/#prs-calculation. 
*default is --missing SET_ZERO, which adds nothing where SNPs are missing. do the equivalent of -this manually.
*/
*to quantify snp-level missingness:
capture gen snp_missingness=.
capture gen name_snp_missingness=""
local k=0
foreach var in $`outcome'_snplist {
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
gen prop_snp_missingness=snp_missingness/93111
summ prop_snp_missingness, det

*UPDATE 12 MAY: NEW VERSION OF PGS adding as if for mean genotype for snp-level missigness.
capture drop `outcome'_pgs*
gen `outcome'_pgs=0 
*if batch_num!=.
*nb: keeping as . for people who weren't genotyped (FIX THIS LATER)
local k=0
foreach var in $`outcome'_snplist {
display "`var'"
fre `var'
local k=`k'+1
local relevantbeta=harm_beta_gwas in `k'
display `relevantbeta'
replace `outcome'_pgs=`outcome'_pgs+(`relevantbeta'*0) if `var'==0
replace `outcome'_pgs=`outcome'_pgs+(`relevantbeta'*1) if `var'==1
replace `outcome'_pgs=`outcome'_pgs+(`relevantbeta'*2) if `var'==2
*for people for whom the snp is missing, ADD BETA*MEAN VALUE
sum `var'
return list
replace `outcome'_pgs=`outcome'_pgs+(`relevantbeta'*r(mean)) if `var'==.
}

*standardize the PGS and rename it for consistency:
egen z_`outcome'_pgs=std(`outcome'_pgs)
summ z_`outcome'_pgs

*Add in the new PCs:
capture drop _merge
merge 1:1 iid using "data/final.1kg-pca_updated.eigenvec.dta"
drop _merge
save "scratch/`outcome'_genetic_pgs.dta", replace
}
*******************************************************************************************


*PART 2:
*PREP AND RESHAPE to merge into phenotype data

cd "N:\data\durable\projects\mh_bmi\"
/*
*start with the full linker file:
*keep everything, but rename:
*get pregnancy id, role etc
use "N:\data\durable\projects\moba_interim_release_post_imp_qc\scratch\cleaned_linkage_file",clear
rename SENTRIXID iid
rename Role role
rename BATCH batch
save "data\2019_10_01_MoBaGenetics_all_2306.dta", replace
*/

foreach outcome in depression anxiety adhd asd {
*start with the linker file:
use "data\2019_10_01_MoBaGenetics_all_2306.dta", clear
*merge to the genetic data including the pgs:
merge 1:1 iid using "scratch/`outcome'_genetic_pgs.dta"
keep if _merge==3
drop _merge
save "scratch/linkable_`outcome'_genetic_pgs.dta", replace

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
*93,208
*now do the reshape 

*first drop the snp-level gwas info (keep stuff from linker file incl batch, individual snp cars, PGSs, and PCs:
keep PREG_ID_2306-sex rs* *pgs* PC* 
*drop CHR-harm_effect_allele_freq_gwas
*then these bits:
drop SampleType PDB SampleType Retrieval_ID
//Need create three files, one for offspring, one for mother and one for fathers
preserve
keep if role=="Child"
*do these renames later to make the mendelian error check easier
*renpfix rs c_rs
*rename `outcome'_pgs c_`outcome'_pgs
*rename z_`outcome'_pgs c_z_`outcome'_pgs
rename batch c_batch
renpfix PC c_PC
save "scratch/offspring",replace
restore
preserve
keep if role=="Father"
renpfix rs f_rs
rename `outcome'_pgs f_`outcome'_pgs
rename z_`outcome'_pgs f_z_`outcome'_pgs
rename batch f_batch
renpfix PC f_PC
save "scratch/father",replace
restore
preserve
keep if role=="Mother"
renpfix rs m_rs
rename `outcome'_pgs m_`outcome'_pgs
rename z_`outcome'_pgs m_z_`outcome'_pgs
rename batch m_batch
renpfix PC m_PC
save "scratch/mother",replace
restore

//Next we need to merge each of these files togeather
use  "scratch/offspring",clear
*31,374
duplicates report PREG_ID_2306
*none
*use  "scratch/father",clear
*30,239 
*use  "scratch/mother",clear
*31,595
//Just keep ids, sex, and genetic data
keep PREG_ID_2306 BARN_NR fid iid pid mid sex *rs* *`outcome'* *batch* *PC*
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
pwcorr z_`outcome'_pgs m_z_`outcome'_pgs  f_z_`outcome'_pgs 

*and the PCs:
forvalues i= 1(1)20 {
    regress c_PC`i' m_PC`i' f_PC`i'
}

*then rename the child's genetic vars with equivalent prefix
renpfix rs c_rs
rename `outcome'_pgs c_`outcome'_pgs
rename z_`outcome'_pgs c_z_`outcome'_pgs

*UPDATE 12 MAY: RESTRICT TO PEOPLE IN UPDATED SV_INFO FILE (i.e. exclude the new consent withdrawls) 
merge m:1 PREG_ID_2306 using "N:\data\durable\data\MoBaPhenoData\PDB2306_MoBa_v12\Statafiles\PDB2306_SV_INFO_v12_UPDATED_April_2021"
*the 11 not matched are the ones that need to be dropped as new exclusions.
drop if _merge!=3
drop _merge
count
*31,363

save "scratch/reshaped_linkable_`outcome'_genetic_pgs.dta", replace
}

*suspect there are space issues, so erase intermediate files:
foreach outcome in depression anxiety adhd asd {
*capture erase "scratch/`outcome'_genetic_pgs.dta"
capture erase "scratch/linkable_`outcome'_genetic_pgs.dta"
}

**************************************************************************************

*SANITY CHECK IN UNIMPUTED DATA - DO THEY PREDICT OUTCOMES?
use  "data\statafiles\phenotypes_bmi_analyses.dta", clear
foreach outcome in depression anxiety adhd asd {
merge 1:1 PREG_ID_2306 BARN using "scratch/reshaped_linkable_`outcome'_genetic_pgs.dta"
*no point in keeping ones without genetic info:
keep if _merge==3
drop _merge
}

*delineate the usable trios:
gen fullgeneticinfo=0
replace fullgeneticinfo=1 if (c_batch!="" & m_batch!="" & f_batch!="")
tab fullgeneticinfo
*26,644 trios

*standardize all the outcomes (in imputed version, done with mi passive):
foreach var in MFQ_8yr SCARED_8yr RSDBD_ADHD_8yr RSDBD_inADHD_8yr RSDBD_hyADHD_8yr SCQ_8yr SCQ_SCI_8yr SCQ_RRB_8yr  {
egen z_`var'=std(`var')
}

regress z_MFQ_8yr c_depression_pgs KJONN
regress z_SCARED_8yr c_anxiety_pgs KJONN
regress z_RSDBD_ADHD_8yr c_adhd_pgs KJONN
regress z_SCQ_8yr c_asd_pgs KJONN
*does something for depression and adhd, for anxiety and asd does f-all


