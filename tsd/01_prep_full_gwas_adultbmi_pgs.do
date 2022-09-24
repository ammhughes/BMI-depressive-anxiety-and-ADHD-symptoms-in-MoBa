

*IMPORT AND RESHAPE THE FULL GENOME PGS FOR ADULT AND CHILD BMI.

clear all
set maxvar 20000
cd "N:\durable\projects\mh_bmi\fullsample\"
sysdir set PLUS "N:\durable\people\Mandy\Statafiles\ado\plus"

*ADULTBMI:

*import the full score:
import delimited "data\genetic\noclump_full_pgs_adultbmi\full_pgs_adultbmi.all_score", delimiter(space, collapse) clear
save "data\genetic\noclump_full_pgs_adultbmi\full_pgs_adultbmi.all_score.dta", replace

*import delimited "data\genetic\full_gwas_pgs\full_pgs_adultbmi.all_score", delimiter(space, collapse) clear
*save "data\genetic\full_pgs_adultbmi.all_score.dta", replace

*like the other PGS, need to reshape via the covariate file.
use "data\genetic\noclump_full_pgs_adultbmi\full_pgs_adultbmi.all_score.dta", clear
*use "data\genetic\full_pgs_adultbmi.all_score.dta", clear

capture drop _merge
merge 1:1 iid using "data\genetic\MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.dta"
drop _merge

preserve
keep if role=="Child"
*sort out linking id:
list id_2306 
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
*then restrict and rename:
keep PREG_ID_2306 BARN_NR fid iid pt_1 pt_5e08
rename fid c_fid
rename iid c_iid 
rename pt_1 c_full_pgs_adultbmi
rename pt_5e08 c_full_pgs_adultbmi_5e08
save "scratch/c_full_pgs_adultbmi.dta", replace
restore

preserve
keep if role=="Mother"
keep id_2306 fid iid role pt_1 pt_5e08
rename id_2306 M_ID_2306
rename fid m_fid
rename iid m_iid 
rename pt_1 m_full_pgs_adultbmi
rename pt_5e08 m_full_pgs_adultbmi_5e08
save "scratch/m_full_pgs_adultbmi.dta", replace
restore

preserve
keep if role=="Father"
keep iid fid id_2306 pt_1 pt_5e08
rename iid f_iid 
rename fid f_fid
rename id_2306 F_ID_2306
rename pt_1 f_full_pgs_adultbmi
rename pt_5e08 f_full_pgs_adultbmi_5e08
save "scratch/f_full_pgs_adultbmi.dta", replace
restore


*now attach all three to the SV_INFO file,
*then merge these all together using the SV_INFO file.
use "N:\durable\projects\mh_bmi\fullsample\PDB2306_SV_INFO_v12.dta", clear
count
*112,198
merge 1:m PREG_ID_2306 using "scratch/c_full_pgs_adultbmi.dta"
rename _merge _merge1
*N increased here because of multiple births, presumably (1:m merge)
count
*112,949
merge m:1 M_ID_2306 using "scratch/m_full_pgs_adultbmi.dta"
rename _merge _merge2
*check duplicates:
duplicates tag PREG_ID_2306 BARN_NR, gen (m_dup)
list if m_dup!=0
*42 not matched because missing PREG_ID_2306 BARN_NR - not inc child's part of the file

merge m:1 F_ID_2306 using "scratch/f_full_pgs_adultbmi.dta"
rename _merge _merge3
*check duplicates:
duplicates tag PREG_ID_2306 BARN_NR, gen (f_dup)
list if f_dup!=0

*for each role, a group of PREG_ID_2306 values not merged from master (i.e., non-genetic SV_INFO file)
*nobody is in the genetic covariate file but not in the SV_INFO file
*of the unmerged, there are the least for mothers, more for children, most for fathers. as expected.

count
**113,024
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
*N=42272

*having checked that, make a single variable for fid, initially equal to that of the child: 
gen fid=c_fid 
count if fid==""
*36,494
*sub in for where child isn't genotyped:
replace fid=m_fid if fid==""
replace fid=f_fid if fid==""
count if fid==""
*8,141
*AS EXPECTED - THESE ARE THE PREGNANCIES FOR WHICH NONE OF THE THREE PEOPLE WERE IN THE GENOTYPE FILE.

drop if _merge1==1 & _merge2==1 & _merge3==1
count
*104,883

*restrict to full trios, otherwise problems merging:
keep if complete_trio==1 
*42,238


*standaridize the scores.
foreach var in c_full_pgs_adultbmi m_full_pgs_adultbmi f_full_pgs_adultbmi c_full_pgs_adultbmi_5e08 m_full_pgs_adultbmi_5e08 f_full_pgs_adultbmi_5e08 {
	egen z_`var'=std(`var')
	summ z_`var'
}

save "scratch/reshaped_linkable_full_pgs_adultbmi.dta", replace


