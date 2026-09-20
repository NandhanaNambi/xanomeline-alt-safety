# Step 14 | Independently check the new TLFs against SOURCE domain data.
# This is focused QC, not formal submission validation / fully independent ADaM re-derivation.
needed<-c('pharmaversesdtm','dplyr','readr')
absent<-needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(absent)) stop('Install packages: ',paste(absent,collapse=', '))
suppressPackageStartupMessages({library(dplyr);library(readr)})
local <- 'outputs/tlf/private'; pub<-'outputs/tlf/public'
files<-file.path(local,c('table_10_safety_denominators.csv','table_11_demographics_age.csv',
 'table_11_demographics_sex.csv','table_12_all_reported_ae_summary.csv',
 'table_13_ae_terms_LOCAL.csv','listing_01_week2_alt_SUBJECT_LEVEL_DO_NOT_UPLOAD.csv'))
if(!all(file.exists(files))) stop('Run R/11_build_tlf_package.R first.')
f<-lapply(files,readr::read_csv,show_col_types=FALSE)
n<-f[[1]]; age<-f[[2]]; sex<-f[[3]]; ae_summary<-f[[4]]; terms<-f[[5]]; lst<-f[[6]]
dm<-pharmaversesdtm::dm; ae<-pharmaversesdtm::ae; lb<-pharmaversesdtm::lb
adsl<-readRDS('outputs/adam/adsl.rds')
arms<-c('Placebo','Xanomeline Low Dose','Xanomeline High Dose')
s<-adsl[adsl$SAFFL=='Y' & adsl$TRT01P %in% arms,,drop=FALSE]
source_s<-dm[dm$USUBJID %in% s$USUBJID,,drop=FALSE]
checks<-data.frame(check=character(),pass=logical(),details=character())
add<-function(name,pass,detail='') {checks <<- rbind(checks,data.frame(check=name,pass=isTRUE(pass),details=detail))}
add('One ADSL subject per safety participant',!anyDuplicated(s$USUBJID))
add('All listing subjects match ADSL safety population',nrow(lst)==nrow(s)&&!anyDuplicated(lst$USUBJID)&&setequal(lst$USUBJID,s$USUBJID))
add('Demographic denominator matches DM subject count by assigned ARM',
 isTRUE(all.equal(unname(n$denominator_n[match(arms,n$ARM)]),
  unname(as.integer(table(factor(source_s$ARM,levels=arms)))),check.attributes=FALSE)))
add('Age availability + missingness equal arm denominator',all(age$age_available+age$age_missing==age$denominator_n))
for(a in arms){
 expected<-s$AGE[s$TRT01P==a]; observed<-age$age_mean[age$ARM==a]
 add(paste('Age mean directly recalculated for',a),
  length(observed)==1L && (if(all(is.na(expected))) is.na(observed) else abs(mean(expected,na.rm=TRUE)-observed)<1e-9))
 add(paste('Sex counts sum to denominator for',a),
  sum(sex$n[sex$ARM==a])==n$denominator_n[n$ARM==a])
}
# Independently recompute sex counts by assigned arm from source DM.
for(a in arms){
  src_sex<-source_s$SEX[source_s$ARM==a]
  src_sex[is.na(src_sex)|trimws(as.character(src_sex))=='']<-'MISSING'
  expected<-table(as.character(src_sex))
  got<-sex[sex$ARM==a,,drop=FALSE]
  add(paste('Sex categories independently match source DM for',a),
    nrow(got)==length(expected) &&
    setequal(got$SEX,names(expected)) &&
    all(got$n[match(names(expected),got$SEX)]==as.integer(expected)))
}
# Independent AE count: base R split / unique against original AE rather than TLF summaries.
ae_src<-ae[ae$USUBJID %in% s$USUBJID,,drop=FALSE]
for(a in arms){
 ids<-s$USUBJID[s$TRT01P==a]
 x<-ae_src[ae_src$USUBJID %in% ids,,drop=FALSE]
 observed<-ae_summary[ae_summary$ARM==a,,drop=FALSE]
 add(paste('AE record count source replay',a),nrow(observed)==1L&&observed$ae_records==nrow(x))
 add(paste('AE subject incidence source replay',a),
  nrow(observed)==1L&&observed$participants_with_ae==length(unique(x$USUBJID)))
 add(paste('Serious AE incidence source replay',a),
  nrow(observed)==1L&&observed$participants_with_serious_ae==length(unique(x$USUBJID[!is.na(x$AESER)&x$AESER=='Y'])))
 add(paste('Term counts recover all AE records',a),sum(terms$records[terms$ARM==a])==nrow(x))
}
# Exact source-value replay using recorded LBSEQ rather than comparing the new
# listing solely to the earlier derived ADLB.
alt<-lb[lb$LBTESTCD=='ALT',c('USUBJID','LBSEQ','LBSTRESN'),drop=FALSE]
stopifnot(!anyDuplicated(alt[c('USUBJID','LBSEQ')]))
key<-function(id,seq) paste(id,seq,sep='|')
base_ix<-match(key(lst$USUBJID,lst$BASE_LBSEQ),key(alt$USUBJID,alt$LBSEQ))
w2_ix<-match(key(lst$USUBJID,lst$W2_LBSEQ),key(alt$USUBJID,alt$LBSEQ))
match_num<-function(a,b) all((is.na(a)&is.na(b))|(!is.na(a)&!is.na(b)&abs(a-b)<1e-9))
add('Listing baseline values independently match LBSEQ source',
 match_num(lst$BASE,alt$LBSTRESN[base_ix]))
add('Listing Week-2 values independently match LBSEQ source',
 match_num(lst$AVAL,alt$LBSTRESN[w2_ix]))
add('Listing change is correctly derived',match_num(lst$CHG,lst$AVAL-lst$BASE))
add('Listing evaluable flag matches paired ALT',
 all((lst$ANL01FL=='Y')==(!is.na(lst$BASE)&!is.na(lst$AVAL))))
# Only compare the public aggregate copies; private rows never exported by QC.
for(pair in list(c('table_11_demographics_age.csv','demographics_age_AGGREGATE.csv'),
 c('table_11_demographics_sex.csv','demographics_sex_AGGREGATE.csv'),
 c('table_12_all_reported_ae_summary.csv','all_reported_ae_AGGREGATE.csv'))){
 a<-readLines(file.path(local,pair[1]),warn=FALSE)
 b<-readLines(file.path(pub,pair[2]),warn=FALSE)
 add(paste('Public aggregate matches local',pair[2]),identical(a,b))
}
write_csv(checks,file.path(local,'tlf_qc_checks.csv'))
print(checks,n=nrow(checks))
if(any(!checks$pass)) stop('TLF QC FAILED: inspect outputs/tlf/private/tlf_qc_checks.csv')
cat('STEP 14 TLF QC PASSED:',nrow(checks),'checks; outputs/tlf/private remains LOCAL.\n')
