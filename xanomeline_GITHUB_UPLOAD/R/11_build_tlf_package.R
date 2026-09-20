# Step 14 | Additional regulatory-style TLFs (EDUCATIONAL; not submission-ready)
# Run from PROJECT ROOT after R/02_build_adam.R and R/03_ancova_tlf.R.
# Input: official pharmaversesdtm AE plus locally derived ADSL / Week-2 ADLB.
# Individual-subject listing written to outputs/tlf/private ONLY; NEVER commit it.
needed <- c('pharmaversesdtm','dplyr','tidyr','readr')
absent <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(absent)) stop('Install packages first: ',paste(absent,collapse=', '))
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr)})
private_dir <- 'outputs/tlf/private'
public_dir <- 'outputs/tlf/public'
dir.create(private_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(public_dir,recursive=TRUE,showWarnings=FALSE)
required_files <- c('outputs/adam/adsl.rds','outputs/adam/adlb_week2.rds')
if(!all(file.exists(required_files))) stop('Run 02_build_adam.R first; required analysis datasets missing.')
adsl <- readRDS(required_files[1]); week2 <- readRDS(required_files[2])
ae <- pharmaversesdtm::ae
lb <- pharmaversesdtm::lb
arms <- c('Placebo','Xanomeline Low Dose','Xanomeline High Dose')
stopifnot(all(c('USUBJID','TRT01P','SAFFL','AGE','SEX') %in% names(adsl)),
          all(c('USUBJID','TRT01P','BASE','AVAL','CHG','BASE_LBSEQ','W2_LBSEQ','ANL01FL','SAFFL') %in% names(week2)),
          all(c('USUBJID','AESEQ') %in% names(ae)),
          !anyDuplicated(adsl$USUBJID),!anyDuplicated(week2$USUBJID),
          !anyDuplicated(ae[c('USUBJID','AESEQ')]))
safety <- adsl %>% filter(SAFFL=='Y',TRT01P %in% arms) %>%
  mutate(ARM=factor(TRT01P,levels=arms))
stopifnot(nrow(safety)==n_distinct(safety$USUBJID),setequal(as.character(unique(safety$ARM)),arms))
N <- safety %>% count(ARM,name='denominator_n') %>% mutate(ARM=as.character(ARM))
write_csv(N,file.path(private_dir,'table_10_safety_denominators.csv'))

# TABLE 1: demographics. One subject is counted once; only treated assigned arms.
# Missing AGE / SEX are documented, not imputed.
age <- safety %>% group_by(ARM) %>% summarise(
  denominator_n=n(),age_available=sum(!is.na(AGE)),age_missing=sum(is.na(AGE)),
  age_mean=if(all(is.na(AGE))) NA_real_ else mean(AGE,na.rm=TRUE),
  age_sd=if(sum(!is.na(AGE))<2) NA_real_ else sd(AGE,na.rm=TRUE),
  age_median=if(all(is.na(AGE))) NA_real_ else median(AGE,na.rm=TRUE),
  age_min=if(all(is.na(AGE))) NA_real_ else min(AGE,na.rm=TRUE),
  age_max=if(all(is.na(AGE))) NA_real_ else max(AGE,na.rm=TRUE),.groups='drop') %>%
  mutate(ARM=as.character(ARM))
sex <- safety %>% mutate(SEX=if_else(is.na(SEX)|trimws(SEX)=='','MISSING',as.character(SEX))) %>%
  count(ARM,SEX,name='n') %>% mutate(ARM=as.character(ARM)) %>%
  left_join(N,by='ARM') %>% mutate(pct=100*n/denominator_n)
write_csv(age,file.path(private_dir,'table_11_demographics_age.csv'))
write_csv(sex,file.path(private_dir,'table_11_demographics_sex.csv'))

# TABLE 2: incidence of ALL REPORTED AE records in safety subjects.
# This is NOT called TEAE: onset window / treatment-emergence rules have not
# been approved; SAE counts reflect the recorded AESER flag only.
ae_safety <- ae %>% inner_join(safety %>% select(USUBJID,ARM),by='USUBJID')
if(nrow(ae_safety)==0) stop('No AE records in safety population; inspect AE input.')
if(!'AESER' %in% names(ae_safety)) stop('AESER not available; SAE summary requires source review.')
if(any(!is.na(ae_safety$AESER) & !ae_safety$AESER %in% c('Y','N',''))) {
  stop('Unexpected AESER codes: inspect before producing SAE incidence.')
}
ae_summary <- N %>% left_join(
  ae_safety %>% group_by(ARM) %>% summarise(
    ae_records=n(),participants_with_ae=n_distinct(USUBJID),
    serious_ae_records=sum(AESER=='Y',na.rm=TRUE),
    participants_with_serious_ae=n_distinct(USUBJID[AESER=='Y' & !is.na(AESER)]),
    missing_aeser=sum(is.na(AESER)|AESER==''),.groups='drop') %>%
    mutate(ARM=as.character(ARM)),by='ARM') %>%
  mutate(across(c(ae_records,participants_with_ae,serious_ae_records,
                  participants_with_serious_ae,missing_aeser),~replace_na(.x,0L)),
    pct_with_ae=100*participants_with_ae/denominator_n,
    pct_with_serious_ae=100*participants_with_serious_ae/denominator_n)
stopifnot(all(ae_summary$participants_with_ae<=ae_summary$denominator_n),
          all(ae_summary$participants_with_serious_ae<=ae_summary$participants_with_ae))
write_csv(ae_summary,file.path(private_dir,'table_12_all_reported_ae_summary.csv'))
# Preserve coding provenance. If no standardized preferred term exists, use
# verbatim text and label it explicitly; never claim MedDRA coding.
term_var <- if('AEDECOD' %in% names(ae_safety)) 'AEDECOD' else if('AETERM' %in% names(ae_safety)) 'AETERM' else NA_character_
if(is.na(term_var)) stop('Neither AEDECOD nor AETERM is present.')
term_label <- if(term_var=='AEDECOD') 'Source decoded term (AEDECOD)' else 'Source verbatim term (AETERM)'
ae_terms <- ae_safety %>% mutate(AE_TERM=as.character(.data[[term_var]]),
          AE_TERM=if_else(is.na(AE_TERM)|trimws(AE_TERM)=='','MISSING',AE_TERM)) %>%
  group_by(ARM,AE_TERM) %>% summarise(records=n(),participants=n_distinct(USUBJID),.groups='drop') %>%
  mutate(ARM=as.character(ARM),source_term_variable=term_var) %>%
  left_join(N,by='ARM') %>% mutate(pct=100*participants/denominator_n) %>%
  arrange(ARM,desc(participants),AE_TERM)
write_csv(ae_terms,file.path(private_dir,'table_13_ae_terms_LOCAL.csv'))
# Date audit: avoid inventing TEAE status from potentially incomplete onset dates.
onset <- if('AESTDTC' %in% names(ae_safety)) {
 ae_safety %>% mutate(date_status=case_when(
   is.na(AESTDTC)|AESTDTC==''~'MISSING',
   grepl('^[0-9]{4}-[0-9]{2}-[0-9]{2}',AESTDTC)~'FULL_DATE_AVAILABLE',
   TRUE~'PARTIAL_OR_OTHER')) %>% count(date_status,name='AE_records')
} else tibble(date_status='AESTDTC_NOT_IN_SOURCE',AE_records=nrow(ae_safety))
write_csv(onset,file.path(private_dir,'ae_onset_date_audit.csv'))

# LISTING 1: traceable Week-2 ALT observations for the ENTIRE safety population.
# Includes exclusions so that missing observations are visible; subject-level
# output stays local, not in GitHub/Pages or Power BI public materials.
listing <- week2 %>% filter(SAFFL=='Y',TRT01P %in% arms) %>%
  select(USUBJID,TRT01P,TRT01A,SAFFL,BASE,BASE_LBSEQ,BASE_LBDTC,
         AVAL,W2_LBSEQ,W2_LBDTC,AVALU,CHG,ANL01FL,EXCLUSION) %>%
  arrange(TRT01P,USUBJID)
stopifnot(nrow(listing)==nrow(safety),!anyDuplicated(listing$USUBJID),
          setequal(listing$USUBJID,safety$USUBJID))
write_csv(listing,file.path(private_dir,'listing_01_week2_alt_SUBJECT_LEVEL_DO_NOT_UPLOAD.csv'))
listing_summary <- listing %>% group_by(TRT01P,EXCLUSION) %>% summarise(n=n(),.groups='drop')
write_csv(listing_summary,file.path(private_dir,'listing_01_exclusions_AGGREGATE.csv'))

# Explicit, limited public exports; never copy the complete listing / AE rows.
write_csv(age,file.path(public_dir,'demographics_age_AGGREGATE.csv'))
write_csv(sex,file.path(public_dir,'demographics_sex_AGGREGATE.csv'))
write_csv(ae_summary,file.path(public_dir,'all_reported_ae_AGGREGATE.csv'))
# AE-by-term counts are intentionally NOT exported publicly: rare events and
# complementary cells can disclose counts even after simple suppression.
writeLines(c('Educational TLFs only; not submission-ready.',
 'Demographic denominator: safety subjects (SAFFL=Y), grouped by assigned arm TRT01P.',
 'AE table: ALL REPORTED adverse events in safety subjects, regardless of onset timing. NOT a TEAE table.',
 paste('AE term display:',term_label),
 'AE subject incidence: distinct subjects per treatment; AE records and subjects are distinct metrics.',
 'SAE flag: AESER=Y as recorded; missing AESER counted separately.',
 'ALT listing: local subject-level output ONLY; one row per safety subject.',
 'AE-by-term counts remain private pending a disclosure-control review.',
 'No individual data, source ADaM, XPT, or unapproved submission metadata may be committed.'),
 file.path(public_dir,'TLF_METHODS_AND_LIMITATIONS.txt'))
cat('STEP 14 TLF BUILD COMPLETE: demographics, reported AE summary, local subject listing.\n')
cat('Safety N:',nrow(safety),' | AE records:',nrow(ae_safety),' | Listing rows:',nrow(listing),'\n')
cat('Next: source("R/12_qc_tlf_package.R") before treating outputs as validated.\n')
