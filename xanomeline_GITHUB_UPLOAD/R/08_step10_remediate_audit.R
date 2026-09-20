# Step 10: non-destructive remediation, source replay, and review exports.
# Run from the SAME project root as R/06_build_full_adlb.R.
# Does not modify adlb_full.rds, adlb_week2.rds, or any ANCOVA output.
required <- c('dplyr','readr','pharmaversesdtm')
miss <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(miss)) stop('Install packages first: ',paste(miss,collapse=', '))
suppressPackageStartupMessages({library(dplyr);library(readr)})
source_dir <- 'outputs/full_adlb'
stopifnot(file.exists(file.path(source_dir,'adlb_full.rds')))
x <- readRDS(file.path(source_dir,'adlb_full.rds'))
out <- 'outputs/step10'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_out <- function(x,name) readr::write_csv(x,file.path(out,paste0(name,'.csv')),na='')
stopifnot(nrow(x)>0, all(c('STUDYID','USUBJID','SRCSEQ','BASE_SRCSEQ','FIRST_DOSE_DATE',
  'PARAMCD','ABLFL','ANL01FL','AVISIT','BASE','AVAL','CHG','LBDTC','LBSTRESC','LBSTRESN') %in% names(x)))

# SOURCE REPLAY: compare independent Pharmaverse LB with mapped source keys and source values.
# Independent source input, but NOT a separately programmed complete ADaM derivation.
lb <- pharmaversesdtm::lb
stopifnot(!anyDuplicated(lb[c('STUDYID','USUBJID','LBSEQ')]),
          !anyDuplicated(x[c('STUDYID','USUBJID','SRCSEQ')]))
source_keys <- lb %>% transmute(STUDYID,USUBJID,SRCSEQ=LBSEQ,
  src_test=LBTESTCD, src_value=LBSTRESN, src_char=LBSTRESC,
  src_unit=LBSTRESU, src_date=LBDTC)
replay <- x %>% select(STUDYID,USUBJID,SRCSEQ,PARAMCD,LBSTRESN,LBSTRESC,LBSTRESU,LBDTC) %>%
 full_join(source_keys,by=c('STUDYID','USUBJID','SRCSEQ')) %>%
 mutate(issue=case_when(is.na(src_test)~'ADLB_KEY_NOT_IN_LB',is.na(PARAMCD)~'LB_KEY_NOT_IN_ADLB',
    PARAMCD!=src_test~'PARAMETER_MISMATCH',
    is.na(LBSTRESN)!=is.na(src_value)~'NUMERIC_MISSINGNESS_MISMATCH',
    !is.na(LBSTRESN)&!is.na(src_value)&abs(LBSTRESN-src_value)>1e-9~'NUMERIC_VALUE_MISMATCH',
    coalesce(LBSTRESC,'<NA>')!=coalesce(src_char,'<NA>')~'CHAR_VALUE_MISMATCH',
    coalesce(LBSTRESU,'<NA>')!=coalesce(src_unit,'<NA>')~'UNIT_MISMATCH',
    coalesce(LBDTC,'<NA>')!=coalesce(src_date,'<NA>')~'DATE_MISMATCH',TRUE~'MATCH'))
write_out(replay %>% count(issue,name='records'),'source_replay_summary')
# Only discrepant rows: local QC only; never upload to public website.
write_out(replay %>% filter(issue!='MATCH'),'source_replay_discrepancies_PRIVATE')
if(any(replay$issue!='MATCH')) stop('Source replay mismatch: review outputs/step10/source_replay_discrepancies_PRIVATE.csv')

# Document missing baseline explicitly rather than assigning/imputing one.
missing_baseline <- x %>% filter(ANL01FL=='Y',is.finite(AVAL),is.na(BASE))
write_out(missing_baseline %>% count(PARAMCD,PARAM,name='selected_rows_without_base') %>%
 arrange(desc(selected_rows_without_base)),'baseline_missing_by_parameter')
write_out(missing_baseline %>% distinct(STUDYID,USUBJID,PARAMCD) %>%
 count(PARAMCD,name='subject_parameter_combinations'),'baseline_missing_subject_parameter')

# Preserve nonassigned analysis visits; count actual source visit descriptions.
unmapped <- x %>% filter(is.na(AVISIT))
write_out(unmapped %>% count(VISIT,VISITNUM,name='records') %>%
 arrange(desc(records)),'unmapped_visits_by_source_visit')

# Non-destructive, 8-character-name-friendly analysis-only CSV candidate.
# Deliberately omit source LB columns duplicated from SDTM, BASE_SRCSEQ and
# FIRST_DOSE_DATE. Trace chosen baseline with ABLFL and SRCSEQ inside dataset;
# original extended source linkage stays in internal adlb_full.rds.
columns <- c('STUDYID','USUBJID','PARAMCD','PARAM','PARAMN','PARCAT1',
 'AVAL','AVALC','AVALU','ADT','ADY','AVISIT','AVISITN','BASETYPE',
 'ABLFL','BASE','CHG','PCHG','ANL01FL','ANRLO','ANRHI','ANRIND',
 'BNRIND','SHIFT1','R2ANRHI','TRTP','TRTA','SAFFL','SRCDOM','SRCSEQ','SRCVAR')
stopifnot(all(columns %in% names(x)),all(nchar(columns)<=8),!anyDuplicated(columns))
candidate <- x %>% select(all_of(columns))
stopifnot(nrow(candidate)==nrow(x),
  !anyDuplicated(candidate[c('STUDYID','USUBJID','SRCSEQ')]))
write_out(candidate,'adlb_exchange_CANDIDATE_NOT_VALIDATED')
# High-level reviewer-facing aggregate evidence.
summary <- tibble::tibble(metric=c('source_records','participants','parameters',
 'selected_rows_without_baseline','unassigned_analysis_visits','source_replay_mismatches',
 'exchange_candidate_columns','exchange_candidate_names_over_8'),
 value=c(nrow(x),n_distinct(x$USUBJID),n_distinct(x$PARAMCD),nrow(missing_baseline),
 nrow(unmapped),sum(replay$issue!='MATCH'),ncol(candidate),sum(nchar(names(candidate))>8)))
write_out(summary,'step10_summary')
cat('STEP10 AUDIT COMPLETE. Full ADLB and ANCOVA unchanged.\n')
print(summary,n=nrow(summary))
cat('The output CSV is a candidate only: metadata approval, derivation review, XPT/Define-XML and formal validator remain OPEN.\n')
