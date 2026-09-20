# 06_build_full_adlb.R
# Educational, multi-parameter, full-visit ADLB BDS candidate.
# Source from project root AFTER running R/02_build_adam.R.
# Never overwrite adlb_week2.rds or the original ANCOVA output.
required <- c('pharmaversesdtm', 'dplyr', 'tidyr', 'readr')
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) install.packages(missing_packages, repos = 'https://cloud.r-project.org')
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr)})

out <- 'outputs/full_adlb'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stop_if <- function(condition, message) if (isTRUE(condition)) stop(message, call. = FALSE)
write_issue <- function(data, name) readr::write_csv(data, file.path(out, paste0(name, '.csv')), na = '')

stop_if(!file.exists('outputs/adam/adsl.rds'), 'ADSL missing: run source("R/02_build_adam.R") first.')
adsl <- readRDS('outputs/adam/adsl.rds')
lb <- pharmaversesdtm::lb
must_lb <- c('STUDYID','USUBJID','LBSEQ','LBTESTCD','LBTEST','LBCAT','LBSTRESN',
             'LBSTRESC','LBSTRESU','LBDTC','LBBLFL','VISIT','VISITNUM',
             'LBSTNRLO','LBSTNRHI')
stop_if(length(setdiff(must_lb,names(lb)))>0,
        paste('LB columns not found:',paste(setdiff(must_lb,names(lb)),collapse=', ')))
stop_if(!all(c('STUDYID','USUBJID','TRT01P','TRT01A','SAFFL','FIRST_DOSE_DATE') %in% names(adsl)),
        'Existing ADSL is missing essential variables; inspect 02_build_adam.R.')
stop_if(anyDuplicated(adsl$USUBJID)>0, 'ADSL does not have one row per USUBJID.')
stop_if(anyNA(lb$USUBJID) || anyNA(lb$LBSEQ), 'LB missing source identifiers: resolve before derivation.')
stop_if(anyDuplicated(lb[c('STUDYID','USUBJID','LBSEQ')])>0, 'Duplicate LB source keys: resolve before derivation.')
stop_if(any(!lb$USUBJID %in% adsl$USUBJID), 'LB has subjects not present in ADSL.')
stop_if(anyNA(lb$LBTESTCD) || any(lb$LBTESTCD==''), 'Unmapped LBTESTCD: resolve before derivation.')

# Validate parameter dictionary and standardized units before generating analysis data.
# NEVER treat values in different units as directly comparable.
unit_audit <- lb %>% filter(is.finite(LBSTRESN)) %>%
  group_by(LBTESTCD) %>% summarise(
    n_units=n_distinct(LBSTRESU,na.rm=TRUE),
    units=paste(sort(unique(LBSTRESU[!is.na(LBSTRESU) & LBSTRESU!=''])),collapse=' | '),
    missing_unit=sum(is.na(LBSTRESU) | LBSTRESU==''),.groups='drop')
write_issue(unit_audit,'unit_audit')
# Study-specific, audited unit exceptions (source audit: 9 LBTESTCDs).
# PH and SPGRAV are dimensionless numeric measurements. Morphology findings and
# urine KETONES/UROBIL 0/1 are treated as coded categories, NOT concentrations.
# The dataset alone does not establish what 0/1 mean clinically.
dimensionless_numeric <- c('PH','SPGRAV')
categorical_codes <- c('ANISO','KETONES','MACROCY','MICROCY','POIKILO','POLYCHR','UROBIL')
allowed_no_unit <- c(dimensionless_numeric,categorical_codes)
unexpected_units <- unit_audit %>%
  filter(n_units > 1 | (missing_unit > 0 & !LBTESTCD %in% allowed_no_unit) |
           (LBTESTCD %in% allowed_no_unit & n_units > 0))
write_issue(unexpected_units,'unit_audit_unresolved')
stop_if(nrow(unexpected_units)>0,
        'Unresolved unit inconsistency: inspect outputs/full_adlb/unit_audit_unresolved.csv.')
# Verify the specific 0/1 source pattern; do not silently generalize this mapping
# if a future version of the sample data contains different codes or actual units.
code_audit <- lb %>% filter(LBTESTCD %in% c('KETONES','UROBIL')) %>%
  count(LBTESTCD,LBSTRESC,LBSTRESN,LBSTRESU,name='records')
write_issue(code_audit,'urinalysis_coded_result_audit')
bad_codes <- code_audit %>% filter(
  is.na(LBSTRESC) | !LBSTRESC %in% c('0','1') |
  is.na(LBSTRESN) | !LBSTRESN %in% c(0,1) |
  (!is.na(LBSTRESU) & LBSTRESU!=''))
stop_if(nrow(bad_codes)>0,
        'Ketones/urobil source codes changed; review urinalysis_coded_result_audit.csv.')
name_audit <- lb %>% distinct(LBTESTCD,LBTEST) %>% count(LBTESTCD,name='n_names')
write_issue(name_audit %>% filter(n_names>1),'parameter_name_conflicts')
stop_if(any(name_audit$n_names>1),
        'One LBTESTCD maps to multiple names: review parameter_name_conflicts.csv.')
stop_if(any(nchar(lb$LBTESTCD)>8),
        'Some LBTESTCD exceed 8 characters: establish documented PARAMCD mapping.')

# Key maps are frozen from available source metadata; parameters with no numeric records
# are included with character analysis values, but have no baseline/change calculation.
param_map <- lb %>% group_by(LBTESTCD) %>%
  summarise(PARAM=first(LBTEST), AVALU=first(LBSTRESU[!is.na(LBSTRESU) & LBSTRESU!=''],default=NA_character_),
            PARCAT1=first(LBCAT,default=NA_character_),.groups='drop') %>%
  arrange(LBTESTCD) %>% mutate(PARAMCD=LBTESTCD, PARAMN=row_number(),
      PARAM=if_else(is.na(AVALU),PARAM,paste0(PARAM,' (',AVALU,')'))) %>%
  select(LBTESTCD,PARAMCD,PARAM,PARAMN,PARCAT1,AVALU)
write_issue(param_map,'parameter_dictionary')

# Preserve every LB source record; do not discard unscheduled, screen failure,
# missing, qualitative, or repeat measurements.
full <- lb %>%
  left_join(adsl %>% select(USUBJID,TRT01P,TRT01A,SAFFL,FIRST_DOSE_DATE),by='USUBJID') %>%
  left_join(param_map,by='LBTESTCD') %>%
  mutate(
    SRCDOM='LB', SRCSEQ=LBSEQ,
    SRCVAR=if_else(LBTESTCD %in% categorical_codes | !is.finite(LBSTRESN),
                   'LBSTRESC','LBSTRESN'),
    # Preserve numeric source codes in the retained LBSTRESN column, but never
    # treat categorical 0/1 or morphology codes as continuous analysis values.
    AVAL=if_else(!LBTESTCD %in% categorical_codes & is.finite(LBSTRESN),
                 as.numeric(LBSTRESN),NA_real_),
    AVALC=if_else(!is.na(LBSTRESC) & LBSTRESC!='',LBSTRESC,NA_character_),
    ADT=suppressWarnings(as.Date(if_else(grepl('^[0-9]{4}-[0-9]{2}-[0-9]{2}',LBDTC),substr(LBDTC,1,10),NA_character_))),
    ADY=if_else(!is.na(ADT) & !is.na(FIRST_DOSE_DATE),
                as.integer(ADT-FIRST_DOSE_DATE)+if_else(ADT>=FIRST_DOSE_DATE,1L,0L),NA_integer_),
    ANRLO=as.numeric(LBSTNRLO),ANRHI=as.numeric(LBSTNRHI),
    ANRIND=case_when(is.na(AVAL)~NA_character_,!is.na(ANRLO)&AVAL<ANRLO~'LOW',
                     !is.na(ANRHI)&AVAL>ANRHI~'HIGH',
                     !is.na(ANRLO)&!is.na(ANRHI)~'NORMAL',TRUE~NA_character_),
    BASETYPE='LAST',
    AVISIT=case_when(LBBLFL=='Y'~'BASELINE',
                     grepl('^WEEK [0-9]+$',VISIT)~VISIT,
                     TRUE~NA_character_),
    AVISITN=case_when(AVISIT=='BASELINE'~0,
                      !is.na(AVISIT)~as.numeric(VISITNUM),TRUE~NA_real_),
    TRTP=TRT01P, TRTA=TRT01A,
    ABLFL=NA_character_,ANL01FL=NA_character_)

# A flagged baseline after first dose is a blocking ambiguity; no invented time ordering.
bad_baseline <- full %>% filter(LBBLFL=='Y',!is.na(ADT),!is.na(FIRST_DOSE_DATE),ADT>FIRST_DOSE_DATE) %>%
  select(STUDYID,USUBJID,LBSEQ,LBTESTCD,LBDTC,FIRST_DOSE_DATE)
write_issue(bad_baseline,'baseline_after_first_dose')
stop_if(nrow(bad_baseline)>0,'Baseline lab samples after first dose: inspect baseline_after_first_dose.csv.')

# Document missing/partial dates instead of inventing analysis dates.
date_issues <- full %>% filter(!is.na(LBDTC),LBDTC!='',is.na(ADT)) %>%
  select(STUDYID,USUBJID,LBSEQ,LBTESTCD,LBDTC)
write_issue(date_issues,'partial_or_invalid_dates')

# Selection rule: among LB source baseline-flagged numeric records, choose latest
# date-time string then highest LBSEQ. Preserve all other LB rows; never invent
# an unflagged baseline. Unique LBSEQ resolves exact timestamp ties reproducibly.
baseline_candidates <- full %>% filter(LBBLFL=='Y',is.finite(AVAL))
chosen_baseline <- baseline_candidates %>%
  arrange(STUDYID,USUBJID,PARAMCD,LBDTC,LBSEQ) %>%
  group_by(STUDYID,USUBJID,PARAMCD) %>% slice_tail(n=1) %>% ungroup() %>%
  select(STUDYID,USUBJID,PARAMCD,SRCSEQ,AVAL,ANRIND) %>%
  rename(BASE_SRCSEQ=SRCSEQ,BASE=AVAL,BNRIND=ANRIND)
write_issue(baseline_candidates %>% count(STUDYID,USUBJID,PARAMCD,name='n_baselines') %>%
              filter(n_baselines>1),'multiple_baseline_candidates')
full <- full %>% left_join(chosen_baseline,by=c('STUDYID','USUBJID','PARAMCD')) %>%
  mutate(ABLFL=if_else(!is.na(BASE_SRCSEQ)&SRCSEQ==BASE_SRCSEQ,'Y',NA_character_),
         CHG=if_else(!is.na(AVAL)&!is.na(BASE)&AVISIT!='BASELINE',AVAL-BASE,NA_real_),
         PCHG=if_else(!is.na(CHG)&!is.na(BASE)&BASE!=0,100*CHG/abs(BASE),NA_real_),
         R2ANRHI=if_else(!is.na(AVAL)&!is.na(ANRHI)&ANRHI>0,AVAL/ANRHI,NA_real_),
         SHIFT1=if_else(!is.na(BNRIND)&!is.na(ANRIND)&AVISIT!='BASELINE',
                        paste(BNRIND,ANRIND,sep=' to '),NA_character_))

# For scheduled post-baseline visits, select latest valid numeric specimen per
# subject/parameter/analysis visit. Unscheduled records remain unflagged;
# baseline source records are indicated by ABLFL, not ANL01FL.
analysis_selected <- full %>%
  filter(!is.na(AVISIT),AVISIT!='BASELINE',is.finite(AVAL),
         is.na(ADT)|is.na(FIRST_DOSE_DATE)|ADT>=FIRST_DOSE_DATE) %>%
  arrange(STUDYID,USUBJID,PARAMCD,AVISITN,LBDTC,LBSEQ) %>%
  group_by(STUDYID,USUBJID,PARAMCD,AVISITN) %>% slice_tail(n=1) %>%
  ungroup() %>% select(STUDYID,USUBJID,PARAMCD,SRCSEQ) %>%
  mutate(ANL01FL='Y')
full <- full %>% left_join(analysis_selected,
    by=c('STUDYID','USUBJID','PARAMCD','SRCSEQ'),suffix=c('','.selected')) %>%
  mutate(ANL01FL=ANL01FL.selected) %>% select(-ANL01FL.selected)

# Standard BDS variables up front, source columns retained at end for inspection.
front <- c('STUDYID','USUBJID','PARAMCD','PARAM','PARAMN','PARCAT1','AVAL','AVALC','AVALU',
           'ADT','ADY','AVISIT','AVISITN','BASETYPE','ABLFL','BASE','CHG','PCHG',
           'ANL01FL','ANRLO','ANRHI','ANRIND','BNRIND','SHIFT1','R2ANRHI',
           'TRTP','TRTA','SAFFL','SRCDOM','SRCSEQ','SRCVAR','BASE_SRCSEQ')
full <- full %>% select(all_of(front),everything())
stop_if(nrow(full)!=nrow(lb),'Row-preservation check failed.')
stop_if(anyDuplicated(full[c('STUDYID','USUBJID','SRCSEQ')])>0,'Source row key became nonunique.')
stop_if(anyDuplicated(full %>% filter(ABLFL=='Y') %>% select(STUDYID,USUBJID,PARAMCD))>0,
        'More than one selected baseline per subject/parameter.')
stop_if(anyDuplicated(full %>% filter(ANL01FL=='Y') %>% select(STUDYID,USUBJID,PARAMCD,AVISITN))>0,
        'More than one analysis record per subject/parameter/visit.')
stop_if(any(full$ABLFL=='Y' & is.na(full$BASE),na.rm=TRUE),'Selected baseline missing BASE.')

write_issue(full %>% count(PARAMCD,PARAM,PARCAT1,AVALU,name='records') %>%
              left_join(full %>% filter(ANL01FL=='Y') %>% count(PARAMCD,name='selected_records'),by='PARAMCD'),
            'parameter_inventory')
write_issue(full %>% count(PARAMCD,AVISIT,AVISITN,name='records'), 'visit_parameter_coverage')
write_issue(full %>% summarise(source_lb_rows=n(),subjects=n_distinct(USUBJID),
              parameters=n_distinct(PARAMCD),numeric_records=sum(!is.na(AVAL)),
              baseline_flags=sum(ABLFL=='Y',na.rm=TRUE),
              analysis_flags=sum(ANL01FL=='Y',na.rm=TRUE)), 'build_summary')
write_issue(full %>% filter(!is.na(AVISIT),AVISIT!='BASELINE',
                            !is.na(ADT),!is.na(FIRST_DOSE_DATE),ADT<FIRST_DOSE_DATE) %>%
              select(STUDYID,USUBJID,SRCSEQ,PARAMCD,LBDTC,AVISIT,ADY),
            'postbaseline_predose_review')
saveRDS(full,file.path(out,'adlb_full.rds'))
write_issue(full,'adlb_full')

# Regression guard: compare the observed ALT Week 2 records independently of
# baseline measurements. Legacy extract has ONE ROW PER SAFETY SUBJECT, including
# subjects without a Week 2 LB record; full ADLB has ONE ROW PER SOURCE LB RECORD.
# Absence of a Week 2 source row must not be interpreted as a missing baseline.
old_path <- 'outputs/adam/adlb_week2.rds'
if (file.exists(old_path)) {
  old <- readRDS(old_path)
  stop_if(anyDuplicated(old$USUBJID)>0, 'Legacy Week 2 extract has duplicate subjects.')
  stop_if(!all(c('USUBJID','AVAL','BASE','ANL01FL') %in% names(old)),
          'Legacy Week 2 extract lacks comparison variables.')
  old <- old %>% filter(SAFFL=='Y')
  new_week2 <- full %>%
    filter(PARAMCD=='ALT', AVISIT=='WEEK 2', ANL01FL=='Y', SAFFL=='Y') %>%
    select(USUBJID, NEW_AVAL=AVAL, NEW_BASE=BASE)
  stop_if(anyDuplicated(new_week2$USUBJID)>0,
          'Full ADLB has duplicate selected Week 2 ALT records.')
  # Chosen baseline exists independently of a Week 2 visit (even when it is absent).
  new_baseline <- chosen_baseline %>% filter(PARAMCD=='ALT') %>%
    select(USUBJID, BASE_FROM_LB=BASE)
  stop_if(anyDuplicated(new_baseline$USUBJID)>0,
          'Full ADLB has duplicate chosen ALT baselines.')
  check <- old %>%
    select(USUBJID, AVAL, BASE, ANL01FL) %>%
    left_join(new_week2,by='USUBJID') %>%
    left_join(new_baseline,by='USUBJID') %>%
    mutate(
      aval_match=(is.na(AVAL)&is.na(NEW_AVAL))|
        (!is.na(AVAL)&!is.na(NEW_AVAL)&abs(AVAL-NEW_AVAL)<1e-9),
      base_match=(is.na(BASE)&is.na(BASE_FROM_LB))|
        (!is.na(BASE)&!is.na(BASE_FROM_LB)&abs(BASE-BASE_FROM_LB)<1e-9),
      comparison_status=case_when(
        !aval_match | !base_match ~ 'TRUE_NUMERIC_DISCREPANCY',
        is.na(AVAL) ~ 'NO_WEEK2_SOURCE_ROW_EXPECTED',
        TRUE ~ 'MATCHED_OBSERVED_WEEK2'))
  # A new selected Week 2 ALT observation without a legacy observation is a
  # separate population/flag issue, not something an old-side left join sees.
  extra_new <- new_week2 %>% anti_join(old,by='USUBJID')
  write_issue(extra_new,'legacy_alt_extra_new_week2')
  write_issue(check %>% filter(comparison_status=='TRUE_NUMERIC_DISCREPANCY'),
              'legacy_alt_discrepancies')
  write_issue(check %>% filter(comparison_status=='NO_WEEK2_SOURCE_ROW_EXPECTED'),
              'legacy_alt_expected_missing_week2')
  write_issue(check %>% count(comparison_status,name='subjects'),
              'legacy_alt_reconciliation_summary')
  stop_if(nrow(check)!=nrow(old) || any(!check$aval_match | !check$base_match) ||
            nrow(extra_new)>0,
          paste('True ALT reconciliation discrepancy: inspect',
                'legacy_alt_discrepancies.csv and legacy_alt_extra_new_week2.csv.'))
  cat('PASS: legacy ALT reconciled across ',nrow(old),' safety subjects; ',
      sum(!is.na(old$AVAL)),' observed Week 2 results and ',
      sum(is.na(old$AVAL)),' correctly missing Week 2 results.\n',sep='')
} else cat('NOTE: legacy Week 2 extract not found; cross-version numerical QC not run.\n')
cat('FULL ADLB BUILD COMPLETE (candidate only). Check outputs/full_adlb/build_summary.csv and QA gate.\n')
