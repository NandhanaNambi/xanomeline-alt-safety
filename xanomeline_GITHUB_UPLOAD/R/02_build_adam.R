# 02_build_adam.R | derived analysis datasets and QC
# Run from the project root. Does not calculate treatment effects.
required <- c('pharmaversesdtm','dplyr','tidyr','readr')
miss <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(miss)) install.packages(miss, repos='https://cloud.r-project.org')
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr)})
dir.create('outputs/adam',recursive=TRUE,showWarnings=FALSE)
dir.create('outputs/qc',recursive=TRUE,showWarnings=FALSE)
DM <- pharmaversesdtm::dm; LB <- pharmaversesdtm::lb; EX <- pharmaversesdtm::ex
arms <- c('Placebo','Xanomeline Low Dose','Xanomeline High Dose')
stopifnot(!anyDuplicated(DM$USUBJID), all(c('USUBJID','ARM','ACTARM','RFXSTDTC') %in% names(DM)))
stopifnot(all(c('USUBJID','LBTESTCD','LBSTRESN','LBSTRESU','LBBLFL','VISITNUM','VISIT','LBSEQ','LBDTC') %in% names(LB)))
# ADSL: one row/subject. SAFFL is provisional evidence of receiving any study treatment;
# check discrepancy with EX before finalizing.
adsl <- DM %>% mutate(TRT01P=ARM,TRT01A=ACTARM,
  SAFFL=if_else(!is.na(RFXSTDTC) & RFXSTDTC!='' & ARM %in% arms,'Y','N'),
  FIRST_DOSE_DATE=as.Date(substr(RFXSTDTC,1,10)))
ex_subjects <- EX %>% distinct(USUBJID) %>% mutate(HAS_EX='Y')
adsl <- adsl %>% left_join(ex_subjects,by='USUBJID') %>%
  mutate(HAS_EX=replace_na(HAS_EX,'N'))
write_csv(adsl %>% count(TRT01P,TRT01A,SAFFL,HAS_EX,name='n'), 'outputs/qc/treatment_population_reconciliation.csv')
write_csv(adsl %>% filter((SAFFL=='Y' & HAS_EX=='N') | (SAFFL=='N' & HAS_EX=='Y')) %>%
  select(USUBJID,TRT01P,TRT01A,SAFFL,HAS_EX,RFXSTDTC), 'outputs/qc/dose_evidence_discrepancies.csv')
# Protocol-defined analysis here: randomized/assigned treatment, treated patients.
safety <- adsl %>% filter(SAFFL=='Y', TRT01P %in% arms)
if(!setequal(unique(safety$TRT01P),arms)) stop('Expected three assigned arms in safety population.')
alt <- LB %>% filter(LBTESTCD=='ALT') %>% inner_join(
  safety %>% select(USUBJID,TRT01P,TRT01A,FIRST_DOSE_DATE),by='USUBJID') %>%
  mutate(LBDATE=as.Date(substr(LBDTC,1,10)))
units <- unique(alt$LBSTRESU[is.finite(alt$LBSTRESN)])
if(length(units)!=1L || is.na(units) || !units %in% c('U/L','IU/L'))
  stop('Inspect ALT units before proceeding: ',paste(units,collapse=', '))
# Flag duplicate baseline and Week 2 values explicitly and stop; no arbitrary picking.
base <- alt %>% filter(LBBLFL=='Y',is.finite(LBSTRESN))
w2 <- alt %>% filter(VISITNUM==4, VISIT=='WEEK 2',is.finite(LBSTRESN),
                     is.na(LBBLFL) | LBBLFL!='Y')
write_csv(base %>% count(USUBJID,name='n') %>% filter(n!=1), 'outputs/qc/duplicate_baseline.csv')
write_csv(w2 %>% count(USUBJID,name='n') %>% filter(n!=1), 'outputs/qc/duplicate_week2.csv')
write_csv(alt %>% filter(VISITNUM==4) %>% count(VISIT,LBSTRESU,LBBLFL,name='n'),
          'outputs/qc/week2_record_dictionary.csv')
if(anyDuplicated(base$USUBJID) || anyDuplicated(w2$USUBJID))
  stop('Duplicate records present: inspect outputs/qc/duplicate_*.csv and specify a selection rule.')
if(any(!is.na(base$LBDATE) & !is.na(base$FIRST_DOSE_DATE) & base$LBDATE>base$FIRST_DOSE_DATE))
  stop('Baseline after first dose: inspect source dates before proceeding.')
if(any(!is.na(w2$LBDATE) & !is.na(w2$FIRST_DOSE_DATE) & w2$LBDATE<w2$FIRST_DOSE_DATE))
  stop('Week 2 sample predates first dose: inspect visit/date coding.')
# Re-evaluate coverage by ASSIGNED arm; earlier audit grouped by actual arm.
pop_n <- safety %>% count(TRT01P,name='denominator_n')
coverage <- pop_n %>% left_join(base %>% count(TRT01P,name='baseline_n'),by='TRT01P') %>%
  left_join(w2 %>% count(TRT01P,name='week2_n'),by='TRT01P') %>%
  mutate(across(c(baseline_n,week2_n),~replace_na(.x,0L)),
         baseline_pct=100*baseline_n/denominator_n,
         week2_pct=100*week2_n/denominator_n)
write_csv(coverage,'outputs/qc/assigned_arm_coverage.csv')
print(coverage)
if(any(coverage$baseline_pct<80 | coverage$week2_pct<80))
  stop('Pre-specified 80% per-arm Week 2/baseline threshold fails by assigned arm. Review coverage before analysis.')
# Source traceability via LBSEQ, original date/visit retained. ANALYSIS flag marks
# evaluable baseline + Week2 pairs only, with no missing-data imputation.
base1 <- base %>% transmute(USUBJID,BASE=LBSTRESN,BASE_LBSEQ=LBSEQ,BASE_LBDTC=LBDTC)
w21 <- w2 %>% transmute(USUBJID,AVAL=LBSTRESN,W2_LBSEQ=LBSEQ,W2_LBDTC=LBDTC)
adlb_week2 <- safety %>% select(USUBJID,TRT01P,TRT01A,SAFFL,AGE,SEX) %>%
  left_join(base1,by='USUBJID') %>% left_join(w21,by='USUBJID') %>%
  mutate(PARAMCD='ALT',PARAM='Alanine aminotransferase',AVISIT='WEEK 2',AVISITN=4,
         AVALU=units[1],CHG=AVAL-BASE,
         ANL01FL=if_else(is.finite(BASE)&is.finite(AVAL),'Y','N'),
         EXCLUSION=case_when(is.na(BASE)&is.na(AVAL)~'Missing baseline and Week 2',
                             is.na(BASE)~'Missing baseline',
                             is.na(AVAL)~'Missing Week 2',TRUE~'Included'))
write_csv(adlb_week2 %>% count(TRT01P,EXCLUSION,name='n'), 'outputs/qc/analysis_exclusions.csv')
stopifnot(!anyDuplicated(adsl$USUBJID),!anyDuplicated(adlb_week2$USUBJID),
          all(adlb_week2$USUBJID %in% adsl$USUBJID))
saveRDS(adsl,'outputs/adam/adsl.rds');saveRDS(adlb_week2,'outputs/adam/adlb_week2.rds')
write_csv(adlb_week2,'outputs/adam/adlb_week2.csv')
cat('\nADSL/ADLB created; see outputs/adam and outputs/qc.\n')
