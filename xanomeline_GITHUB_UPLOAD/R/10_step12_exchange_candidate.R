# STEP 12: non-destructive XPT v5 exchange candidate + strict read-back.
# Run from the existing project root, after Steps 10 and 11.
# This DOES NOT establish submission readiness or create Define-XML.
needed <- c('dplyr','readr','haven','tibble')
missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_pkgs)) stop('Install missing packages: install.packages(c(',paste(sprintf('"%s"',missing_pkgs),collapse=','),'))')
suppressPackageStartupMessages({library(dplyr);library(readr)})
input <- 'outputs/full_adlb/adlb_full.rds'
stopifnot(file.exists(input))
x <- readRDS(input)
out <- 'outputs/step12'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
cols <- c('STUDYID','USUBJID','PARAMCD','PARAM','PARAMN','PARCAT1',
 'AVAL','AVALC','AVALU','ADT','ADY','AVISIT','AVISITN','BASETYPE',
 'ABLFL','BASE','CHG','PCHG','ANL01FL','ANRLO','ANRHI','ANRIND',
 'BNRIND','SHIFT1','R2ANRHI','TRTP','TRTA','SAFFL','SRCDOM','SRCSEQ','SRCVAR')
stopifnot(all(cols %in% names(x)), all(nchar(cols)<=8L),nrow(x)==59580L,
          dplyr::n_distinct(x$USUBJID)==254L,dplyr::n_distinct(x$PARAMCD)==47L)
candidate <- x %>% select(all_of(cols))
# ADT is stored as source date text; convert only known ISO 8601 YYYY-MM-DD (never fabricate).
if (is.character(candidate$ADT)) {
  bad_date <- !is.na(candidate$ADT) & !grepl('^\\d{4}-\\d{2}-\\d{2}$',candidate$ADT)
  if(any(bad_date)) stop('ADT has non-ISO dates: review before exporting.')
  candidate$ADT <- as.Date(candidate$ADT)
}
# fail instead of silently truncating text in SAS transport v5.
char_cols <- names(candidate)[vapply(candidate,is.character,logical(1))]
char_lengths <- vapply(candidate[char_cols],function(v) if(all(is.na(v))) 0L else max(nchar(v[!is.na(v)],type='bytes')),integer(1))
readr::write_csv(tibble(variable=char_cols,max_bytes=as.integer(char_lengths)),file.path(out,'character_lengths.csv'))
if(any(char_lengths>200L)) stop('Character length over 200 bytes. Review outputs/step12/character_lengths.csv')
if(anyDuplicated(candidate[c('STUDYID','USUBJID','SRCSEQ')])) stop('Duplicate source keys')
# Reserve exact baseline and visit policies as review decisions, not automated fixes.
missing_base <- x %>% filter(ANL01FL=='Y',is.finite(AVAL),is.na(BASE))
unmapped <- x %>% filter(is.na(AVISIT))
ambiguous <- unmapped %>% filter(toupper(trimws(VISIT))=='BASELINE')
readr::write_csv(missing_base %>% count(PARAMCD,name='selected_without_baseline') %>% arrange(desc(selected_without_baseline)),file.path(out,'baseline_review_AGGREGATE.csv'))
readr::write_csv(unmapped %>% count(VISIT,VISITNUM,name='records') %>% arrange(desc(records)),file.path(out,'visit_review_AGGREGATE.csv'))
readr::write_csv(ambiguous %>% count(VISIT,VISITNUM,name='records'),file.path(out,'same_day_baseline_review_AGGREGATE.csv'))
stopifnot(nrow(missing_base)==917L,nrow(unmapped)==1664L,nrow(ambiguous)==12L,
          sum(unmapped$ANL01FL=='Y',na.rm=TRUE)==0L)
# Actual SAS transport v5; appropriate metadata approval and Define-XML still open.
# haven uses the XPT filename stem as the default SAS member name (maximum 8 characters).
# Submission status is documented separately, not encoded in the member name.
xpt_path <- file.path(out,'ADLB.xpt')
haven::write_xpt(candidate,path=xpt_path,version=5)
if (!file.exists(xpt_path)) stop('XPT export did not create the expected file: ', xpt_path)
reloaded <- haven::read_xpt(xpt_path)
stopifnot(nrow(reloaded)==nrow(candidate),identical(names(reloaded),names(candidate)),
          !anyDuplicated(reloaded[c('STUDYID','USUBJID','SRCSEQ')]))
# All-column equality: permit only floating-point transport precision tolerance.
comparison <- lapply(names(candidate),function(nm){
 a<-candidate[[nm]];b<-reloaded[[nm]]
 if(inherits(a,'Date')) b <- as.Date(b)
 if(is.numeric(a)) {
  bad <- xor(is.na(a),is.na(b)) | (!is.na(a)&!is.na(b)&abs(a-b)>1e-8*pmax(1,abs(a)))
 } else {
  a <- as.character(a); b <- as.character(b)
  a[!is.na(a) & trimws(a)==''] <- NA_character_
  b[!is.na(b) & trimws(b)==''] <- NA_character_
  bad <- xor(is.na(a),is.na(b)) | (!is.na(a)&!is.na(b)&a!=b)
 }
 tibble(variable=nm,mismatches=sum(bad),type_before=paste(class(a),collapse='/'),type_after=paste(class(b),collapse='/'))
}) %>% bind_rows()
readr::write_csv(comparison,file.path(out,'xpt_readback_all_column_checks.csv'))
if(any(comparison$mismatches)) stop('XPT read-back differences: review outputs/step12/xpt_readback_all_column_checks.csv')
summary <- tibble(metric=c('rows','subjects','parameters','candidate_columns','xpt_bytes','readback_mismatch_columns',
 'selected_without_baseline_OPEN','unassigned_visits_OPEN','baseline_labeled_unmapped_OPEN'),
 value=c(nrow(candidate),n_distinct(candidate$USUBJID),n_distinct(candidate$PARAMCD),ncol(candidate),
 unname(file.info(xpt_path)$size),sum(comparison$mismatches>0),nrow(missing_base),nrow(unmapped),nrow(ambiguous)))
readr::write_csv(summary,file.path(out,'step12_summary.csv'))
cat('STEP12 XPT v5 WRITTEN AND READ BACK; source ADLB and ANCOVA unchanged.\n')
print(summary,n=nrow(summary))
cat('NOT SUBMISSION-READY: variable and parameter metadata, Define-XML, formal validation, independent derivations and baseline/visit decisions remain OPEN.\n')
