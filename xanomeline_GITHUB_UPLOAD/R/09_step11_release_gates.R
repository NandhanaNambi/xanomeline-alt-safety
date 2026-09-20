# Step 11: read-only release-gate audit; run from xanomeline_alt_project root.
# Does not edit full ADLB or the pre-existing Week 2 ALT analysis.
if (!requireNamespace('dplyr',quietly=TRUE) || !requireNamespace('readr',quietly=TRUE))
  stop('Install dplyr and readr before running.')
stopifnot(file.exists('outputs/full_adlb/adlb_full.rds'))
x <- readRDS('outputs/full_adlb/adlb_full.rds')
need <- c('USUBJID','PARAMCD','ANL01FL','ABLFL','AVAL','BASE','AVISIT','VISIT','SRCSEQ')
stopifnot(all(need %in% names(x)))
missing <- x |> dplyr::filter(ANL01FL=='Y',is.finite(AVAL),is.na(BASE))
unmapped <- x |> dplyr::filter(is.na(AVISIT))
stopifnot(nrow(x)==59580L,dplyr::n_distinct(x$PARAMCD)==47L)
dir.create('outputs/step11',recursive=TRUE,showWarnings=FALSE)
readr::write_csv(missing |> dplyr::count(PARAMCD,name='rows_without_baseline'),
 'outputs/step11/missing_baseline_by_parameter_AGGREGATE.csv')
readr::write_csv(unmapped |> dplyr::count(VISIT,name='records'),
 'outputs/step11/unassigned_visits_AGGREGATE.csv')
counts <- data.frame(metric=c('rows','subjects','parameters','selected_rows_without_baseline',
 'unassigned_analysis_visit','unassigned_selected_for_analysis'),
 value=c(nrow(x),dplyr::n_distinct(x$USUBJID),dplyr::n_distinct(x$PARAMCD),
 nrow(missing),nrow(unmapped),sum(unmapped$ANL01FL=='Y',na.rm=TRUE)))
readr::write_csv(counts,'outputs/step11/step11_reconciliation_AGGREGATE.csv')
print(counts,row.names=FALSE)
if(nrow(missing)!=917L || nrow(unmapped)!=1664L) {
 stop('Expected historical counts differ: investigate data version before reviewing Step 11 metadata.')
}
cat('STEP 11 READ-ONLY RECONCILIATION PASSED; OPEN GATES ARE NOT CLEARED. No XPT / Define-XML / formal conformance validation performed.\n')
