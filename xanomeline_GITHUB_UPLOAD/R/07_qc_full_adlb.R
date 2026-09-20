# 07_qc_full_adlb.R | independent post-build checks, reporting not a conformance certification.
required <- c('dplyr','readr','pharmaversesdtm')
miss <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(miss)) install.packages(miss,repos='https://cloud.r-project.org')
suppressPackageStartupMessages({library(dplyr);library(readr)})
p <- 'outputs/full_adlb'
if(!file.exists(file.path(p,'adlb_full.rds'))) stop('Build ADLB first: source("R/06_build_full_adlb.R")')
x <- readRDS(file.path(p,'adlb_full.rds'))
s <- pharmaversesdtm::lb
b <- x %>% filter(ABLFL=='Y') %>% select(USUBJID,PARAMCD,BASE,AVAL)
y <- x %>% filter(ANL01FL=='Y')
results <- tibble::tibble(
  check=c('All LB records retained','LB source keys unique','Each baseline unique','Each visit analysis row unique',
          'All analysis rows numeric','Baseline flag equals measured value',
          'Change equals analysis minus baseline','Source-to-result values match',
          'Every source has study subject','PARAMCD <= 8 characters'),
  pass=c(nrow(x)==nrow(s),!anyDuplicated(x[c('STUDYID','USUBJID','SRCSEQ')]),
         !anyDuplicated(b[c('USUBJID','PARAMCD')]),
         !anyDuplicated(y[c('USUBJID','PARAMCD','AVISITN')]),
         all(is.finite(y$AVAL)),all(abs(b$BASE-b$AVAL)<1e-9),
         all(abs(x$CHG[!is.na(x$CHG)]-(x$AVAL-x$BASE)[!is.na(x$CHG)])<1e-9),
         isTRUE(all.equal(x$AVAL,ifelse(x$LBTESTCD %in%
           c('ANISO','KETONES','MACROCY','MICROCY','POIKILO','POLYCHR','UROBIL'),
           NA_real_,ifelse(is.finite(x$LBSTRESN),x$LBSTRESN,NA_real_)),
           check.attributes=FALSE)),
         all(!is.na(x$USUBJID)),all(nchar(x$PARAMCD)<=8)),
  meaning=c('full source coverage','traceability','baseline uniqueness','visit selection uniqueness',
            'valid numerical analysis records','baseline consistency','arithmetic',
            'LBSTRESN traceability with categorical code exclusions','subject linkage','parameter code length'))
write_csv(results,file.path(p,'qc_checks.csv'))
print(results,n=nrow(results))
if(any(!results$pass)) stop('QC FAILED: inspect outputs/full_adlb/qc_checks.csv')
cat('ALL 10 PROGRAMMATIC CHECKS PASSED. This does NOT certify CDISC conformance or submission readiness.\n')
