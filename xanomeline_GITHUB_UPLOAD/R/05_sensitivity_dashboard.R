# 05_sensitivity_dashboard.R | sensitivity checks and Power BI-ready CSVs
# Run from project root AFTER 01-04. Retains assigned arm for the primary estimand.
required <- c('dplyr','readr','ggplot2','emmeans')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) install.packages(missing,repos='https://cloud.r-project.org')
suppressPackageStartupMessages({library(dplyr);library(readr);library(ggplot2);library(emmeans)})
dir.create('outputs/sensitivity',recursive=TRUE,showWarnings=FALSE)
dir.create('outputs/dashboard',recursive=TRUE,showWarnings=FALSE)
adsl <- readRDS('outputs/adam/adsl.rds')
adlb <- readRDS('outputs/adam/adlb_week2.rds')
arms <- c('Placebo','Xanomeline Low Dose','Xanomeline High Dose')
dat <- adlb %>% filter(SAFFL=='Y',ANL01FL=='Y') %>%
  mutate(TRT01P=factor(TRT01P,levels=arms))
stopifnot(!anyDuplicated(adlb$USUBJID), all(table(dat$TRT01P)>0),
          all(is.finite(dat$AVAL)),all(is.finite(dat$BASE)))
# CI and p-values are model-based and nominal; no causal inference for actual-arm analyses.
get_contrasts <- function(dataset, formula, method){
 fit <- lm(formula,data=dataset)
 if(anyNA(coef(fit))) stop('Singular sensitivity model: ',method)
 emm <- emmeans(fit,~TRT01P)
 z <- as.data.frame(summary(contrast(emm,method=list(
  'Low dose - Placebo'=c(-1,1,0),
  'High dose - Placebo'=c(-1,0,1))),infer=c(TRUE,TRUE),adjust='none'))
 z$method <- method
 z$n <- nobs(fit)
 z
}
base <- get_contrasts(dat,AVAL~TRT01P+BASE,'Primary ANCOVA (complete cases)')
change <- get_contrasts(dat,CHG~TRT01P,'Unadjusted change-score comparison')
# Log ANCOVA estimates exponentiated ratios, NOT differences in U/L.
log_results <- NULL
if(all(dat$AVAL>0 & dat$BASE>0)) {
 log_dat <- dat %>% mutate(LOG_AVAL=log(AVAL),LOG_BASE=log(BASE))
 log_results <- get_contrasts(log_dat,LOG_AVAL~TRT01P+LOG_BASE,'Log-scale ANCOVA') %>%
   mutate(ratio=exp(estimate),ratio_lower=exp(lower.CL),ratio_upper=exp(upper.CL))
 write_csv(log_results,'outputs/sensitivity/table_07_log_alt_ratios.csv')
} else {
 writeLines('Log-scale sensitivity not performed: one or more nonpositive ALT values.',
            'outputs/sensitivity/log_sensitivity_status.txt')
}
# Fixed robust-outlier sensitivity: trim by OUTCOME-POOLED 99th percentile of AVAL;
# exploratory, excludes genuine high values, not a preferred alternative analysis.
cutoff <- unname(quantile(dat$AVAL,probs=.99,na.rm=TRUE,type=7))
trimmed <- dat %>% filter(AVAL<=cutoff)
trim_result <- if(all(table(trimmed$TRT01P)>0) && nrow(trimmed)>10) {
 get_contrasts(trimmed,AVAL~TRT01P+BASE,'Exploratory: exclude ALT above pooled 99th percentile')
} else NULL
write_csv(bind_rows(base,change,trim_result),
          'outputs/sensitivity/table_06_sensitivity_comparisons.csv')
write_csv(data.frame(method='Pooled observed Week 2 ALT 99th percentile',
  cutoff_U_L=cutoff, excluded_n=nrow(dat)-nrow(trimmed),
  note='Exploratory; trimming observed outcomes can bias treatment contrasts.'),
  'outputs/sensitivity/outlier_sensitivity_metadata.csv')
# Missingness is descriptive; do NOT interpret complete-case estimates as missing-data robust.
miss <- adlb %>% group_by(TRT01P) %>% summarise(
 eligible_n=n(),paired_n=sum(ANL01FL=='Y',na.rm=TRUE),
 baseline_missing_n=sum(!is.finite(BASE)),week2_missing_n=sum(!is.finite(AVAL)),
 paired_pct=round(100*paired_n/eligible_n,1),.groups='drop')
write_csv(miss,'outputs/sensitivity/table_08_missingness.csv')
# Exact subject-level analysis inputs (deidentified study IDs); share only if study license permits.
patient <- adlb %>% transmute(USUBJID,assigned_arm=TRT01P,actual_arm=TRT01A,
 baseline_alt=BASE,week2_alt=AVAL,change_alt=CHG,
 paired_alt=ANL01FL,missing_reason=EXCLUSION,
 baseline_source_seq=BASE_LBSEQ,week2_source_seq=W2_LBSEQ)
write_csv(patient,'outputs/dashboard/patient_week2_alt.csv')
means <- read_csv('outputs/results/table_01_adjusted_means.csv',show_col_types=FALSE)
contr <- read_csv('outputs/results/table_02_placebo_comparisons.csv',show_col_types=FALSE)
write_csv(means,'outputs/dashboard/adjusted_means.csv')
write_csv(contr,'outputs/dashboard/placebo_comparisons.csv')
write_csv(miss,'outputs/dashboard/analysis_population.csv')
write_csv(read_csv('outputs/results/table_00_descriptives.csv',show_col_types=FALSE),
          'outputs/dashboard/arm_descriptives.csv')
# Subject-level two-visit trajectories: no spurious multi-visit interpolation.
trajectory <- bind_rows(
 patient %>% transmute(USUBJID,assigned_arm,actual_arm,visit='Baseline',visit_order=0L,alt_U_L=baseline_alt),
 patient %>% transmute(USUBJID,assigned_arm,actual_arm,visit='Week 2',visit_order=2L,alt_U_L=week2_alt)) %>%
 filter(is.finite(alt_U_L))
write_csv(trajectory,'outputs/dashboard/patient_alt_trajectories.csv')
# Figure of adjusted means and model-based CI; y-axis in U/L.
p <- ggplot(means,aes(x=TRT01P,y=emmean)) +
 geom_point(size=3)+geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),width=.15)+
 labs(title='Week 2 ALT: baseline-adjusted means',subtitle='Model-based 95% confidence intervals; assigned treatment',
 x='Assigned treatment',y='Adjusted Week 2 ALT (U/L)')+theme_minimal()+
 theme(axis.text.x=element_text(angle=15,hjust=1))
ggsave('outputs/sensitivity/figure_03_adjusted_means.png',p,width=9,height=5,dpi=180)
# Readable limitations and instructions for an accurate dashboard.
notes <- c('# Step 4: sensitivity results and dashboard specification','',
 'The primary model remains assigned-treatment Week 2 ANCOVA adjusted for baseline ALT.',
 'The unadjusted change-score analysis asks a related but different question.',
 'Log ANCOVA (only when all ALT values are positive) reports ratios; never label these U/L differences.',
 'Pooled outcome-based 99th percentile trimming is exploratory and potentially biased. Compare direction and size; do not select the most favorable estimate.',
 'Missing Week 2 ALT may be related to early discontinuation or adverse events; this workflow does not establish a missing-at-random mechanism.',
 'The 12 assigned-versus-actual mismatches must be retained. This script intentionally does not relabel the main analysis by treatment received.',
 'No ALT >3x ULN safety classification is created because a validated subject-/lab-specific ULN derivation is not present in the Week 2 ADLB.',
 'No patient-level responder classification or machine learning is presented as validated.', '',
 '## Recommended Power BI pages',
 '1. Primary results: placebo_comparisons.csv (differences, CIs, p-values), adjusted_means.csv, analysis_population.csv.',
 '2. Patient ALT: patient_alt_trajectories.csv, filter assigned_arm, show subject trajectories and number observed.',
 '3. Data QC: analysis_population.csv, patient_week2_alt.csv; display missingness and treatment mismatches.',
 'Connect CSVs via Get Data > Text/CSV; do not recompute ANCOVA in Power BI.',
 'Use USUBJID only to relate subject-level tables. Adjusted-means/contrast tables are aggregate and should NOT join to subject-level data.',
 'Comparisons are nominal and do not establish a monotonic dose response.', '',
 '## Execution', 'source("R/05_sensitivity_dashboard.R")',
 'Review sensitivity CSVs, figure_03, and residual diagnostics produced by 04_report_qc.R.')
writeLines(notes,'outputs/sensitivity/step4_readme.md')
cat('\nSTEP 4 COMPLETE: sensitivity checks and dashboard CSVs created.\n')
cat('Primary high-vs-placebo: ',round(base$estimate[base$contrast=='High dose - Placebo'],3),' U/L.\n',sep='')
cat('Paired sample N: ',nrow(dat),'; trimmed N: ',nrow(trimmed),'; excluded: ',nrow(dat)-nrow(trimmed),'\n',sep='')
cat('Dashboard files: outputs/dashboard/ | Review: outputs/sensitivity/step4_readme.md\n')
