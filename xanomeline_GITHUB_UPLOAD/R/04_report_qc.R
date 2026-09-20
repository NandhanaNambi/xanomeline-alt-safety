# 04_report_qc.R | descriptive and model QC plus a shareable report
# Run from project root after scripts 02 and 03 complete.
required <- c('dplyr','readr','emmeans','ggplot2')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) install.packages(missing,repos='https://cloud.r-project.org')
suppressPackageStartupMessages({library(dplyr);library(readr);library(emmeans);library(ggplot2)})
dir.create('outputs/report',recursive=TRUE,showWarnings=FALSE)
adsl <- readRDS('outputs/adam/adsl.rds')
adlb <- readRDS('outputs/adam/adlb_week2.rds')
expected <- c('Placebo','Xanomeline Low Dose','Xanomeline High Dose')
stopifnot(!anyDuplicated(adsl$USUBJID),!anyDuplicated(adlb$USUBJID),
          all(adlb$USUBJID %in% adsl$USUBJID),
          all(adlb$PARAMCD=='ALT'),all(adlb$AVISITN==4),
          all(is.na(adlb$CHG) | abs(adlb$CHG-(adlb$AVAL-adlb$BASE))<1e-9))
complete <- adlb %>% filter(SAFFL=='Y',ANL01FL=='Y') %>%
  mutate(TRT01P=factor(TRT01P,levels=expected))
stopifnot(all(table(complete$TRT01P)>0), all(is.finite(complete$BASE)),
          all(is.finite(complete$AVAL)), nrow(complete)>10)
fit <- lm(AVAL~TRT01P+BASE,data=complete)
means <- emmeans(fit,~TRT01P)
new_contr <- as.data.frame(summary(contrast(means,method=list(
  'Low dose - Placebo'=c(-1,1,0),'High dose - Placebo'=c(-1,0,1))),
  infer=c(TRUE,TRUE),adjust='none'))
saved_contr <- read_csv('outputs/results/table_02_placebo_comparisons.csv',show_col_types=FALSE)
check <- new_contr %>% select(contrast,estimate,SE,df,lower.CL,upper.CL,p.value) %>%
  left_join(saved_contr %>% select(contrast,estimate,SE,df,lower.CL,upper.CL,p.value),
            by='contrast',suffix=c('_recalc','_saved'))
stopifnot(nrow(check)==2L)
for(nm in c('estimate','SE','df','lower.CL','upper.CL','p.value')){
  stopifnot(all(abs(check[[paste0(nm,'_recalc')]]-check[[paste0(nm,'_saved')]])<1e-6))
}
# Independent algebraic contrast verification: regression coefficient for high dose
# equals adjusted high-minus-placebo under this additive common-slope model.
primary <- new_contr %>% filter(contrast=='High dose - Placebo')
stopifnot(nrow(primary)==1L,
          abs(primary$estimate-unname(coef(fit)['TRT01PXanomeline High Dose']))<1e-8)
write_csv(check,'outputs/report/independent_contrast_reconciliation.csv')

coverage <- adlb %>% group_by(TRT01P) %>% summarise(
  treated_n=n(),baseline_n=sum(is.finite(BASE)),week2_n=sum(is.finite(AVAL)),
  complete_case_n=sum(ANL01FL=='Y'),
  complete_case_pct=round(100*complete_case_n/treated_n,1),
  missing_baseline_n=sum(!is.finite(BASE)),
  missing_week2_n=sum(!is.finite(AVAL)),.groups='drop')
write_csv(coverage,'outputs/report/table_04_analysis_population.csv')

# Planned-versus-actual mismatches and EX evidence need review before final safety claims.
recon <- adsl %>% filter(SAFFL=='Y',TRT01P %in% expected) %>%
  count(TRT01P,TRT01A,HAS_EX,name='n')
write_csv(recon,'outputs/report/table_05_treatment_reconciliation.csv')
mismatch_n <- sum(adsl$SAFFL=='Y' & adsl$TRT01P %in% expected &
                    !is.na(adsl$TRT01A) & adsl$TRT01P!=adsl$TRT01A)
ex_discrepancy_n <- sum((adsl$SAFFL=='Y' & adsl$HAS_EX=='N') |
                        (adsl$SAFFL=='N' & adsl$HAS_EX=='Y'))

# Diagnostics for model assumptions, not significance cherry-picking.
diag <- data.frame(USUBJID=complete$USUBJID,assigned_arm=as.character(complete$TRT01P),
                   fitted=fitted(fit),residual=residuals(fit),
                   standardized_residual=rstandard(fit),cooks_distance=cooks.distance(fit))
write_csv(diag,'outputs/report/model_diagnostics.csv')
plot_diag <- ggplot(diag,aes(fitted,residual)) +
  geom_point(alpha=.55) + geom_hline(yintercept=0,linetype=2) +
  labs(title='ANCOVA residuals versus fitted ALT',x='Fitted Week 2 ALT (U/L)',y='Residual (U/L)') +
  theme_minimal()
ggsave('outputs/report/figure_02_residuals.png',plot_diag,width=8,height=5,dpi=180)

# Descriptive comparison of change, explicitly not the primary adjusted contrast.
desc <- read_csv('outputs/results/table_00_descriptives.csv',show_col_types=FALSE)
trend <- read_csv('outputs/results/table_03_exploratory_assignment_trend.csv',show_col_types=FALSE)
trend_rank <- trend %>% filter(term=='DOSE_RANK')
fmt <- function(x,d=2) ifelse(is.na(x),'NA',formatC(x,format='f',digits=d))
fmt_p <- function(x) ifelse(x<.0001,'<0.0001',fmt(x,4))
render_table <- function(x){
  paste(c(paste0('| ',paste(names(x),collapse=' | '),' |'),
          paste0('|',paste(rep('---',ncol(x)),collapse='|'),'|'),
          apply(x,1,function(r) paste0('| ',paste(r,collapse=' | '),' |'))),collapse='\n')
}
primary_line <- sprintf('High dose minus placebo: %s U/L (95%% CI %s to %s), p = %s.',
                        fmt(primary$estimate),fmt(primary$lower.CL),fmt(primary$upper.CL),fmt_p(primary$p.value))
low <- new_contr %>% filter(contrast=='Low dose - Placebo')
low_line <- sprintf('Low dose minus placebo: %s U/L (95%% CI %s to %s), p = %s.',
                    fmt(low$estimate),fmt(low$lower.CL),fmt(low$upper.CL),fmt_p(low$p.value))
interpret <- if(primary$lower.CL>0) {
  'The observed baseline-adjusted high-dose minus placebo difference is positive and its confidence interval excludes zero.'
} else if(primary$upper.CL<0) {
  'The observed baseline-adjusted high-dose minus placebo difference is negative and its confidence interval excludes zero.'
} else {
  'The 95% confidence interval includes zero; these data do not establish a nonzero high-dose versus placebo difference at Week 2. The interval also permits both a small decrease and an increase.'
}
trend_sentence <- if(nrow(trend_rank)==1) sprintf(
  'Ordered assigned-arm rank (placebo=0, low=1, high=2): slope %s U/L per rank step (95%% CI %s to %s), exploratory, assumes equal spacing of ranks, not actual mg exposure.',
  fmt(trend_rank$Estimate),fmt(trend_rank$lower_95),fmt(trend_rank$upper_95)) else
  'Trend coefficient could not be extracted; inspect the exploratory trend CSV.'
md <- c('# Xanomeline and early ALT: Week 2 analysis report','',
  '**Project type:** retrospective educational analysis of the CDISC Pilot Study example data; not a prospective regulatory analysis.','',
  '## Research question','Does assigned high-dose xanomeline differ from placebo in Week 2 ALT after adjustment for baseline ALT? Low dose and an ordered assigned-dose trend are secondary/exploratory.','',
  '## Study population and derivation','Treated participants are identified provisionally from DM first-exposure dates. Assigned ARM is the primary treatment classification. ALT baseline is the SDTM LBBLFL=Y record; Week 2 is VISITNUM=4 and VISIT=WEEK 2; only paired finite measurements enter the ANCOVA. No imputation is performed.','',
  '## Primary statistical result',primary_line,interpret,'',
  '## Secondary and exploratory results',low_line,trend_sentence,
  'Secondary/exploratory p-values are unadjusted for multiplicity. A randomized dose group trend cannot establish individual exposure-response.','',
  '## Analysis population and missingness',render_table(coverage),'',
  '## Descriptive ALT values (unadjusted)',render_table(desc),'',
  '## Treatment and source-data QC',
  sprintf('Participants with different assigned and actual arm: %d. DM-versus-EX treated-evidence discrepancies: %d. Review these before making definitive on-treatment safety claims.',mismatch_n,ex_discrepancy_n),
  'Separate assigned-arm and actual-treatment reconciliation is saved in `table_05_treatment_reconciliation.csv`.','',
  '## Model and independent checks',
  sprintf('ANCOVA: `Week2_ALT ~ assigned_arm + baseline_ALT`; complete-case N = %d; residual degrees of freedom = %s.',nrow(complete),fmt(df.residual(fit),0)),
  'Recomputed model contrasts agree with stored outputs to < 0.000001, and the high-dose contrast equals the corresponding regression coefficient. See `independent_contrast_reconciliation.csv`.','',
  '## Limitations',
  'Single early visit; complete-case analysis; no multiple-comparison adjustment for secondary/exploratory analyses; study dataset is educational; treatment misclassification/mismatch requires review; residual assumptions and outliers require inspection; group mean ALT does not measure incidence of clinically important ALT elevations.','',
  '## Outputs','Primary and secondary results: `outputs/results/table_02_placebo_comparisons.csv`; adjusted means: `outputs/results/table_01_adjusted_means.csv`; trend: `outputs/results/table_03_exploratory_assignment_trend.csv`; QC: `outputs/report/`; ALT change figure: `outputs/results/figure_01_alt_change.png`.','',
  '## Reproduction','From project root, run scripts `01_sdtm_audit.R`, `02_build_adam.R`, `03_ancova_tlf.R`, and `04_report_qc.R` in order. Review QC flags rather than bypassing script stops.','')
writeLines(md,'outputs/report/analysis_report.md')
# Minimal self-contained HTML report; the markdown text remains the authoritative editable report.
esc <- function(x){x <- gsub('&','&amp;',x,fixed=TRUE); x <- gsub('<','&lt;',x,fixed=TRUE);x <- gsub('>','&gt;',x,fixed=TRUE);x}
html_lines <- vapply(md,function(line){
  if(startsWith(line,'# ')) paste0('<h1>',esc(sub('^# ','',line)),'</h1>')
  else if(startsWith(line,'## ')) paste0('<h2>',esc(sub('^## ','',line)),'</h2>')
  else if(grepl('^\\|',line)) paste0('<pre>',esc(line),'</pre>')
  else if(line=='') '' else paste0('<p>',esc(line),'</p>')
},character(1))
html <- c('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Xanomeline ALT | Analysis report</title>',
 '<style>body{font:17px/1.6 system-ui, sans-serif;max-width:920px;margin:2rem auto;padding:0 1.3rem;color:#183047;background:#fafcfd}h1,h2{line-height:1.2}h2{margin-top:2rem;border-bottom:1px solid #ccd;padding-bottom:.4rem}p{margin:.65rem 0}pre{overflow-x:auto;white-space:pre-wrap;font-size:12px;background:#edf3f6;padding:.3rem}img{max-width:100%;border-radius:9px}</style></head><body>',html_lines,
 '<h2>Figures</h2><p>Week 2 ALT change (descriptive, unadjusted):</p><img src="../results/figure_01_alt_change.png" alt="ALT change by treatment arm">',
 '<p>ANCOVA residuals versus fitted ALT:</p><img src="figure_02_residuals.png" alt="Model residual diagnostic">',
 '</body></html>')
writeLines(html,'outputs/report/analysis_report.html')
cat('\nREPORT COMPLETE. Primary: ',primary_line,'\n',interpret,'\n',low_line,
    '\nTreatment mismatches: ',mismatch_n,'; EX evidence discrepancies: ',ex_discrepancy_n,
    '\nOpen outputs/report/analysis_report.html and review diagnostics.\n',sep='')
