# 03_ancova_tlf.R | run only AFTER 02_build_adam.R passes QC.
required <- c('dplyr','readr','emmeans','ggplot2')
miss <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(miss)) install.packages(miss,repos='https://cloud.r-project.org')
suppressPackageStartupMessages({library(dplyr);library(readr);library(emmeans);library(ggplot2)})
dir.create('outputs/results',recursive=TRUE,showWarnings=FALSE)
adlb <- readRDS('outputs/adam/adlb_week2.rds')
arms <- c('Placebo','Xanomeline Low Dose','Xanomeline High Dose')
dat <- adlb %>% filter(SAFFL=='Y',ANL01FL=='Y') %>%
  mutate(TRT01P=factor(TRT01P,levels=arms))
stopifnot(all(table(dat$TRT01P)>0),nrow(dat)>10,all(is.finite(dat$BASE)),all(is.finite(dat$AVAL)))
# Main model: randomized assigned treatment, common baseline ALT covariate.
fit <- lm(AVAL ~ TRT01P + BASE,data=dat)
lsm <- emmeans(fit,~TRT01P)
lsm_t <- as.data.frame(summary(lsm,infer=c(TRUE,FALSE)))
write_csv(lsm_t,'outputs/results/table_01_adjusted_means.csv')
contrasts <- contrast(lsm,method=list('Low dose - Placebo'=c(-1,1,0),
                                     'High dose - Placebo'=c(-1,0,1)))
contr_t <- as.data.frame(summary(contrasts,infer=c(TRUE,TRUE),adjust='none'))
write_csv(contr_t,'outputs/results/table_02_placebo_comparisons.csv')
write_csv(dat %>% group_by(TRT01P) %>% summarise(n=n(),
  baseline_mean=mean(BASE),week2_mean=mean(AVAL),mean_change=mean(CHG),
  median_change=median(CHG),.groups='drop'), 'outputs/results/table_00_descriptives.csv')
# Exploratory ordered assignment trend (0/1/2 ranks, NOT actual drug exposure).
dat$DOSE_RANK <- as.integer(dat$TRT01P)-1L
trend <- lm(AVAL~DOSE_RANK+BASE,data=dat)
trend_out <- as.data.frame(coef(summary(trend)))
trend_ci <- confint(trend)
trend_out$term <- rownames(trend_out)
trend_out$lower_95 <- trend_ci[,1]; trend_out$upper_95 <- trend_ci[,2]
write_csv(trend_out,'outputs/results/table_03_exploratory_assignment_trend.csv')
# Plot is unadjusted patient-level descriptive ALT change, not ANCOVA adjusted.
p <- ggplot(dat,aes(x=TRT01P,y=CHG)) + geom_boxplot(outlier.alpha=.4) +
  labs(title='ALT change from baseline to Week 2',subtitle='Descriptive, unadjusted; treated participants with paired ALT',
       x='Assigned treatment',y='ALT change (U/L)') + theme_minimal() +
  theme(axis.text.x=element_text(angle=15,hjust=1))
ggsave('outputs/results/figure_01_alt_change.png',p,width=9,height=5,dpi=180)
# model diagnostics and provenance
capture.output(summary(fit),file='outputs/results/ancova_model_summary.txt')
writeLines(c(paste('R:',R.version.string),paste('pharmaversesdtm:',as.character(packageVersion('pharmaversesdtm'))),
             paste('emmeans:',as.character(packageVersion('emmeans'))),
             paste('Complete-case N:',nrow(dat)),
             'Exploratory educational retrospective analysis; p-values unadjusted for multiplicity.'),
           'outputs/results/run_metadata.txt')
cat('\nPRIMARY RESULT (high dose minus placebo; U/L):\n')
print(contr_t %>% filter(contrast=='High dose - Placebo'))
cat('\nRESULTS saved in outputs/results/.\n')
