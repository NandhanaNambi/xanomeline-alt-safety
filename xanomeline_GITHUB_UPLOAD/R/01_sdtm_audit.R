# 01_sdtm_audit.R — CDISC Pilot Study ALT audit
# Run from project root: source("R/01_sdtm_audit.R")
# This script does NOT calculate treatment-effect estimates.

required <- c("pharmaversesdtm", "dplyr", "tidyr", "readr", "tibble")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
})
dir.create("outputs/audit", recursive = TRUE, showWarnings = FALSE)

dm <- pharmaversesdtm::dm
lb <- pharmaversesdtm::lb
ex <- pharmaversesdtm::ex
ds <- pharmaversesdtm::ds
ae <- pharmaversesdtm::ae
stopifnot(!anyDuplicated(dm$USUBJID), all(c("USUBJID", "ARM", "ACTARM") %in% names(dm)))
stopifnot(all(c("USUBJID", "LBTESTCD", "LBSTRESN", "LBSTRESU", "VISITNUM", "VISIT", "LBBLFL", "LBSEQ", "LBDTC") %in% names(lb)))

write_csv(tibble(domain = c("DM", "LB", "EX", "DS", "AE"),
                 rows = c(nrow(dm), nrow(lb), nrow(ex), nrow(ds), nrow(ae)),
                 subjects = c(n_distinct(dm$USUBJID), n_distinct(lb$USUBJID),
                              n_distinct(ex$USUBJID), n_distinct(ds$USUBJID),
                              n_distinct(ae$USUBJID))), "outputs/audit/domain_inventory.csv")

# Actual treatment labels are reported, not inferred from assumed dose codes.
arm_counts <- dm %>% count(ARM, ACTARM, name = "subjects")
write_csv(arm_counts, "outputs/audit/arm_counts.csv")
write_csv(ex %>% count(EXTRT, EXDOSE, EXDOSU, name = "records"),
          "outputs/audit/exposure_codes.csv")

alt <- lb %>% filter(LBTESTCD == "ALT") %>%
  left_join(dm %>% select(USUBJID, ARM, ACTARM, RFXSTDTC), by = "USUBJID") %>%
  mutate(ARM_ANALYSIS = if_else(!is.na(ACTARM) & ACTARM != "", ACTARM, ARM),
         BASELINE_FLAG = !is.na(LBBLFL) & LBBLFL == "Y",
         HAS_VALUE = is.finite(LBSTRESN),
         LBDATE = suppressWarnings(as.Date(substr(LBDTC, 1, 10))),
         FIRST_DOSE_DATE = suppressWarnings(as.Date(substr(RFXSTDTC, 1, 10))))
stopifnot(nrow(alt) > 0L)
write_csv(alt %>% count(LBSTRESU, name = "records"), "outputs/audit/alt_units.csv")
write_csv(alt %>% count(VISITNUM, VISIT, LBBLFL, name = "records"),
          "outputs/audit/alt_visit_dictionary.csv")

# Denominator: DM participants in each arm with a recorded first-treatment date.
# This is a provisional treated-population denominator, NOT a finalized SAFFL.
pop <- dm %>% mutate(ARM_ANALYSIS = if_else(!is.na(ACTARM) & ACTARM != "", ACTARM, ARM),
                     TREATED_PROXY = !is.na(RFXSTDTC) & RFXSTDTC != "") %>%
  filter(TREATED_PROXY, !is.na(ARM_ANALYSIS), ARM_ANALYSIS != "") %>%
  count(ARM_ANALYSIS, name = "denominator_n")
write_csv(pop, "outputs/audit/provisional_population.csv")

baseline <- alt %>% filter(BASELINE_FLAG, HAS_VALUE) %>%
  distinct(USUBJID, ARM_ANALYSIS) %>% count(ARM_ANALYSIS, name = "baseline_n")
write_csv(pop %>% left_join(baseline, by = "ARM_ANALYSIS") %>%
            mutate(baseline_n = replace_na(baseline_n, 0L),
                   baseline_pct = 100 * baseline_n / denominator_n),
          "outputs/audit/baseline_coverage.csv")

# Exclude records marked baseline; visits with non-positive planned days are
# excluded if VISITDY is available. Keep true visit IDs; do not invent 'Week 2'.
post <- alt %>% filter(!BASELINE_FLAG, HAS_VALUE,
                       is.finite(VISITNUM),
                       is.na(VISITDY) | VISITDY > 1,
                       !is.na(ARM_ANALYSIS), ARM_ANALYSIS != "")
coverage <- post %>% group_by(VISITNUM, VISIT, ARM_ANALYSIS) %>%
  summarise(participants = n_distinct(USUBJID), records = n(),
            .groups = "drop") %>%
  right_join(tidyr::crossing(
    post %>% distinct(VISITNUM, VISIT), pop %>% select(ARM_ANALYSIS, denominator_n)),
    by = c("VISITNUM", "VISIT", "ARM_ANALYSIS")) %>%
  mutate(participants = replace_na(participants, 0L),
         records = replace_na(records, 0L),
         coverage_pct = round(100 * participants / denominator_n, 1),
         missing_n = denominator_n - participants) %>%
  arrange(VISITNUM, ARM_ANALYSIS)
write_csv(coverage, "outputs/audit/alt_visit_coverage.csv")

# Inspect repeats rather than silently selecting an arbitrary measurement.
duplicates <- alt %>% filter(HAS_VALUE) %>%
  count(USUBJID, VISITNUM, VISIT, name = "records") %>% filter(records > 1L)
write_csv(duplicates, "outputs/audit/alt_repeated_visit_records.csv")
write_csv(alt %>% filter(is.na(VISITNUM) | grepl("UNSCHED", VISIT, ignore.case = TRUE)),
          "outputs/audit/alt_unscheduled_records.csv")
write_csv(alt %>% filter(BASELINE_FLAG, HAS_VALUE) %>%
            count(USUBJID, name = "baseline_records") %>%
            filter(baseline_records > 1L),
          "outputs/audit/multiple_baselines.csv")
write_csv(alt %>% filter(HAS_VALUE, is.na(LBSTRESU) | LBSTRESU == ""),
          "outputs/audit/missing_units.csv")
write_csv(alt %>% filter(HAS_VALUE, !is.na(LBDATE), !is.na(FIRST_DOSE_DATE),
                         !BASELINE_FLAG, LBDATE < FIRST_DOSE_DATE),
          "outputs/audit/postbaseline_before_first_dose.csv")

# Decision support ONLY: >=80% in EACH arm, baseline completeness >=80% in
# EACH arm, exactly three arms. Review duplicates, dates, treatment coding and
# unscheduled visits before manually freezing the selected visit in the spec.
visit_screen <- coverage %>% group_by(VISITNUM, VISIT) %>%
  summarise(arms = n(), minimum_arm_coverage_pct = min(coverage_pct),
            passes_80pct_visit_coverage = arms == 3L && minimum_arm_coverage_pct >= 80,
            .groups = "drop") %>% arrange(VISITNUM)
write_csv(visit_screen, "outputs/audit/early_visit_screen.csv")

cat("\nAUDIT COMPLETE. Review outputs/audit/*.csv, especially:\n",
    "arm_counts.csv, alt_units.csv, baseline_coverage.csv,\n",
    "alt_visit_coverage.csv, early_visit_screen.csv and duplicates.\n",
    "Do not run ANCOVA until the visit, populations and repeated-record rules are frozen.\n")
print(arm_counts)
print(coverage, n = 100)
print(visit_screen, n = 100)
