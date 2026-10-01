# Xanomeline ALT Safety Analysis

An end-to-end clinical statistical programming case study using the Pharmaverse CDISC Pilot Study example data.

This project evaluates whether xanomeline exposure is associated with an early increase in alanine aminotransferase (ALT), with emphasis on reproducible data derivation, analysis traceability, quality control, sensitivity analysis, regulatory-style reporting, and downstream Power BI visualization.

**Live case study:**  
https://nandhananambi.github.io/xanomeline-alt-safety/index.html?v=2

---

## Research Question

**Does xanomeline exposure lead to a greater early increase in ALT than placebo, and is the magnitude of that difference related to dose?**

### Primary comparison

High-dose xanomeline vs placebo.

### Secondary comparison

Low-dose xanomeline vs placebo.

### Exploratory objective

Evaluate whether the pattern of treatment effects is consistent with an increasing dose-response relationship.

### Primary endpoint

Baseline-adjusted ALT at Week 2.

---

## Key Result

The primary ANCOVA estimated the following treatment difference:

**High-dose xanomeline − placebo: +1.68 U/L**

- SE: 1.18 U/L
- 95% CI: −0.65 to 4.01 U/L
- p = 0.1569

The confidence interval includes zero, so the primary analysis does not establish a nonzero Week 2 ALT difference between the high-dose and placebo groups.

The secondary low-dose comparison was:

**Low-dose xanomeline − placebo: +2.64 U/L**

- 95% CI: 0.31 to 4.96 U/L
- nominal p = 0.0266

Because the low-dose estimate exceeded the high-dose estimate, the results did not show a clear monotonic increasing dose-response pattern.

Secondary and exploratory p-values were not adjusted for multiplicity.

---

## Analysis Population

The source DM domain contained **306 participants**.

The treatment/safety population used for the analysis contained **254 participants**:

| Assigned treatment | N |
|---|---:|
| Placebo | 86 |
| Xanomeline Low Dose | 84 |
| Xanomeline High Dose | 84 |

The complete-case Week 2 ALT ANCOVA population contained **239 participants** with both baseline and Week 2 ALT.

No missing-value imputation was performed.

---

## Week 2 Coverage

Week 2 was selected after auditing scheduled ALT coverage across treatment arms.

| Assigned treatment | Treated | Paired baseline + Week 2 | Coverage |
|---|---:|---:|---:|
| Placebo | 86 | 83 | 96.5% |
| Low Dose | 84 | 78 | 92.9% |
| High Dose | 84 | 78 | 92.9% |

Week 2 represented an early follow-up visit with high ALT completeness across all three assigned treatment groups.

---

## Statistical Model

The primary analysis used baseline-adjusted ANCOVA:

```r
lm(AVAL ~ TRT01P + BASE)
```

where:

- `AVAL` = Week 2 ALT
- `TRT01P` = assigned treatment
- `BASE` = baseline ALT
- placebo = reference treatment

The primary estimand was the adjusted mean difference between assigned high-dose xanomeline and placebo.

Adjusted treatment means and contrasts were obtained from the fitted model.

---

## End-to-End Workflow

```text
CDISC SDTM source data
        ↓
Source audit
        ↓
Treatment and visit reconciliation
        ↓
ADSL derivation
        ↓
Week 2 ALT ADLB derivation
        ↓
Baseline-adjusted ANCOVA
        ↓
Independent QC
        ↓
Sensitivity analyses
        ↓
Expanded multi-parameter ADLB
        ↓
XPT transport validation
        ↓
Regulatory-style TLFs
        ↓
Aggregate reporting
        ↓
Power BI dashboard
        ↓
Public case-study website
```

---

## Source Data

The project uses example clinical trial data from the Pharmaverse CDISC Pilot Study.

Key SDTM domains used include:

- **DM** — Demographics
- **LB** — Laboratory Tests
- **EX** — Exposure
- **DS** — Disposition
- **AE** — Adverse Events

The workflow demonstrates the transition from standardized source data to analysis-ready datasets and reporting outputs.

---

## ADSL and ADLB

### ADSL

The subject-level analysis dataset includes one row per participant and treatment/population variables used across downstream analyses.

Examples include:

- planned treatment
- actual treatment
- safety flag
- demographics
- treatment dates

### Week 2 ALT ADLB

A focused ALT analysis dataset was derived containing:

- baseline ALT
- Week 2 ALT
- change from baseline
- analysis flags
- assigned treatment
- source record sequence
- source date
- exclusion or missingness information

The derivation retained source-level traceability through variables linking analysis records back to the laboratory domain.

---

## Treatment Assignment Reconciliation

During source review, **12 participants assigned to the high-dose arm were recorded with low-dose actual treatment**.

The project preserves both treatment concepts rather than silently recoding them.

The primary analysis uses **assigned treatment**.

This distinction is explicitly documented to maintain reproducibility and transparency.

---

## Sensitivity Analyses

Several analyses were used to assess whether conclusions were dependent on the primary model specification.

### Unadjusted change-score analysis

High dose vs placebo:

**+1.40 U/L**

95% CI: −0.98 to 3.78 U/L

p = 0.2485

### Log-scale ANCOVA

High dose vs placebo adjusted ratio:

**1.157**

95% CI: 1.055 to 1.268

nominal p = 0.0020

This is a multiplicative treatment ratio and should not be interpreted as a U/L difference.

### Exploratory 99th-percentile outcome trim

Three observations above the pooled Week 2 ALT 99th percentile were excluded.

High dose vs placebo:

**+2.90 U/L**

95% CI: 1.01 to 4.79 U/L

p = 0.0028

This analysis is explicitly exploratory because exclusions were based on observed outcome values.

It does not replace the prespecified primary analysis.

---

## Model Diagnostics

The primary ANCOVA was reviewed using residual and influence diagnostics.

Review flags included:

- 4 observations with absolute standardized residual > 3
- 10 observations exceeding the heuristic Cook's distance threshold of `4/N`
- maximum Cook's distance approximately 1.08

These flags identify observations requiring review.

They were not used as automatic deletion criteria.

---

## Expanded Multi-Parameter ADLB

The original ALT analysis was later extended into a broader laboratory analysis dataset.

The expanded ADLB candidate contains:

- **59,580 laboratory records**
- **254 participants**
- **47 laboratory parameters**

The dataset preserves source records while deriving analysis variables including:

- `PARAMCD`
- `PARAM`
- `PARCAT1`
- `AVAL`
- `AVALC`
- `AVALU`
- `ADT`
- `ADY`
- `AVISIT`
- `AVISITN`
- `ABLFL`
- `BASE`
- `CHG`
- `PCHG`
- `ANL01FL`
- reference-range variables
- treatment variables
- source-traceability variables

Numeric and categorical laboratory tests were handled separately rather than forcing all source measurements into a numeric analysis schema.

---

## Expanded ADLB QC

The expanded laboratory dataset underwent programmatic source-to-analysis checks.

Core QC included:

- source-row retention
- unique record keys
- baseline uniqueness
- analysis-visit uniqueness
- numeric-value consistency
- baseline consistency
- change-from-baseline arithmetic
- source-to-analysis value reconciliation
- participant linkage
- parameter naming constraints

The core QC suite passed **10/10 programmed checks**.

Additional release-review outputs were also generated.

Important documented review items include:

- 917 selected analysis rows without baseline across subject-parameter combinations
- 1,664 source rows without assigned analysis visit
- zero unassigned-visit records selected for analysis
- baseline timing cases requiring additional review

These records were retained and documented rather than silently removed or imputed.

---

## XPT Transport Validation

A SAS Transport v5 analysis candidate was created:

```text
ADLB.xpt
```

The exchange candidate contained **31 exported variables**.

After normalizing the expected representation difference between R missing character values and blank XPT character values, all exported variables matched on read-back.

**31/31 variables showed zero read-back differences.**

This demonstrates transport consistency under the documented comparison rules.

It does **not** constitute formal CDISC conformance certification or regulatory submission validation.

---

## Tables, Listings and Figures

Regulatory-style reporting outputs were generated from the analysis datasets.

Outputs include:

- demographic age summaries
- demographic sex summaries
- all-reported adverse-event summaries
- treatment-comparison tables
- sensitivity-analysis outputs
- Week 2 ALT subject listing
- diagnostic figures

The subject-level listing remains local and is not distributed publicly.

The adverse-event table summarizes **all reported adverse events**.

It should not be interpreted as a validated TEAE table because a formal treatment-emergent onset-window derivation was not implemented.

---

## TLF Quality Control

Focused source-based QC was performed on the reporting outputs.

The TLF QC workflow passed:

**32 programmed checks**

Public-facing aggregate tables were then separated from private participant-level outputs.

Small serious-event counts were excluded from public aggregate reporting as part of the disclosure review.

---

## Power BI Dashboard

A three-page Power BI dashboard was built using aggregate outputs generated by the R analysis workflow.

### Page 1 — ALT Results

Includes:

- eligible safety population
- paired ALT analysis population
- baseline-adjusted placebo comparisons
- adjusted mean Week 2 ALT by assigned treatment

### Page 2 — Safety & Demographics

Includes:

- participants with at least one reported adverse event
- sex distribution by treatment arm
- age summaries by treatment arm

### Page 3 — Data Quality & Analysis Population

Includes:

- paired ALT completeness
- missing Week 2 ALT counts
- eligible vs analysis population summaries

The Power BI report is used as a reporting and visualization layer.

Statistical derivations and model fitting remain in R rather than being recomputed in Power BI.

### Power BI Files

See:

```text
powerbi/
```

The folder contains:

- `Xanomeline_ALT_Safety_Dashboard.pbix`
- `Xanomeline_ALT_Safety_Dashboard.pdf`

---

## Public Aggregate Data

Public aggregate TLF outputs are stored under:

```text
data/tlf/
```

Including:

- `all_reported_ae_AGGREGATE.csv`
- `demographics_age_AGGREGATE.csv`
- `demographics_sex_AGGREGATE.csv`
- `TLF_METHODS_AND_LIMITATIONS.txt`

Participant-level laboratory and listing data are not published.

---

## Reproducible R Workflow

The repository contains modular R programs covering the full workflow.

```text
R/
```

Key scripts include:

```text
01_sdtm_audit.R
02_build_adam.R
03_ancova_tlf.R
04_report_qc.R
05_sensitivity_dashboard.R
06_build_full_adlb.R
07_qc_full_adlb.R
08_step10_remediate_audit.R
09_step11_release_gates.R
10_step12_exchange_candidate.R
11_build_tlf_package.R
12_qc_tlf_package.R
13_run_tlf_only.R
14_publish_aggregate_tlfs.R
```

The scripts separate:

- source auditing
- derivation
- statistical analysis
- quality control
- sensitivity analysis
- expanded ADLB generation
- exchange-format checks
- reporting
- public aggregate publication

---

## Specifications

Dataset and reporting specifications are stored under:

```text
specifications/
```

Examples include:

- `ADLB_31_variable_spec_DRAFT.csv`
- `ADLB_47_parameter_spec_DRAFT.csv`
- `ADLB_SPECIFICATION.md`
- `TLF_SPECIFICATION.md`

These documents describe analysis-variable intent, scope and limitations.

---

## Repository Structure

```text
xanomeline-alt-safety/
│
├── R/
│   └── reproducible R analysis and QC programs
│
├── specifications/
│   └── ADLB and TLF specifications
│
├── data/
│   └── public aggregate outputs
│
├── powerbi/
│   ├── Xanomeline_ALT_Safety_Dashboard.pbix
│   └── Xanomeline_ALT_Safety_Dashboard.pdf
│
├── website/
│   └── public case-study website and figures
│
├── REPRODUCIBILITY.md
├── package_versions.csv
├── README.md
└── index.html
```

---

## Data Governance and Disclosure

The public repository intentionally excludes participant-level analysis data.

Publicly distributed outputs are restricted to aggregate results, figures, documentation and reproducible code.

Participant-level laboratory listings, source keys and internal QC artifacts remain local.

---

## Important Limitations

This project is an educational clinical-data and statistical-programming case study.

It should not be interpreted as a regulatory analysis or clinical safety determination.

Key limitations include:

- retrospective analysis of example clinical trial data
- complete-case Week 2 analysis
- no missing-value imputation
- influential observations in the original-scale ANCOVA
- secondary and exploratory analyses without multiplicity adjustment
- no validated treatment-emergent AE derivation
- provisional expanded ADLB rules for some laboratory parameters
- no formal CDISC conformance validator
- no approved Define-XML
- no independent production-programming validation
- no regulatory submission or agency review

---

## Interpretation

The primary high-dose-versus-placebo estimate was positive but imprecise:

**+1.68 U/L, 95% CI −0.65 to 4.01**

The low-dose estimate was larger than the high-dose estimate, so the results do not demonstrate a clear monotonic dose-response pattern.

Sensitivity analyses produced different estimates under different modeling assumptions, emphasizing the importance of model choice, data provenance, influential observations and transparent reporting.

---

## Skills Demonstrated

This project demonstrates practical experience with:

- clinical trial data
- CDISC SDTM concepts
- ADaM-style derivations
- ADSL
- ADLB
- treatment reconciliation
- visit and baseline derivation
- ANCOVA
- estimated marginal means
- sensitivity analysis
- residual diagnostics
- R programming
- programmatic QC
- source-to-analysis traceability
- laboratory-data harmonization
- SAS XPT transport
- regulatory-style TLF development
- aggregate disclosure control
- Git/GitHub
- Power BI
- reproducible analytical workflows

---

## Project Links

**Live case study**

https://nandhananambi.github.io/xanomeline-alt-safety/index.html?v=2

**GitHub repository**

https://github.com/NandhanaNambi/xanomeline-alt-safety

---

## Disclaimer

This repository is a portfolio and educational project based on public/example CDISC Pilot Study data.

The analysis, ADaM-style datasets, TLFs, XPT files and validation outputs are educational development artifacts.

They are not intended for clinical decision-making, regulatory submission or formal CDISC certification.
