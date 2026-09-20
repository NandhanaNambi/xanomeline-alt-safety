# Xanomeline ALT safety — reproducible programming portfolio

**Educational retrospective CDISC Pilot Study analysis; not a regulatory submission.**

[View the case-study website](https://nandhananambi.github.io/xanomeline-alt-safety/) · Data source: `pharmaversesdtm` R example datasets. Code includes AI-assisted development; source-derived outputs and programmatic checks require independent review before use.

## Scope

- Source audit, provisional ADSL, baseline/Week-2 ALT ADLB, ANCOVA, LS means, sensitivity analyses and report (scripts 01–05).
- Expanded 47-parameter all-visit ADLB and focused QC (06–07). Source replay, metadata/release gate audits and XPT v5 transport read-back (08–10).
- Supplemental TLFs: demographics, **all-reported adverse-event** incidence and private traceable Week-2 ALT listing, with source-based QC (11–13). **Not a treatment-emergent adverse-event analysis:** onset-window derivation and SAP rules not finalized.

## Reproduce from clean checkout

1. Install R and use the repository root as your working directory (RStudio Project or `setwd(...)`). Run `source("requirements.R")` to install missing CRAN packages and print versions. Internet is required for initial package installation. `pharmaversesdtm` provides the educational DM/EX/LB/AE/DS domains. No trial source records are committed here.
2. Run in order:

```r
source("R/01_sdtm_audit.R")
source("R/02_build_adam.R")
source("R/03_ancova_tlf.R")
source("R/04_report_qc.R")
source("R/05_sensitivity_dashboard.R")
source("R/06_build_full_adlb.R")
source("R/07_qc_full_adlb.R")
source("R/08_step10_remediate_audit.R")
source("R/09_step11_release_gates.R")
source("R/10_step12_exchange_candidate.R")
source("R/13_run_tlf_only.R") # calls 11_build and 12_qc
```

3. Review console outputs and `outputs/` for failures. The first five scripts preserve the original Week-2 results. The extra TLF build does **not** modify the ANCOVA or ADLB. Do not claim TLF QC until `R/12_qc_tlf_package.R` passes locally.

## File map

- `R/`: numbered analysis, derivation and QC scripts.
- `specifications/`: ADLB and supplemental TLF derivation specification drafts.
- `validation/`: supplementary independent-check scripts; these are additional focused checks, not formal CDISC conformance certification.
- `requirements.R`: package installation and version report; **no renv lockfile supplied**. Record versions before attempting exact reproduction.
- `website/`: pre-existing GitHub Pages website (leave it and `.github/workflows/deploy.yml` untouched). The present upload contains only code/docs additions, not the website itself.

## Exclusions and limitations

- `.gitignore` excludes `outputs/`, `.rds`, `.xpt`, etc. Do not commit `outputs/full_adlb`, the XPT candidate, patient-level dashboard CSVs, the private listing or source-key-level audits. The educational source package can be installed by the reviewer directly.
- ADSL safety flag is provisional. Treatment mismatches, baseline timing, 917 selected observations without baseline and 1,664 unassigned visits remain documented review items.
- The XPT read-back passed in the project after NA/blank comparison normalization, **not** formal CDISC validation. Define-XML, approved metadata and conformance signoff remain open.
- AE summary is **all reported** events in the provisional safety population, not a TEAE or incidence-adjusted causal safety inference.
- Expanded TLF scripts were packaged, but cannot be confirmed as executed until the project owner runs them locally and reviews the QC report.
