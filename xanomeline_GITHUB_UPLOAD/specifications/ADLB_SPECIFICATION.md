# ADLB derivation specification — educational candidate v0.1

## Scope

Source: CDISC Pilot Study `pharmaversesdtm::lb`; merge with project ADSL by `USUBJID`. Retain all original LB rows, **every available test and visit**, qualitative and unscheduled observations, missing numeric results, and screen-failure rows. No fabricated measurements or visits. Existing Week 2 ALT extract and its ANCOVA stay immutable.

## Key variables and rules

| Variable | Source / derivation |
|---|---|
| `STUDYID`, `USUBJID` | LB source identifiers; ADSL joins on USUBJID. |
| `SRCDOM`, `SRCSEQ`, `SRCVAR` | `LB`, `LBSEQ`, numeric result source `LBSTRESN` or character source `LBSTRESC`. Keep original LB columns. |
| `PARAMCD`, `PARAM`, `PARAMN`, `PARCAT1` | Source LBTESTCD, LBTEST with standardized unit, sorted PARAMN, LBCAT. Block inconsistent names, long codes, mixed numerical units. |
| `AVAL`, `AVALC`, `AVALU` | Numeric LBSTRESN; original LBSTRESC; LBSTRESU. Do not convert, impute, or mix units silently. |
| `ADT`, `ADY` | Complete date from LBDTC; relative to project ADSL FIRST_DOSE_DATE. No treatment day zero. Partial dates produce missing ADT and a review file. |
| `AVISIT`, `AVISITN` | LBBLFL=Y => BASELINE/0; scheduled WEEK n => source VISIT and VISITNUM; unscheduled/non-week visits retained without an analysis visit. This mapping is study-specific and needs sponsor confirmation. |
| `BASETYPE`, `ABLFL` | LAST; last flagged, numeric LB baseline by LBDTC then LBSEQ; first-dose-date conflict stops build. Multiple candidates logged. |
| `BASE`, `BNRIND` | Selected baseline measurement and reference-range indicator by subject and parameter. BASE_SRCSEQ links baseline to LBSEQ. |
| `CHG`, `PCHG` | Postbaseline AVAL-BASE; percent change 100*CHG/abs(BASE) if BASE nonzero. |
| `ANRLO`, `ANRHI`, `ANRIND`, `SHIFT1`, `R2ANRHI` | LB standard range, derive LOW/NORMAL/HIGH only with sufficient reference bounds; baseline-to-postbaseline category and result/ULN where ULN positive. **No toxicity grade or Hy's-law classification inferred.** |
| `ANL01FL` | Y for last valid numeric sample per scheduled postbaseline subject/parameter/visit by LBDTC, then LBSEQ. Retain other rows with blank flag. Requires review of repeat sample rules. |
| `TRTP`, `TRTA`, `SAFFL` | From original project ADSL. `SAFFL` is still provisional: receipt inferred from RFXSTDTC with EX presence reconciled, not a final sponsor safety flag. |

## Blocking review items before any 'submission-ready' claim

1. Confirm appropriate ADaMIG version, controlled terminology, required/permissible BDS variables and variable order/type/length against the study standards. Source LB standardization itself is not independently validated.
2. Independently review the baseline definition and time-of-dose details. Same-day sample times may matter; original code compares dates only. Verify actual exposure start/end from EX; `ONTRTFL`, treatment periods and last-on-treatment flags deliberately NOT asserted here.
3. Confirm analysis visit windowing, scheduled/unscheduled mapping, handling of repeated tests, baseline duplicate preference, nonnumeric results and missing/partial dates in an approved analysis spec.
4. Review reference-range applicability, lab units and lab-specific ULNs before reporting clinically meaningful multiples of ULN. Do not call this a liver-injury adjudication.
5. Validate one-to-one source trace, analysis flags, missingness and independent derivations; reconcile full ADLB to the original Week 2 ALT analysis without changing its estimates.
6. Produce formal variable metadata / Define-XML, controlled terminology, a version-specific conformance report (e.g., appropriate Pinnacle 21 configuration), independent QC signoff and a submission package if an actual sponsor submission is required. None is automatically produced by these scripts.

## Outputs

`outputs/full_adlb/adlb_full.rds`, `adlb_full.csv`, `parameter_dictionary.csv`, `unit_audit.csv`, `parameter_inventory.csv`, `visit_parameter_coverage.csv`, `build_summary.csv`, issue files, `legacy_alt_discrepancies.csv`, and `qc_checks.csv`. These are **local internal files** with row-level patient records; do NOT publish to the portfolio website or GitHub Pages.

## v0.2 cross-dataset reconciliation
The legacy ALT extract is subject-level with placeholders for missing Week 2 data, whereas full ADLB is source-record-level. For safety participants, compare selected Week 2 numeric observations separately from selected baseline values. Expected absent Week 2 source records are logged in `legacy_alt_expected_missing_week2.csv`, not counted as numeric discrepancies. Any actual mismatch or extra new Week 2 subject still blocks the build.


## v0.3 audited unit and result-type mapping

The numeric-source unit audit identified 9 LBTESTCDs with no LBSTRESU. Explicit
study-specific exceptions: `PH` and `SPGRAV` are dimensionless numeric parameters,
retaining numeric `AVAL` with missing `AVALU`; `ANISO`, `MACROCY`, `MICROCY`,
`POIKILO`, `POLYCHR`, `KETONES`, and `UROBIL` are handled as categorical findings.
For `KETONES`, source results are 0 (850 rows) and 1 (24 rows); for `UROBIL`,
0 (873 rows) and 1 (1 row), source `LBORRESU` is `NO UNITS`, and standardized
`LBSTRESU` is missing. The clinical meaning of 0 vs 1 is **not established**:
retain raw `LBSTRESC` in `AVALC`, original numeric code in source `LBSTRESN`,
and set numeric `AVAL` missing to avoid inappropriate means, baseline/change and
continuous modelling. Confirm code definitions from study metadata before deriving
clinical categories (negative/positive) or safety interpretations. All other
missing units and all mixed standard units remain blocking. The script saves
`unit_audit_unresolved.csv` and `urinalysis_coded_result_audit.csv` for audit.
