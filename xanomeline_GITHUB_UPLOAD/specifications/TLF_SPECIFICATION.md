# Supplemental TLF specification (educational study)

Population: treated participants (`SAFFL == "Y"`) in three **assigned** arms (`TRT01P`), excluding screen failures. Assignment is not substituted with actual dose; known assignment/actual-treatment mismatches are disclosed in the existing report. Analysis inputs: `pharmaversesdtm::dm`, `pharmaversesdtm::ae`, `pharmaversesdtm::lb`, locally created `outputs/adam/adsl.rds`, and `outputs/adam/adlb_week2.rds`.

| Output | Definition / denominator | Review status |
|---|---|---|
| Demographics Table 1 | Each safety subject once. AGE: nonmissing n, missing n, mean, SD, median, min, max. SEX: subject count and % of total arm safety N including explicit missing category. | Focused source QC script required. |
| AE Table 2 | **All reported AEs** attached to treated subjects, irrespective of onset date. Total AE *records* and distinct participants with ≥1 AE; serious AE record and distinct participant counts use source `AESER == "Y"`, and missing AESER flags counted separately. Denominator safety N. | Not a TEAE table. No first-dose to treatment-end window or study SAP TEAE rule confirmed. |
| AE terms Table 3 | Source `AEDECOD`, if available, else `AETERM`; show actual origin, do not claim MedDRA version. Count distinct participants and records by assigned arm and term. Public release only cells with ≥5 participants; further disclosure review may be necessary. | Exploratory supplemental table. |
| Listing 1 | One row per safety subject, retained locally. Source LBSEQ pointers for baseline and Week 2, observed ALT, units, change, inclusion/exclusion, assigned and actual arm. Includes missing-visit subjects. | Contains study subject IDs; **never push to GitHub**. |

Output `outputs/tlf/public/` is aggregate-only; outputs under `outputs/tlf/private/` are local and excluded by `.gitignore`. This is a study training project, not an approved regulatory TLF package; metadata, controlled terminology, SAP, Define-XML and formal validator remain open.
