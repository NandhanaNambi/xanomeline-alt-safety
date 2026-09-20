# From project root; requires existing outputs/adam/adsl.rds and adlb_week2.rds.
stopifnot(file.exists('outputs/adam/adsl.rds'),file.exists('outputs/adam/adlb_week2.rds'))
source('R/11_build_tlf_package.R')
source('R/12_qc_tlf_package.R')
source('R/14_publish_aggregate_tlfs.R')
