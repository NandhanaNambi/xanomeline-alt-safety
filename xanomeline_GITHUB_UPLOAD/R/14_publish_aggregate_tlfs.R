# Run only AFTER 11_build_tlf_package.R and 12_qc_tlf_package.R pass.
# Copy ONLY arm-level aggregates into GitHub Pages data folder.
stopifnot(requireNamespace('readr',quietly=TRUE))
qc_path <- 'outputs/tlf/private/tlf_qc_checks.csv'
if (!file.exists(qc_path)) stop('Run TLF QC first: missing ',qc_path)
qc <- readr::read_csv(qc_path,show_col_types=FALSE)
if(nrow(qc)==0L || any(is.na(qc$pass)) || !all(qc$pass)) stop('TLF QC did not pass; public export blocked')
from <- 'outputs/tlf/public'
to <- 'data/tlf'
files <- c('demographics_age_AGGREGATE.csv','demographics_sex_AGGREGATE.csv',
           'all_reported_ae_AGGREGATE.csv','TLF_METHODS_AND_LIMITATIONS.txt')
if (!all(file.exists(file.path(from,files)))) stop('Missing an expected aggregate export; run build again')
dir.create(to,recursive=TRUE,showWarnings=FALSE)
for(f in files) file.copy(file.path(from,f),file.path(to,f),overwrite=TRUE)
stopifnot(all(file.exists(file.path(to,files))))
cat('PUBLIC ARM-LEVEL AGGREGATES COPIED TO data/tlf/ AFTER QC.\n')
cat('No subject-level listing or AE term breakdown copied. Review disclosure policy before publishing.\n')
