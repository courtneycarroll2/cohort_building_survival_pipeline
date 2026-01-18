library(UCSCXenaTools)
library(jsonlite)

# Deserialize JSON file to get cohort specifications
get_cohort_specs <- function(file) {
  cohort_specs <- fromJSON(file)
}

# If multiple cancer types, put them into correct search format string
make_cancer_filter_term <- function(cancer_types = cancer_types) {
  return(paste(cancer_types, collapse = "|"))
}

# download TCGA LUAD clinical data from UCSC Xena
get_clinical_data <- function(cohort_specs) {
  cancer_filter_term <- make_cancer_filter_term(cohort_specs$cancer_type)
  
  XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
    XenaFilter(filterDatasets = "clinical") %>% 
    XenaFilter(filterDatasets = cancer_filter_term) -> clin_data
  if (clin_data == "No valid cohorts or datasets find! Please check your input."){
    stop("Error: cancer_types must be a vector of one or more valid cancer abbreviations available in TCGA. e.g., c('LUAD'), c('UCS', 'UCEC'), etc.")
  }
  
  XenaQuery(clin_data) %>%
    XenaDownload() %>%
    XenaPrepare() -> clin
  return(clin)
}

# Download TCGA LUAD gene-level mutation data from UCSC Xena
get_mutation_data <- function(cohort_specs) {
  cancer_filter_term <- make_cancer_filter_term(cohort_specs$cancer_type)
  
  XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
    XenaFilter(filterDatasets = "mc3_gene_level") %>% 
    XenaFilter(filterDatasets = cancer_filter_term) -> mutation_data
  
  XenaQuery(mutation_data) %>%
    XenaDownload() %>%
    XenaPrepare() -> mut
  return(mut)
}




