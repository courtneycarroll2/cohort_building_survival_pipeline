# get_data.R provides functions to deserialize a JSON file of cohort specifications and
# retrieve the clinical (and, if applicable, genomic) TCGA data requested in the JSON file.

library(UCSCXenaTools)
library(jsonlite)

# Deserialize JSON file to get cohort specifications
# requires path to JSON file as input.
get_cohort_specs <- function(file) {
  cohort_specs <- fromJSON(file)
}

# Helper function: If multiple cancer types are listed in the cohort specifications, put them into correct search format string
make_cancer_filter_term <- function(cancer_types = cancer_types) {
  return(paste(cancer_types, collapse = "|"))
}

# Helper function: download TCGA LUAD clinical data from UCSC Xena
get_clinical_data <- function(cohort_specs) {
  cancer_filter_term <- make_cancer_filter_term(cohort_specs$cancer_type)
  
  XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
    XenaFilter(filterDatasets = "clinical") %>% 
    XenaFilter(filterDatasets = cancer_filter_term) -> clin_data
  if (typeof(clin_data) == "character"){
    stop("Error: cancer_types must be a vector of one or more valid cancer abbreviations available in TCGA. e.g., c('LUAD'), c('UCS', 'UCEC'), etc.")
  }
  
  XenaQuery(clin_data) %>%
    XenaDownload() %>%
    XenaPrepare() -> clin
  return(clin)
}

# Helper function: Download TCGA LUAD gene-level mutation data from UCSC Xena
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


# Use this function to retrieve the clinical and, if desired, genomic data, based on the criteria and variables specified in the JSON. Input is a deserialized JSON file (use get_cohort_specs())
# Will return the clinical data, plus any genomic criteria or genomic survival variables specified.
get_data <- function(cohort_specs) {
  clin <- get_clinical_data(cohort_specs)
  
  # get number of genes to keep (criteria and variables)
  covars <- cohort_specs$model$covars
  genomic_vars <- c()
  genomic_vars <- append(genomic_vars, cohort_specs$genomic_criteria$field)
  genomic_vars <- append(genomic_vars, covars$field[which(covars$ds == "genomic")])
  
  # If there are no genomic criteria, and no genes are listed as covariates, don't get mutation data; just return clinical.
  if (is.null(cohort_specs$genomic_criteria) & (length(genomic_vars) == 0)) {
    return(clin)
  }
  
  else {
    mut <- get_mutation_data(cohort_specs)
    sample_names <- mut$sample
    mut2 <- as.data.frame(t(mut[,2:ncol(mut)]))
    colnames(mut2) <- sample_names
    mut2 <- data.frame("sampleID" = rownames(mut2), mut2)
    rownames(mut2) <- NULL
    
    mut3 <- mut2[,which(colnames(mut2) %in% c("sampleID", genomic_vars))]
    
    output_df <- merge(clin, mut3, by="sampleID")
    return(output_df)
  }
}




