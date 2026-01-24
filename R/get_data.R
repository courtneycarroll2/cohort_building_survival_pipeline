# get_data.R provides functions to deserialize a JSON file of cohort specifications and
# retrieve the clinical (and, if applicable, genomic) TCGA data requested in the JSON file.

library(UCSCXenaTools)
library(jsonlite)

source("R/testing.R")

# Deserialize JSON file to get cohort and model specifications
# requires path to JSON file as input.
get_cohort_model_specs <- function(file) {
  cohort_model_specs <- fromJSON(file)
}

# Helper function: If multiple cancer types are listed in the cohort specifications, put them into correct search format string.
# Input is the "cancer_type" object from the deserialized JSON file.
make_cancer_filter_term <- function(cancer_types = cancer_types) {
  return(paste(cancer_types, collapse = "|"))
}

# Helper function: download TCGA LUAD clinical data from UCSC Xena.
# Input is deserialized JSON file with cohort and model specifications.
# Throws an exception if the "cancer_type" object is specified in the JSON file, but is not a vector of one or more valid cancer abbreviations available in TCGA.
get_clinical_data <- function(cohort_model_specs) {
  cancer_filter_term <- make_cancer_filter_term(cohort_model_specs$cancer_type)
  
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

# Download TCGA LUAD gene-level mutation/genomic data from UCSC Xena if mutation data is needed for criteria or analysis.
# Input is deserialized JSON file with cohort and model specifications.
# Output is a dataframe gene-level mutation data, with genes in the first column and subsequent columns representing sample IDs; or NULL, if mutation data is not needed.
get_genomic_data <- function(cohort_model_specs) {
  
  # Get mutation/genomic data
  cancer_filter_term <- make_cancer_filter_term(cohort_model_specs$cancer_type)
  
  XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
    XenaFilter(filterDatasets = "mc3_gene_level") %>% 
    XenaFilter(filterDatasets = cancer_filter_term) -> genomic_data
  
  XenaQuery(genomic_data) %>%
    XenaDownload() %>%
    XenaPrepare() -> gen
  
  # get list of needed genes (either for criteria or analysis)
  needed_genes <- get_needed_genes(cohort_model_specs, gen)
  
  # if no needed genes, return NULL
  if (length(needed_genes) == 0) {
    return(NULL)
  }
  
  return(gen)
}


# Use this function to retrieve the clinical and, if desired, genomic data, based on the criteria and variables specified in the JSON. Input is a deserialized JSON file (use get_cohort_model_specs())
# Will return the clinical data, plus any genomic criteria or genomic survival variables specified.
get_data <- function(cohort_model_specs, clin, gen) {
  
  
  # If there are no genomic criteria, and no genes are listed as covariates, don't get genomic data; just return clinical.
  if (is.null(cohort_model_specs$genomic_criteria) & (length(needed_genes) == 0)) {
    return(clin)
  }
  
  else {
    gen <- get_genomic_data(cohort_model_specs)
    sample_names <- gen$sample
    gen2 <- as.data.frame(t(gen[,2:ncol(gen)]))
    colnames(gen2) <- sample_names
    gen2 <- data.frame("sampleID" = rownames(gen2), gen2)
    rownames(gen2) <- NULL
    
    gen3 <- gen2[,which(colnames(gen2) %in% c("sampleID", genomic_vars))]
    
    output_df <- merge(clin, gen3, by="sampleID")
    return(output_df)
  }
}

# Helper function to get a vector of gene columns needed, based on 1) cohort criteria, 2) survival subset specification, and 3) survival model variables.
# Input is the deserialized JSON cohort and model specifications; and genomic dataset.
# Output is a vector of gene columns.
get_needed_genes <- function(cohort_model_specs, gen) {
  # get number of genes to keep (criteria and variables)
  covars <- cohort_model_specs$model$covars
  subset_var <- get_subset(cohort_model_specs)
  needed_genes <- c()
  
  # add criteria genes (if they exist)
  if (!is.null(cohort_model_specs$genomic_criteria)){
    needed_genes <- append(needed_genes, cohort_model_specs$genomic_criteria$field)
  }
  
  # If survival model will be subset by presence/absence of a gene mutation, add the gene
  if (!is.null(subset_var)) {
    if (col_exists_in_dataset(subset_var, gen)) {
      needed_genes <- append(needed_genes, subset_var)
    }
  }
  
  # Add any survival model gene covariates
  needed_genes <- append(needed_genes, covars$field[which(covars$ds == "genomic")])
  
  return(needed_genes)
}

# Helper function to get the variable to subset the data on for survival analysis, if it exists.
# Input is the deserialized JSON cohort and model specifications.
# Output is a string containing the variable name if it exists, or NULL if no subset is specified in the JSON. 
get_subset <- function(cohort_model_specs) {

  if (!is.null(cohort_model_specs$model$subset)) {
    subset_var <- cohort_model_specs$model$subset$field
    return(subset_var)
  }
  
  return(NULL)
}




