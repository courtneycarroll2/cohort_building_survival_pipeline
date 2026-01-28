# helper.R contains helper functions for the other scripts.


# Checks whether a specified column name exists in the specified dataset.
# Returns a boolean.
col_exists_in_dataset <- function(col_name, ds){
  if (is.null(ds)) {
    return(FALSE)
  }
  else {
    return(col_name %in% colnames(ds))
  }
}

# Checks whether a specified column name exists in either the clinical or genomic/mutation dataset
# Returns a boolean.
col_exists_any_dataset <- function(col_name, clin, gen) {
  if (is.null(gen)) {
    return(col_exists_in_dataset(col_name, clin))
  }
  else {
    return(col_exists_in_dataset(col_name, clin) | col_exists_in_dataset(col_name, gen))
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