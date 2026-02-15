# qc.R contains code to do quality control on the data and report missing rates of variables.
# Also plan to add multiple imputation


# Identifies samples with >X% genomic/mutation data missing; percentage specified in deserialized JSON specs file under cohort_model_specs$qc$proportion_missing_rm_variable, if 1) genomic_data is not NULL, and 2) cohort_model_specs$qc$proportion_missing_rm_variable exists.
# If the two criteria above are met, this function returns a list containing 1) a vector of sample IDs to retain and 2) stats to add to a qc file. Otherwise, it returns NULL.
qc_genomic_ids <- function(genomic_data, cohort_model_specs) {
  # Find samples missing < specified proportion of genomic data
  if (!is.null(cohort_model_specs$qc$proportion_missing_rm_variable) & !is.null(genomic_data)) {
    pct <- cohort_model_specs$qc$proportion_missing_rm_variable
    keep_indices <- (rowSums(is.na(genomic_data[,2:ncol(genomic_data)]))/(ncol(genomic_data) - 1)) <= pct
    keep_samples <- genomic_data$sampleID[keep_indices]
    
    # qc output:
    stage <- paste0("rm_missing_genomic_gt_", pct, "_pct")
    n_in <- nrow(genomic_data)
    n_keep <- length(keep_samples)
    n_rm <- n_in - n_keep
    pct_rm <- round((n_rm / n_in) * 100, 2)
    description <- paste0("removed samples missing > ", pct*100 , "% genomic data")
    qc_row <- c(stage, n_in, n_keep, n_rm, pct_rm, description)
    return(
      list(
        keep_samples = keep_samples,
        qc_row = qc_row
      )
    )
  } 
  
  # otherwise, return NULL
  else {
    return(NULL)
  }
}




