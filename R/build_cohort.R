# The functions in build_cohort.R filter the dataset to consist of individuals meeting all specified cohort criteria.

library(dplyr)

# Filters the dataset to meet cohort criteria.
# Input: Takes in clinical data (output of get_clinical_data()), genomic data (output of get_genomic_data()) and the deserialized JSON cohort criteria
# Output: returns merged clinical and genomic data (if genomic data is not null), with samples filtered to meet cohort specifications.
# Throws: Error thrown if cohort_model_specs$clinical_criteria$logic or cohort_model_specs$genomic_criteria$logic terms are not valid. Valid terms are "equals", "in","greater_than", "less_than", "gte", "lte".
build_cohort <- function(cohort_model_specs, clinical_data, genomic_data) {
  criteria_data <- get_criteria_data(cohort_model_specs)
  output_df <- merge(clinical_data, genomic_data, by="sampleID")
  
  for (i in 1:nrow(criteria_data)) {
    logic_term <- criteria_data$logic[i]
    
    # check that logic term is supported
    if (!(logic_term %in% c("equals", "in","greater_than", "less_than", "gte", "lte"))) {
      stop(sprintf("Error: does not support filter logic of %s '%s' %s; please use a supported logic term. Supported terms are 'equals', 'in', 'greater_than', 'less_than', 'gte', and 'lte'", criteria_data$field[i], logic_term, as.character(unlist(criteria_data$value[i]))))
    }
    
    # Filter on "equals"
    if (logic_term == "equals") {
      output_df <- output_df %>%
        dplyr::filter(.data[[criteria_data$field[i]]] == criteria_data$value[[i]])
    }
    
    # Filter on "in"
    if (logic_term == "in") {
      output_df <- output_df %>%
        dplyr::filter(.data[[criteria_data$field[i]]] %in% criteria_data$value[[i]])
    }
    
    # Filter on "greater_than"
    if (logic_term == "greater_than") {
      output_df <- output_df %>%
        dplyr::filter(.data[[criteria_data$field[i]]] > criteria_data$value[[i]])
    }
    
    # Filter on "less_than"
    if (logic_term == "lte") {
      output_df <- output_df %>%
        dplyr::filter(.data[[criteria_data$field[i]]] < criteria_data$value[[i]])
    }
    
    # Filter on "gte" (greater than or equal)
    if (logic_term == "gte") {
      output_df <- output_df %>%
        dplyr::filter(.data[[criteria_data$field[i]]] >= criteria_data$value[[i]])
    }
    
    # Filter on "lte" (less than or equal)
    if (logic_term == "lte") {
      output_df <- output_df %>%
        dplyr::filter(.data[[criteria_data$field[i]]] <= criteria_data$value[[i]])
    }
  }
  
  return(output_df)
  
}


# Helper function: Combine the clinical and genomic criteria. Takes in the deserialized JSON cohort specifications.
get_criteria_data <- function(cohort_model_specs) {
  clin_criteria_exists <- !is.null(cohort_model_specs$clinical_criteria)
  gen_criteria_exists <- !is.null(cohort_model_specs$genomic_criteria)
  
  if (clin_criteria_exists & gen_criteria_exists) {
    criteria_data <- cohort_model_specs$clinical_criteria
    criteria_data <- rbind(criteria_data, cohort_model_specs$genomic_criteria)
    return(criteria_data)
  }
  
  else {
    if (clin_criteria_exists) {
      return(cohort_model_specs$clinical_criteria)
    }
    return(cohort_model_specs$genomic_criteria)
  }
  
}

