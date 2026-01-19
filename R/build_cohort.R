# The functions in build_cohort.R filter the dataset to consist of individuals meeting all specified cohort criteria.

library(dplyr)

# Filters the dataset to meet cohort criteria Takes in dataframe (output of get_data()) and the deserialized JSON cohort criteria and returns the appropriate cohort.
build_cohort <- function(cohort_specs, combined_data) {
  criteria_data <- get_criteria_data(cohort_specs)
  output_df <- combined_data
  
  for (i in 1:nrow(criteria_data)) {
    logic_term <- criteria_data$logic[i]
    
    # check that logic term is supported
    if (!(logic_term %in% c("equals", "in"))) {
      stop(sprintf("Error: does not support filter logic of %s '%s' %s; please use a supported logic term. Supported terms are 'equals' and 'in'.", criteria_data$field[i], logic_term, as.character(unlist(criteria_data$value[i]))))
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
    
  }
  
}


# Helper function: Combine the clinical and genomic criteria. Takes in the deserialized JSON cohort specifications.
get_criteria_data <- function(cohort_specs) {
  clin_criteria_exists <- !is.null(cohort_specs$clinical_criteria)
  gen_criteria_exists <- !is.null(cohort_specs$genomic_criteria)
  
  if (clin_criteria_exists & gen_criteria_exists) {
    criteria_data <- cohort_specs$clinical_criteria
    criteria_data <- rbind(criteria_data, cohort_specs$genomic_criteria)
    return(criteria_data)
  }
  
  else {
    if (clin_criteria_exists) {
      return(cohort_specs$clinical_criteria)
    }
    return(cohort_specs$genomic_criteria)
  }
  
}

