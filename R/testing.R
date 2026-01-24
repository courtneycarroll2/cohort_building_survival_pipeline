# testing.R contains code to check if columns exist in a dataframe.


col_exists_in_dataset <- function(col_name, ds){
  if (is.null(ds)) {
    return(FALSE)
  }
  else {
    return(col_name %in% colnames(ds))
  }
}

col_exists_any_dataset <- function(col_name, clin, gen) {
  if (is.null(gen)) {
    return(col_exists_in_dataset(col_name, clin))
  }
  else {
    return(col_exists_in_dataset(col_name, clin) | col_exists_in_dataset(col_name, gen))
  }
}