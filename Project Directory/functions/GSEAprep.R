prepare_data <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  condition <- tryCatch({
    read.csv(file_path)
  }, error = function(e) {
    stop("Error reading file: ", e$message)
  })
  
  conditionRefined <- condition[c(1,3)]
  
  if (ncol(conditionRefined) != 2) {
    stop("conditionRefined is not a two-column data frame.")
  }
  
  conditionDeframed <- deframe(conditionRefined)
  
  if (length(conditionDeframed) > 0) {
    return(conditionDeframed)
  } else {
    stop("Input data is empty.")
  }
}
