metaGrouping <- function(data) {
  
  data %>%
    group_by(Sample_Name) %>%
    count()
  
}

