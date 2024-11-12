# Recursive merging function for a list of data frames
merge_all_samples <- function(file_paths, 
                              merge_by){

  # Load the first file as the starting data frame
  merged_data <- loadRData(file_paths[1])
  
  # Loop through remaining files and merge iteratively
  for (i in 2:length(file_paths)) {
    next_data <- loadRData(file_paths[i])
    merged_data <- merge(merged_data, next_data, by = merge_by, all = FALSE)  # Only common sites
  }
  
  return(merged_data)
}
