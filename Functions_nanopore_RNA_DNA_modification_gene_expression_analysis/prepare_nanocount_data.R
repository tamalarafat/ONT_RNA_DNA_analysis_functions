# Function to prepare and process nanocount output table
prepare_nanocount_data <- function(file_paths, sample_names){
  
  # Let's create a directory to store the GO annotation files
  if (!dir.exists("temporary_data")){
    dir.create("temporary_data", showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = "temporary_data/"
  
  # Loop through remaining files and merge iteratively
  for (i in 1:length(file_paths)) {
    
    file_name = basename(file_paths[i])
    
    # Find matches using grep
    matching_sample <- sample_names[sapply(sample_names, function(x) grepl(x, file_name))]
    
    # read the file, in tsv format
    temp = read_tsv(file_paths[i], show_col_types = FALSE)
    
    # rearrange the columns to a desired format
    temp = temp[, c("transcript_name", "transcript_length", "transcript_length", "est_count", "tpm")]
    
    # assign column names to the count table
    colnames(temp) = c("target_id", "length", "eff_length", "est_counts", "tpm")
    
    temp$target_id <- sub(pattern = "\\..*", replacement = "", temp$target_id)
    
    # Sort the data frames using numeric extraction from target_id
    temp <- temp[order(as.numeric(sub("ENST", "", temp$target_id))), ]
    
    # write the table in tsv format
    write.table(temp, file = paste0(temp_dir, matching_sample, ".tsv"), row.names = FALSE, sep = "\t")
  }
}