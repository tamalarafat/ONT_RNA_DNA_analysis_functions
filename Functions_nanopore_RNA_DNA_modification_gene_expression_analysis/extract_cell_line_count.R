extract_cell_line_count <- function(input_file,
                                    cell_line_name,
                                    write_output = TRUE,
                                    store_dir = NULL,
                                    store_folder = "07_cell_line_count_table"){
  
  if (write_output == TRUE){
    if (missing(store_dir)){
      store_dir = getwd()
    }
    
    if (!dir.exists(str_c(store_dir, "/", store_folder))){
      dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # storing the directory information in a temporary variable
    temp_dir = str_c(store_dir, "/", store_folder, "/")
  }
  
  # Check if the input_file is a character string containing the file path or directory path
  if (class(input_file) == "character") {
    
    # Check if the input_file exists; if yes, load the input_file from the character vector containing the file
    if (file_test("-f", input_file) == TRUE){
      # Load the gene count table
      temp_all_samples <- loadRData(input_file)
    }
    else {
      paste0("File name is missing. Provide the file path containing the file name.")
    }
  }
  
  else if (class(input_file) %in% c("list") == TRUE){
    temp_all_samples = input_file
  }
  
  else {
    paste0("Error: The input_file can take a character vector containing the file path of the gene counts list object, or directly the list object generated using generate_gene_count_table(). Please check the provided input_file.")
  }
  
  # Create an empty list to store results
  temp_list <- list()

  # Extract the columns corresponding to each cell line
  temp_columns <- grep(cell_line_name, colnames(temp_all_samples$counts))
  
  if (length(temp_columns) == 0){
    paste0("Invalid cell line name.")
  }
  
  else {
    # Subset each element in the original list based on identified columns
    temp_list$abundance <- temp_all_samples$abundance[, temp_columns]
    temp_list$counts <- temp_all_samples$counts[, temp_columns]
    temp_list$length <- temp_all_samples$length[, temp_columns]
    temp_list$countsFromAbundance <- temp_all_samples$countsFromAbundance
    
    if (write_output == TRUE){
      save(temp_list, file = str_c(temp_dir, paste0(cell_line_name, "samples_gene_count_file", collapse = "_"), ".RData"))
    }
    
    else {
      return(temp_list)
    }
  }
}
