generate_simulated_replicates <- function(input_file,
                                          cell_line = NULL,
                                          proportion_to_permutate = 0.05,
                                          variable_std = 0.05,
                                          num_replicates = 2,
                                          write_output = TRUE,
                                          store_dir = NULL,
                                          set_seed = c(1, 42, 5, 23, 17, 49, 52, 33, 11, 24),
                                          store_folder = "08_simulated_counts"){
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
  
  
  if (!missing(cell_line)){
    temp_all_samples = extract_cell_line_count(input_file = temp_all_samples, cell_line_name = cell_line, write_output = FALSE)
  }
  
  # Filter the data- remove rows containing only 0's
  keep = (!rowSums(temp_all_samples$counts != 0) == 0)
  
  # Colnames in the original count matrix
  temp_colnames = colnames(temp_all_samples$counts)
  
  # Proportion to permutate or proportion of null genes
  proportion_genes_to_permutate = (0.1 - proportion_to_permutate)
  
  # Subset the counts
  temp_all_samples$counts = temp_all_samples$counts[keep, ]
  
  # Assign the counts to a variable
  temp_count = temp_all_samples$counts
  
  # Convert the column values to integer
  temp_count <- apply(temp_count, c(1, 2), as.integer)
  
  # Subset the abundances
  temp_all_samples$abundance = temp_all_samples$abundance[keep, ]
  
  # Assign the abundances to a variable
  temp_abundance = temp_all_samples$abundance
  
  # Convert the column values to integer
  temp_abundance <- apply(temp_abundance, c(1, 2), as.integer)
  
  # Subset the length
  temp_all_samples$length = temp_all_samples$length[keep, ]
  
  # Assign the abundances to a variable
  temp_length = temp_all_samples$length
  
  # Convert the column values to integer
  temp_length <- apply(temp_length, c(1, 2), as.integer)
  
  # Create lists to store the tables
  temp_list_counts = list()
  temp_list_abundances = list()
  temp_list_length = list()
  
  
  for (i in c(1:num_replicates)){
    
    # Set the seed to make the results reproducible
    set.seed(set_seed[i])
    
    ## Counts 
  
    # Generate the thinning object of counts
    thout_count <- thin_2group(mat = temp_count, 
                         prop_null = proportion_genes_to_permutate, 
                         signal_fun = stats::rnorm,
                         signal_params = list(mean = 0, sd = variable_std))
    
    # Let's check the original counts and simulated counts
    temp_permutated_count = thout_count[["mat"]]
    
    colnames(temp_permutated_count) <- str_c("rep", (i + 1), colnames(temp_count), sep = "_")
    
    rownames(temp_permutated_count) <- rownames(temp_count)
    
    temp_list_counts[[i]] <- temp_permutated_count
    
    names(temp_list_counts)[i] <- str_c("rep", (i + 1), sep = "_")
    
    ## Abundances
    thout_abundance <- thin_2group(mat = temp_abundance, 
                                   prop_null = proportion_genes_to_permutate, 
                                   signal_fun = stats::rnorm,
                                   signal_params = list(mean = 0, sd = variable_std))
    
    # Let's check the original counts and simulated counts
    temp_permutated_abundance = thout_abundance[["mat"]]
    
    colnames(temp_permutated_abundance) <- str_c("rep", (i + 1), colnames(temp_count), sep = "_")
    
    rownames(temp_permutated_abundance) <- rownames(temp_count)
    
    temp_list_abundances[[i]] <- temp_permutated_abundance
    
    names(temp_list_abundances)[i] <- str_c("rep", (i + 1), sep = "_")
    
    ## Length
    
    thout_length <- thin_2group(mat = temp_length, 
                                prop_null = proportion_genes_to_permutate, 
                                signal_fun = stats::rnorm,
                                signal_params = list(mean = 0, sd = variable_std))
    
    # Let's check the original counts and simulated counts
    temp_permutated_length = thout_length[["mat"]]
    
    colnames(temp_permutated_length) <- str_c("rep", (i + 1), colnames(temp_count), sep = "_")
    
    rownames(temp_permutated_length) <- rownames(temp_count)
    
    temp_list_length[[i]] <- temp_permutated_length
    
    names(temp_list_length)[i] <- str_c("rep", (i + 1), sep = "_")

  }
  
  # Original coutn matrix: change the colnames to add replicate information 
  colnames(temp_count) <- str_c("rep", 1, colnames(temp_count), sep = "_")
  
  temp_counts = cbind(temp_count, do.call(cbind, temp_list_counts))
  
  temp_counts <- apply(temp_counts, c(1, 2), as.numeric)
  
  temp_all_samples$counts = temp_counts
  
  # Original coutn matrix: change the colnames to add replicate information 
  colnames(temp_abundance) <- str_c("rep", 1, colnames(temp_abundance), sep = "_")
  
  temp_abundances = cbind(temp_abundance, do.call(cbind, temp_list_abundances))
  
  temp_abundances <- apply(temp_abundances, c(1, 2), as.numeric)
  
  temp_all_samples$abundance = temp_abundances
  
  # Original coutn matrix: change the colnames to add replicate information 
  colnames(temp_length) <- str_c("rep", 1, colnames(temp_length), sep = "_")
  
  temp_lengths = cbind(temp_length, do.call(cbind, temp_list_length))
  
  temp_lengths <- apply(temp_lengths, c(1, 2), as.numeric)
  
  temp_all_samples$length = temp_lengths
  
  if (write_output == TRUE){
    save(temp_all_samples, file = str_c(temp_dir, str_c(temp_colnames, collapse = "_"), str_c("_simulated", num_replicates, "replicates_gene_counts", sep = "_"), ".RData"))
  }
  
  else {
    return(temp_all_samples)
  }

}
