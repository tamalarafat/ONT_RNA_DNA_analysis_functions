extract_simulated_deg_deseq <- function(input_file,
                              reference_group = NULL,
                              control_for_cell_line = TRUE,
                              filter_gene_exp_threshold = 10,
                              filter_exp_across_n_sample = NULL,
                              write_output = FALSE,
                              write_deseq_object = FALSE,
                              store_dir = NULL,
                              store_folder = "05_DEG_files"){
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  # Check if the input_file is a character string containing the file path or directory path
  if (class(input_file) == "character") {
    # Check if the input_file exists; if yes, load the input_file from the character vector containing the file
    if (file_test("-f", input_file) == TRUE){
      # Load the gene count table
      temp_gct <- loadRData(input_file)
      
      # Generate a sample table to track the coldata
      coldata <- data.frame(names = factor(colnames(temp_gct[["counts"]])), 
                            cell = factor(sub(pattern = "_.*", replacement = "", colnames(temp_gct[["counts"]]))),
                            condition = factor(sub(pattern = ".*_", replacement = "", colnames(temp_gct[["counts"]]))))
      
      # Assign rownames to the coldata
      rownames(coldata) <- coldata$names
      
      # Assign the reference group: the group against which comparison will be performed
      if (!missing(reference_group)){
        if (levels(coldata$condition)[2] != reference_group){
          coldata$condition %<>% relevel(levels(coldata$condition)[2])
        }
      }
      
      # Initialize the DESeq object - controlling for cell line effect
      if (control_for_cell_line == TRUE){
        dds <- DESeqDataSetFromTximport(txi = temp_gct, 
                                        colData = coldata, 
                                        design = ~ cell + condition)
      }
      # Initialize the DESeq object - without controlling for cell line effect
      else {
        dds <- DESeqDataSetFromTximport(txi = temp_gct, 
                                        colData = coldata, 
                                        design = ~ condition)
      }
      
    }
    else {
      paste0("File name is missing. Provide the file path containing the file name.")
    }
  }
  
  # Check if the input file is a list object that was generated using the "generate_gene_count_table()" from the transcript expression read counts
  else if (class(input_file) %in% c("list") == TRUE){
    # Assign the gene count object to a variable
    temp_gct = input_file
    
    # Generate a sample table to track the coldata
    coldata <- data.frame(names = factor(colnames(temp_gct[["counts"]])), 
                          cell = factor(sub(pattern = "_.*", replacement = "", colnames(temp_gct[["counts"]]))),
                          condition = factor(sub(pattern = ".*_", replacement = "", colnames(temp_gct[["counts"]]))))
    
    # Assign rownames to the coldata
    rownames(coldata) <- coldata$names
    
    # Assign the reference group: the group against which comparison will be performed
    if (!missing(reference_group)){
      if (levels(coldata$condition)[2] != reference_group){
        coldata$condition %<>% relevel(levels(coldata$condition)[2])
      }
    }
    
    # Initialize the DESeq object - controlling for cell line effect
    if (control_for_cell_line == TRUE){
      dds <- DESeqDataSetFromTximport(txi = temp_gct, 
                                      colData = coldata, 
                                      design = ~ cell + condition)
    }
    # Initialize the DESeq object - without controlling for cell line effect
    else {
      dds <- DESeqDataSetFromTximport(txi = temp_gct, 
                                      colData = coldata, 
                                      design = ~ condition)
    }
  }
  
  # Check if the input file is matrix or array containing the count table
  else if (class(input_file) %in% str_c(c("matrix", "array"), collapse = ", ") == TRUE){
    # Assign the count matrix to a variable
    temp = input_file
    
    # Generate a sample table to track the coldata
    coldata <- data.frame(names = factor(colnames(temp)), 
                          cell = factor(sub(pattern = "_.*", replacement = "", colnames(temp))),
                          condition = factor(sub(pattern = ".*_", replacement = "", colnames(temp))))
    
    # Assign rownames to the coldata
    rownames(coldata) <- coldata$names
    
    # Assign the reference group: the group against which comparison will be performed
    if (!missing(reference_group)){
      if (levels(coldata$condition)[2] != reference_group){
        coldata$condition %<>% relevel(levels(coldata$condition)[2])
      }
    }
    
    # Initialize the DESeq object - controlling for cell line effect
    if (control_for_cell_line == TRUE){
      dds <- DESeqDataSetFromMatrix(countData = round(temp),
                                    colData = coldata,
                                    design = ~ cell + condition)
    }
    # Initialize the DESeq object - without controlling for cell line effect
    else {
      dds <- DESeqDataSetFromMatrix(countData = round(temp),
                                    colData = coldata,
                                    design = ~ condition)
    }
  }
  
  else {
    paste0("Error: The input_file can take a character vector containing the file path, or a list object generated usign the generate_gene_count_table() function, or a count matrix. Please check the provided input_file.")
  }
  
  # Filter the data
  # criteria 1: remove rows containing only 0's
  temp_cr_1 = (!rowSums(counts(dds) != 0) == 0)
  
  # criteria 2:remove rows if the gene expression is less than 10 (default) or user-specified in "filter_exp_across_n_sample" args
  if (missing(filter_exp_across_n_sample)){
    temp_cr_2 = rowSums(counts(dds) >= 10) == ncol(counts(dds))
  }
  
  else{
    temp_cr_2 = (rowSums(counts(dds) >= 10) == filter_exp_across_n_sample)
  }
  
  # combine the two criteria
  keep = (temp_cr_1 & temp_cr_2)
  
  # Subset the deseq object to keep only those genes that satisfy the filtering criteria
  dds <- dds[keep, ]
  
  # Estimate the size factor of the dds object
  dds <- estimateSizeFactors(dds)
  
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  
  # Running the differential expression pipeline
  dds <- DESeq(dds)
  
  # Extract and store the differentially expressed genes result in a variable
  res <- results(dds)
  
  # Store the result file
  if (write_output == TRUE){
    if (!dir.exists(str_c(store_dir, "/", store_folder))){
      dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # storing the directory information in a temporary variable
    temp_dir = str_c(store_dir, "/", store_folder, "/")
    save(res, file = str_c(temp_dir, "DEG_DESeq_", paste0(rownames(coldata), collapse = "_"), ".RData"))
    
    # Write DESeq object file- useful for plotting the figures
    if (write_deseq_object == TRUE){
      save(dds, file = str_c(temp_dir, "DESeq_object_", paste0(rownames(coldata), collapse = "_"), ".RData"))
    }
  }
  # Return the output file - (default)
  else {
    # Create an empty list to store the DESeq object and DEG results file
    temp_list = list()
    
    # Add objects to the list
    temp_list[[1]] <- res # results; differentially expressed genes
    temp_list[[2]] <- dds # DESeq object; 
    
    # Assign names to the list items
    names(temp_list)[c(1, 2)] <- c("res", "dds")
    
    return(temp_list)
  }
}
