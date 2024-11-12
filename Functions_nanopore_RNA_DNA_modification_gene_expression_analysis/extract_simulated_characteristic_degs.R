extract_simulated_characteristic_degs <- function(input_file,
                                               filter_by_pval = TRUE,
                                               filter_by_fc = TRUE,
                                               return_sig_by_pval = TRUE,
                                               return_sig_by_padj = FALSE,
                                               log2fc_threshold = 1,
                                               alpha_threshold = 0.05,
                                               padj_threshold = 0.05,
                                               return_result_per_condition = TRUE,
                                               return_geneIDs = FALSE,
                                               write_output = FALSE,
                                               tag_output_name = NULL,
                                               store_dir = NULL,
                                               store_folder = "06_Sample_characteristic_DEG"){
  
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
  
  # If the input is a character vector
  if (class(input_file) == "character") {
    # Check if provided path contains the file - DESeq2 results file
    
    if (all(file_test("-f", input_file) == TRUE)){
      # Load the gene count table
      res <- loadRData(input_file)
    }
    else {
      paste0("File name is missing. Provide the file path containing the file name.")
    }
  }
  
  else if (class(input_file) %in% c("DESeqResults") == TRUE){
    res = input_file
  }
  
  else {
    paste0("Error: The input_file can take a character vector containing the file path, or a DESeqResults object generated either using DESeq() or extract_deg_deseq(). Please check the provided input_file.")
  }
  
  # Criteria 1: Identify top differentially expressed genes based on pval
  temp_pval = which(res$pvalue < alpha_threshold)
  
  # Criteria 2: Identify top differentially expressed genes based on log2FC
  temp_fc = which(abs(res$log2FoldChange) > log2fc_threshold)
  
  # Combine the filtering criteria
  if (all(filter_by_pval, filter_by_fc) == TRUE){
    temp_ind = union(temp_pval, temp_fc)
  }
  
  else if (filter_by_pval == TRUE & filter_by_fc == FALSE){
    temp_ind = temp_pval
  }
  
  else if (filter_by_fc == TRUE & filter_by_pval == TRUE){
    temp_ind = temp_fc
  }
  
  # Create a dataframe with the genes that passed the filtering criteria
  temp_df = data.frame(ID = rownames(res)[temp_ind], log2FC = res$log2FoldChange[temp_ind], pvalue = res$pvalue[temp_ind], padj = res$padj[temp_ind], lfcSE = res$lfcSE[temp_ind], stat = res$stat[temp_ind])
  
  # Assign rownames to gene IDs
  rownames(temp_df) <- temp_df$ID
  
  # Assign group or condition tags to the foldchange direction
  
  # Get group/condition information
  temp_string <- mcols(res)[2, 2]
  
  # Extract the strings before and after "vs"
  temp_matches <- regmatches(temp_string, regexec("([^ ]+) vs ([^ ]+)", temp_string))
  
  # Numerator or the reference group
  temp_numerator = temp_matches[[1]][2]
  
  # Denominator or the effect group where we want to observe differences
  temp_denominator = temp_matches[[1]][3]
  
  # Add group information to the data table
  temp_df$condition = ifelse(temp_df$log2FC > 0, temp_numerator, temp_denominator)
  
  # Add direction of the foldchange information to the data table
  temp_df[str_c(c("exp_observed_in", temp_numerator, "group"), collapse = "_")] = ifelse(temp_df$log2FC > 0, "up", "down")
  
  # Rearrange the data table based on FC
  temp_df = temp_df[order(temp_df$log2FC, decreasing = TRUE), ]
  
  # If to return only significant degs
  if (all(return_sig_by_pval, return_sig_by_padj) == TRUE){
    temp_df = temp_df[temp_df$pvalue < alpha_threshold, ]
    temp_df = temp_df[temp_df$padj < padj_threshold, ]
  }
  
  else if(return_sig_by_pval == TRUE & return_sig_by_padj == FALSE){
    temp_df = temp_df[temp_df$pvalue < alpha_threshold, ]
  }
  
  else if(return_sig_by_padj == TRUE & return_sig_by_pval == FALSE){
    temp_df = temp_df[temp_df$padj < padj_threshold, ]
  }
  
  else{
    temp_df = temp_df
  }

  # If to return by groups
  temp_list = list()
  
  # Add the numerator group
  temp_list[[1]] = temp_df[temp_df$condition == temp_numerator, ]
  
  # Add the denominator group
  temp_list[[2]] = temp_df[temp_df$condition == temp_denominator, ]
  
  # Assign names to the list items
  names(temp_list) = c(temp_numerator, temp_denominator)
  
  # Check if results should be stored or not
  if (write_output == TRUE){

    # Save the complete table
    write.table(temp_df, str_c(temp_dir, str_c(c("DEG_comparison", temp_numerator, temp_denominator, tag_output_name), collapse = "_"), ".csv"))
    
    fileGenerator(markerList = rownames(temp_df), fileName = str_c(temp_dir, str_c("ID_DEG_comparison", temp_numerator, temp_denominator, tag_output_name, sep = "_"), ".txt"))

    # Check if results should be stored per condition
    if (return_result_per_condition == TRUE){
      
      for (i in c(1:length(temp_list))){
        write.table(temp_list[[i]], str_c(temp_dir, str_c(c("DEG_up_in", names(temp_list)[i], "group", tag_output_name), collapse = "_"), ".csv"))
        
        # Save gene IDs if specified
        if (return_geneIDs == TRUE){
          fileGenerator(markerList = rownames(temp_list[[i]]), fileName = str_c(temp_dir, str_c("ID_DEG_up_in", names(temp_list)[i], "group", tag_output_name, sep = "_"), ".txt"))
        }
      }
    }
  }
    
    else {
      
      # Add the denominator group
      temp_list[[3]] = temp_df
      
      # Assign names to the list items
      names(temp_list)[3] = str_c(temp_numerator, temp_denominator, sep = "_")
      
      return(temp_list)
    }
  
} 
