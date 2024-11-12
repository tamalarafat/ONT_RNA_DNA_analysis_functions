perform_DEA_deseq <- function(input_file,
                              config_file_path){
  
  # Load conditions and effect direction from config.yaml
  config <- yaml::read_yaml(config_file_path)
  
  # Get the conditions specified in configuration file - samples belonging to each condition
  conditions <- config$conditions
  
  # Access the ordered effect sequence - direction to calculate the effect of treatment
  effect_order <- config$effect
  
  # Reference condition to which expression will be compared
  reference_group <- effect_order[1]
  
  # Get the cell line information
  cell_line <- config$cell_line
  
  # Controlling for cell line
  control_for_cell_line <- config$control_for_cell_line_effect
  
  # Gene filtering threshold - filter out genes if at least this amount of expression is detected in each of the samples
  gene_filter_threshold <- as.numeric(config$gene_exp_threshold_filter)
  
  # Load the gene count table
  temp_gct <- loadRData(input_file)
  
  # Get the sample names
  sample_names = colnames(temp_gct[["counts"]])
  
  # Get the cell line names
  sample_cell_line <- sapply(sample_names, function(sample) {
    name <- names(cell_line)[sapply(cell_line, function(x) sample %in% x)]
    if (length(name) > 0) name else NA  # Return NA if not found
  })
  
  # Get the condition names
  sample_condition <- sapply(sample_names, function(sample) {
    name <- names(conditions)[sapply(conditions, function(x) sample %in% x)]
    if (length(name) > 0) name else NA  # Return NA if not found
  })
  
  # Generate a sample table to track the coldata
  coldata <- data.frame(names = factor(sample_names), 
                        cell = factor(sample_cell_line),
                        condition = factor(sample_condition))
  
  # Assign rownames to the coldata
  rownames(coldata) <- coldata$names
  
  # Assign the reference group: the group against which comparison will be performed
  if (levels(coldata$condition)[2] != reference_group){
    coldata$condition %<>% relevel(levels(coldata$condition)[2])
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
  
  # Filter the data
  # criteria 1: remove rows containing only 0's
  temp_cr_1 = (!rowSums(counts(dds) != 0) == 0)
  
  # criteria 2:remove rows if the gene expression is less than 10 (default) or user-specified in "filter_exp_across_n_sample" args
  temp_cr_2 = rowSums(counts(dds) >= gene_filter_threshold) == ncol(counts(dds))
  
  # combine the two criteria
  keep = (temp_cr_1 & temp_cr_2)
  
  # Subset the deseq object to keep only those genes that satisfy the filtering criteria
  dds <- dds[keep, ]
  
  # Estimate the size factor of the dds object
  dds <- estimateSizeFactors(dds)
  
  # Running the differential expression pipeline
  dds <- DESeq(dds)
  
  # Extract and store the differentially expressed genes result in a variable
  res <- results(dds)
  
  # Create an list to store the DESeq object and DEG results file
  temp_list = list(res = res, dds = dds)
  
  return(temp_list)
}