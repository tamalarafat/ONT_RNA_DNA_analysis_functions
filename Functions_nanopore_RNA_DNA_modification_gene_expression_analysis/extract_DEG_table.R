# Function to extract differential gene expression table
extract_DEG_table <- function(res_DEseq){
  
  # Create a dataframe with the genes that passed the filtering criteria
  temp_df = data.frame(ID = rownames(res_DEseq), log2FC = res_DEseq$log2FoldChange, pvalue = res_DEseq$pvalue, padj = res_DEseq$padj, lfcSE = res_DEseq$lfcSE, stat = res_DEseq$stat)
  
  # Assign rownames to gene IDs
  rownames(temp_df) <- temp_df$ID
  
  # Assign group or condition tags to the foldchange direction
  
  # Get group/condition information
  temp_string <- mcols(res_DEseq)[2, 2]
  
  # Extract the strings before and after "vs"
  temp_matches <- regmatches(temp_string, regexec("([^ ]+) vs ([^ ]+)", temp_string))
  
  # Numerator or the reference group
  temp_numerator = temp_matches[[1]][2]
  
  # Denominator or the effect group where we want to observe differences
  temp_denominator = temp_matches[[1]][3]
  
  # Column name
  temp_col <- paste0("observed_expression_in_", temp_numerator)
  
  # Add group information to the data table
  temp_df[[temp_col]] = ifelse(temp_df$log2FC > 0, "UP", "DOWN")
  
  return(temp_df)
}