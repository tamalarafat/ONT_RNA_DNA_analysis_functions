# Function to generate fasta subset and generate gene transcript table
detect_signature_motifs <- function(file_paths, config_file_path){
  # Load conditions and effect direction from config.yaml
  config <- yaml::read_yaml(config_file_path)
  
  # Load the merged table
  temp_df = loadRData(file_paths)
  
  # Get motif search parameters
  analyte_type <- config$analyte_type
  
  # Get the motif pattern from the cofig file
  signature_motif_pattern = config$signature_motif_pattern
  
  # Motif column name
  temp_motif_col = paste0("motif_", analyte_type)
  
  # Loop through each signature motif and assign "YES"/"NO" based on matching rows
  for (i in seq_along(signature_motif_pattern)) {
    # Get the signature name and motif pattern
    signature_name <- names(signature_motif_pattern)[i]
    signature_pattern <- paste0("^", signature_motif_pattern[[i]], "$")
    
    # Create a logical vector indicating matches
    matching_rows <- grepl(signature_pattern, temp_df[[temp_motif_col]])
    
    # Directly assign "YES" and "NO" using ifelse (vectorized)
    temp_df[[signature_name]] <- ifelse(matching_rows, "YES", "NO")
  }
  
  return(temp_df)
}