# Function to annotate the mod sites
signature_motif_data_per_condition <- function(file_paths, 
                                               config_file_path){
  # Load conditions and effect direction from config.yaml
  config <- yaml::read_yaml(config_file_path)
  
  # Load configuration again to access signature motifs (already loaded above, but repeated for modularity)
  temp_analyte_id <- config$analyte_type
  
  # Access the ordered effect sequence - direction to calculate the effect of treatment
  effect_order <- config$effect
  
  # Get the names of the signature motif sequences - DRACH and/or any other
  temp_signature_names <- names(config$signature_motif_pattern)
  
  # Load the merged table
  temp_df = loadRData(file_paths)
  
  # Set the first condition in effect_order is the reference (e.g., "control")
  reference_group <- effect_order[1]
  
  temp_diff_mod_cols = colnames(temp_df)[grep("diff_mod_rate_", colnames(temp_df))]
  
  temp_motif_cols = paste0("motif_", temp_analyte_id)
  
  temp_group_list = list()
  
  for (i in 2:length(effect_order)){
    
    ## inside the loop
    temp_treat_group <- effect_order[i]
    
    # Define the new column name for the differential calculation
    diff_col_name <- paste0("diff_mod_rate_", reference_group, "_to_", temp_treat_group)
    
    # Motif sequences in the reference condition
    temp_ref_df <- temp_df[temp_df[[diff_col_name]] > 0, ]
    
    # how many k-mers were detected
    temp_ref_df = temp_ref_df[!is.na(temp_ref_df[[temp_motif_cols]]), ]
    
    # Motif sequences in the treatment condition
    temp_trt_df <- temp_df[temp_df[[diff_col_name]] < 0, ]
    
    # how many k-mers were detected
    temp_trt_df = temp_trt_df[!is.na(temp_trt_df[[temp_motif_cols]]), ]
    
    # Create an empty list to store the data frames
    temp_list <- list()
    
    # Loop through each signature name
    for (j in seq_along(temp_signature_names)) {
      
      # Get the rows containing signature motif for reference and treatment groups
      ref_df <- temp_ref_df[temp_ref_df[[temp_signature_names[j]]] == "YES", ]
      trt_df <- temp_trt_df[temp_trt_df[[temp_signature_names[j]]] == "YES", ]
      
      # Define names for each subset
      temp_ref_name <- paste0(temp_signature_names[j], "_", reference_group, "_motif")
      temp_trt_name <- paste0(temp_signature_names[j], "_", temp_treat_group, "_motif")
      
      # Add non-empty data frames to temp_list
      if (nrow(ref_df) > 0) temp_list[[temp_ref_name]] <- ref_df
      if (nrow(trt_df) > 0) temp_list[[temp_trt_name]] <- trt_df
    }
    
    # Add the list to the group list; needed when you have more than two conditions
    temp_ind = i - 1
    
    temp_group_list[[temp_ind]] <- temp_list
    
    names(temp_group_list)[temp_ind] <- paste0("conditon_comparison_", reference_group, "_to_", temp_treat_group)
  }
  
  temp_group_list[["col_id"]] <- temp_motif_cols
  
  temp_group_list[["analyte_type"]] <- temp_analyte_id
  
  return(temp_group_list)
}