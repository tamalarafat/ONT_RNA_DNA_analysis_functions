# Function to annotate the mod sites
annotate_genes <- function(file_paths, 
                           config_file_path, 
                           gene_transcript_table){
  # Load conditions and effect direction from config.yaml
  config <- yaml::read_yaml(config_file_path)
  
  # Load configuration again to access signature motifs (already loaded above, but repeated for modularity)
  temp_seq_id <- config$assign_gene_ID_type
  
  # Load the merged table
  temp_df = read_tsv(file_paths, show_col_types = FALSE)
  
  # Load symbols and ids table
  temp_symbol_df = loadRData(gene_transcript_table)
  
  # Add Ensembl gene, entrez, and hgnc symbol to the modification table
  temp_mod_ids = temp_df$ID
  
  temp_df$hgnc_symbol = temp_symbol_df$hgnc_symbol[match(temp_df$ID, temp_symbol_df[[temp_seq_id]])]
  
  temp_df$entrezgene_id = temp_symbol_df$entrezgene_id[match(temp_df$ID, temp_symbol_df[[temp_seq_id]])]
  
  temp_df$ensembl_gene_id = temp_symbol_df$ensembl_gene_id[match(temp_df$ID, temp_symbol_df[[temp_seq_id]])]
  
  # Condition names in the DEG table
  temp_conditions <- unique(temp_df$condition)
  
  temp_group_list = list()
  
  for (i in 1:length(temp_conditions)){
    
    # Genes with higher modification rate in the reference condition
    temp_genes <- unique(str_sort(temp_df$entrezgene_id[temp_df$condition == temp_conditions[i]], numeric = TRUE))
    
    temp_group_list[[temp_conditions[i]]] <- temp_genes
  }
  
  # temp_name
  temp_name <- paste0("GSEA_", paste0(names(temp_group_list), collapse = "_and_"))
  
  temp_list <- list()
  temp_list[[temp_name]] <- temp_group_list
  
  temp_list[["genes_table"]] <- temp_df
  
  # Load all entrez ids of the hg
  hg_entrz = as.character(str_sort(temp_symbol_df$entrezgene_id, numeric = TRUE))
  
  hg_entrz = unique(hg_entrz[!is.na(hg_entrz)])
  
  temp_list[["entrezgene_id"]] <- hg_entrz
  
  return(temp_list)
}