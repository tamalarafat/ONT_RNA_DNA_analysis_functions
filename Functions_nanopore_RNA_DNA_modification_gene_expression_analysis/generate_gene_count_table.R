# Function to generate gene count table
generate_gene_count_table <- function(input_file, config_file_path, gene_transcript_symbols){
  
  # List the files
  input_files = list.files(path = input_file, pattern = ".tsv", full.names = TRUE)
  names(input_files) <- sub(pattern = "\\..*", "", basename(input_files))
  
  # Load conditions and effect direction from config.yaml
  config <- yaml::read_yaml(config_file_path)
  
  # Get the name of the gene ID conversion type from the configuration file
  conversion_ID_type <- config[["assign_gene_ID_type"]]
  
  # Load symbols and ids table
  temp_symbol_df = loadRData(gene_transcript_symbols)
  
  # Column to retreive from gene id symbol data
  temp_col_names <- c("ensembl_transcript_id", conversion_ID_type)
  
  tx2gene = temp_symbol_df[, temp_col_names]
  
  # Assign colnames to the id transform
  colnames(tx2gene) = c("TXNAME", "GENEID")
  
  tx2gene = tx2gene[!(is.na(tx2gene$GENEID) | tx2gene$GENEID == ""), ]
  
  tx2gene <- tx2gene[match(unique(tx2gene$TXNAME), tx2gene$TXNAME), ]
  
  # reading read count files to generate the gene-count table 
  txi <- tximport(input_files, type = "none", txOut = TRUE, txIdCol = "target_id", abundanceCol = "tpm",
                  countsCol = "est_counts", lengthCol = "length", importer = read_tsv)
  
  # Summarize the transcript abundance file
  txi.sum <- summarizeToGene(txi, tx2gene)
  
  return(txi.sum)
}