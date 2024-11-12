# Function to subset fasta file
subset_fasta <- function(fasta_file, transcript_id){
  
  if (grepl(pattern = "cdna", names(fasta_file)[1])){
    # Changing the name to keep only the transcript IDs in "ENST....." suffix
    names(fasta_file) <- sub(" cdna.*", "", names(fasta_file))
  }
  
  # Remove version suffixes from Ensembl transcript IDs
  names(fasta_file) <- sub("\\..*$", "", names(fasta_file))
  
  # Find sequences with header that contains gene name (subsetting the fasta file)
  fasta_sub <- fasta_file[str_detect(names(fasta_file), str_c(transcript_id, collapse = "|"))]
  
  return(fasta_sub)
}

# Function to generate transcript gene table
generate_transcript_gene_table <- function(search_ids, 
                                           attributes_to_retrieve,
                                           search_filter){
  
  # Clean the query cache of biomaRt
  biomaRt::biomartCacheClear()
  
  # Generate the mart object for HS 
  hsmart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  # information retrieval from database
  attributes_to_retrieve = unique(c("ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id", "ensembl_transcript_id_version", "entrezgene_id", "hgnc_symbol", "external_gene_name", attributes_to_retrieve))
  
  temp_mart <- getBM(attributes = attributes_to_retrieve,
                     filters = search_filter,
                     values = search_ids, 
                     mart = hsmart)
  
  return(temp_mart)
}


# Function to generate fasta subset and generate gene transcript table
fastaSubset_and_geneTransript_table <- function(file_paths, config_file_path, fasta_path){
  # Load conditions and effect direction from config.yaml
  config <- yaml::read_yaml(config_file_path)
  
  # Get the search filter for biomart
  biomart_search_filter <- config$biomart_search_filter
  
  # Get the attributes to retreive from biomart ensembl database
  biomart_attribute_search <- config$biomart_attribute_search
  
  # Read the fasta file
  temp_fasta <- readDNAStringSet(fasta_path)
  
  if (grepl(pattern = "cdna", names(temp_fasta)[1])){
    # Changing the name to keep only the transcript IDs in "ENST....." suffix
    temp_symbol_retreival <- sub(" cdna.*", "", names(temp_fasta))
    # Remove version suffixes from Ensembl transcript IDs
    temp_symbol_retreival <- sub("\\..*$", "", temp_symbol_retreival)
  } else {
    temp_symbol_retreival <- names(temp_fasta)
    # Remove version suffixes from Ensembl transcript IDs
    temp_symbol_retreival <- sub("\\..*$", "", temp_symbol_retreival)
  }
  
  # Load the merged table
  temp_df = loadRData(file_paths)
  
  # List of transcripts
  transcripts_ids = as.character(unique(temp_df$ref_seq_id))
  
  # get the fasta subset
  temp_subset <- subset_fasta(fasta_file = temp_fasta, transcript_id = transcripts_ids)
  
  # generate gene transcript table
  temp_gene_table <- generate_transcript_gene_table(search_ids = temp_symbol_retreival, attributes_to_retrieve = biomart_attribute_search, search_filter = biomart_search_filter)
  
  # Add to a list and return it
  temp_list <- list(fasta_subset = temp_subset, gene_transcript_table = temp_gene_table)
  
  return(temp_list)
}