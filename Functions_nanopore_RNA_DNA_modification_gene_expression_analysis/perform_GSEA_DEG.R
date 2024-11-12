# Function to perform GSEA analysis
perform_GSEA_DEG <- function(entrezID_list,
                             hg_gene_set_entrezID,
                             DEG_table,
                             store_dir){
  
  # storing the directory information in a temporary variable
  temp_dir = str_c(store_dir, "/")
  
  # Get the item names of the list to generate folder with the name
  temp_item_names = names(entrezID_list)
  
  for (i in 1:length(temp_item_names)){
    
    # Let's create a directory to store the GO annotation files
    if (!dir.exists(str_c(temp_dir, temp_item_names[i]))){
      dir.create(str_c(temp_dir, temp_item_names[i]), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    temp_grp_dir <- str_c(temp_dir, temp_item_names[i], "/")
    
    gene_list <- entrezID_list[[temp_item_names[i]]]
    
    temp_conditions <- names(gene_list)
    
    for (j in 1:length(temp_conditions)) {
      
      # Let's create a directory to store the GO annotation files
      if (!dir.exists(str_c(temp_grp_dir, temp_conditions[j]))){
        dir.create(str_c(temp_grp_dir, temp_conditions[j]), showWarnings = TRUE, recursive = FALSE, mode = "0777")
      }
      
      # Assigning the directory path to a local variable
      temp_go_dir = str_c(temp_grp_dir, temp_conditions[j], "/")
      
      # Let's get the gene ids from the list
      temp_marker_genes = gene_list[[j]]
      
      # Store transcript, gene id table for the condition
      temp_match_ind = sort(match(temp_marker_genes, DEG_table$entrezgene_id))
      temp_match_df = DEG_table[temp_match_ind, ]
      
      # All GO terms (significant, non-significant)
      write.csv(temp_match_df, file = str_c(temp_go_dir, temp_conditions[j], "_condition_gene_transcript_ids.csv"))
      
      # performing the over representation analysis (ORA) for Gene Ontology class - Biological processes
      GeneSet_ORA_BP <- enrichGO(gene = temp_marker_genes,
                                 universe = hg_gene_set_entrezID,
                                 OrgDb = org.Hs.eg.db,
                                 keyType = "ENTREZID", 
                                 ont = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05,
                                 readable = TRUE,
                                 pool = FALSE)
      
      tryCatch({
        # Lets save the GO annotation results for BP class
        GeneSet_all_BP = GeneSet_ORA_BP@result
        
        # All GO terms (significant, non-significant)
        write.csv(GeneSet_all_BP, file = str_c(temp_go_dir, temp_conditions[j], "_condition_genes_all_biological_processes.csv"))
        
      }, 
      error = function(e){str_c("No biological processe annotation was found for ", temp_conditions[j], " genes.")}
      )
      
      # Figure 1 :: Generate dotplot for all GO terms - BP
      if (dim(GeneSet_ORA_BP)[1] < 1){
        print("No module found")
      } 
      else {
        p1 <- dotplot(GeneSet_ORA_BP, showCategory = 10) + 
          labs(colour = "p.adjust") + 
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(color = "black"),
                axis.title.x = element_text(size = 22, face = "bold", color = "black"), 
                axis.ticks.length = unit(.30, "cm"), 
                axis.text.x = element_text(size = 22, color = "black", face = "bold"),
                axis.text.y = element_text(size = 22, color = "black", face = "bold"),
                legend.key = element_rect(size = 22),
                legend.text = element_text(size = 22),
                legend.title = element_text(size = 22, face = "bold"),
                legend.spacing = unit(2.0, 'cm'),
                legend.key.size = unit(3,"line"),
                legend.position = "none")
        
        ggsave(filename = str_c(temp_go_dir, "GO_bp_",temp_conditions[j], "_condition_genes.png"), plot = p1, width = 14, height = 14, dpi = 300)
      }
      
      
      # Figure 2 :: Generate network of GO terms using all GO terms - BP
      if (dim(GeneSet_ORA_BP)[1] > 0){
        # You can also create an enrichment map that connects GO terms with edges between overlapping gene sets. This makes it easier to identify functional modules.
        GeneSet_ORA_BP <- pairwise_termsim(GeneSet_ORA_BP, method = "JC")
        print(dim(GeneSet_ORA_BP@termsim))
        
        if (dim(GeneSet_ORA_BP@termsim)[1] == 1){
          print("No module found")
        } 
        else {
          tryCatch({
            p2 <- emapplot(GeneSet_ORA_BP, color = "qvalue")  + theme_bw() + 
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    axis.title = element_blank(), 
                    text = element_text(size = 24, face = "bold"),
                    axis.line = element_blank(),
                    axis.ticks.length = unit(0, "cm"), 
                    axis.text = element_blank(),
                    legend.key = element_rect(size = 22),
                    legend.spacing = unit(2.0, 'cm'),
                    legend.key.size = unit(3,"line"),
                    legend.position = "none")
            ggsave(filename = str_c(temp_go_dir, "GO_bp_network",temp_conditions[j], "_condition_genes.png"), plot = p2, width = 14, height = 14, dpi = 300)
          }, 
          error = function(e) {
            message <- paste("Error in plotting network for GO terms - BP:", conditionMessage(e))
            cat(message, file = str_c(temp_go_dir, "error_log.txt"), append = TRUE)
            print(message)
          }
          )
        }
      }
      
      else {
        print("No module found")
      }
    }
  }
}
