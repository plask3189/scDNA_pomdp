make_sce<- function(sample_file){
  
  file_name <- basename(sample_file) 
  sample_name <- sub("\\..*", "", file_name)
  data_kp = paste0("data_kp")
  
  if (!dir.exists(data_kp)) { 
    dir.create(data_kp)
  }
  save_sample_data_path = paste0(data_kp,"/",sample_name)
  if (!dir.exists(save_sample_data_path)) {
    dir.create(save_sample_data_path)
  }
  suppressWarnings(
  variant_output<-variant_ID(file=sample_file,
                             panel="MSK_RL", # "UCSC" can be used for other panels
                             GT_cutoff=0,  # mimimum percent of cells where a successful genotyping call was made
                             VAF_cutoff=0) # mimimum variant allele frequency 
  )
  genes_of_interest <- c("IDH2","NRAS","NPM1","TET2","FLT3","IDH1")
  
  variants_of_interest <- variant_output %>%
    dplyr::filter(Class != 'Intronic' & !is.na(Class)) %>%
    dplyr::filter(VAF > 0.01) %>%
    dplyr::filter(genotyping_rate > 85) %>%
    dplyr::filter(!is.na(CONSEQUENCE) & CONSEQUENCE != 'synonymous') %>%
    dplyr::filter(SYMBOL %in% genes_of_interest) %>%
    dplyr::filter(WT != 0) %>% # because if it is then that's fishy
    dplyr::arrange(desc(VAF)) %>%
    dplyr::slice(1:3)
  
  
  sce <- tapestri_h5_to_sce(file = sample_file, variant_set = variants_of_interest)
  sce@metadata[["sample_name"]]<-sample_name 
  sce <- enumerate_clones(sce)
  sce <- compute_clone_statistics(sce, skip_ploidy = FALSE)
  clono<-clonograph(sce, complete_only = TRUE, num_bars_to_keep = 5, title=sample_name); clono
  sce <- trajectory_analysis(sce, use_ADO = FALSE)
  final_vis<-visualize_tree(sce, variants_of_interest, remove_low_reward_edges = TRUE); final_vis
  return(sce)
}