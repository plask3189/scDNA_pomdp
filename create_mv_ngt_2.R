
make_mv_sce <- function(sample_file, variants_of_interest_for_dropout_examining) {
  cat(" For",sample_file, ":", "\n","Creating SCE on variants_of_interest_for_dropout_examining ", "\n")
  more_variants_sce <- tapestri_h5_to_sce(file = sample_file, variant_set = variants_of_interest_for_dropout_examining, GQ_cutoff=50)
  rownames(more_variants_sce)<- make.unique(rownames(more_variants_sce))
  more_variants_sce <- enumerate_clones(more_variants_sce)
  more_variants_sce <- suppressMessages(suppressWarnings(compute_clone_statistics(more_variants_sce, skip_ploidy = FALSE)))
  
  cat("Created SCE on variants_of_interest_for_dropout_examining", "\n" )
  return(more_variants_sce)
}

# add clone_codes from sce to more_variants_ngt based on barcode. 
add_clone_codes_to_mv_ngt <- function(more_variants_ngt, sce) { 
  if ("Clone" %in% names(more_variants_ngt)){
    more_variants_ngt<- more_variants_ngt %>% dplyr::select(-Clone)# remove more_variants_sce NGT's codes.
  }
  if ("Group" %in% names(more_variants_ngt)){
    more_variants_ngt<- more_variants_ngt%>% dplyr::select(-Group)
  } 
  if(!"Cell" %in% names(more_variants_ngt)){
    more_variants_ngt$Cell<- rownames(more_variants_ngt)
  }
  cellbarcodes_with_corresponding_code <- as.data.frame(metadata(sce)[["NGT"]][, c("Cell", "Clone")])
  #Left join to add the 'Clone' column to 'more_variants_ngt_df'
  more_variants_ngt_and_codes <- merge(more_variants_ngt, cellbarcodes_with_corresponding_code[, c("Cell", "Clone")], by = "Cell", all.x = TRUE)
  # Reorder the columns to make "Clone" the third column
  more_variants_ngt_and_codes <- more_variants_ngt_and_codes[, c(1, ncol(more_variants_ngt_and_codes), 2:(ncol(more_variants_ngt_and_codes) - 1))]
  print("Added clone codes from tree to the NGT with more variants.")
  #print(more_variants_ngt_and_codes[1:5,1:5])
  return(more_variants_ngt_and_codes)
}



select_variants_2<- function(variant_output, original_variants_of_interest, just_on_amplicon_param){
  
  if (just_on_amplicon_param == TRUE){
    amplicon_length<- 200# idk
    cat(red("- Including variants within 200bp of a `variants_of_interest` & on same chr:\n", paste(original_variants_of_interest$SYMBOL, collapse = ", "), "\n"))
    
    variants_of_interest_for_dropout_examining <- variant_output %>%
      filter(map_lgl(start, ~ any(abs(.x - original_variants_of_interest$start) <= amplicon_length)))
  } else{
    genes_of_interest <- c("BRAF", "TET2", "ASXL1", "NRAS", "KRAS", "RUNX1", "DNMT3A", "FLT3", "NPM1","SRSF2", "IDH1", "PTPN11", "TP53", "BCOR", "EZH2", "GATA2", "KIT", "PHF6", "STAG2", "U2AF1", "ZRSR2")
    cat("- Filtering for genes of interest: ", genes_of_interest,"\n")
    variants_of_interest_for_dropout_examining <- variant_output %>% # try removing GOIs
      dplyr::filter(SYMBOL %in% genes_of_interest)
  }
  
  cat("- Filtering variants with VAF > 0.01\n")
  variants_of_interest_for_dropout_examining <- variants_of_interest_for_dropout_examining %>%
    dplyr::filter(VAF >0.01)
  
  cat("- Filtering variants with genotyping rate > 80\n")
  variants_of_interest_for_dropout_examining <- variants_of_interest_for_dropout_examining %>%
    dplyr::filter(genotyping_rate > 70)
  

    cat(red("- Excluding variants with id already present in original sce's `variants_of_interest`:", paste(original_variants_of_interest$AA_change, collapse = ", ")), "\n")
    variants_of_interest_for_dropout_examining <- variants_of_interest_for_dropout_examining %>%
      dplyr::filter(!(id %in% original_variants_of_interest$id))#%>%
    #dplyr::filter((WT != 0))
    cat("- Ensuring uniqueness by keeping distinct `final_annot` entries\n")
    variants_of_interest_for_dropout_examining <- variants_of_interest_for_dropout_examining %>%
      dplyr::group_by(final_annot) %>%
      dplyr::slice_max(VAF, n = 1) %>%  # Keep only the row with the highest VAF per AA_change
      dplyr::ungroup() %>%
      dplyr::group_by(AA_change) %>%
      dplyr::slice_max(VAF, n = 1) %>%  # Keep only the row with the highest VAF per AA_change
      dplyr::ungroup()%>% 
      dplyr::arrange(desc(VAF))%>% 
      dplyr::filter(!is.na(AA_change)) %>% 
      dplyr::slice(1:800,)
 
  total_numvars<- nrow(variants_of_interest_for_dropout_examining)
  cat("Final number of variants after filtering:", nrow(variants_of_interest_for_dropout_examining), "\n")
  cat("------------------------------------------", "\n")
  return(variants_of_interest_for_dropout_examining)
}



update_sce<- function(more_variants_ngt, mv_sce){
  # Ensure that more_variants_ngt has row and column names
  if (is.null(rownames(more_variants_ngt)) || is.null(colnames(more_variants_ngt))) {
    stop("The more_variants_ngt matrix must have row and column names.")
  }
  genes_to_keep <- colnames(more_variants_ngt)
  genes_to_keep
  cells_to_keep <- (more_variants_ngt$Cell)
  cells_to_keep
  # Check if these genes and cells exist in mv_sce
  missing_genes <- setdiff(genes_to_keep, rownames(mv_sce))
  missing_cells <- setdiff(cells_to_keep, colnames(mv_sce))
  if (length(missing_genes) > 0) {
    warning("The following genes are not present in mv_sce and will be ignored: ", paste(missing_genes, collapse = ", "))
    genes_to_keep <- setdiff(genes_to_keep, missing_genes)
  }
  if (length(missing_cells) > 0) {
    warning("The following cells are not present in mv_sce and will be ignored: ", paste(missing_cells, collapse = ", "))
    cells_to_keep <- setdiff(cells_to_keep, missing_cells)
  }
  # Subset mv_sce to keep only the specified genes and cells
  mv_sce_filtered <- mv_sce[genes_to_keep, cells_to_keep]
  # Ensure that 'Cell' and 'Clone' columns exist in 'more_variants_ngt'
  if(!all(c("Cell", "Clone") %in% colnames(more_variants_ngt))) {
    stop("The 'more_variants_ngt' data frame must contain 'Cell' and 'Clone' columns.")
  }
  clone_mapping <- setNames(more_variants_ngt$Clone, more_variants_ngt$Cell)
  # chech 'clone_code' column exists in colData; if not, create it
  if(!"clone_code" %in% colnames(colData(mv_sce))) {
    colData(mv_sce)$clone_code <- NA
  }
  matching_cells <- intersect(rownames(colData(mv_sce)), names(clone_mapping))  # Update 'clone_code' in colData by matching row names with 'Cell' identifiers
  colData(mv_sce_filtered)$clone_code[matching_cells] <- clone_mapping[matching_cells]
  dim(mv_sce_filtered)
  return(mv_sce_filtered)
}

create_mv_ngt_2<- function(variant_output, variants_of_interest, sce, final_vis, just_on_amplicon_param = TRUE){
  original_variants_of_interest<- variants_of_interest
  variants_of_interest_for_dropout_examining<- select_variants_2(variant_output, original_variants_of_interest, just_on_amplicon=just_on_amplicon_param)
  mv_sce<- make_mv_sce(sample_file, variants_of_interest_for_dropout_examining)
  
  more_variants_ngt<- as.data.frame(metadata(mv_sce)[["NGT_with_missing"]])
  cat("Created NGT on variants_of_interest_for_dropout_examining.", "\n" )
  
  more_variants_ngt<-add_clone_codes_to_mv_ngt(more_variants_ngt, sce) # matches the barcodes from more_variants_ngt and sce, copies sce's Clone to more_variants_ngt 
  more_variants_ngt<-remove_cells_that_are_not_in_tree(final_vis, more_variants_ngt)
  
  more_variants_ngt <- more_variants_ngt %>%
    mutate(across(everything(), ~ ifelse(. == 3, 0, .)))
  
  mv_sce<- update_sce(more_variants_ngt, mv_sce)
  
  metadata(mv_sce)[["NGT"]]<- more_variants_ngt
  cat("Finished processing more variants NGT for sample:", sce@metadata[["sample_name"]], "\n")
  
  
  return(list(more_variants_ngt= more_variants_ngt, mv_sce=mv_sce ))
}

