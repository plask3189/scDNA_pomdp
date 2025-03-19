library(crayon)
library(ggridges)
library(SingleCellExperiment)
source("./R/pwm_verification/distribution_comparison.R")
source("./R/pwm_verification/visualize_pwm_results.R")
get_lineage <- function(clone_id, ancestry_table) { # Create a lineage path for each clone
  path <- clone_id
  if (!(clone_id %in% ancestry_table$Clone)) {
    return(clone_id)  # Return itself if missing
  }
  ancestor <- ancestry_table %>%
    filter(Clone == clone_id) %>%
    pull(Ancestor)
  while (!is.na(ancestor) && ancestor != "0_0_0") {
    if (!(ancestor %in% ancestry_table$Clone)) {
      break  # Stop if the ancestor is missing
    }
    path <- c(ancestor, path)
    ancestor <- ancestry_table %>%
      filter(Clone == ancestor) %>%
      pull(Ancestor)
  }
  return(paste(rev(path), collapse = " → "))
}


# ----------------------------------------------------------------------------------------------------
# --------Computes likelihood of daughter cell zygosities given the ancestors probability matrix. ----
# ----------------------------------------------------------------------------------------------------
compute_lineage_score <- function(cell_data, cell_lineages, pwm_log_prob = TRUE) { #  to compute lineage probability score
  print(cell_data)
  cat("\n Looking at a cell. This cell has info like ancestry, and the mutation(s) that made it be part of a clone \n")
  this_cell_clone <- unique(cell_data$Clone)
  cat("This Cell's Clone: ", this_cell_clone, "\n")
  variants_by_clone <- cell_data %>%
    group_by(Clone, Variant) %>%
    summarize(Value = sum(Value, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = Clone, values_from = Value, values_fill = list(Value = 0))
  variants_by_clone
  score_list <- c()
  ancestor_matrix_list <- list()
  
  # for each ancestor: like for each "0_0_0" "0_0_1" "1_0_1" "1_1_1" 
  # get all cells with this ancestor clone_code
  ancestor_clones <- unname(cell_data$Lineage_Path[[1]])
  ancestor_cells<- cell_lineages %>% filter(Clone %in% ancestor_clones)
  cat(blue("Our cell has", length(unique(ancestor_cells$Cell)) , "ancestor cells! "))
  # I want to make a matrix with Cells as rows and Variants as columns with corresponding Value. 
  ancestor_matrix <- ancestor_cells %>% # like ngt matrix. cells as rows. variants as cols w/ zygosity data
    dplyr::select(Cell, Variant, Value) %>% 
    pivot_wider(names_from = Variant, values_from = Value, values_fill = list(Value = 0)) %>%  # Reshape into matrix format
    column_to_rownames(var = "Cell") %>%  #  Cell names row names
    as.matrix()
  ancestral_zyg_count_matrix <- rbind( # count how many 0s, how many 1s, and how many 2s in the column (variant)
    Z0 = colSums(ancestor_matrix == 0, na.rm = TRUE),
    Z1 = colSums(ancestor_matrix == 1, na.rm = TRUE),
    Z2 = colSums(ancestor_matrix == 2, na.rm = TRUE)
  )
  # --------- Format the count matrix (num ancestor cells with wt het or hom for each variant)----------------
  ancestral_zyg_count_matrix<- as.data.frame(ancestral_zyg_count_matrix)
  colnames(ancestral_zyg_count_matrix) <- colnames(ancestor_matrix)
  rownames(ancestral_zyg_count_matrix) <- c("Z0", "Z1", "Z2")
  cat(blue("Number of ancestor cells with each wt(Z0), het(Z1), or hom(Z2) variant: \n" ))
  print(ancestral_zyg_count_matrix)
  # ------------------------------------------------------------------------
  total_num_cells<- nrow(ancestor_matrix)
  ancestral_zyg_probability_matrix<- ancestral_zyg_count_matrix /total_num_cells
  cat(blue("Percent of ancestor cells with each wt(Z0), het(Z1), or hom(Z2) variant:", "\n" ))
  print(ancestral_zyg_probability_matrix)
  our_cell_variant_zygosities <-cell_data%>% # our current cell
    dplyr::select(Variant, Value) %>%  # Keep only relevant columns
    pivot_wider(names_from = Variant, values_from = Value, values_fill = list(Value = 0))
  # for pwm_log_prob = TRUE: if PWM score is closer to 0, the daughter_cell follows the ancestral probabilities well
  pwm_score <- compute_pwm_score(ancestral_zyg_probability_matrix, our_cell_variant_zygosities, log_prob = pwm_log_prob) 
  cat(green("PWM Score: ", pwm_score, "\n")) #  (this cell's variant zygosities*prob of ancestor zygosity : 
  cat("--------------------------------- \n")
  return(pwm_score)
}


determine_clone_marker_zyg<- function(sce){
  cat(blue("Adding more variants to the sce... \n"))
  sce_ngt_long<- main_sce_ngt_prep(sce)
  # get new variants: 
  more_variant_data <- create_mv_ngt_2(variant_output, variants_of_interest, sce, final_vis, just_on_amplicon_param = TRUE) # final vis so can get codes in tree
  rownames(more_variant_data$more_variants_ngt) <- more_variant_data$more_variants_ngt$Cell
  more_variants_ngt<- more_variant_data$more_variants_ngt 
  #-----------------------------
  clone_marker_zygosity <- more_variants_ngt %>%
    pivot_longer(cols = -c(Clone, Cell), names_to = "Variant", values_to = "Value")
  cat(blue("Marker zygosities: \n"))
  print(head(clone_marker_zygosity))
  return(clone_marker_zygosity)
}

compute_kruskal<- function(clone_marker_zygosity){
  cat(blue("Running kruskal test...\n"))
  # do mutation distributions for different clones significantly vary for each variant? returns list of variants where these variants differ significantly between codes. 
  kruskal_results <- clone_marker_zygosity %>%
    group_by(Variant) %>%
    summarize(kruskal_pvalue = kruskal.test(Value ~ Clone)$p.value) %>%
    arrange(kruskal_pvalue)
  # gets the names of significantly different variant distribution
  kruskal_variants<- kruskal_analysis(kruskal_results, clone_marker_zygosity, select_significant_variants =TRUE)
  top_sig_kruskal<- (unique(kruskal_variants$Variant))#[1:15]
  cat(green("Identified top markers using the kruskal test: \n"))
  print((top_sig_kruskal))
  return(top_sig_kruskal)
}


make_ancestry_table<- function(adj_linklist_formatted){
  ancestry_table <- adj_linklist_formatted %>% dplyr::select(c("current_state", "next_state", "reward")) %>% dplyr::rename(Clone = current_state)
  cat(blue("Delineated the ancestry using the adj list \n"))
  return(ancestry_table)
}

get_lineage_path <- function(clone, ancestry_df) {
  # Convert the data frame into a named vector for quick lookup
  ancestry_dict <- setNames(ancestry_df$Clone, ancestry_df$next_state)
  # Initialize path with the given clone
  path <- c(clone)
  # Traverse ancestry until we find a root clone (i.e., a clone not present in `next_state`)
  while (!is.na(ancestry_dict[clone]) && ancestry_dict[clone] != clone) {
    clone <- ancestry_dict[clone]  # Move to the parent
    path <- c(path, clone)  # Add to lineage path
  }
  return(rev(path))  # Reverse to get from root to input clone
}

identify_lineages_for_clones <- function(ancestry_table){
  ancestry_dict <- as.data.frame(ancestry_table) %>% # Convert ancestry_table into a lookup dictionary
    dplyr::select(Clone, next_state) #%>% deframe()
  clone_lineages <- ancestry_table %>%
    rowwise() %>%
    mutate(Lineage_Path = list(get_lineage_path(Clone, ancestry_dict))) %>%  # Generate lineage path
    ungroup() %>%
    mutate(Lineage_Path = sapply(Lineage_Path, paste, collapse = " -> "))  %>%  # Convert list to character string
    as.data.frame() %>%
    dplyr::select(c("Clone", "Lineage_Path"))
  cat(blue("Identified lineages for clones (as specified in the tree): \n"))
  print(clone_lineages)
  
  # choose the longest path to represent this clone's lineage. 
  clone_lineages <- clone_lineages %>%
    mutate(Path_Length = str_count(Lineage_Path, "->")) %>%  # Count the number of transitions
    group_by(Clone) %>%  # Group by Clone
    slice_max(Path_Length, with_ties = FALSE) %>%  # Select the row with the longest path per Clone
    dplyr::select(-Path_Length) 
  return(clone_lineages)
}

if_each_tree_clone_was_dropout_what_would_be_its_real_clone <- function(Clones) { # input is all the clones in the tree. 
  dropout_table <- data.frame(
    Clone = character(),
    Clones_that_it_could_actually_be_if_tree_clone_is_dropout = character(),
    stringsAsFactors = FALSE
  )
  for(clone_code in Clones){
    elements <- unlist(strsplit(clone_code, "_"))
    dropout_clone_elements <- elements # init the same
    all_dropout_clones<- list()
    for (i in seq_along(elements)){
      if (elements[i] == "1"){
        dropout_clone_elements[i] <-"2"
        formatted_dropout_clone <- paste(dropout_clone_elements, collapse = "_")
        all_dropout_clones<-c(all_dropout_clones, formatted_dropout_clone) # append to list of dropout clones. 
        dropout_clone_elements <- elements # re-initialize 
      }
      else{ if (elements[i] == "0"){
        dropout_clone_elements[i] <-"1"
        formatted_dropout_clone <- paste(dropout_clone_elements, collapse = "_")
        all_dropout_clones<-c(all_dropout_clones, formatted_dropout_clone) # append to list of dropout clones. 
        dropout_clone_elements <- elements # re-initialize 
      }} 
    }
    new_row <- data.frame(
      Clone= clone_code,
      Clones_that_it_could_actually_be_if_tree_clone_is_dropout = paste(all_dropout_clones, collapse = ", "),
      stringsAsFactors = FALSE
    )
    dropout_table <- rbind(dropout_table, new_row) 
  }
  cat(blue("Determined what each tree clone would actually be if the tree clone is dropout.   \n"))
  print(dropout_table)
  return(dropout_table)
}

kl_divergence <- function(p, q) {
  p <- ifelse(p == 0, 1e-10, p)  # Avoid log(0) errors
  q <- ifelse(q == 0, 1e-10, q)
  res<- sum(p * log(p / q))
  return(res)
}
js_divergence <- function(p, q) {
  m <- (p + q) / 2
  (kl_divergence(p, m) + kl_divergence(q, m)) / 2
}
cosine_similarity <- function(p, q) {
  sum(p * q) / (sqrt(sum(p^2)) * sqrt(sum(q^2)))
}

compute_kl <- function(current_clone_count_matrix, hypothesis_lineage_zyg_count_matrix, log_prob = TRUE) { #zyg_matrix<- ancestral_zyg_prosition_matrix
  #cat(blue("Our current Clone's zygosity probability matrix: \n"))
  current_clone_probability_matrix<- current_clone_count_matrix/colSums(current_clone_count_matrix)
  #print(current_clone_probability_matrix)
  # make the ancestor probability matrix:
  ancestor_probability_matrix<- hypothesis_lineage_zyg_count_matrix/colSums(hypothesis_lineage_zyg_count_matrix)
  #cat(blue("Potential lineage zygosity probability matrix: \n"))
  #print(ancestor_probability_matrix)
  # P is actual distribution. Q is reference or expected.
  kl_divergence <-kl_divergence(current_clone_probability_matrix, ancestor_probability_matrix)
  return(kl_divergence)
}

compute_cosine <- function(current_clone_count_matrix, hypothesis_lineage_zyg_count_matrix, log_prob = TRUE) { #zyg_matrix<- ancestral_zyg_prosition_matrix
  #cat(blue("Our current Clone's zygosity probability matrix: \n"))
  current_clone_probability_matrix<- current_clone_count_matrix/colSums(current_clone_count_matrix)
  #print(current_clone_probability_matrix)
  # make the ancestor probability matrix:
  ancestor_probability_matrix<- hypothesis_lineage_zyg_count_matrix/colSums(hypothesis_lineage_zyg_count_matrix)
  #cat(blue("Potential lineage zygosity probability matrix: \n"))
  #print(ancestor_probability_matrix)
  # P is actual distribution. Q is reference or expected.
  cosine_similarity <-cosine_similarity(current_clone_probability_matrix, ancestor_probability_matrix)
  return(cosine_similarity)
}

compute_js <- function(current_clone_count_matrix, hypothesis_lineage_zyg_count_matrix, log_prob = TRUE) { #zyg_matrix<- ancestral_zyg_prosition_matrix
  cat(blue("Our current Clone's zygosity probability matrix: \n"))
  current_clone_probability_matrix<- current_clone_count_matrix/colSums(current_clone_count_matrix)
  print(current_clone_probability_matrix)
  # make the ancestor probability matrix:
  ancestor_probability_matrix<- hypothesis_lineage_zyg_count_matrix/colSums(hypothesis_lineage_zyg_count_matrix)
  cat(blue("Potential lineage zygosity probability matrix: \n"))
  print(ancestor_probability_matrix)
  # P is actual distribution. Q is reference or expected.
  js_divergence <- js_divergence(current_clone_probability_matrix, ancestor_probability_matrix)
  return(js_divergence)
}

construct_count_matrix<- function(zygosity_matrix){
  count_matrix <- apply(zygosity_matrix, 2, function(col) {
    table(factor(col, levels = c(0, 1, 2)))  # Ensure all levels (0,1,2) exist
  })
  rownames(count_matrix) <- c("Z0", "Z1", "Z2")
  return(count_matrix)
}



identify_alternative_lineages<- function(clone, clone_info){
  #cat(magenta("---------------------- Analyzing", (clone), "---------------------- \n"))
  the_current_cells_clone<- (clone)
  current_lineage <- clone_info %>% filter(Clone == the_current_cells_clone) %>% pull(Lineage_Path)
  #cat(blue("Current Tree Clone:", the_current_cells_clone, "\n"))
  #cat(blue("Current Tree Lineage:", current_lineage, "\n"))
  # get the potential actual clones for this cell. 
  the_current_cells_potential_actual_clones <- clone_info %>% filter(Clone == the_current_cells_clone) %>% pull(Clones_that_it_could_actually_be_if_tree_clone_is_dropout); cat(blue("Cells in this clone might actually belong to these clones (if a victim of dropout):", the_current_cells_potential_actual_clones, "\n"))
  individual_potential_clones <- str_split(the_current_cells_potential_actual_clones, ",\\s*")[[1]]
  # only examine the potential of the clone  if it is already in our tree ie in clone_info
  #individual_potential_clones <- individual_potential_clones[individual_potential_clones %in% clone_info$Clone]
  individual_potential_clones <- intersect(individual_potential_clones, clone_info$Clone)
  # add the tree clone to also compute the PWM on the tree lineage
  individual_potential_clones<- append(individual_potential_clones, the_current_cells_clone)
  return(individual_potential_clones)
}

kl<- function(current_clone, zygosity_matrix, clone_info_with_gene, individual_potential_clones){
  cat(magenta("Analyzing Clone", current_clone, "\n"))
  gene_to_look_for_dropout_on <- as.data.frame(clone_info_with_gene) %>% filter(Clone == current_clone) %>% dplyr::select(Gene) %>% pull()
  current_clone_ancestry <- clone_info_with_gene %>%
    filter(Clone == current_clone) %>%   # Filter for the hypothesis clone
    mutate(Path_Length = str_count(Lineage_Path, "->")) %>%  # Count transitions
    arrange(desc(Path_Length)) %>%   # Sort by longest path
    dplyr::slice(1) %>%   # Select the longest one
    pull(Lineage_Path)
  # current_clone_ancestors <-  current_clone_ancestry%>% str_split(" -> ")%>% unlist() %>% unique() %>% sort()
  # current_clone_ancestors <- setdiff(current_clone_ancestors, current_clone)
  # current_clone_zygosities<- zygosity_matrix %>% dplyr::filter(Clone %in% current_clone_ancestors)%>% dplyr::select(-c("Lineage_Path", "Clone", "Cell"))
  current_clone_zygosities<- zygosity_matrix %>% dplyr::filter(Clone %in% current_clone)%>% dplyr::select(-c("Lineage_Path", "Clone", "Cell"))
  kls_scores_list <- list()
  js_scores_list <- list()
  cosine_scores_list <- list()
  potential_actual_lineages_list<- list()
  potential_new_clone_list<- list()
  compute_pwm_only_on_specific_gene_of_interest = TRUE
  # ACTUALLY COMPARING TO ALL OTHER CLONES NOW:
  individual_potential_clones<-  (unique(zygosity_matrix$Clone)) # TEMP REPLACE WITH ALL CLONES
  # ------- construct zygosity matrices for each hypothesis clone using all the cells in each of the hypothesis clone's ancestors.  ----------------------
  for(hypothesis_clone in individual_potential_clones){# examine lineage of each hypothesis clone. (so this is a hypothesis lineage)
    cat(blue("Is current clone",current_clone, "more similar to the lineage of", hypothesis_clone, "or to its current lineage", current_clone_ancestry,"\n" ))
    # now get the lineages for these potential actual clones. 
    potential_actual_lineage <- clone_info_with_gene %>%
      filter(Clone == hypothesis_clone) %>%   # Filter for the hypothesis clone
      mutate(Path_Length = str_count(Lineage_Path, "->")) %>%  # Count transitions
      arrange(desc(Path_Length)) %>%   # Sort by longest path
      dplyr::slice(1) %>%   # Select the longest one
      pull(Lineage_Path)  # Extract the lineage path
    
    cat(blue("A potential lineage: ", paste(potential_actual_lineage, collapse = ", "), "\n"))
    potential_actual_lineages_list<- append(potential_actual_lineages_list, potential_actual_lineage)
    # get the ancestor clones of the hypothesis_clone
    ancestor_clones_of_our_hypothesis_clone <- potential_actual_lineage%>% str_split(" -> ")%>% unlist() %>% unique() %>% sort()
    # remove our clone to examine 
    ancestor_clones_of_our_hypothesis_clone <- setdiff(ancestor_clones_of_our_hypothesis_clone, current_clone)
    # all the cells in the hypothesis lineage
    hypothesis_lineage_zyg_matrix<- zygosity_matrix %>% filter(Clone %in% c(ancestor_clones_of_our_hypothesis_clone)) %>% dplyr::select(-c("Clone", "Cell", "Lineage_Path"))
    if (compute_pwm_only_on_specific_gene_of_interest== TRUE){
      # filter the zygosity matrix with all the variants of interest to only include ones on the gene of interest.
      hypothesis_lineage_zyg_matrix <- hypothesis_lineage_zyg_matrix %>% dplyr::select(starts_with(gene_to_look_for_dropout_on))
      current_clone_zygosities <-current_clone_zygosities %>% dplyr::select(starts_with(gene_to_look_for_dropout_on))
    }
    
    hypothesis_lineage_zyg_count_matrix<- construct_count_matrix(hypothesis_lineage_zyg_matrix)
    current_clone_zyg_count_matrix<- construct_count_matrix(current_clone_zygosities)
    kl_score = compute_kl(current_clone_zyg_count_matrix,hypothesis_lineage_zyg_count_matrix , log_prob = TRUE)
    kls_scores_list<- append(kls_scores_list, kl_score)
    cosine_score = compute_cosine(current_clone_zyg_count_matrix,hypothesis_lineage_zyg_count_matrix , log_prob = TRUE)
    cosine_scores_list<- append(cosine_scores_list, cosine_score)
    js_score = compute_js(current_clone_zyg_count_matrix,hypothesis_lineage_zyg_count_matrix , log_prob = TRUE)
    js_scores_list<- append(js_scores_list, js_score)
    potential_new_clone_list<- append(potential_new_clone_list, hypothesis_clone)
  }
  # make a df of the potential_actual_lineage and the resulting pwm score. 
  this_clones_pwm_scores_for_potential_lineages <- data.frame(
    #Lineage_Path = unlist(potential_actual_lineages_list), 
    Potential_New_Lineage =  unlist(potential_new_clone_list), 
    KL_Divergence = unlist(kls_scores_list),
    JS = unlist(js_scores_list),
    Cosine = unlist(cosine_scores_list)
  )
}

add_genes_to_clone_info <- function(clone_info, adj_linklist_formatted){
  df_with_clone_gene_info<- adj_linklist_formatted
  df_with_clone_gene_info$Clone <-  df_with_clone_gene_info$next_state
  df_with_clone_gene_info$label_prefix <- sub("\\..*", "", df_with_clone_gene_info$action_type)
  df_with_clone_gene_info<- df_with_clone_gene_info %>% dplyr::select(c("Clone", "label_prefix")) %>% dplyr::rename(Gene = label_prefix) %>% dplyr::filter(!Gene == "none") %>% distinct()
  clone_info_with_gene<- left_join(clone_info, df_with_clone_gene_info, by= "Clone")
  return(clone_info_with_gene)
}





find_lineage_specific_variants<- function(zygosity_matrix){
  lineage_comparisons <- list()
  lineages<- unique(zygosity_matrix$Lineage_Path)
  variants<-zygosity_matrix %>% dplyr::select(-c("Cell", "Clone", "Lineage_Path")) %>% colnames()
  
  lineage_summary <- data.frame(matrix(ncol = length(variants), nrow = 3)) 
  colnames(lineage_summary) <- variants  
  rownames(lineage_summary) <- c("Z0", "Z1", "Z2")
  
  for (lin in lineages){
    clones_in_this_lineage<- lin%>% str_split(" -> ")%>% unlist()
    lineage_zygosities<- zygosity_matrix %>% dplyr::filter(Clone %in% clones_in_this_lineage) %>% dplyr::select(-c("Cell", "Clone", "Lineage_Path"))
    lineage_summary <- apply(lineage_zygosities, 2, function(col) {
      c(sum(col == 0, na.rm = TRUE), # Count 0s
        sum(col == 1, na.rm = TRUE), # Count 1s
        sum(col == 2, na.rm = TRUE)) # Count 2s
    }) %>% as.data.frame()
    rownames(lineage_summary) <- c("Z0", "Z1", "Z2")
    # Reorder columns based on the values in the "Z1" row (second row)
    lineage_summary <- lineage_summary[, order(-as.numeric(lineage_summary["Z1", ]))]
    lineage_comparisons[[lin]] <- lineage_summary
    
    
  }
  lineage_comparison_df <- bind_rows(lineage_comparisons, .id = "Lineage")
}
