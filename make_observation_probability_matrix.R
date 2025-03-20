# use DP or AF or GQ in obs matrix? or just prob of dropout? like 22 is 0% likely to be victim of  dropout
# DP = DP = Reference_reads + Mutant reads.
# A cell could have 40 reads of reference aligned reads and 70 of mutant reads. then gq = 110

# Allele Frequency (AF) or Variant Allele Frequency (VAF) or GQ 
make_observation_prob_matrix<- function(sce, legal_actions){

  legal_actions2 <- lapply(legal_actions, function(x) {
    if (is.list(x) && "state_trans_matrix" %in% names(x)) {
      return(x$state_trans_matrix)  # Move matrix up one level
    } else {
      return(x)  # Leave as is if no nested matrix
    }
  })
  # get clone and data info
  mutation_names<-sub("_.*", "", names(legal_actions2)) ; mutation_names
  clone_states<- colnames(legal_actions2[[1]]); clone_states
  state_matrix <- matrix(0, nrow = length(clone_states), ncol = length(clone_states),
                         dimnames = list(clone_states, clone_states))
  # subset the quality data by mutation
  clone_quality_data<- sce@metadata[["Clones"]] %>% dplyr::select(-c("Group"))

  which_quality_metrix<- "GQ_med"
  qual_obs_matrix<- clone_quality_data %>% dplyr::select(c("Clone", "variants", which_quality_metrix)) %>% arrange(Clone)
  
  qual_variant_subsets_list<- list()
  for(variant in mutation_names){
    qual_variant_subset <- qual_obs_matrix %>% filter(variants == variant)
    qual_variant_subset<- qual_variant_subset %>%
      group_by(Clone) %>%
      summarise(qual = mean(.data[[which_quality_metrix]]) / 100, .groups = 'drop')
    qual_variant_subsets_list[[variant]] <- qual_variant_subset
  }
  
  qual_variant_matrices_list <- list()
  
  for (variant in names(qual_variant_subsets_list)) {
    # Extract the tibble for the current variant
    qual_variant_subset <- qual_variant_subsets_list[[variant]]
    
    # Create an empty matrix with state names as rows and columns
    state_matrix <- matrix(0, nrow = length(clone_states), ncol = length(clone_states),
                           dimnames = list(clone_states, clone_states))
    
    # Fill in the diagonal with DP_mean values (assuming self-transition storage)
    for (i in seq_along(clone_states)) {
      state <- clone_states[i]
      if (state %in% qual_variant_subset$Clone) {
        state_matrix[state, state] <- qual_variant_subset$qual[qual_variant_subset$Clone == state]
      }
    }
    qual_variant_matrices_list[[variant]] <- state_matrix
  }
  
  # Apply normalization to each matrix in qual_variant_matrices_list
  qual_variant_matrices_list_normal <- lapply(qual_variant_matrices_list, normalize_matrix_rows)
  
  # add the obs matrix like format of transition_probs
  main_observation_matrix<- legal_actions2
  for (mutation_action in names(main_observation_matrix)) {
    for (variant_name in names(qual_variant_matrices_list_normal)) {
      if (grepl(paste0("^", variant_name), mutation_action)) {  # Check if mutation_action starts with variant_name
        #cat("Mutation action:", mutation_action, "| Variant name:", variant_name, "\n") 
        # get the variant associated obs matrix and add it to the mutation action
        main_observation_matrix[[mutation_action]]<-qual_variant_matrices_list_normal[[variant_name]]
      }
    }
  }
  cat(blue("Made observation probability matrix using", which_quality_metrix, " which looks like:\n"))
  print(main_observation_matrix[1])
  return(main_observation_matrix)
}
