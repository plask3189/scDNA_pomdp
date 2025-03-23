make_observation_prob_matrix3<- function(sce, transition_matrices){
  states<- colnames(transition_matrices[[1]])
  # -------------------------------- Observation -------------------------------- 
  observation_matrices <-transition_matrices# make observation matrices structure same as transition matrices
  quality_metric<- "DP_med"
  # DP med: (total # reads * read length) / (total # bp in target region)
  
  quality<- sce@metadata[["Clones"]] %>%filter(Group == "Complete") %>% dplyr::select(c("Clone", "variants", quality_metric))
  quality$Clone <- gsub("_", "", quality$Clone ) %>% as.numeric()
  quality[[quality_metric]] <- quality[[quality_metric]] /100 ; cat(blue("Qualiy df: \n")); head(quality)
  
  # joint_quality_metrics<- quality %>%
  #   group_by(Clone) %>%
  #   mutate(join_qual = prod(GQ_med)) # probability of mut 1 and mut2 and mut 3 and mut n
  # apply the quality metrics to the listen action 
  listen_matrix<- observation_matrices[["listen"]] 
  # the prob that the state is 1 and we hear 0 is x%. this is probability that a call is 1 but is mistakenly heard as 0.
  # the prob that the state is 0 and we hear 1 is 0% bc dropout cant happen to state 0. 
  # dummy probs for now. 
  clones_with_potential_dropout<- make_potential_dropout_clones(states)
  probabilities <- clones_with_potential_dropout %>%
    mutate(prob = 1 / sapply(could_be_heard_as, length))
  
  # Assume that listen_matrix has row and column names corresponding to your states.
  # For each row in probabilities, update the cells in listen_matrix.
  for (i in seq_len(nrow(probabilities))) {
    # Get the current row's og_mutation and convert to character (matching row names)
    row_id <- as.character(probabilities$og_mutation[i])
    # Get the vector of potential heard states for that mutation.
    # Convert them to character so they match the column names.
    cols_to_update <- as.character(probabilities$could_be_heard_as[[i]])
    # Get the probability value for this row
    prob_val <- probabilities$prob[i]
    # Substitute the probability value into listen_matrix at the specified row and columns
    listen_matrix[row_id, cols_to_update] <- prob_val
  }
  observation_matrices[["listen"]]  <- listen_matrix
  # Function to normalize rows so they sum to 1
  normalize_rows <- function(mat) {
    row_sums <- rowSums(mat)
    row_sums[row_sums == 0] <- 1  # Avoid division by zero
    return(mat / row_sums)
  }
  
  # Normalize all matrices in the list
  observation_matrices <- lapply(observation_matrices, normalize_rows)
  return(observation_matrices)
}





make_observation_prob_matrix2<- function(sce, transition_matrices){
  states<- colnames(transition_matrices[[1]])
  # -------------------------------- Observation -------------------------------- 
  observation_matrices <-transition_matrices# make observation matrices structure same as transition matrices
  quality_metric<- "DP_med"
  quality<- sce@metadata[["Clones"]] %>%filter(Group == "Complete") %>% dplyr::select(c("Clone", "variants", quality_metric))
  quality$Clone <- gsub("_", "", quality$Clone ) %>% as.numeric()
  quality[[quality_metric]] <- quality[[quality_metric]] /100 ; cat(blue("Qualiy df: \n")); head(quality)
  
  # make a quality matrix for just valid transitions for now. 
  variants<- unique(quality$variants)
  for (variant in variants){
    quality_for_var<- quality %>% filter(variants == variant)
    for(state_index in 1:nrow(observation_matrices[[variant]])){
      clone_state<- rownames(observation_matrices[[variant]])[state_index]
      cat("Clone state:",clone_state , "\n")
      quality_score<- quality_for_var %>% filter(Clone== clone_state) %>% pull(quality_metric)
      if(length(quality_score)==0){
        quality_score<-1
      }
      cat("quality score:", quality_score, "\n")
      row_vals<- observation_matrices[[variant]][state_index,]
      observation_matrices[[variant]][ state_index,] <- row_vals*quality_score
    }
    cat("observation matrix for this variant:", variant, "\n")
    print(observation_matrices[[variant]])
  }
  # HANDLE DROPOUT HERE:
  # if the row clone is victim of dropout, it would look like (be observed as) as set of column clones
  # identify the column clones that our clone might look like if dropout. 
  clone_dropout_table<-identify_clones_that_this_clone_could_actually_be(states)
  
  # distribute leftover probabilities to the clones that it could actually be if dropout happened.
  for (variant in names(observation_matrices)){
    for(state_index in 1:nrow(observation_matrices[[variant]])){
      clone_state<- rownames(observation_matrices[[variant]])[state_index]
      sum_of_probs<- sum(observation_matrices[[variant]][state_index,])
      leftover_prob <- (1 - sum_of_probs)  
      clone_if_dropout<- clone_dropout_table %>% filter(Clone == clone_state) %>% pull(Neighbor)
      cat("Actual clone", clone_state, "but might incorrectly observe", clone_if_dropout, "with probability", leftover_prob*100,"%\n")
      colindex <- which(colnames(observation_matrices[[variant]]) == clone_if_dropout)
      observation_matrices[[variant]][ state_index,colindex] <-leftover_prob
    }
    cat("observation matrix for this variant:", variant, ":\n")
    print(format(observation_matrices[[variant]], digits = 4, nsmall = 4))
  }
  # Function to normalize rows so they sum to 1
  normalize_rows <- function(mat) {
    row_sums <- rowSums(mat)
    row_sums[row_sums == 0] <- 1  # Avoid division by zero
    return(mat / row_sums)
  }
  
  # Normalize all matrices in the list
  observation_matrices <- lapply(observation_matrices, normalize_rows)
  return(observation_matrices)
}









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




