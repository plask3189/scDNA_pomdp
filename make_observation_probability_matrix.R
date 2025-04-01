

get_variant_positions_in_clone<- function(transition_matrices){
  variant_clone_pos <- data.frame(variant = character(), clone_position = numeric(), example = character(), stringsAsFactors = FALSE)
  for(action_name in names(transition_matrices)){
    if (!action_name %in% c("listen","none")){
      #cat(action_name, "\n")
      action_data<- transition_matrices[[action_name]]
      longest_clone_length<-nchar(colnames(action_data)[length(colnames(action_data))-1])
      first_row<- action_data[1,]
      keeper<- names(first_row[first_row != 0])
      pos<- nchar(keeper)
      #pos<- (longest_clone_length -position)+1
      variant_clone_pos <- rbind(variant_clone_pos, data.frame(variant = action_name, clone_position = pos, example = keeper, stringsAsFactors = FALSE))
    }
  }
  print(variant_clone_pos)
  return(variant_clone_pos)
}
# this one actually uses quality but dont think it works  
make_observation_prob_matrix3<- function(sce, transition_matrices){
  
  #' update the none action's matrix. while the action of doing nothing leaves the state 
  #' exactly as it was (identity transition), it also offers no opportunity to gather 
  #' information about the state (observation is always “done”)
  observation_matrices[["none"]][,] <- 0
  observation_matrices[["none"]][,"absorption_state"] <-1
  
  variant_clone_pos<- get_variant_positions_in_clone(transition_matrices)
  states<- colnames(transition_matrices[[1]])
  # -------------------------------- Observation -------------------------------- 
  observation_matrices <-transition_matrices# make observation matrices structure same as transition matrices
  quality_metric<- "DP_med"
  # DP med: (total # reads * read length) / (total # bp in target region)
  quality<- sce@metadata[["Clones"]] %>%filter(Group == "Complete") %>% dplyr::select(c("Clone", "variants", quality_metric))
  quality$Clone <- gsub("_", "", quality$Clone ) %>% as.numeric()
  quality[[quality_metric]] <- quality[[quality_metric]] /100 ; cat(blue("Qualiy df: \n")); head(quality)
  
 #every_combo<- expand.grid(unique(quality$Clone), unique(quality$variants)) %>% dplyr::rename(Clone =Var1, variants=Var2) %>% arrange(Clone)
  
  # the prob that the state is 1 and we hear 0 is x%. this is probability that a call is 1 but is mistakenly heard as 0.
  # the prob that the state is 0 and we hear 1 is 0% bc dropout cant happen to state 0. 
  clones_with_potential_dropout<- make_potential_dropout_clones(states)

  # apply the quality metrics to the listen action 
  listen_matrix<- observation_matrices[["listen"]] 
  probabilities <- clones_with_potential_dropout %>%
    mutate(prob = 1 / sapply(could_be_heard_as, length))
  probabilities<- rbind(probabilities, data.frame(og_mutation= "absorption_state", could_be_heard_as= "absorption_state", prob=1))
  # use variant_clone_pos to determine which quality metric to use
  for(i in seq_along(rownames(listen_matrix))){
    actual_state <- rownames(listen_matrix)[i]
    observation_states <- as.character(probabilities$could_be_heard_as[[i]])
    cat("actual state:", actual_state, "observation states:", observation_states, "\n")
    for(obs_state in observation_states){
      obs_state_index <- which(observation_states == obs_state)
      actual_chars <- strsplit(actual_state, "")[[1]]
      obs_chars <- strsplit(obs_state, "")[[1]]
      if(length(actual_chars) == length(obs_chars)){
        position_of_var <- which(actual_chars != obs_chars)
        print(position_of_var)
        if(length(position_of_var !=0)){
          gene_for_index<- variant_clone_pos%>% filter(clone_position == position_of_var) %>% pull(variant)
          print(gene_for_index)
          # get row data
          row_data<- listen_matrix[i,]
          quality_for_pos<- quality%>% filter(Clone == actual_state) %>% filter(variants == gene_for_index) %>% pull(DP_med)
          print(quality_for_pos) # put this quality at colbindex
          listen_matrix[i,obs_state_index] 
          if(length(quality_for_pos) != 0){
            listen_matrix[i,obs_state_index]<- quality_for_pos
          } else{ 
            listen_matrix[i,obs_state_index]<-probabilities$prob[i]
          }
        }
      } else {
        cat("Cannot compare states ",
            actual_state, "and", obs_state, "\n")
      }
    }
  }

  observation_matrices[["listen"]]  <- listen_matrix
  #  to normalize rows so they sum to 1 just in case
  normalize_rows <- function(mat) {
    row_sums <- rowSums(mat)
    row_sums[row_sums == 0] <- 1  # Avoid division by zero
    normalized_mat <- sweep(mat, 1, row_sums, FUN = "/")
    cat("Original row sums: ", row_sums, "\n")
    cat("Normalized row sums: ", rowSums(normalized_mat), "\n")
    cat("--------------------------\n")
    return(normalized_mat)
  }
  # Normalize all matrices in the list
  observation_matrices <- lapply(observation_matrices, normalize_rows)
  return(observation_matrices)
}



make_observation_prob_matrix4<- function(sce, transition_matrices){
  states<- colnames(transition_matrices[[1]])
  observation_matrices <-transition_matrices# make observation matrices structure same as transition matrices
  listen_matrix<- observation_matrices[["listen"]] 

  clones_with_potential_dropout<- make_potential_dropout_clones(states)
  probabilities <- clones_with_potential_dropout %>%
    mutate(prob = 1 / sapply(could_be_heard_as, length))
  probabilities<- rbind(probabilities, data.frame(og_mutation= "absorbing_state", could_be_heard_as= "absorbing_state", prob=1))
  
  states
  # For each row in probabilities, update the cells in listen_matrix.
  for (i in seq_len(nrow(probabilities))) {
    cat(blue("looking at", probabilities$og_mutation[i], "\n"))
    action_state <- as.character(probabilities$og_mutation[i]); action_state
    observation_states <- as.character(probabilities$could_be_heard_as[[i]]); observation_states
    observation_states_indices<- which(states == observation_states); observation_states_indices
    prob_val <- probabilities$prob[i]
    cat(" action_state:", action_state,  "\n")
    cat("observation_states",observation_states, "prob:", prob_val, "\n")

    for(os in observation_states){
      observation_state_i<- which(states == os)
      listen_matrix[i, observation_state_i] <- prob_val
      cat("os:", observation_state_i, "i:", i,"prob", prob_val, "\n")
      if( (observation_state_i == i)){
        listen_matrix[i, observation_state_i] <- prob_val*2
      }
      
    }
    print("----------------")
  }
  observation_matrices[["listen"]]  <- listen_matrix
  observation_matrices[["listen"]] 
  normalize_rows <- function(mat) {
    row_sums <- rowSums(mat)
    row_sums[row_sums == 0] <- 1  # Avoid division by zero
    return(mat / row_sums)
  }
  
  # Normalize all matrices in the list
  observation_matrices <- lapply(observation_matrices, normalize_rows)
  observation_matrices[["listen"]]
  return(observation_matrices)
}








#-------


make_observation_prob_matrix6<- function(sce, transition_matrices, data_for_transition_matrix){
  states<- colnames(transition_matrices[[1]])
  actions<- names(transition_matrices)
  observation_matrices <-transition_matrices# make observation matrices structure same as transition matrices
  
  
  # ----------- LISTEN -----------------
  listen_matrix<- observation_matrices[["listen"]] 
  clones_with_potential_dropout<- make_potential_dropout_clones(states)
  probabilities <- clones_with_potential_dropout %>%
    mutate(prob = 1 / sapply(could_be_heard_as, length))
  # For each row in probabilities, update the cells in listen_matrix.
  for (i in seq_len(nrow(probabilities))) {
    action_state <- as.character(probabilities$og_mutation[i])
    cat(blue("looking at", action_state, "\n"))
    observation_states <- as.character(probabilities$could_be_heard_as[[i]]); observation_states
    observation_states_indices <- which(states %in% observation_states); observation_states_indices
    prob_val <- probabilities$prob[i]
    cat(" action_state:", action_state,  "\n")
    cat("observation_states",observation_states, "prob:", prob_val, "\n")
    for(os in observation_states){
      observation_state_i<- which(states == os)
      listen_matrix[i, observation_state_i] <- prob_val
      cat("os:", observation_state_i, "i:", i,"prob", prob_val, "\n")
      if( (observation_state_i == i)){
        listen_matrix[i, observation_state_i] <- prob_val
      }
    }
    print("----------------")
  }
  observation_matrices[["listen"]]  <- listen_matrix
  observation_matrices[["listen"]] 
  
  # data_for_transition_matrix[sapply(data_for_transition_matrix, is.factor)] <- 
  #   lapply(data_for_transition_matrix[sapply(data_for_transition_matrix, is.factor)], as.character)
  # # ---------- other probs ---------
  # actions<- setdiff(actions, "listen")
  # for(action in actions){
  #   print(action)
  #   data_for_action<- data_for_transition_matrix %>% filter(action_type == action)
  #   for(i in seq_len(nrow(data_for_action))){
  #     row<-data_for_action[i,]
  #     print(row)
  #     observation_state<- row[["next_state"]]
  #     prob <- row[["Probability_transition"]]
  #     observation_matrices[[action]][,observation_state] <- prob
  #   }
  # }
  
  
  #--------------------------------
  normalize_rows <- function(mat) {
    row_sums <- rowSums(mat)
    row_sums[row_sums == 0] <- 1  # Avoid division by zero
    return(mat / row_sums)
  }
  
  # Normalize all matrices in the list
  observation_matrices <- lapply(observation_matrices, normalize_rows)
  observation_matrices[["listen"]]
  return(observation_matrices)
}
