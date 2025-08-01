
library(scDNA)
library(dplyr)
library(ggplot2)
library(magick)
library(crayon)
library(stringr)
library(HDF5Array)
library(pomdp)
source("./R/RL_dev/mdp_Q_learning_with_linklist_kp.R")
source("./R/RL_dev/attach_weights_kp.R")
source("./R/RL_dev/reformat_actions.R")
source("./R/RL_dev/make_observation_probability_matrix.R")
source("./R/RL_dev/make_rewards.R")
source("./R/new_dev/visualize_clonal_evolution_kp.R")
source("./R/clonograph.R")

sce<- make_sce(sample_file= "BRAF/A5330braf.dna+protein.h5") # normal run. Makes sce on 2 variants. 
#---------POMDP-----------
mutation_states<-length(unique(sce@metadata$Architecture$final_annot)) ; mutation_states
paste("Building MDP initially with all possible mutation combinations = ", mutation_states) # different mutations
adj_list<-BuildMDP(mutation_states,use_ADO= FALSE)
adj_list<-attach_weights_kp(sce,adj_list)
adj_list <-adj_list %>% dplyr::select(-c("observed_states"))
# format adjacency list
adj_list_formatted<- adj_list
adj_list_formatted$current_state <- as.character(adj_list_formatted$current_state)
# Find the maximum length of the numbers in the column
max_length <- max(nchar(adj_list_formatted$current_state))
# Format the numbers with leading zeros
adj_list_formatted$current_state <- sprintf(paste0("%0", max_length, "d"), as.integer(adj_list_formatted$current_state))
adj_list_formatted$next_state <- as.character(adj_list_formatted$next_state)
max_length <- max(nchar(adj_list_formatted$next_state))
adj_list_formatted$next_state <- sprintf(paste0("%0", max_length, "d"), as.integer(adj_list_formatted$next_state))

# a global key sorry
adj_list_formatted <- adj_list_formatted %>%
  mutate(digit_position_changed = mapply(identify_digit_change, current_state, next_state))
position_of_variants_in_clone<- as.data.frame(adj_list_formatted) %>% dplyr::select(c("action_type", "digit_position_changed")) %>%  filter(action_type != "none") %>% distinct()
assign("position_of_variants_in_clone", position_of_variants_in_clone, envir = .GlobalEnv)
reformatted_adj_list2<- reformatted_actions_adj_list(adj_list_formatted)
reformatted_adj_list2 <- reformatted_adj_list2 %>%
  mutate(
    current_state = add_underscore(current_state),
    next_state = add_underscore(next_state)
  )
# a global key sorry. need to do this again to deal with the _
reformatted_adj_list2 <- reformatted_adj_list2 %>%
  mutate(digit_position_changed = mapply(identify_digit_change, current_state, next_state))
position_of_variants_in_clone<- as.data.frame(reformatted_adj_list2) %>% dplyr::select(c("action_type", "digit_position_changed")) %>%  filter(action_type != "none") %>% distinct()
assign("position_of_variants_in_clone", position_of_variants_in_clone, envir = .GlobalEnv)

states<-unique(reformatted_adj_list2$current_state) ; cat(blue("States:", paste(states, collapse = ", ")))

actions<-make_actions(reformatted_adj_list2) 

length(actions)
legal_actions<- make_legal_action_matrix(actions, reformatted_adj_list2)
length(legal_actions)
# make and format the trasnition probs
transition_probs<- make_transition_matrix(legal_actions)
length(transition_probs)
for (action in names(transition_probs)) {
  if (!is.null(transition_probs[[action]]$state_prob_matrix)) {
    transition_probs[[action]] <- transition_probs[[action]]$state_prob_matrix
    dimnames(transition_probs[[action]]) <- list(states, states)
    transition_probs[[action]]<- as.matrix(transition_probs[[action]])
  }
}


main_observation_matrix<- make_observation_prob_matrix(sce, legal_actions)
transition_probs
#OVERWRITE OTHER OBS MATRIX TEMP.
# main_observation_matrix <- lapply(main_observation_matrix, function(mat) {
#   diag(nrow(mat))  # Create an identity matrix of the same size
# })


rewards<- make_reward_df2(states, actions, reformatted_adj_list2) #%>% distinct()
g <- transition_graph(sc, action = NULL, episode = NULL, 
                      epoch = NULL, state_col = NULL, simplify_transitions = TRUE)

plot(g, edge.label= NA,
                      edge.label.cex = 0.5,
                      vertex.label.cex = 0.5,
                      vertex.size = 30, edge.arrow.size = .3, margin = .5)




start <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)
horizon_val <- Inf
sc <- POMDP(
  name = "sc",
  discount = 0.3,# if 1,  values future rewards as much as immediate rewards
  states = c(states),
  actions = c(actions),
  horizon=horizon_val, # number of epochs. 
  observations = c(states),
  #start = start,
  transition_prob = transition_probs,
  observation_prob = main_observation_matrix,
  reward = rewards
)

# c("grid", "enum", "twopass", "witness", "incprune")
solution <- solve_POMDP(sc, method = "incprune", verbose = TRUE)





add_underscore <- function(state) {
  gsub("(\\d)(\\d)", "\\1_\\2", state) # Uses regex to separate digits with "_"
}

# helper for observation matrix construction
normalize_matrix_rows <- function(mat) {
  for (i in 1:nrow(mat)) {
    row_values <- mat[i, ]  # Extract row
    row_sum <- sum(row_values[row_values != 0])  # Sum of non-zero values
    vals_to_fill <- 1 - row_sum  # Remaining proportion
    num_zeros <- sum(row_values == 0)  # Count zeros
    if (num_zeros > 0) {# Only update zero values
      mat[i, row_values == 0] <- vals_to_fill / num_zeros
    }
  }
  return(mat)
}


# Function to identify digit position change
identify_digit_change <- function(current, next_state) {
  # Convert strings to character vectors
  split_current <- unlist(strsplit(current, ""))
  split_next <- unlist(strsplit(next_state, ""))
  # Find the position where the digits differ
  diff_pos <- which(split_current != split_next)
  # Return the position of change (if multiple, concatenate)
  if (length(diff_pos) == 0) {
    return(NA)  # No change
  } else {
    return(paste(diff_pos, collapse = ","))
  }
}



#check_if_action_is_legal_for_state("NRAS.G12R_1-->2", "10", "next")
check_if_action_is_legal_for_state<- function(mutation, clone_state, is_start_or_next = 'start'){
  mutations_starting_zygosity<-  sub(".*_(\\d+)-->.*", "\\1", mutation)
  mutations_resulting_zygosity<- sub(".*-->\\s*(\\d+).*", "\\1", mutation)
  cat(blue("Checking if",mutation, "is possible if Clone", clone_state, "is",is_start_or_next, "clone and the mutation follows.\n"))
  variant_name<- sub("_.*", "", mutation)
  index_of_variant_in_clone<- position_of_variants_in_clone%>% filter(action_type == variant_name) %>% pull(digit_position_changed) 
  state_zygosity <- substr(clone_state, index_of_variant_in_clone, index_of_variant_in_clone)
  if(is_start_or_next == 'start'){
    if(state_zygosity == mutations_starting_zygosity){
      cat(green(mutation, "is possible if Clone", clone_state, "is",is_start_or_next, "clone and the mutation follows.\n"))
      return(TRUE)
    } else {
      cat(red(mutation, "is not possible if Clone", clone_state, "is",is_start_or_next, "clone and the mutation follows.\n"))
      return(FALSE)}
  }
  if(is_start_or_next == 'next'){
    if(state_zygosity == mutations_resulting_zygosity){
      cat(green(mutation, "is possible if Clone", clone_state, "is",is_start_or_next, "clone and the mutation follows.\n"))
      
      return(TRUE)
    } else{ 
      cat(red(mutation, "is not possible if Clone", clone_state, "is",is_start_or_next, "clone and the mutation follows.\n"))
      return(FALSE)
      }
  }
}

make_actions<- function(adj_linklist_formatted){
  original_actions<- unique(adj_list_formatted$action_type)
  original_actions <- original_actions[original_actions != "none"]
  transitions <- list(c(0, 1, 2))  # Possible states
  # Generate all possible combinations of state transitions
  legal_actions <- expand.grid(state_variable1 = transitions[[1]], state_variable2 = transitions[[1]])
  # Remove invalid transitions (e.g., direct jump from 0 to 2)
  valid_transitions <- function(state1, state2) {
    return((state1 == 0 & state2 %in% c(0, 1)) |  # 0 → 0, 0 → 1
             (state1 == 1 & state2 %in% c(1, 2)) |  # 1 → 1, 1 → 2
             (state1 == 2 & state2 == 2))           # 2 → 2
  }
  legal_actions <- subset(legal_actions, valid_transitions(state_variable1, state_variable2))
  
  actions_list <- character()  
  # Loop through original actions and legal transitions
  for (action in original_actions) {
    for (i in seq_along(legal_actions$state_variable1)) {
      formatted_action <- paste0(action, "_", legal_actions$state_variable1[i], "-->", legal_actions$state_variable2[i])
      actions_list <- c(actions_list, formatted_action)
    }
  }
  print(actions_list)
  return(actions_list)
}

# what states could happen for an action
make_legal_action_matrix<- function(actions, reformatted_adj_list2){
  transition_matricies<- list()
  unique_states <- unique(reformatted_adj_list2$current_state)
  state_adj_matrix <- as.data.frame(matrix(0, nrow = length(unique_states), ncol = length(unique_states),
                             dimnames = list(unique_states, unique_states)))
  action_data<- parse_action(actions, position_of_variants_in_clone)
  
  # Initialize a named list where each action has its own list of valid indices
  state_transitions <- list()
  state_transition_matrices<- list()
  for (i in seq_len(nrow(action_data))) { 
    #reset adj matrix
    state_adj_matrix <- as.data.frame(matrix(0, nrow = length(unique_states), ncol = length(unique_states),
                               dimnames = list(unique_states, unique_states)))
    action_name <- action_data$action[i]  # Get the action name
    state_transition_matrices[[action_name]] <- list(state_trans_matrix = data.frame())
    index_of_variant_in_clone <- action_data$position_in_clone[i]
    # Find valid NEXT states in the columns
    for (next_state in colnames(state_adj_matrix)) { # these are the clones we potentially transition to for the action. 
      # get the clone's zygosity by looking at what num is at index_of_variant_in_clone
      next_state_zygosity <- substr(next_state, index_of_variant_in_clone, index_of_variant_in_clone)
      if (next_state_zygosity == action_data$next_zygosity[i]) {
        state_transitions[[action_name]]$valid_next_states <- append(
          state_transitions[[action_name]]$valid_next_states, next_state)
      }
    }
    # Find valid current states
    for (current_state in rownames(state_adj_matrix)) { 
      current_state_zygosity <- substr(current_state, index_of_variant_in_clone, index_of_variant_in_clone)
      if (current_state_zygosity == action_data$curent_zygosity[i]) {
        state_transitions[[action_name]]$valid_current_states <- append(
          state_transitions[[action_name]]$valid_current_states, current_state)
      }
    }
     state_adj_matrix[unlist(state_transitions[[action_name]]$valid_current_states), 
                     unlist(state_transitions[[action_name]]$valid_next_states)] <- 1
      cat(blue(action_name, "at clone index", index_of_variant_in_clone,"\n"))
      print(state_adj_matrix)
      state_transition_matrices[[action_name]]$state_trans_matrix <- rbind(
        state_transition_matrices[[action_name]]$state_trans_matrix, 
        state_adj_matrix)
  }
  return(state_transition_matrices)
}

# transitino matrix helper functions
normalize_rows <- function(matrix) {
  row_sums <- rowSums(matrix)  # Compute sum of each row
  row_sums[row_sums == 0] <- 1  # Avoid division by zero
  return(matrix / row_sums)  # Element-wise division
}
# ensure the probs sum to 1. 
make_self_loop<- function(normalized_action_matrix){
  for (i in 1:nrow(normalized_action_matrix)) {
    row_sum <- sum(normalized_action_matrix[i, ])
    if (row_sum < 1) {
      normalized_action_matrix[i, i] <- 1 - row_sum  # Add self-loop
    }
  }
  return(normalized_action_matrix)
}

make_transition_matrix<- function(state_transition_matrices){
  state_probabilities<-list()
  # get the probability transition from adj_list_formatted
  for (action_name in names(state_transition_matrices)) {
    state_probabilities[[action_name]] <- list(state_prob_matrix = data.frame())
    action_matrix <- state_transition_matrices[[action_name]][["state_trans_matrix"]]
    normalized_action_matrix <- normalize_rows(action_matrix) # make probability of all possible transitions
    # add self loop or absorbing state for non valid transitions. 
    # make absorbing state for illegal transitions.
    #transition_matrix<-normalized_action_matrix # make_absorbing_state(normalized_action_matrix)
    # for now self loop:
    transition_matrix<- make_self_loop(normalized_action_matrix)
    cat(blue("When", action_name, "Probabilities of transitions: \n"))
    print(transition_matrix)
    # add to list
    state_probabilities[[action_name]]$state_prob_matrix <- data.frame(
      transition_matrix, 
      row.names = rownames(action_matrix)  # Preserve row names
    )
  }
  return(state_probabilities)
}

parse_action<- function(actions, position_of_variants_in_clone){
  action_df_all <- data.frame(action_formatted = character(),
                              variant = character(),
                              position_in_clone = character(),
                              curent_zygosity = character(),
                              next_zygosity = character(),
                              stringsAsFactors = FALSE)
  for (action in actions){
    action_df<- data.frame(action_formatted = NA, variant= NA, position_in_clone=NA, curent_zygosity = NA, next_zygosity = NA)
    variant <-  sub("_.*", "", action) 
    position_in_clone<- position_of_variants_in_clone %>% dplyr::filter(action_type == variant) %>% pull(digit_position_changed)
    curent_zygosity <- sub(".*_(\\d+)-->.*", "\\1", action)
    next_zygosity <- sub(".*-->(\\d+)", "\\1", action)
    action_df <- data.frame(action, variant, position_in_clone, curent_zygosity, next_zygosity, stringsAsFactors = FALSE)
    action_df_all <- rbind(action_df_all, action_df)
  }
  cat(blue("Constructed a matrix for probability of state transitions for the mutations \n"))
  return(action_df_all)
}

