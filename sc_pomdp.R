
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
source("./R/RL_dev/make_sce.R")
source("./R/new_dev/visualize_clonal_evolution_kp.R")
source("./R/clonograph.R")
sample_file= "BRAF/A5330braf.dna+protein.h5"
sce<- make_sce(sample_file) # normal run. Makes sce on 2 variants. 
#---------POMDP-----------
mutation_states<-length(unique(sce@metadata$Architecture$final_annot)) ; mutation_states
paste("Building MDP initially with all possible mutation combinations = ", mutation_states) # different mutations
adj_list<-BuildMDP(mutation_states,use_ADO= FALSE)
adj_list<-attach_weights_kp(sce,adj_list)

states<-unique(levels(adj_list$current_state))
absorbing_state <- "absorbing_state"
states <- c(states, absorbing_state)  ; cat(blue("States:", paste(states, collapse = ", ")))
actions<-unique(adj_list$action_type) ; actions

data_for_transition_matrix<- adj_list %>% dplyr::select(c("current_state", "action_type", "next_state", "Probability_transition"))
transition_matrices <- list()
# Generate a transition matrix for each action type
for (action in actions) {
  transition_matrix <- matrix(0, nrow = length(states), ncol = length(states),
                              dimnames = list(states, states))
  if(action == "none"){ #make identity matrix bc if no change then current_state == next_state.
    diag(transition_matrix) <- 1
  } else {
    action_data <- data_for_transition_matrix %>% 
      filter(action_type == action) %>%
      dplyr::select(current_state, next_state, Probability_transition)
    for (i in seq_len(nrow(action_data))) {
      from <- as.character(action_data$current_state[i])
      to <- as.character(action_data$next_state[i])
      prob <- 1 #action_data$Probability_transition[i]
      transition_matrix[from, to] <- prob
    }
    for (state in states) {# add absorbing state if the action is not valid for any resulting state. 
      row_sum <- sum(transition_matrix[state, ])
      if (row_sum < 1) {
        transition_matrix[state, absorbing_state] <- 1 - row_sum
      }
    }
  }
  transition_matrices[[action]] <- transition_matrix
}
transition_matrices

# -------------------------------- Rewards -------------------------------- 
rewards<- adj_list %>% dplyr::select(c("current_state", "action_type", "next_state", "reward")) %>%
  arrange(next_state)
rewards2 <- as.data.frame(rewards) %>% dplyr::rename(start.state=next_state) %>% dplyr::rename(action=action_type) %>% dplyr::rename(value=reward) %>% dplyr::select(-c("current_state"))
rewards2 <- rewards2 %>%
  mutate(value = ifelse(action == "none", -1, value))
absorbing_state_rewards <- merge(actions, absorbing_state, by = NULL) %>% dplyr:: rename(action = x) %>% dplyr:: rename(start.state = y) 
absorbing_state_rewards$value <- -Inf
rewards3<- rbind(rewards2, absorbing_state_rewards)
rewards3_fixed <- rewards3 %>%
  mutate(
    end.state = "*", 
    observation = "*"
  ) %>%
  dplyr::select(action, start.state, end.state, observation, value) # Ensure correct order

# -------------------------------- Observation -------------------------------- 
# make observation matrices structure
observation_matrices <-transition_matrices
quality<- sce@metadata[["Clones"]] %>%filter(Group == "Complete") %>% dplyr::select(c("Clone", "variants", "GQ_med"))
quality$Clone <- gsub("_", "", quality$Clone ) %>% as.numeric()
quality$DP_med <- quality$GQ_med /100

# make a quality matrix for just valid transitions for now. 
variants<- unique(quality$variants)
for (variant in variants){
  quality_for_var<- quality %>% filter(variants == variant)
  for(state_index in 1:nrow(observation_matrices[[variant]])){
    clone_state<- rownames(observation_matrices[[variant]])[state_index]
    cat("Clone state:",clone_state , "\n")
    quality_score<- quality_for_var %>% filter(Clone== clone_state) %>% pull("DP_med")
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
    leftover_prob<- 1- sum_of_probs
    clone_if_dropout<- clone_dropout_table %>% filter(Clone == clone_state) %>% pull(Neighbor)
    cat("potential clones that", clone_state, "could actually be if dropout", clone_if_dropout, "\n")
    colindex <- which(colnames(observation_matrices[[variant]]) == clone_if_dropout)
    observation_matrices[[variant]][ state_index,colindex] <-leftover_prob
  }
  cat("observation matrix for this variant:", variant, ":\n")
  print(observation_matrices[[variant]])
}

# -------------------------------- POMDP -------------------------------- 
start <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)
horizon_val <- Inf
sc <- POMDP(
  name = "sc",
  discount = 0.5,# if 1,  values future rewards as much as immediate rewards
  states = c(states),
  actions = c(actions),
  horizon=horizon_val, # number of epochs. 
  observations = c(states),
  #start = start,
  transition_prob = transition_matrices,
  observation_prob = observation_matrices,
  reward = rewards3_fixed
)

g <- transition_graph(sc, action = NULL, episode = NULL, 
                      epoch = NULL, state_col = NULL, simplify_transitions = TRUE)
plot(g, edge.label= NA,
     edge.label.cex = 0.5,
     vertex.label.cex = 0.5,
     vertex.size = 30, edge.arrow.size = .3, margin = .5)


# c("grid", "enum", "twopass", "witness", "incprune")
solution <- solve_POMDP(sc, method = "incprune", verbose = TRUE)
g_vis<- plot_policy_graph(solution, engine = "visNetwork", show_belief =TRUE)
g_vis

sim<-simulate_POMDP(solution,
                    n = 1000, 
                    horizon = 10, 
                    return_beliefs = TRUE,
                    return_trajectories = TRUE)

identify_clones_that_this_clone_could_actually_be<- function(clones){
  clones <- clones[clones != "absorbing_state"]
  clones_num <- as.numeric(clones)
  clone_neighbors <- list()
  for (clone in clones_num) {
    neighbors <- clones_num[clones_num == clone - 1]
    clone_neighbors[[as.character(clone)]] <- as.character(neighbors)
  }
  clone_df <- do.call(rbind, lapply(names(clone_neighbors), function(clone) {
    if (length(clone_neighbors[[clone]]) == 0) {
      data.frame(Clone = clone, Neighbor = NA)
    } else {
      data.frame(Clone = clone, Neighbor = clone_neighbors[[clone]])
    }
  }))
  print(clone_df)
  return(clone_df)
}
