
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
source("./R/RL_dev/generate_potential_dropout_clones.R")
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
adj_list<-attach_weights_kp(sce,adj_list) %>% arrange(current_state)

states<-unique(levels(adj_list$current_state))
absorbing_state <- "absorbing_state"
states <- c(states, absorbing_state)  ; cat(blue("States:", paste(states, collapse = ", ")))
actions<-c("listen", unique(adj_list$action_type)); actions

data_for_transition_matrix<- adj_list %>% dplyr::select(c("current_state", "action_type", "next_state", "Probability_transition")) ;data_for_transition_matrix
transition_matrices <- list()
# Generate a transition matrix for each action type
for (action in actions) {
  transition_matrix <- matrix(0, nrow = length(states), ncol = length(states), dimnames = list(states, states))
  if(action == "none" |action == "listen" ){ #make identity matrix bc if no change then current_state == next_state.
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

#---------------------------------- Observations ---------------------------------------
observation_matrices<- make_observation_prob_matrix3(sce, transition_matrices)
# -------------------------------- Rewards -------------------------------- 
rewards3_fixed<- make_rewards3(adj_list)
# -------------------------------- POMDP -------------------------------- 
vector_of_0s <- rep(0, length(states))     # create a vector of 0s
# start 
start<- vector_of_0s
start[1] <- 1  # 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

#terminal
terminal <- vector_of_0s
terminal[length(vector_of_0s)] <- 1 


horizon_val <-Inf
sc <- POMDP(
  name = "sc",
  discount = 0.95,# if 1,  values future rewards as much as immediate rewards
  states = c(states),
  actions = c(actions),
  horizon=horizon_val, # number of epochs. 
  observations = c(states),
  #start = vec,
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
solution <- solve_POMDP(sc, method = "grid", verbose = TRUE)
g_vis<- plot_policy_graph(solution, engine = "visNetwork", show_belief =TRUE)
g_vis
g<- plot_policy_graph(solution,
                      edge.label.cex = 0.5,
                      vertex.label.cex = 0.5,
                      vertex.size = 30,
                      layout = layout_in_circle
                      ); g


absorbing_states(sc)
simulate(solution)
