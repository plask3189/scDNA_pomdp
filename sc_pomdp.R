
library(scDNA)
source("./R/RL_dev/make_sce.R")
source("./R/clonograph.R")
library(dplyr)
library(ggplot2)
source("./R/new_dev/visualize_clonal_evolution_kp.R")

sample_file= "BRAF/A5330braf.dna+protein.h5"
sce<- make_sce(sample_file) # normal run. Makes sce on 2 variants. 

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
    diag(transition_matrix) <- 1 # identity
  } else {
    action_data <- data_for_transition_matrix %>%  
      filter(action_type == action) %>%
      dplyr::select(current_state, next_state, Probability_transition)
    for (i in seq_len(nrow(action_data))) {
      from <- as.character(action_data$current_state[i])
      to <- as.character(action_data$next_state[i])
      prob <- 1 
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

# ------------------------Start info-------------------------------
vector_of_0s <- rep(0, length(states))     # create a vector of 0s
start<- vector_of_0s
#start[1] <- 1 # 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
start[1:(length(vector_of_0s)-1)] <-1/(length(vector_of_0s)-1)

#---------------------------------- Observations ---------------------------------------
observation_matrices<- make_observation_prob_matrix4(sce, transition_matrices)
# if do nothing go to absorbign state
# observation_matrices[["none"]][] <- 0
# observation_matrices[["none"]][,"absorbing_state"] <- 1
# 
# observation_matrices[["NRAS.G12R"]][] <- 0
# observation_matrices[["NRAS.G12R"]][,"absorbing_state"] <- 1
# 
# observation_matrices[["listen"]]["absorbing_state",] <- 0
# observation_matrices[["listen"]]["absorbing_state",1:(length(colnames(observation_matrices[["listen"]]))-1)]<- 1/(length(colnames(observation_matrices[["listen"]]))-1)
# -------------------------------- Rewards -------------------------------- 
rewards3_fixed<- make_rewards3(adj_list); rewards3_fixed

# rewards3_fixed = rbind(R_("listen", start.state = "*", end.state = "*", observation = "*",v = -10),
#                   R_("nothing", start.state = "*", end.state = "*", observation = "*",v = -Inf),
#                   R_("*", start.state = "absorbing_state", end.state = "*", observation = "*",v = -Inf),
#                   R_("nothing", start.state = "absorbing_state", end.state = "*", observation = "*",v = 0),
#                   R_("NRAS.G12R", start.state = "*", end.state = "*", observation = "*",v = 0),
#                   #R_("NRAS.G12R", start.state = "0", end.state = "1", observation = "*",v = 1),
#                   R_("NRAS.G12R", start.state = "0", end.state = "2", observation = "*",v = -Inf))
                  

# rewards3_fixed = rbind(R_(action, start.state = "*", end.state = "*", observation = "*",v = -Inf),
#                        rewards3_fixed)

# -------------------------------- POMDP -------------------------------- 

horizon_val <-Inf
sc <- POMDP(
  name = "sc",
  discount = 0.95,# if 1,  values future rewards as much as immediate rewards
  states = c(states),
  actions = c(actions),
  horizon=horizon_val, # number of timesteps in episode
  observations = c(states),
  #start = start,
  transition_prob = transition_matrices,
  observation_prob = observation_matrices,
  reward = rewards3_fixed
)

g <- transition_graph(sc, action = NULL, episode = NULL, 
                      epoch = NULL, state_col = NULL, simplify_transitions = TRUE)
unique_labels <- unique(E(g)$label)
label_colors <- setNames(rainbow(length(unique_labels)), unique_labels)
E(g)$color <- label_colors[E(g)$label]
plot(g, edge.label =NA,
     edge.label.cex = 0.5,  # you can remove this if you don't want edge labels at all
     vertex.label.cex = 0.5,
     vertex.size = 20, 
     edge.arrow.size = 0.1, 
     layout = layout_in_circle)
legend("topright", 
       legend = names(label_colors), 
       col = label_colors, 
       lty = 1,cex = 0.4)

# c("grid", "enum", "twopass", "witness", "incprune")
solution <- solve_POMDP(sc, method = "grid", initial_belief = start ,verbose = TRUE)
g_vis<- plot_policy_graph(solution, engine = "visNetwork", show_belief =TRUE)
g_vis
g<- plot_policy_graph(solution,
                      edge.label.cex = 0.5,
                      vertex.label.cex = 0.5,
                      vertex.size = 20,
                      layout = layout_in_circle); g


absorbing_states(sc)
#simulate(solution)

sim<-simulate_POMDP(solution,
                    belief = start,
                    n = 100, # number of times to run the simulation
                    horizon = 5, # The simulation horizon (maximum number of time steps per episode)
                    method= "trajectories",
                    epsilon= NULL,
                    return_beliefs = TRUE,
                    return_trajectories = TRUE,
                    verbose = TRUE); sim

