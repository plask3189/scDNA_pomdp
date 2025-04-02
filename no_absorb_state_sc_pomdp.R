
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
final_vis<-visualize_tree(sce, variants_of_interest, remove_low_reward_edges = TRUE); final_vis
#---------POMDP-----------
mutation_states<-length(unique(sce@metadata$Architecture$final_annot)) ; mutation_states
paste("Building MDP initially with all possible mutation combinations = ", mutation_states) # different mutations
adj_list<-BuildMDP(mutation_states,use_ADO= FALSE)
adj_list<-attach_weights_kp(sce,adj_list) %>% arrange(current_state)
max_len <- max(nchar(c(as.character(adj_list$current_state),
                       as.character(adj_list$next_state))))

clones_with_potential_dropout<- make_potential_dropout_clones(unique(levels(adj_list$current_state)))
probabilities<- compute_probability_of_dropout(clones_with_potential_dropout)

adj_list <- adj_list %>%
  mutate(
    current_state = sprintf(paste0("%0", max_len, "d"), as.numeric(as.character(current_state))),
    next_state    = sprintf(paste0("%0", max_len, "d"), as.numeric(as.character(next_state)))
  ) %>%
  mutate(
    current_state = gsub("(?<=.)(?=.)", "_", current_state, perl = TRUE),
    next_state    = gsub("(?<=.)(?=.)", "_", next_state, perl = TRUE)
  )

states<-unique((adj_list$current_state)); cat(blue("States:", paste(states, collapse = ", ")))
actions<- setdiff(unique(adj_list$action_type), "none"); actions

# probabilities calculated by quantity of current_state
data_for_transition_matrix<- adj_list %>% dplyr::select(c("current_state", "action_type", "next_state", "Probability_transition")) ;data_for_transition_matrix
transition_matrices <- list()
# Generate a transition matrix for each action type
for (action in actions) {
  transition_matrix <- matrix(0, nrow = length(states), ncol = length(states), dimnames = list(states, states))
  action_data <- data_for_transition_matrix %>%  
    filter(action_type == action) %>%
    dplyr::select(current_state, next_state, Probability_transition)
  for (i in seq_len(nrow(action_data))) {
    from <- as.character(action_data$current_state[i])
    to <- as.character(action_data$next_state[i])
    prob <- 1 
    transition_matrix[from, to] <- prob
  }
  transition_matrices[[action]] <- transition_matrix
}
for (action in seq_along(transition_matrices)) {
  matrix <- transition_matrices[[action]]
  for (i in seq_len(nrow(matrix))) {
    if (sum(matrix[i, ]) == 0)
      matrix[i, i] <- 1
  }
  transition_matrices[[action]] <- matrix 
}
transition_matrices

# ----------------------start -----------------------

vector_of_0s <- rep(0, length(states))     # create a vector of 0s
uniform_start<- vector_of_0s
uniform_start[1:(length(vector_of_0s)-1)] <-1/(length(vector_of_0s)-1);uniform_start

wt_start<- vector_of_0s
wt_start[1:1]<-1;wt_start
#-----------------------------------------

for (action in seq_along(transition_matrices)) {
  matrix <- transition_matrices[[action]]
  for (i in seq_len(nrow(matrix))) {
    if (i == nrow(matrix)){# if the last clone
      restart_probs<- 1/(ncol(matrix)-1)
      matrix[i,]<-uniform_start
      transition_matrices[[action]]<- matrix
    }
  }
}

transition_matrices
#---------------------------------- Observations ---------------------------------------
# For a POMDP, there is uncertainty about current state

# observations<- c("nothing", "goal")
# observation_matrices<- list()
# for (action in actions) {
#   observation_matrix <- matrix(0, nrow = length(states), ncol = length(observations), dimnames = list(states, observations))
#   observation_matrix[,"nothing"]<-1
#   observation_matrix[nrow(observation_matrix),]<-1-observation_matrix[nrow(observation_matrix),]
#   observation_matrices[[action]]<- observation_matrix
# }
# observation_matrices
observation_matrices<- list()
observable <- TRUE
if(observable == TRUE){
  for (action in actions) { # 1 for sames
    observation_matrices[[action]]<- diag(1, length(states), length(states))
  }
} else{
  quality<- sce@metadata[["Clones"]] %>% filter(Group == "Complete")
  all_candidates <- unique(unlist(probabilities$could_be_heard_as))
  obs_result_matrix <- matrix(0, 
                              nrow = nrow(probabilities), 
                              ncol = length(all_candidates),
                              dimnames = list(probabilities$og_mutation, all_candidates))
  for (i in seq_len(nrow(probabilities))) {
    for (candidate in probabilities$could_be_heard_as[[i]]) {
      obs_result_matrix[i, candidate] <- probabilities$prob[i]
    }
  }
  
  observations<- states
  for (action in actions) { # instead of just probs of dropout, double the same so more likely legit
    observation_matrices[[action]]<- obs_result_matrix
    obs_matrix<- observation_matrices[[action]]
    for (i in seq_len(nrow(obs_matrix))){
      obs_matrix[i,i]<- obs_matrix[i,i]*2 #double the same sames
    }
    observation_matrices[[action]]<- obs_matrix
  }
}

normalize_rows <- function(mat) {
  row_sums <- rowSums(mat)
  row_sums[row_sums == 0] <- 1  # Avoid division by zero
  return(mat / row_sums)
}
observation_matrices <- lapply(observation_matrices, normalize_rows)
observation_matrices



# -------------------------------- Rewards -------------------------------- 
max_clone<- colnames(transition_matrices[[1]])[length(colnames(transition_matrices[[1]]))];max_clone
wt_clone<- colnames(transition_matrices[[1]])[1];wt_clone

adj_list_reward <- adj_list %>% filter(action_type !="none") %>% arrange(next_state)
rewards3_fixed <- data.frame(action = character(),start.state = character(),end.state = character(),observation = character(),v = numeric(),stringsAsFactors = FALSE)
for(i in seq_len(nrow(adj_list_reward))){
  rewards3_fixed <- rbind(rewards3_fixed,R_("*", start.state = as.character(adj_list_reward$next_state[i]), end.state = "*", observation = "*", v = adj_list_reward$reward[i]))
}
max_reward<-max(rewards3_fixed$value);max_reward
rewards3_fixed <- rbind(rewards3_fixed, R_("*", start.state = "*", end.state = max_clone, observation = max_clone, v = (max_reward+2)))
rewards3_fixed
# -------------------------------- POMDP -------------------------------- 
#horizon_val<-length(states)


horizon_val <-Inf
sc <- POMDP(
  name = "sc",
  discount = 0.85,# if 1,  values future rewards as much as immediate rewards
  states = c(states),
  actions = c(actions),
  start= "uniform",
  horizon=horizon_val, # number of timesteps in episode
  #start = wt_start,
  observations = observations,
  transition_prob = transition_matrices,
  observation_prob = observation_matrices,
  reward = rewards3_fixed
)


g <- transition_graph(sc, action = NULL, episode = NULL,epoch = NULL, state_col = NULL, simplify_transitions = TRUE)
unique_labels <- unique(E(g)$label)
label_colors <- setNames(rainbow(length(unique_labels)), unique_labels)
E(g)$color <- label_colors[E(g)$label]
plot(g, edge.label =NA,
     edge.label.cex = 0.5,  
     vertex.label.cex = 0.5,
     vertex.size = 20, 
     edge.arrow.size = 0.1, 
     layout = layout_in_circle)
legend("topright", 
       legend = names(label_colors), 
       col = label_colors, 
       lty = 1,cex = 0.4)

mdp_sc<- MDP(c(states), c(actions), transition_matrices, rewards3_fixed, discount = 0.85,
             horizon = Inf, start = "uniform", name = "mdp_sc")
# c("grid", "enum", "twopass", "witness", "incprune")
solution <- solve_POMDP(sc, method = "grid", verbose = TRUE, parameter = list(timeout = 10))

mdp_solution<- solve_MDP(mdp_sc, method = "value_iteration", verbose = TRUE)



plot_transition_graph(solution, vertex.size = 20, edge.arrow.size = .3, margin = .5, edge.label.cex = 0.4,vertex.label.cex = 0.4)

g_vis<- plot_policy_graph(solution, engine = "visNetwork", show_belief =TRUE)
g_vis
g<- plot_policy_graph(solution,
                      edge.label.cex = 0.5,
                      vertex.label.cex = 0.5,
                      vertex.size = 20,
                      layout = layout_in_circle); g



sim<-simulate_POMDP(solution,
                      n = 10, # number of times to run the simulation
                      #horizon =10,
                      #initial_belief = wt_start,
                      return_beliefs = TRUE,
                      method= "trajectories",
                      return_trajectories = TRUE)

trajectories <- sim$trajectories
cumulative_rewards <- trajectories %>%
  group_by(episode) %>%
  arrange(time) %>%
  mutate(cumulative_reward = cumsum(r))

ggplot(cumulative_rewards, aes(x = time, y = cumulative_reward, color = factor(episode))) +
  geom_line() +
  theme_minimal() +
  labs(title = "Cumulative Rewards per Episode",
       x = "Time Step",
       y = "Cumulative Reward",
       color = "Episode")


make_transition_tree(trajectories)







# --------------- Helper function --------------------
compute_probability_of_dropout<- function(clones_with_potential_dropout){
  
  probabilities <- clones_with_potential_dropout %>%
    mutate(prob = 1 / sapply(could_be_heard_as, length))
  probabilities <- probabilities %>%
    mutate(
      og_mutation = sprintf(paste0("%0", max_len, "d"), as.numeric(as.character(og_mutation)))
    ) %>%
    mutate(
      og_mutation = gsub("(?<=.)(?=.)", "_", og_mutation, perl = TRUE),
    )
  for (row in seq_len(nrow(probabilities))) {
    formatted <- sprintf(
      paste0("%0", max_len, "d"),
      as.numeric(as.character(probabilities$could_be_heard_as[[row]]))
    )
    formatted <- gsub("(?<=.)(?=.)", "_", formatted, perl = TRUE)
    probabilities$could_be_heard_as[[row]] <- formatted
  }; print(probabilities)
  return(probabilities)
}


make_transition_tree<- function(trajectories){
  sim$trajectories <- sim$trajectories %>%
    mutate(simulation_state = as.character(simulation_state))
  
  transitions <- sim$trajectories %>%
    group_by(episode) %>%
    arrange(time) %>%
    mutate(next_state = lead(simulation_state),
           action_taken = a) %>%
    filter(!is.na(next_state)) %>%
    ungroup() %>%
    dplyr::select(from = simulation_state, to = next_state, action = action_taken) %>%
    filter(from != max_clone) %>% 
    filter(from != to)%>% 
    arrange(from)
  transitions
  transitions <- transitions %>% 
    mutate(across(where(is.factor), as.character))
  
  g <- graph_from_data_frame(transitions, directed = TRUE)
  paths <- all_simple_paths(g, from = wt_clone, mode = "out")
  paths_list <- lapply(paths, function(path) V(g)[path]$name)
  
  print(paths_list)
  
  nodes_list <- list()   # will store node id and label
  edges_list <- list()   # will store edges as data frames
  
  for(traj in paths_list) {
    for(i in seq_along(traj)) {
      action<- transitions %>% filter(from ==traj[i] & to ==traj[i+1]) %>% pull(action)
      action<- unique(action)
      cat("i:", traj[i], "takes action:",action , "\n")
      node_id <- paste(traj[1:i], collapse = "-")
      node_label <- traj[i]
      nodes_list[[node_id]] <- node_label
      
      # For every step beyond the first, create an edge from the previous node to the current one
      if(i > 1) {
        from_id <- paste(traj[1:(i-1)], collapse = "-")
        to_id <- node_id
        edges_list[[length(edges_list) + 1]] <- data.frame(from = from_id, to = to_id, stringsAsFactors = FALSE)
      }
    }
  }
  nodes_df <- data.frame(
    id = names(nodes_list),
    label = unlist(nodes_list),
    stringsAsFactors = FALSE
  )
  
  edges_df <- do.call(rbind, edges_list)
  edges_df <- unique(edges_df)
  
  visNetwork(nodes_df, edges_df)%>% 
    visEdges(arrows = "to") %>% 
    visHierarchicalLayout()
  
  
}
