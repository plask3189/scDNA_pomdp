simulate<- function(solution){

  #'Simulate trajectories through a POMDP.
  #' The start state for each trajectory is randomly chosen using the specified belief. 
  #' The belief is used to choose actions from the the epsilon-greedy policy and then 
  #' updated using observations.
  sim<-simulate_POMDP(solution,
                      n = 500, # number of times to run the simulation
                      horizon = 5, # The simulation horizon (maximum number of time steps per episode)
                      method= "trajectories",
                      epsilon= NULL,
                      initial_belief = start,
                      return_beliefs = TRUE,
                      return_trajectories = TRUE,
                      verbose = TRUE); sim
  # calculate the percentage that each action is used in the simulation
  round_stochastic(sim$action_cnt / sum(sim$action_cnt), 2)
  # reward distribution
  hist(sim$reward)
  head(sim$belief_states)
  head(sim$trajectories)
  
  trajectories <- sim$trajectories # each row of this df is a time step. 
  
  cat("head of trajectories:\n"); print(head(trajectories))
  
  
  # Process trajectories to summarize state transitions per episode
  transitions <- trajectories %>%
    group_by(episode) %>%
    mutate(next_state = lead(o)) %>%  # Use observed state for branching
    filter(!is.na(next_state)) %>%
    ungroup() %>%
    count(simulation_state, next_state, a)
  transitions
  
  transitions<- transitions %>% filter(a != "listen") %>%
    filter(next_state != "absorbing_state")
  
  edges <- transitions %>%
    mutate(label = paste(a, "(", n, ")", sep = "")) %>%
    dplyr::select(simulation_state, next_state, label, n) # n is number of times the state transition occurs. 

 
    
  # edges <- edges %>%
  #   filter(simulation_state != "absorbing_state" & next_state != "absorbing_state")
  # edges <- edges %>%
  #   filter(!grepl("listen", label))
  
  g <- graph_from_data_frame(edges, directed = TRUE)
  plot(g, edge.label = E(g)$label,
       vertex.label = V(g)$name,
       edge.arrow.size = 0.3,
       edge.label.cex = 0.5,
       vertex.label.cex = 0.5,
       vertex.size = 10,
       layout = layout_in_circle)
  x
  
  # Calculate cumulative reward per episode
  cumulative_rewards <- trajectories %>%
    group_by(episode) %>%
    arrange(time) %>%
    mutate(cumulative_reward = cumsum(r))
  
  ggplot(cumulative_rewards, aes(x = time, y = cumulative_reward, color = factor(episode))) +
    geom_line() +
    theme_minimal()
  
  ggplot(cumulative_rewards, aes(x = time, y = cumulative_reward, color = factor(episode))) +
    geom_line() +
    theme_minimal() +
    labs(title = "Cumulative Rewards per Episode",
         x = "Time Step",
         y = "Cumulative Reward",
         color = "Episode")
  
  
}
