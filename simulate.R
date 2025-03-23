
sim<-simulate_POMDP(solution,
                    n = 20, # number of times to run the simulation
                    horizon = Inf, 
                    return_beliefs = TRUE,
                    method= "trajectories",
                    return_trajectories = TRUE)



sim



trajectories <- sim$trajectories # each row of this df is a time step. 
# a: The action taken at the current time step. 
# o: The observation received after taking the action.

# Process trajectories to summarize state transitions per episode
transitions <- trajectories %>%
  # Group by each episode to ensure operations are episode-specific
  group_by(episode) %>%
  # Create a new column 'next_state' which is the state at the next time step
  dplyr::mutate(next_state = lead(simulation_state)) %>%
  # Remove rows where there is no next state (e.g., the last row in each episode)
  filter(!is.na(next_state)) %>%
  ungroup() %>%
  # Count the occurrences of each unique combination of current state, next state, and action
  dplyr::count(simulation_state, next_state, a)

edges <- transitions %>%
  mutate(label = paste(a, "(", n, ")", sep = "")) %>%
  dplyr::select(simulation_state, next_state, label, n) # n is number of times the state transition occurs. 

g <- graph_from_data_frame(edges, directed = TRUE)
plot(g, edge.label = E(g)$label,
     vertex.label = V(g)$name,
     edge.arrow.size = 0.5,
     
     vertex.size = 20,
     layout = layout_with_kk)


# Calculate cumulative reward per episode
cum_rewards <- trajectories %>%
  group_by(episode) %>%
  arrange(time) %>%
  mutate(cum_reward = cumsum(r))

ggplot(cum_rewards, aes(x = time, y = cum_reward, color = factor(episode))) +
  geom_line() +
  theme_minimal() +
  labs(title = "Cumulative Rewards per Episode",
       x = "Time Step",
       y = "Cumulative Reward",
       color = "Episode")
