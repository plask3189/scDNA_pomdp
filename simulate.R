sim<-simulate_POMDP(solution,n = 20, horizon=10)
# calculate the percentage that each action is used in the simulation
round_stochastic(sim$action_cnt / sum(sim$action_cnt), 2)
# reward distribution
hist(sim$reward)

sim<-simulate_POMDP(solution,
                    belief = start,
                    n = 20, # number of times to run the simulation
                    horizon = 10, # The simulation horizon (maximum number of time steps per episode)
                    return_beliefs = TRUE,
                    method= "trajectories",
                    return_trajectories = TRUE,
                    verbose = TRUE); sim
# calculate the percentage that each action is used in the simulation
round_stochastic(sim$action_cnt / sum(sim$action_cnt), 2)
# reward distribution
hist(sim$reward)
head(sim$belief_states)
head(sim$trajectories)

trajectories <- sim$trajectories # each row of this df is a time step. 
(trajectories)



# Process trajectories to summarize state transitions per episode
transitions <- trajectories %>%
  # Group by each episode to ensure operations are episode-specific
  group_by(episode) %>%
  dplyr::mutate(next_state = lead(simulation_state)) %>%
  # Remove rows where there is no next state (e.g., the last row in each episode)
  filter(!is.na(next_state)) %>%
  ungroup() %>%
  # Count the occurrences of each unique combination of current state, next state, and action
  dplyr::count(simulation_state, next_state, a)
transitions
edges <- transitions %>%
  mutate(label = paste(a, "(", n, ")", sep = "")) %>%
  dplyr::select(simulation_state, next_state, label, n) # n is number of times the state transition occurs. 

g <- graph_from_data_frame(edges, directed = TRUE)
plot(g, edge.label = E(g)$label,
     vertex.label = V(g)$name,
     edge.arrow.size = 0.1,
     edge.label.cex = 0.5,
     vertex.label.cex = 0.5,
     vertex.size = 10,
     layout = layout_with_kk)


# Calculate cumulative reward per episode
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

