
library(pomdp)
library(dplyr)
library(ggplot2)
library(igraph)
#----------------------
problem <-"https://www.pomdp.org/examples/mini-hall2.POMDP"
problem<- "https://www.pomdp.org/examples/tiger.95.POMDP"
problem<- "/Users/kateplas/Downloads/tiger-grid.POMDP 2"
problem <- "/Users/kateplas/Downloads/sunysb.POMDP 2"
problem<-"https://www.pomdp.org/examples/1d.POMDP"


sol<- solve_POMDP(problem,method = "grid", parameter = list(timeout = 10), verbose = TRUE)
sol <- POMDP(
  name = "sc",
  discount = 0.95,# if 1,  values future rewards as much as immediate rewards
  states = sol$states,
  actions = sol$actions,
  start= "uniform",
  horizon=Inf, # number of timesteps in episode
  observations = sol[["observations"]],
  transition_prob = sol[["transition_prob"]],
  observation_prob = sol[["observation_prob"]],
  reward =  sol[["reward"]]
)
sol <- solve_POMDP(sol, method = "grid", verbose = TRUE, parameter = list(timeout = 10))
       
g_vis_sol<- plot_policy_graph(sol, engine = "igraph", show_belief =TRUE)
g_vis_sol
g_vis_sol<- plot_policy_graph(sol, engine = "visNetwork", show_belief =TRUE)
g_vis_sol

sol[["transition_prob"]]
plot_transition_graph(sol, vertex.size = 20, edge.arrow.size = .3, margin = .5, edge.label.cex = 0.4)



sim_t<-simulate_POMDP(sol,
                      n = 20, # number of times to run the simulation
                      horizon =10,
                      return_beliefs = TRUE,
                      method= "trajectories",
                      return_trajectories = TRUE)

# calculate the percentage that each action is used in the simulation
round_stochastic(sim_t$action_cnt / sum(sim_t$action_cnt), 2)
# reward distribution
hist(sim_t$reward)
trajectories_t <- sim_t$trajectories 
transitions_t <- trajectories_t %>%
  # Group by each episode to ensure operations are episode-specific
  group_by(episode) %>%
  dplyr::mutate(next_state = lead(simulation_state)) %>%
  # Remove rows where there is no next state (e.g., the last row in each episode)
  filter(!is.na(next_state)) %>%
  ungroup() %>%
  dplyr::count(simulation_state, next_state, a) # Count the occurrences of each unique combination of current state, next state, and action
transitions_t
edges <- transitions_t %>%
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
cumulative_rewards <- trajectories_t %>%
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
