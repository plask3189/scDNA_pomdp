suny<- solve_SARSOP("/Users/kateplas/Downloads/sunysb.POMDP 2", parameter = list(timeout = 10))


library(pomdp)
library(dplyr)
library(ggplot2)
library(igraph)
#----------------------
sol <- solve_SARSOP("http://www.pomdp.org/examples/cheese.95.POMDP", parameter = list(timeout = 10))
sol <- solve_SARSOP("https://www.pomdp.org/examples/shuttle.95.POMDP", parameter = list(timeout = 10))
sol<- solve_SARSOP("https://www.pomdp.org/examples/stand-tiger.95.POMDP", parameter = list(timeout = 10))
sol[["transition_prob"]]
plot_transition_graph(sol, vertex.size = 20, edge.arrow.size = .3, margin = .5, edge.label.cex = 0.4)

absorbing_states(sol)

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

# ------------------
data("RussianTiger")
print(RussianTiger)

RussianTiger$states  
RussianTiger$actions 
RussianTiger[["reward"]] # reward (-Inf indicates unavailable actions)
RussianTiger[["transition_prob"]]
RussianTiger[["observation_prob"]]


plot_transition_graph(RussianTiger, vertex.size = 30, edge.arrow.size = .3, margin = .5, vertex.label.cex = 0.5, edge.label.cex = 0.5)

# absorbing states
absorbing_states(RussianTiger)

solution_RussianTiger <- solve_POMDP(RussianTiger, 
                                     method = "grid")


# Check the updated observation probabilities
print(RussianTiger$transition_prob)

solution_RussianTiger <- solve_POMDP(RussianTiger, 
                              method = "grid")
                             # initial_belief = c(0.8,0.2)) # if we take an action, we get to a different belief state. (probabilities of states)
print(solution_RussianTiger)

g_vis_RussianTiger<- plot_policy_graph(solution_RussianTiger, engine = "igraph",
                                vertex.label.cex = 0.5,
                                edge.label.cex = 0.5,)

g_vis_RussianTiger
belief_space <- plot_belief_space(solution_RussianTiger, epoch = 9)

g_vis_RussianTiger<- plot_policy_graph(solution_RussianTiger, engine = "visNetwork"); g_vis_RussianTiger

sim_rt<-simulate_POMDP(solution_RussianTiger,
                    n = 20, # number of times to run the simulation
                    horizon =10,
                    return_beliefs = TRUE,
                    method= "trajectories",
                    return_trajectories = TRUE)

sim_rt
trajectories_rt <- sim_rt$trajectories 

# ---------------------- 


library(pomdp)
data("Tiger")

Tiger[["states"]]
Tiger[["actions"]]
Tiger[["transition_prob"]]
Tiger[["observation_prob"]]
Tiger[["reward"]]

plot_transition_graph(Tiger, vertex.size = 30, edge.arrow.size = .3, margin = .5, edge.label.cex = 0.5,vertex.label.cex = 0.5)

solution_Tiger <- solve_POMDP(Tiger, method = "grid", verbose = TRUE)

g_vis_Tiger<- plot_policy_graph(solution_Tiger, engine = "igraph",
                                       vertex.label.cex = 0.5,
                                       edge.label.cex = 0.5,)

g_vis_Tiger
belief_space <- plot_belief_space(solution_Tiger, epoch = 9)

g_vis_Tiger<- plot_policy_graph(solution_Tiger, engine = "visNetwork"); g_vis_Tiger

sim_t<-simulate_POMDP(solution_Tiger,
                       n = 20, # number of times to run the simulation
                       horizon =10,
                       return_beliefs = TRUE,
                       method= "trajectories",
                       return_trajectories = TRUE)

trajectories_t <- sim_t$trajectories 
transitions_t <- trajectories_t %>%
  # Group by each episode to ensure operations are episode-specific
  group_by(episode) %>%
  dplyr::mutate(next_state = lead(simulation_state)) %>%
  # Remove rows where there is no next state (e.g., the last row in each episode)
  filter(!is.na(next_state)) %>%
  ungroup() %>%
  # Count the occurrences of each unique combination of current state, next state, and action
  dplyr::count(simulation_state, next_state, a)
transitions
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








