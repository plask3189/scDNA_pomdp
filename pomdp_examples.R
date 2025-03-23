
library(pomdp)

data("Tiger")
print(Tiger)

Tiger$states  
Tiger$actions 
Tiger$observations
# reward (-Inf indicates unavailable actions)
Tiger$reward
Tiger[["reward"]]
plot_transition_graph(Tiger, vertex.size = 30, edge.arrow.size = .3, margin = .5, edge.label.cex = 0.5)

# absorbing states
absorbing_states(Tiger)

solution_Tiger <- solve_POMDP(Tiger, 
                                     method = "grid")
Tiger$transition_prob <- list(
  listen = diag(length(Tiger[["states"]])),  # Identity matrix (listening does not change the state)
  "open-left" = diag(length(Tiger[["states"]])),
  "open-right" = diag(length(Tiger[["states"]]))
)
# 
# Tiger$observation_prob <- list(
#   # Improved accuracy for "listen" action
#   listen = matrix(c(0.90, 0.10,  # If the actual state is Tiger-left
#                     0.10, 0.90), # If the actual state is Tiger-right
#                   nrow = length(Tiger[["states"]]), byrow = TRUE,
#                   dimnames = list(Tiger$states, Tiger$observations)),
# 
#   # Biased observation when opening left
#   "open-left" = matrix(c(0.30, 0.70,  # If the actual state is Tiger-left
#                          0.70, 0.30), # If the actual state is Tiger-right
#                        nrow = length(Tiger[["states"]]), byrow = TRUE,
#                        dimnames = list(Tiger$states, Tiger$observations)),
# 
#   # Biased observation when opening right
#   "open-right" = matrix(c(0.70, 0.30,  # If the actual state is Tiger-left
#                           0.30, 0.70), # If the actual state is Tiger-right
#                         nrow = length(Tiger[["states"]]), byrow = TRUE,
#                         dimnames = list(Tiger$states, Tiger$observations))
# )

# Check the updated observation probabilities
print(Tiger$transition_prob)
# 
# Tiger$reward <- Tiger$reward %>%
#   mutate(end.state = case_when(
#     action == "listen" ~ "listening",
#     action == "open-left" ~ "game-over",
#     action == "open-right" ~ "game-over",
#     TRUE ~ NA_character_
#   ))
# print(Tiger$reward)
solution_Tiger <- solve_POMDP(Tiger, 
                              method = "grid")
                             # initial_belief = c(0.8,0.2)) # if we take an action, we get to a different belief state. (probabilities of states)
print(solution_Tiger)

g_vis_Tiger<- plot_policy_graph(solution_Tiger, engine = "igraph",
                                vertex.label.cex = 0.5,
                                edge.label.cex = 0.5,)

g_vis_Tiger
belief_space <- plot_belief_space(solution_Tiger, epoch = 9)

g_vis_Tiger<- plot_policy_graph(solution_Tiger, engine = "visNetwork"); g_vis_Tiger

