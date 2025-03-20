
library(pomdp)

data("RussianTiger")
print(RussianTiger)

RussianTiger$states  
RussianTiger$actions 
RussianTiger$observations
# reward (-Inf indicates unavailable actions)
RussianTiger$reward
RussianTiger[["reward"]]
plot_transition_graph(RussianTiger, vertex.size = 30, edge.arrow.size = .3, margin = .5, edge.label.cex = 0.5)

# absorbing states
absorbing_states(RussianTiger)

solution_RussianTiger <- solve_POMDP(RussianTiger, 
                                     method = "grid")
# RussianTiger$transition_prob <- list(
#   listen = diag(length(RussianTiger[["states"]])),  # Identity matrix (listening does not change the state)
#   "open-left" = diag(length(RussianTiger[["states"]])),
#   "open-right" = diag(length(RussianTiger[["states"]]))
# )
# 
# RussianTiger$observation_prob <- list(
#   # Improved accuracy for "listen" action
#   listen = matrix(c(0.90, 0.10,  # If the actual state is RussianTiger-left
#                     0.10, 0.90), # If the actual state is RussianTiger-right
#                   nrow = length(RussianTiger[["states"]]), byrow = TRUE,
#                   dimnames = list(RussianTiger$states, RussianTiger$observations)),
# 
#   # Biased observation when opening left
#   "open-left" = matrix(c(0.30, 0.70,  # If the actual state is RussianTiger-left
#                          0.70, 0.30), # If the actual state is RussianTiger-right
#                        nrow = length(RussianTiger[["states"]]), byrow = TRUE,
#                        dimnames = list(RussianTiger$states, RussianTiger$observations)),
# 
#   # Biased observation when opening right
#   "open-right" = matrix(c(0.70, 0.30,  # If the actual state is RussianTiger-left
#                           0.30, 0.70), # If the actual state is RussianTiger-right
#                         nrow = length(RussianTiger[["states"]]), byrow = TRUE,
#                         dimnames = list(RussianTiger$states, RussianTiger$observations))
# )

# Check the updated observation probabilities
print(RussianTiger$transition_prob)
# 
# RussianTiger$reward <- RussianTiger$reward %>%
#   mutate(end.state = case_when(
#     action == "listen" ~ "listening",
#     action == "open-left" ~ "game-over",
#     action == "open-right" ~ "game-over",
#     TRUE ~ NA_character_
#   ))
# print(RussianTiger$reward)
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

