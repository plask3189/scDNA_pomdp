# the columns (clone states) are observations.
# policy graph shows actions for each belief state(probabilities of states).
# belief first take action to get to new belief state. 

g_vis<- plot_policy_graph(solution, engine = "visNetwork", show_belief =TRUE)
g_vis


g_vis
b<- estimate_belief_for_nodes(solution, method = "trajectories")
p<- policy_graph(solution)
plot(p, 
     layout = layout_as_tree(p, root = 3, mode = "out"),
     edge.curved = curve_multiple(p, .2))


g<- plot_policy_graph(solution, belief,
                      vertex.label.cex = 0.5,
                      edge.label.cex = 0.5,
                      edge.curved = 0.3,
                      edge.arrow.size = 0.5,
                      vertex.size = 20,
                      vertex.label.dist=0.5,
                      layout = layout.circle)
g
sim<-simulate_POMDP(solution,
                    n = 1000, 
                    belief = start,
                    horizon = 10, 
                    return_beliefs = TRUE,
                    return_trajectories = TRUE)


trajectories<- sim[["trajectories"]]
trajectories
# Assuming simulate_POMDP is a list containing reward and action sequences

# Find the index of the maximum reward
max_reward_index <- which.max(sim$reward)

# Extract the corresponding actions and transitions (if available)
most_likely_trajectory <- sim$action_cnt[max_reward_index]

# Print the most likely trajectory
cat("Most likely trajectory (highest reward):\n")
print(most_likely_trajectory)

hmm <- initHMM(states, observations, transProbs, obsProbs, initProbs)

# t <- transition_graph(sc)
# plot_transition_graph(sc, vertex.size = 20, edge.label.cex = 0.8, layout = layout.circle)
# head(solution[["solution"]][["belief_points_solver"]] )
# head(solution[["solution"]][["alpha"]] )

# pg <- policy_graph(solution, show_belief = TRUE,
#                    simplify_observations = TRUE, remove_unreachable_nodes = TRUE)
# plot(pg, layout = layout_as_tree(pg, root = 3, mode = "out"))
# 

# want to get most likely belief state (distribution of states) for each action



  
  
  
explain_policy_graph <- function(sample_policy_graph) {
  for (i in 1:nrow(sample_policy_graph)){
    row<- (sample_policy_graph[i,])
    cat(blue("After taking action", row$action, "\n"))
    observed_clones<- (colnames(sample_policy_graph)[3:length(colnames(sample_policy_graph))])
    for (c in seq_along(observed_clones)){
      new_node<- sample_policy_graph[i,c+3]
      next_mutation<- sample_policy_graph %>% filter(node == new_node) %>% pull(action)
      cat("the agent will transition to Node", new_node, "when", observed_clones[c], "is observed. \n")
    }
  }
}


# Load necessary library
library(rvest)

# Sample data: g_vis[["x"]][["nodes"]][["title"]]
titles <- g_vis[["x"]][["nodes"]][["title"]]
extract_node_info <- function(html_title) {
  print(html_title)
  print("----------")
  parsed_html <- read_html(html_title)
  node_id <- parsed_html %>%
    html_node(xpath = "//b[contains(text(), 'node id:')]/following-sibling::text()[1]") %>%
    html_text() %>%
    str_trim()  # Trim leading/trailing spaces
  # Convert to numeric (if extraction is successful)
  node_id <- as.integer(node_id)
  # Extract action
  action <- parsed_html %>%
    html_node(xpath = "//b[contains(text(), 'action:')]/following-sibling::text()[1]") %>%
    html_text() %>%
    str_trim()
  # Extract belief table
  belief_table <- parsed_html %>%
    html_node("table") %>%
    html_table(header = TRUE) # Ensure headers are considered
  # Check for missing or empty column names and assign defaults if necessary
  if (any(is.na(names(belief_table)) | names(belief_table) == "")) {
    names(belief_table) <- paste0("V", seq_along(belief_table))
  }
  # Add node_id and action to each row of the belief_table
  belief_table <- belief_table %>%
    mutate(node_id = node_id, action = action)
  return(belief_table)
}

# Apply the function to each title
node_info_list <- lapply(titles, extract_node_info)

# Combine all data frames into one
node_info_df <- bind_rows(node_info_list)

# Display the combined data frame
print(node_info_df)
