
add_absorbing_states<- function(states_modified, transition_probs, main_observation_matrix,rewards){
  absorbing_state_name<- "absorbing_state"
  states_modified_abs<- c(states_modified, absorbing_state_name) 
  # FIX TRANSITION MATRIX
  transition_probs_with_abs<- transition_probs
  for (action_name in names(transition_probs_with_abs)){
    print(action_name)
    individual_trans_prob_matrix<- transition_probs_with_abs[[action_name]]
    prob_matrix_with_absorbing_state<- make_absorbing_state_in_transition_matrix(individual_trans_prob_matrix)
    print(prob_matrix_with_absorbing_state)
    # replace the original trans matrix.
    transition_probs_with_abs[[action_name]]<- prob_matrix_with_absorbing_state
  }
  # FIX OBSERVATIONS PROBS.
  main_observation_matrix_with_abs<- main_observation_matrix
  for (action_name in names(main_observation_matrix_with_abs)){
    individual_main_observation_matrix<- main_observation_matrix_with_abs[[action_name]]
    observation_prob_matrix_with_absorbing_state<- make_absorbing_state_in_observation_matrix(individual_main_observation_matrix)
    print(observation_prob_matrix_with_absorbing_state)
    main_observation_matrix_with_abs[[action_name]]<- observation_prob_matrix_with_absorbing_state
  }
  # FIX REWARD MATRIX
  new_rows <- rewards %>%
    dplyr::select(action) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      start.state = "absorbing_state", 
      end.state = "absorbing_state", 
      observation = "*", 
      value = -100
    )
  updated_rewards <- bind_rows(rewards, new_rows)

  return(list(states = states_modified_abs, 
              transition_probs= transition_probs_with_abs, 
              observation_probs = main_observation_matrix_with_abs, 
              rewards= updated_rewards))
}

make_absorbing_state_in_observation_matrix<- function(main_observation_matrix){
  main_observation_matrix<- as.data.frame(main_observation_matrix)
  main_observation_matrix$absorbing_state <- 0
  absorbing_row <- setNames(rep(0, ncol(main_observation_matrix)), colnames(main_observation_matrix))
  absorbing_row["absorbing_state"] <- 1
  main_observation_matrix <- rbind(main_observation_matrix, absorbing_row)
  rownames(main_observation_matrix)[nrow(main_observation_matrix)] <- "absorbing_state"
  return(main_observation_matrix)
}

make_absorbing_state_in_transition_matrix<- function(normalized_action_matrix){
  normalized_action_matrix<- as.data.frame(normalized_action_matrix)
  normalized_action_matrix$absorbing_state <- 0
  for (i in 1:nrow(normalized_action_matrix)) {
    row_sum <- sum(normalized_action_matrix[i, ])
    if (row_sum < 1) {
      normalized_action_matrix[i, "absorbing_state"] <- 1 
    }
  }
  absorbing_row <- setNames(rep(0, ncol(normalized_action_matrix)), colnames(normalized_action_matrix))
  absorbing_row["absorbing_state"] <- 1
  normalized_action_matrix <- rbind(normalized_action_matrix, absorbing_row)
  rownames(normalized_action_matrix)[nrow(normalized_action_matrix)] <- "absorbing_state"
  normalized_action_matrix_with_absorbing_state <- normalized_action_matrix
  return(normalized_action_matrix_with_absorbing_state)
}
