reformatted_actions_adj_list <- function(adj_list_formatted) {

  # Initialize formatted_action column
  adj_list_formatted$formatted_action <- NA_character_
  
  formatted_actions <- character()  # Store formatted actions
  
  for (i in 1:nrow(adj_list_formatted)){
    mutation <- adj_list_formatted$action_type[i]
    if (mutation != "none"){
      index_of_variant_in_clone <- adj_list_formatted$digit_position_changed[i]
      current_state_zygosity <- substr(adj_list_formatted$current_state[i], index_of_variant_in_clone, index_of_variant_in_clone)
      next_state_zygosity <- substr(adj_list_formatted$next_state[i], index_of_variant_in_clone, index_of_variant_in_clone)
      formatted_action <- paste0(mutation, "_", current_state_zygosity, "-->", next_state_zygosity)
    }
    else {
      formatted_action<- "none"
    }
    adj_list_formatted$formatted_action[i] <- formatted_action
  }
  
  return(adj_list_formatted)
}
