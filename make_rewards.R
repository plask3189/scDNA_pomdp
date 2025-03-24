
make_rewards3<- function(adj_list){
  # every combo of state and action
  # states <- unique(as.character(adj_list$current_state))
  # actions<- unique(adj_list$action_type)
  # every_combo <- expand.grid(current_state = states, action_type = actions) %>% arrange(current_state)

  # reward obtained when action a is executed in state s. 
  rewards<- adj_list %>% dplyr::select(c("current_state", "action_type", "next_state", "reward")) %>%
    arrange(current_state)
  
  # rewards<- left_join(every_combo,rewards, by= c("current_state", "action_type"))
  # rewards$reward[is.na(rewards$reward)] <- -Inf # illegal. 
  
  rewards2 <- as.data.frame(rewards) %>% 
    dplyr::rename(start.state=current_state) %>% 
    #dplyr::rename(observation=next_state) %>%
    dplyr::rename(action=action_type) %>% 
    dplyr::rename(value=reward) 
  rewards2 <- rewards2 %>% # replace 0s with -1 to really discourage staying in sample place.
    mutate(value = ifelse(action == "none", -1, value))
  # absorbing_state_rewards <- merge(actions, absorbing_state, by = NULL) %>% dplyr:: rename(action = x) %>% 
  #   dplyr:: rename(start.state = y) 
  rewards2 <- rewards2 %>% filter(action !="none")
  # absorbing_state_rewards$value <- -Inf
  # rewards3<- rbind(rewards2, absorbing_state_rewards)
  rewards3<-rewards2
  rewards3_fixed <- rewards3 %>%
    mutate(
      end.state = "*",
      observation = "*"
    ) %>%
    dplyr::select(action, start.state, end.state, observation, value) 
  
  new_row1 <- data.frame( # Doing nothing forever is forbidden
    action = "none",
    start.state = "*",
    end.state = "*",
    observation = "*",    
    value = -Inf     
  )
  rewards3_fixed <- rbind( new_row1,rewards3_fixed)

  new_row3 <- data.frame( #If the game is already done (in absorbing_state), taking no action has no cost
    action = "*", # we use absorbing state to represent illegal transitions 
    start.state = "absorbing_state",
    end.state = "*",
    observation = "*",    
    value = -Inf  
  )
  rewards3_fixed <- rbind(new_row3, rewards3_fixed)
  
  new_row4 <- data.frame( 
    action = "listen",
    start.state = "*",
    end.state = "*",
    observation = "*",    
    value = -1
  )
  rewards3_fixed <- rbind(new_row4, rewards3_fixed)
  
  new_row5 <- data.frame( 
    action = "none",
    start.state = "absorbing_state",
    end.state = "*",
    observation = "*",    
    value = 0
  )
  rewards3_fixed <- rbind(new_row5, rewards3_fixed)
  
  
  
  # actions_set <- setdiff(unique(adj_list[, "action_type"]) %>% pull(), "none")
  # for (action in actions_set){
  #   new_row <- data.frame( 
  #     action = action,
  #     start.state = "*",
  #     end.state = "*",
  #     observation = "*",    
  #     value = -100
  #   )
  #   rewards3_fixed <- rbind(new_row, rewards3_fixed)
  # }
  

  return(rewards3_fixed)
}

# generate_illegal_state_transitions<- function(adj_list){
#   legal_transiitons<- adj_list %>% dplyr::select(c("current_state", "next_state"))
#   all_states<- unique(as.character(legal_transiitons$current_state))
#   all_possible_combos <- expand.grid(current_state = all_states, next_state = all_states)
#   illegal_transitions<- setdiff(all_possible_combos, legal_transiitons)
#   illegal_transitions<- illegal_transitions %>% dplyr::rename(start.state= current_state, observation = next_state )
#   unique_illegal_transitions <- unique(illegal_transitions)
#   return(illegal_transitions)
# }

make_reward_df2<- function(states, actions, reformatted_adj_list2){
  actions <- factor(actions); actions
  states <- factor(states); states
  reward_matrix <- data.frame(
    action = rep(actions, each = length(states)),  
    start.state = as.character(rep(states, times = length(actions))),
    end.state = "*",   
    observation = "*",   
    stringsAsFactors = FALSE 
  )
  
  state_action_combinations<- as.data.frame(reformatted_adj_list2) %>% 
    dplyr:: select(c("next_state","reward", "formatted_action"))
    
  state_action_combinations<- state_action_combinations %>% 
    dplyr::rename(start.state= next_state) %>%
    dplyr::rename(value= reward) %>% 
    dplyr::rename(action= formatted_action) %>% 
    arrange(start.state)
  res<- left_join(reward_matrix, state_action_combinations, by = c("action", "start.state"))
  res$value[is.na(res$value)]<- -Inf
  print(res)
  
  res <- res %>%
    mutate(is_legal = mapply(check_if_action_is_legal_for_state_and_change_happens, 
                             mutation = action, 
                             clone_state = start.state, 
                             MoreArgs = list(is_start_or_next = 'next')))
  res
  for (x in 1:nrow(res)) {
    if (res$is_legal[x] == TRUE & res$value[x] == -Inf) {
      res$value[x] <- 0  
    } 
    # no change gets small negative reward. 
    if ((res$is_legal[x]) == "No change") {
      res$value[x] <- -1 
    } 
  }
  res
  return(res)
}


#check_if_action_is_legal_for_state("NRAS.G12R_1-->2", "10", "next")
check_if_action_is_legal_for_state_and_change_happens<- function(mutation, clone_state, is_start_or_next = 'start'){
  mutations_starting_zygosity<-  sub(".*_(\\d+)-->.*", "\\1", mutation)
  mutations_resulting_zygosity<- sub(".*-->\\s*(\\d+).*", "\\1", mutation)
  cat(blue("Checking if",mutation, "is possible if Clone", clone_state, "is",is_start_or_next, "clone and the mutation follows.\n"))
  variant_name<- sub("_.*", "", mutation)
  index_of_variant_in_clone<- position_of_variants_in_clone%>% filter(action_type == variant_name) %>% pull(digit_position_changed) 
  state_zygosity <- substr(clone_state, index_of_variant_in_clone, index_of_variant_in_clone)
  if(is_start_or_next == 'start'){
    if(state_zygosity == mutations_starting_zygosity){
      cat(green(mutation, "is possible if Clone", clone_state, "is",is_start_or_next, "clone and the mutation follows.\n"))
      return(TRUE)
    } else {
      cat(red(mutation, "is not possible if Clone", clone_state, "is",is_start_or_next, "clone and the mutation follows.\n"))
      return(FALSE)}
  }
  if(is_start_or_next == 'next'){
    # if no change from starting zygosity. 
    if (mutations_starting_zygosity == mutations_resulting_zygosity & state_zygosity == mutations_resulting_zygosity){
      return("No change")
    } else {
    if(state_zygosity == mutations_resulting_zygosity){
      cat(green(mutation, "is possible if Clone", clone_state, "is",is_start_or_next, "clone and the mutation follows.\n"))
      return(TRUE)
      
    } else{ 
      cat(red(mutation, "is not possible if Clone", clone_state, "is",is_start_or_next, "clone and the mutation follows.\n"))
      return(FALSE)
    }
    }
  }
}
