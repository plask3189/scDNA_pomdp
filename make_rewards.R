

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
