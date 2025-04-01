#' the agent and environment interact at each of a sequence of discrete time steps, t = 0, 1, 2, 3,...
#' At each time step t, the agent receives some representation of the environment’s state, and on that basis selects
#' an action.  One time step later, in part as a consequence of its action, the agent receives a numerical reward, Rt+1
make_rewards3<- function(adj_list){
  states_1<- unique(as.character(adj_list$current_state))
  actions_1<- unique((adj_list$action_type))
  # every_combo_state_action <- expand.grid(states_1 ,actions_1, stringsAsFactors = FALSE)
  # every_combo_state_action_formatted<- every_combo_state_action %>%
  #   dplyr::rename(start.state=Var1) %>% 
  #   dplyr::rename(action=Var2) %>% 
  #   mutate(
  #     end.state = "*",
  #     observation = "*",
  #     value = -Inf
  #   ) 
  # pair action and obs. * for current state
  
  # reward obtained when action a is executed in state s. 
  rewards<- adj_list %>% dplyr::select(c("current_state", "action_type", "next_state", "reward")) %>%
    arrange(current_state) %>% 
    dplyr::rename(start.state=current_state) %>% 
    dplyr::rename(end.state=next_state) %>% 
    mutate(
      end.state = "*", # overwrite next state
      observation = "*"
    ) %>%
    dplyr::rename(action=action_type) %>% 
    dplyr::rename(value=reward) %>% 
    filter(action !="none")
  rewards
  
  rewards$start.state<- as.character(rewards$start.state, stringsAsFactors = FALSE)
  rewards$end.state<- as.character(rewards$end.state, stringsAsFactors = FALSE)
  rewards$action<- as.character(rewards$action, stringsAsFactors = FALSE)
  
  rewards_for_mutations <- rewards
  rewards3_fixed <- rewards[0, ]
  
  new_row4 <- data.frame( action = "listen",start.state = "*",end.state = "*", observation = "*",value = -0.5)
  rewards3_fixed <- rbind(rewards3_fixed, new_row4)
  
  rewards3_fixed = rbind(rewards3_fixed,
                         R_(action = "none", start.state = "*", end.state = "*", observation = "*",v = -Inf))

  new_row3 <- data.frame( #If the game is already done (in absorbing_state)
    action = "*", # we use absorbing state to represent illegal transitions 
    start.state = "absorbing_state",
    end.state = "*", observation = "*",value = -Inf)
  rewards3_fixed <- rbind(rewards3_fixed, new_row3)
  
  
  new_row5 <- data.frame( 
    action = "none",
    start.state = "absorbing_state",
    end.state = "*",
    observation = "*",    
    value = 0
  )
  rewards3_fixed <- rbind(rewards3_fixed, new_row5)

  
  rewards3_fixed<- rbind(rewards3_fixed, rewards_for_mutations)
  return(rewards3_fixed)
}


make_rewards4<- function(adj_list){
  rewards<- adj_list %>% dplyr::select(c("current_state", "action_type", "next_state", "reward")) %>%
    arrange(current_state);rewards
  rewards2 <- as.data.frame(rewards) %>% 
    dplyr::rename(start.state=current_state) %>% 
    dplyr::rename(end.state=next_state) %>% 
    dplyr::rename(action=action_type) %>% 
    dplyr::rename(value=reward) %>% 
    filter(action !="none") ;head(rewards2)
  
  
  states_1<- unique(as.character(adj_list$current_state))
  actions_1<- unique((adj_list$action_type))
  every_combo_state_action <- expand.grid(states_1 ,actions_1, stringsAsFactors = FALSE)
  every_combo_state_action_formatted<- every_combo_state_action %>%
    dplyr::rename(start.state=Var1) %>%
    dplyr::rename(action=Var2) %>%
    mutate(
     # end.state = "*",
      observation = "*"#,
      #value = -Inf
    )
  every_combo_state_action_formatted <- every_combo_state_action_formatted%>% filter(action != "none")
  rew<- left_join(every_combo_state_action_formatted,rewards2, by = c("start.state", "action"))
  rew
  rew$end.state<- as.character(rew$end.state, stringsAsFactors = FALSE)
  rew$action<- as.character(rew$action, stringsAsFactors = FALSE)
  rew <- rew %>% mutate(across(everything(), ~ifelse(is.na(.), "*", .)))
  rew$value[rew$value == "*"] <- 0
  
  # Doing nothing forever is forbidden
  new_row1<-data.frame(action = "none",start.state = "*",end.state = "*",observation = "*",value = -Inf)
  rew <- rbind( new_row1,rew)
  
  
  new_row3 <- data.frame( #If the game is already done (in absorbing_state)
    action = "*", # we use absorbing state to represent illegal transitions 
    start.state = "absorbing_state",
    end.state = "*",
    observation = "*",    
    value = -Inf  
  )
  rew <- rbind(new_row3, rew)
  
  new_row4 <- data.frame( 
    action = "listen",
    start.state = "*",
    end.state = "*",
    observation = "*",    
    value = -10
  )
  rew <- rbind(new_row4, rew)
  
  new_row5 <- data.frame( 
    action = "none",
    start.state = "absorbing_state",
    end.state = "*",
    observation = "*",    
    value = 0
  )
  rew <- rbind(new_row5, rew)
  rew
  
  return(rew)
}

