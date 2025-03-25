
make_rewards3<- function(adj_list){
  # reward obtained when action a is executed in state s. 
  rewards<- adj_list %>% dplyr::select(c("current_state", "action_type", "next_state", "reward")) %>%
    arrange(current_state);rewards
  
  rewards2 <- as.data.frame(rewards) %>% 
    dplyr::rename(start.state=current_state) %>% 
    dplyr::rename(action=action_type) %>% 
    dplyr::rename(value=reward) %>% 
    filter(action !="none") ;head(rewards2)

  rewards3_fixed <- rewards2 %>%
    mutate(
      end.state = "*",
      observation = "*"
    ) %>%
    dplyr::select(action, start.state, end.state, observation, value) ; head(rewards3_fixed)
  
  # Doing nothing forever is forbidden
  new_row1 <- data.frame( action = "none", start.state = "*", end.state = "*",observation = "*",    
    value = -Inf     
  )
  rewards3_fixed <- rbind( new_row1,rewards3_fixed)

  new_row3 <- data.frame( #If the game is already done (in absorbing_state)
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


  last_possible_clone<- max(unique(as.character(adj_list$current_state)))
  # no purpose in listening any further in last possible state.
  new_row7 <- data.frame( 
    action = "listen",
    start.state =as.integer(last_possible_clone),
    end.state = "*",
    observation = "*",    
    value = -Inf
  )
  rewards3_fixed <- rbind(new_row7, rewards3_fixed)
  
  new_row6 <- data.frame( 
    action = "*",
    start.state = "*",
    end.state = as.integer(last_possible_clone),
    observation = "*",    
    value = 10
  )
  rewards3_fixed <- rbind(new_row6, rewards3_fixed) 
  
  new_row8 <- data.frame(  
    action = "*",
    start.state = "*",
    end.state = "*",
    observation = "absorbing_state",    
    value = -Inf
  )
  rewards3_fixed <- rbind(new_row8, rewards3_fixed) 
  
  return(rewards3_fixed)
}

