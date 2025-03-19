#' This file run Reinforcement Learning (model-free Q-learning) to evaluate the most likely mutation paths from the data
#'
#' @param adj_list this is the adjacency link link to
#' @param discount the factor determines how much future mutations impact current actions. 0 only looks at next step, close to 1 will have a very long horizon time
#' @param N Number of iterations, also called episodes to explore the MDP
#' @export
#' @return A dataframe that contains the MDP adjacency list, normalized Q-values and mutations taken
mdp_Q_learning_with_linklist_kp<-function(adj_list, discount, N) {
  #N<- 10000
  #discount<-0.7
  # initialization of optional argument
  if (nargs() < 3) {
    N <- 5000
  }
  # total list of potential states
  total_state_list <-adj_list$current_state
  #sample initial random state
  state <- sample(total_state_list,1,replace=T)
  next_state <- state
  # These are measures used to grade each episode in the RL.
  mean_discrepancy =NULL
  discrepancy = NULL
  
  # This code consists of 2 loops. The outer loop is the iterations(or episodes)
  # the system is currently on. The
  # max number is empirically chosen by the user. Future releases will have an 
  # automated hueristic option based on the number of states/mutuations selected.
  
  
  # The next while loop is for achieving the goal state through an action. 
  # Note:The next_step MUST be included into the conditional statement of the 
  # while loop or else it will be terminated immediately upon randomly selecting
  # the last state! The next step is to find a random VALID action to another state. 
  # This was accomplished by using a while loop were a random sample is only 
  # valid if the selected state is in the list. The while loop ensures it will keep
  # picking an action until it picks a valid one. Once this action to another
  # state is selected the value for the Q_matrix can be updated.Then the action to the new
  # state becomes the current state.
  for(n in 1:N){
    if(n %% 5000 == 0){
      print(paste0("Finished Iteration: ",n))
    }
    # Reinitialisation of trajectories for a few different reasons:
    # 1) every terminal node hit (when we break while loop below)
    if(n >1){
      state <- sample(total_state_list,1,replace=T)
      next_state<-state
    }
    # epsilon = 1 means total exploration. equal prob of exploring each mutation. 
    
    epsilon <- 0.99  # decrease to focus more on exploitation. 
    # explore until we reach terminal state (last one in  total_state_list)
    # and we have already stayed at this location for 1 iteration.(needed so we build up and explore)
    # there may be another terminating step in the future, but need to think on this.
    
    while((state!=total_state_list[length(total_state_list)])&(next_state!=total_state_list[length(total_state_list)])){   
      # get list of potential actions to take    
      Action_list <-adj_list$next_state[adj_list$current_state==state]
      # get random action possible action (this currently assumes equal probability to take actions)
      if (runif(1) < epsilon) {
        #print("exploration")
        # Exploration: choose a random action
        next_state <- sample(Action_list, 1, replace = TRUE)
      } else {
        #print(paste("exploitation."))
        # Exploitation: choose the action with the highest Q-value
        Q_values_for_actions <- adj_list$Q_values[adj_list$current_state == state]
        best_action_index <- which.max(Q_values_for_actions)
        next_state <- Action_list[best_action_index]
      }
      # get reward with associated state-action transition
      r <- adj_list$reward[(adj_list$current_state==state)&(adj_list$next_state==next_state)] + 0.1 # ******* Added 0.1 --kp. # Add a small "bonus reward" to rare transitions
      # Updating the value of Q   
      # Decaying update coefficient (1/sqrt(n+2)) can be changed
      # Greedy search right now. 
      delta <- r + discount*max(adj_list$Q_values[adj_list$current_state==next_state]) - adj_list$Q_values[adj_list$current_state==state&adj_list$next_state==next_state]
      dQ <- (1/sqrt(n+2))*delta
      adj_list$Q_values[adj_list$current_state==state&adj_list$next_state==next_state] <- adj_list$Q_values[adj_list$current_state==state&adj_list$next_state==next_state] + dQ
      # Current state is updated
      state <- next_state
    }
  }
  # compute the normalized Q_values to treat them as percentages
  adj_list$Q_values_normalized <- 100*(adj_list$Q_values/max(adj_list$Q_values))
  
  # getting Value function and optimal policies are in a different function this way we don't need to retrain the model
  
  return(adj_list)
  
}