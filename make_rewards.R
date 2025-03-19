# prep_data_frame_clarity<- function(){
#   data_frame_clarity<-dplyr::inner_join(sce@metadata$Architecture,
#                                         as.data.frame(sce@metadata$Clones)%>%
#                                           dplyr::mutate(Count=ifelse(Group=="Other",n_Other,n_Complete))%>%
#                                           dplyr::filter(Group=="Complete")%>%
#                                           dplyr::select(Clone,Count)%>%dplyr::distinct(),
#                                         by="Clone")%>% dplyr::group_by(Clone)
#   weights_unique_list <-dplyr::distinct(data.frame(reward=data_frame_clarity$Count/sum(data_frame_clarity$Count)*100,
#                                                    Clone = data_frame_clarity$Clone ))
#   
#   
# }


make_reward_df2<- function(states, actions, reformatted_adj_list2){
  actions <- factor(actions)
  states <- factor(states)
  reward_matrix <- data.frame(
    action = rep(actions, each = length(states)),   # Repeat each action for all states
    start.state = as.character(rep(states, times = length(actions))),  # Assign each state per action
    end.state = "*",   
    observation = "*",   
    stringsAsFactors = FALSE  # Ensure character columns remain character
  )
  
  state_action_combinations<- as.data.frame(reformatted_adj_list2) %>% 
    dplyr:: select(c("current_state","reward", "formatted_action"))%>%
    dplyr::rename(start.state= current_state) %>%
    dplyr::rename(value= reward) %>% 
    dplyr::rename(action= formatted_action)
  res<- left_join(reward_matrix, state_action_combinations, by = c("action", "start.state"))
  res$value[is.na(res$value)]<- -1
  print(res)
  return(res)
}
