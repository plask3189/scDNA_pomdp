#' Assigns observed counts as weights to the Markov Decision Process
#'
#' @param sce The single cell experiment object that carries all data
#' @param adj_linklist this is an adjacency list for the theoretical states describing the MDP
#' @export
attach_weights_kp <-function(sce,adj_linklist){
  # first we obtain weights from our final_summary dataframe (this was taken from the myeloid one)
  # TODO: Change the filtering, not sure if we always want "complete" or "other" or both?

  data_frame_clarity<-dplyr::inner_join(sce@metadata$Architecture,
                                       as.data.frame(sce@metadata$Clones)%>%
                                         dplyr::mutate(Count=ifelse(Group=="Other",n_Other,n_Complete))%>%
                                         dplyr::filter(Group=="Complete")%>%
                                                      dplyr::select(Clone,Count)%>%dplyr::distinct(),
                                      by="Clone")%>% 
    dplyr::group_by(Clone)

  check_sub_var<-data.frame(current_state=gsub('^.|.$', '',gsub("", "_",stringr::str_pad(adj_linklist$current_state, length(unique(sce@metadata$Architecture$final_annot)), pad = "0"))),
                            next_state=gsub('^.|.$', '',gsub("", "_",stringr::str_pad(adj_linklist$next_state, length(unique(sce@metadata$Architecture$final_annot)), pad = "0"))))
  check2<-array2DF(t(rbind(mapply(function(x,y) unique(sce@metadata$Architecture$final_annot)[(which(bitwXor(utf8ToInt(x),utf8ToInt(y))>0)+1)/2], 
         check_sub_var$current_state,check_sub_var$next_state,SIMPLIFY = FALSE))))
  check_sub_var$mutation_taken<-"none"
  check_sub_var$mutation_taken[as.numeric(rownames(check2))]<-check2$Value
  if (!is.null(sce@metadata$FalsePositive)) {
    renamed_false_positive <- sce@metadata$FalsePositive %>%
      dplyr::rename(mutation_taken = final_annot)
  } else {
    renamed_false_positive <- data.frame(mutation_taken = character(0)) # Create an empty placeholder
  }
  check_sub_var <- dplyr::full_join(
    check_sub_var, 
    renamed_false_positive, 
    by = "mutation_taken"
  )
  
  
  weights_unique_list <-dplyr::distinct(data.frame(reward=data_frame_clarity$Count/sum(data_frame_clarity$Count)*100,
                                            names=as.numeric(gsub("_", "",data_frame_clarity$Clone))))

  colnames(weights_unique_list)<-c("reward","next_state")
  weights_unique_list$next_state<-as.factor(weights_unique_list$next_state)
  
  # here, we do three things:
  #   1) attach the weights as a reward,
  #   2) use these as observed clonograph states for searching later
  #   3) assign staying in place as a reward of 0 
  adj_linklist <- dplyr::left_join(adj_linklist,weights_unique_list, by="next_state")%>%
    dplyr::mutate(
      dplyr::across(where(is.numeric), ~tidyr::replace_na(.x,0)))
  
  # for search later
  adj_linklist$observed_states <- ifelse(adj_linklist$reward>0,1,0)
  
    adj_linklist$reward[which(adj_linklist$rank==2)-1] <- (1-check_sub_var$false_positiveHom[which(adj_linklist$rank==2)-1])*adj_linklist$reward[which(adj_linklist$rank==2)-1]
  adj_linklist<-adj_linklist%>%
    dplyr::mutate(action_type=ifelse(action_type=="mutation",check_sub_var$mutation_taken,action_type))
  adj_linklist$reward[which(adj_linklist$action_type=="ADO")]<-(-1)*check_sub_var$false_positiveWT[which(adj_linklist$action_type=="ADO")]*ifelse((adj_linklist$reward[which(adj_linklist$action_type=="ADO")])>1e-4, (adj_linklist$reward[which(adj_linklist$action_type=="ADO")]),1e-4)
  adj_linklist$reward[which(adj_linklist$action_type=="forward_ADO")]<-adj_linklist$reward[which(adj_linklist$action_type=="forward_ADO")]*check_sub_var$false_positiveHom[which(adj_linklist$action_type=="forward_ADO")]
  
  # for ensuring we move forward through the MDP
  #If a state points back to itself (no change), the reward is set to 0. This ensures the MDP always progresses to new states.
  adj_linklist$reward[which(adj_linklist$current_state==adj_linklist$next_state)]=0

  adj_linklist<-adj_linklist%>%dplyr::select("current_state","action_type","next_state","reward","Q_values","observed_states","legal_action","total_actions","Probability_transition")
  # Now we can run RL in the next function
  return(adj_linklist)
}

# 
# weight_reward_by_kl<- function(adj_linklist, sce){
#   adj_linklist_formatted <- adj_linklist %>%
#     mutate(
#       current_state = format_state(current_state),
#       next_state = format_state(next_state)
#     )
# 
#   cells_and_all_their_variants_zygosities<- determine_clone_marker_zyg(sce)
#   #top_significant_markers_using_kruskal<- compute_kruskal(cells_and_all_their_variants_zygosities)
#   # ALL VARIANTS:
#   top_significant_markers_using_kruskal<- unique(cells_and_all_their_variants_zygosities$Variant)
#   clone_marker_zygosity <- cells_and_all_their_variants_zygosities %>% dplyr::filter(Variant %in% top_significant_markers_using_kruskal) %>% as.data.frame() %>% dplyr::rename(Zygosity = Value); cat(blue("Filtered the dataset with all cells and all variants to now just include the top significant variants for each cell. \n"))
#   ancestry_table<- make_ancestry_table(adj_linklist_formatted)
#   clone_lineages_lookup_table <- identify_lineages_for_clones(ancestry_table)
#   all_cell_lineages <- clone_marker_zygosity %>% left_join(clone_lineages_lookup_table, by = c("Clone"))%>% as.data.frame()
#   all_cell_lineages
#   #Determine what each tree clone would actually be if the tree clone is dropout. (generate table of hypothesis clones)
#   clone_dropout_table<- if_each_tree_clone_was_dropout_what_would_be_its_real_clone(unique(clone_lineages_lookup_table$Clone))
#   clone_info<- inner_join(clone_lineages_lookup_table, clone_dropout_table, by = "Clone")
#   clone_info_with_gene<- add_genes_to_clone_info(clone_info, adj_linklist_formatted)
#   clone_info_with_gene
#   wt_id<- ancestry_table$Clone[1]
#   
#   # just the unique cells 
#   cell_lins <- all_cell_lineages %>%
#     dplyr::filter(Clone != wt_id) %>%
#     group_by(Cell) %>%  
#     dplyr::slice(1) %>% 
#     ungroup()   %>%
#     dplyr::select(-c("Variant")) %>% as.data.frame()  %>% dplyr::select(-c("Zygosity"))
#   
#   zygosity_matrix<- all_cell_lineages %>% distinct() %>%
#     tidyr::pivot_wider(names_from = Variant, values_from = Zygosity, values_fill = 0)
#   
#   kl_for_each_Clone <- data.frame( clone = character(), pwm = list(), stringsAsFactors = FALSE)
#   unique_clones <- unique(cell_lins$Clone)
#   for (clone in unique_clones){
#     individual_potential_clones <- identify_alternative_lineages(clone, clone_info)
#     pwm_result<- kl(clone, zygosity_matrix, clone_info_with_gene, individual_potential_clones)
#     kl_for_each_Clone <- rbind(kl_for_each_Clone, data.frame(clone = clone, pwm = pwm_result))
#   }
#   
#   # NORMALIZE WITHIN THE CLONE RATHER THAN NORMALIZE ACROSS CLONES?
#   kl_for_each_Clone <- kl_for_each_Clone %>%
#     group_by(clone) %>%
#     mutate(Normalized_JS = (pwm.JS- min(pwm.JS)) /
#              (max(pwm.JS) - min(pwm.JS))) %>% ungroup()
#   vis_kl(kl_for_each_Clone)
#   adj_linklist_weighted <- adj_linklist_formatted %>%
#     left_join(kl_for_each_Clone, by = c("current_state" = "clone", "next_state" = "pwm.Potential_New_Lineage")) %>%
#     mutate(New_Reward = reward.x * (1 - Normalized_JS)) %>%
#     dplyr::select(current_state, action_type, next_state, reward.x, Normalized_JS, New_Reward)
#   
# 
# }
# 
# 
# 
# #  to format state values like in data_frame_clarity
# format_state <- function(state) {
#   state <- as.character(state)  # Convert factors to characters first
#   state <- as.numeric(state)  # Convert characters to numeric
#   state <- sprintf("%03d", state)  # Ensure 3-digit format
#   state <- gsub("(.)", "\\1_", state)  # Insert "_"
#   formatted_state <- substr(state, 1, nchar(state) - 1)  # Remove last "_"
#   return(formatted_state)
# }
# 
# 
# vis_kl<- function(kl_for_each_Clone){
#   ggplot(kl_for_each_Clone, aes(x = clone, y = pwm.KL_Divergence, color = pwm.Potential_New_Lineage)) +
#     geom_boxplot(outlier.shape = NA) +  # Avoid overlap of outliers with jittered points
#     geom_jitter(width = 0.2, alpha = 0.7) +  # Add jittered points with some transparency
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     labs(title = "KL Divergence Distribution by Clone", y = "KL Divergence", color = "Potential New lineage") +
#     theme(legend.position = "right")
# 
#   which_metric<- "Normalized_JS"
#   heatmap_data <- kl_for_each_Clone %>%
#     dplyr::select(clone, pwm.Potential_New_Lineage, which_metric) %>%
#     pivot_wider(names_from = pwm.Potential_New_Lineage, values_from = which_metric, values_fill = 0)
#   
#   heatmap_long <- heatmap_data %>%
#     pivot_longer(cols = -clone, names_to = "Potential_New_Lineage", values_to = which_metric)
#   heatmap_long$which_metric <- as.numeric(heatmap_long[[which_metric]])
#   ggplot(heatmap_long, aes(x = clone, y = Potential_New_Lineage, fill = which_metric)) +
#     geom_tile() +
#     scale_fill_gradient(low = "white", high = "red") +
#     labs(title = "Heatmap of Divergence Between Clones",
#          x = "Clone",
#          y = "Potential New Lineage",
#          fill = which_metric) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
# }
# 
