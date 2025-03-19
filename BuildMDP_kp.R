
# to build adjacency list. 
# only thing dependent on your data is num_mutations. 
# this is just set up. Calculating posibilities . 
BuildMDP<-function(num_mutations,use_ADO=FALSE){

  getidx <-function(n) {ifelse( (n==0|n==1),TRUE,FALSE)}
  dsum <- function(n) {ifelse(n < 10, n, n %% 10 + dsum(floor(n / 10)))}
  getTernaryFromState<-function(x,base=3){
    ifelse(x<base,x,getTernaryFromState(x%/% base)*10+x%%base)
  }
  backwardCheck<-function(x){(ifelse(dsum((-1)*x)==1,TRUE,FALSE))
  }

  
  if(num_mutations>=9){
    load(system.file(paste0('data/AdjacencyList_',num_mutations,'mutations.rDa'), package = 'scDNA'))
    return(adj_list)
  }
  num_type = 3
  num_states =num_type^num_mutations
  # this builds our masking function, we loop over upper triangle.
  # check if the summed digits are different by 1, if so it is an ok state we accept
  # EXAMPLE, check if StateA subtracted by stateC and stateB
  # stateC 121, stateB 110
  # stateA 010
  # subtract = 111 , 100
  # now sum digits 1+1+1 = 3, so reject, vs. 1+0+0 = 1 so accept!
  temp_vals <-(data.frame(x=(0:(num_states-1)))%>%apply(.,1,getTernaryFromState))
  # 0   1   2  10  11  12  20  21  22 100 101 102 110 111 112 120 121 122 200 201 202 210 211 212
  #220 221 222
  #length(temp_vals)  is 27 bc number of unique elements 3^3 = 27
  
  
  j<-NULL
  i<-NULL
  state_tern = temp_vals
  
  f.compfun3<-function(state_tern,num_muts){
    f<-compiler::cmpfun(function(state_tern=temp_vals,num_muts=num_mutations){
      j<-NULL
      i<-NULL
      state_tern = temp_vals
      #print(length(state_tern))
      r<-foreach::foreach(state_to_check=1:length(state_tern), .combine='cbind', .multicombine=TRUE,.export=c("dsum","getidx")) %dopar%{
        temp<-matrix(NA_real_,nrow=2)
        state_val <-(dsum(state_tern-state_tern[state_to_check]))
        vec<-which(getidx(state_val))
        ADO <-which(ifelse(dsum((-1)*state_val)==1,TRUE,FALSE))
        check1 <- unlist(lapply(ADO,
                                function(x) ifelse((intToUtf8(
                                  utf8ToInt(
                                    stringr::str_pad(state_tern[state_to_check], 
                                                     num_muts, 
                                                     pad = "0"))[(which((utf8ToInt(stringr::str_pad(state_tern[x], 
                                                                                                    num_muts, 
                                                                                                    pad = "0"))-utf8ToInt(stringr::str_pad(state_tern[state_to_check], 
                                                                                                                                           num_muts, 
                                                                                                                                           pad = "0"))) !=0))]))=="1",x,NA)))
        forwardADO <- unlist(lapply(vec,
                                    function(x) ifelse((intToUtf8(
                                      utf8ToInt(
                                        stringr::str_pad(state_tern[state_to_check], 
                                                         num_muts, 
                                                         pad = "0"))[(which((utf8ToInt(stringr::str_pad(state_tern[x], 
                                                                                                        num_muts, 
                                                                                                        pad = "0"))-utf8ToInt(stringr::str_pad(state_tern[state_to_check], 
                                                                                                                                               num_muts, 
                                                                                                                                               pad = "0"))) !=0))]))=="1",x,NA)))
        
        
        vec<-append(vec,check1[!is.na(check1)])
        vec<-append(vec,forwardADO[!is.na(forwardADO)])
        #print(state_to_check)
        j<-t(vec)
        i<-t(rep(state_to_check,length(vec)))
        temp <-rbind(i,j)
        return(temp)
      }
      r
    })
    f(f(0))
  }
  
  f.compfun2<-function(state_tern){
    f<-compiler::cmpfun(function(state_tern=temp_vals){
      j<-NULL
      i<-NULL
      state_tern = temp_vals
      #print(length(state_tern))
      r<-foreach::foreach(state_to_check=1:length(state_tern), .combine='cbind', .multicombine=TRUE,.export=c("dsum","getidx")) %dopar%{
        temp<-matrix(NA_real_,nrow=2)
        vec<-which(getidx(dsum(state_tern-state_tern[state_to_check])))
        #print(state_to_check)
        j<-t(vec)
        i<-t(rep(state_to_check,length(vec)))
        temp <-rbind(i,j)
        return(temp)
      }
      r
    })
    f(f(0))
  }
  library(doParallel)
  # some parallelization clusters, and records time.
  if(use_ADO){
    cl <- parallel::makeCluster(7)
    doParallel::registerDoParallel(cl)
    #start<-Sys.time()
    output_mat <-f.compfun3(temp_vals,num_mutations)
    #print( Sys.time() - start )
    parallel::stopCluster(cl)
  }else{
    cl <- parallel::makeCluster(7)
    doParallel::registerDoParallel(cl)
    #start<-Sys.time()
    output_mat <-f.compfun2(temp_vals)
    #print( Sys.time() - start )
    parallel::stopCluster(cl)
  }
  
  # The first row (output_mat[1, ]): Contains the indices (or IDs) of the current states.
  #The second row (output_mat[2, ]): Contains the indices (or IDs) of the next states to which transitions occur.

  #output_mat = each state is like 102 or 010 or 221 or something i think. the number are numbers of legal transitions? 
  #[,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17]
  #[1,]    1    1    1    1    2    2    2    2    3     3     3     4     4     4     4     5     5
  #[2,]    1    2    4   10    2    3    5   11    3     6    12     4     5     7    13     5     6
  
  # Column 1: Transition from state 1 to state 1.
  # Column 2: Transition from state 1 to state 2.
  # Column 3: Transition from state 1 to state 4.
  # Column 4: Transition from state 1 to state 10.
  # Column 5: Transition from state 2 to state 2.
  # Column 6: Transition from state 2 to state 3.
  # Column 7: Transition from state 2 to state 5.
  # Column 8: Transition from state 2 to state 11.
  # And so on.
  
  # Now we have a 2 row matrix we want to make a sparse representation
  reward_mat<-Matrix::sparseMatrix(output_mat[1,],output_mat[2,],,1)
  
  reward_mat
  rownames(reward_mat) <-temp_vals
  colnames(reward_mat) <-temp_vals
  # this next line is basically a pivot longer but ditches all 0s
  adj_list<-mefa4::Melt(reward_mat)
  colnames(adj_list) <-c("current_state","next_state","legal_action")
  # this unsets the duplicate rows so we can set forward ADO and mutation
  library(tidyr)
  library(dplyr)
  adj_list<-adj_list %>%
    tidyr::uncount(legal_action)%>%
    dplyr::mutate(legal_action=1)
  
  # labeling forward ADO, backward ADO, mutation, or none
  adj_list$action_type <-"mutation"
  adj_list$action_type <-ifelse(dsum(as.numeric(adj_list$next_state)-as.numeric(adj_list$current_state))<0,"ADO","mutation")
  adj_list$action_type <-ifelse(dsum(as.numeric(adj_list$next_state)-as.numeric(adj_list$current_state))==0,"none",adj_list$action_type)
  adj_list<-adj_list %>%
    dplyr::group_by(current_state,next_state,legal_action,action_type)%>%
    dplyr::mutate(rank = rank(action_type,ties.method="first"))%>%
    dplyr::mutate(action_type=ifelse(rank==2,"forward_ADO",action_type))#%>%


  # For now we will say equal transition probability
  adj_list%>%
    dplyr::group_by(current_state)%>% # that subsequent operations will be performed within each group of rows that share the same current_state.
    dplyr::mutate(total_actions=sum(legal_action))%>% #  calculates how many total legal actions are possible from each current_state.
    dplyr::ungroup()%>%
    dplyr::mutate(Probability_transition=legal_action/total_actions)->adj_list
  
  #assign Q-values to be 0 for each one right now, these get updated in the mdp_Q_learning_with_linklist.R function
  adj_list$Q_values<-rep(0,length(adj_list$current_state))


  return(adj_list)
}
