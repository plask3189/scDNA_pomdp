
# get potential dropout clones
make_potential_dropout_clones<- function(states){
  states<-states[!is.na(as.numeric(states))]
  max_length <- max(nchar(as.character(states)))
  states_padded <- sprintf(paste0("%0", max_length, "d"), as.numeric(states))
  legal <- data.frame(
    initial_zyg = c(0, 1, 2),
    potential_zygs = I(list(c(0), c(1, 0), c(2, 1)))
  )
  clones_heard_as <- data.frame(
    og_mutation = states_padded,
    could_be_heard_as = I(as.list(states_padded))
  )
  
  for (mut_zyg in states_padded){ # loop through each state
    cat("----------- \n")
    cat("Clone:", mut_zyg, "\n")
    digits <- unlist(strsplit(mut_zyg, split = ""))
    for (i in seq_along(digits)){# loop through each mutation (digit) in a state 
      cat("  Digit:", digits[i], "\n")# get potential zygs for that.
      potential_zygs <- legal %>% filter(initial_zyg == digits[i]) %>% pull(potential_zygs) %>% unlist()
      for (potential_zyg in potential_zygs){
        cat("     Potential zyg", potential_zyg, "\n")
        if (potential_zyg != digits[i]){
          new_digits<- digits
          # replace with new zygosity.
          cat("digit", new_digits[i], "to be converted to", potential_zyg, "\n")
          new_digits[i]<- potential_zyg
          new_clone<- paste(new_digits, collapse ="")
          cat(green("Clone", mut_zyg, "could be heard as",  new_clone, "\n"))
          # append to the could_be_heard_as column of appropriate og_mutation
          # Find the row where og_mutation equals the current mut_zyg
          row_idx <- which(clones_heard_as$og_mutation == mut_zyg)
          # Append new_clone to that list element
          clones_heard_as$could_be_heard_as[[row_idx]] <- c(clones_heard_as$could_be_heard_as[[row_idx]], new_clone)
          
        }
      }
    }
  }
  clones_heard_as$og_mutation <- as.numeric(as.character(clones_heard_as$og_mutation))
  clones_heard_as$could_be_heard_as <- lapply(clones_heard_as$could_be_heard_as, as.numeric)
  print(clones_heard_as)
  return(clones_heard_as)
}
