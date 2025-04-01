
transition_graph_kp<- function (x, action = NULL, episode = NULL, epoch = NULL, state_col = NULL, 
          simplify_transitions = TRUE, remove_unavailable_actions = TRUE) 
{
  x<- sc
  action = NULL
  episode = NULL
  epoch = NULL
  state_col = NULL
  simplify_transitions = TRUE
  remove_unavailable_actions = TRUE
  state_col <- colors_discrete(length(x$states), state_col)
  m <- transition_matrix(x, action = NULL, episode = episode, 
                         epoch = epoch, sparse = FALSE)
  if (is.null(action)) 
    action <- x$actions
  gs <- sapply(action, FUN = function(a) {
    cat("action:", a, "\n")
    g <- graph_from_adjacency_matrix(m[[a]], mode = "directed", 
                                     weighted = TRUE)
    E(g)$label <- a
    df <- as_data_frame(g)
    cat("df:\n"); print(df)
    if (remove_unavailable_actions) {
      available <- sapply(seq_len(nrow(df)), FUN = function(i) all(reward_val(x, 
                                                                              action = df$label[i], start.state = df$from[i], 
                                                                              end.state = df$to[i]) != -Inf))
      df <- df[available, , drop = FALSE]
    }
    df
  }, simplify = FALSE)
  g <- graph_from_data_frame(do.call(rbind, gs))
  g <- igraph::permute(g, match(V(g)$name, x$states))
  E(g)$label <- paste0(E(g)$label, ifelse(E(g)$weight != 1, 
                                          paste0(" (", round(E(g)$weight, 2), ")"), ""))
  if (simplify_transitions) 
    g <- igraph::simplify(g, edge.attr.comb = list(label = function(x) paste(x, 
                                                                             collapse = "/\n"), "ignore"), remove.loops = FALSE)
  if (!any(is.na(state_col))) 
    V(g)$color <- state_col
  g
}


.get_alpha_file<- function (file_prefix, model, number = "") 
{
  alpha <- pomdpSolve::read_alpha_file(paste0(file_prefix, 
                                              "-0.alpha", number))
  colnames(alpha) <- model$states
  alpha
}
.get_pg_file<- function (file_prefix, model, number = "") 
{
  filename <- paste0(file_prefix, "-0.pg", number)
  pg <- pomdpSolve::read_pg_file(filename)
  colnames(pg) <- c("node", "action", as.character(model$observations))
  pg[, 2] <- factor(pg[, 2], levels = seq(length(model$actions)), 
                    labels = model$actions)
  pg
}
.get_belief_file<- function (file_prefix, model) 
{
  belief <- pomdpSolve::read_belief_file(paste0(file_prefix, 
                                                "-0.belief"))
  if (!is.null(belief)) 
    colnames(belief) <- as.character(model$states)
  belief
}

.max_abs_reward<-function (x, episode = NULL, epoch = NULL) 
{
  r <- reward_matrix(x, episode = episode, epoch = epoch, sparse = NULL)
  if (is.function(r)) 
    r <- reward_matrix(x, episode = episode, epoch = epoch, 
                       sparse = FALSE)
  if (is.data.frame(r)) {
    rew <- r$value
  }
  else {
    rew <- unlist(r, recursive = TRUE)
  }
  rew[rew == -Inf] <- NA
  max(abs(rew), na.rm = TRUE)
}

.format_number_fixed<-function (x, digits = 7, debug = "unknown") 
{
  if (is.null(x)) 
    stop("missing field ", debug)
  if (is.vector(x)) {
    if (!is.numeric(x)) 
      stop("Write_POMDP expects numbers, but got: ", dQuote(paste(x, 
                                                                  collapse = ", ")))
    paste(sprintf(paste0("%.", digits, "f"), x), collapse = " ")
  }
  else if (is.matrix(x)) {
    paste(apply(x, MARGIN = 1, .format_number_fixed, digits = digits), 
          collapse = "\n")
  }
  else stop("formating not implemented for ", class(x), " in field ", 
            debug)
}
  
write_POMDP_kp<- function (x, file, digits = 7, labels = FALSE) {
  labels = FALSE
  digits = 7
  if (!inherits(x, "POMDP")) 
    stop("model needs to be a POMDP model use POMDP()!")
  mr <- .max_abs_reward(x)
  reward_max <- 1e+10
  if (mr > reward_max) 
    stop("Reward values supported by the solver need to be in [", 
         -reward_max, ", ", +reward_max, "].")
  reward_neg_Inf <- -2 * mr
  pomdp_solve_OK_chars <- "[^a-zA-Z0-9_-]"
  x <- check_and_fix_MDP_kp(x)
  if (is.character(file)) {
    file <- file(file, "w")
    on.exit(close(file))
  }
  with(x, {
    number_of_states <- length(states)
    number_of_observations <- length(observations)
    number_of_actions <- length(actions)
    values <- "reward"
    if (is.function(transition_prob)) 
      transition_prob <- transition_matrix(x, sparse = FALSE)
    if (is.function(observation_prob)) 
      observation_prob <- observation_matrix(x, sparse = FALSE)
    if (is.function(reward)) 
      reward <- reward_matrix(x, sparse = FALSE)
    if (labels) {
      if (is.character(states) && any(grepl(pomdp_solve_OK_chars, 
                                            states))) 
        stop("Some state labels use characters unsupported by the solver! Only use ", 
             pomdp_solve_OK_chars)
      if (is.character(actions) && any(grepl(pomdp_solve_OK_chars, 
                                             actions))) 
        stop("Some action labels use characters unsupported by the solver! Only use ", 
             pomdp_solve_OK_chars)
      if (is.character(observations) && any(grepl(pomdp_solve_OK_chars, 
                                                  observations))) 
        stop("Some observation labels use characters unsupported by the solver! Only use ", 
             pomdp_solve_OK_chars)
    }
    preamble <- paste0("# POMDP File: ", name, "\n", "# Produced with R package pomdp (created: ", 
                       date(), ")\n", "\n", "discount: ", format(discount, 
                                                                 digits = digits), "\n", "values: ", values, "\n", 
                       "states: ", ifelse(!labels, number_of_states, paste(states, 
                                                                           collapse = " ")), "\n", "actions: ", ifelse(!labels, 
                                                                                                                       number_of_actions, paste(actions, collapse = " ")), 
                       "\n", "observations: ", ifelse(!labels, number_of_observations, 
                                                      paste(observations, collapse = " ")), "\n")
    if (!is.null(start)) {
      if (is.numeric(start)) {
        if (length(start) == length(states) && sum1(start)) {
          preamble_start <- paste0("start: ", .format_number_fixed(round_stochastic(start, 
                                                                                    digits), digits, "start"), "\n")
        }
        else {
          start_ids <- as.integer(abs(start)) - 1L
          if (all(start < 0)) 
            preamble_start <- paste0("start include: ", 
                                     paste(start_ids, collapse = " "), "\n")
          if (all(start > 0)) 
            preamble_start <- paste0("start exclude: ", 
                                     paste(start_ids, collapse = " "), "\n")
          else stop("Illegal specification of start. ", 
                    "State ids start with 1 need to be all positive or negative (for exclusion).")
        }
      }
      else if (is.character(start)) 
        if (length(start) == 1 && start[1] == "uniform") {
          preamble_start <- paste0("start: ", paste(start, 
                                                    collapse = " "), "\n")
        }
      else {
        if (start[1] == "-") {
          pre <- "start exclude: "
          start_states <- start[-1]
        }
        else {
          pre <- "start include: "
          start_states <- start
        }
        if (!labels) {
          start_states <- factor(start_states, levels = states)
          start_states <- as.integer(start_states) - 
            1L
        }
        preamble_start <- paste0(pre, paste(start_states, 
                                            collapse = " "), "\n")
      }
      else stop("Not suported!")
    }
    cat(preamble, "\n", file = file)
    cat(preamble_start, "\n", file = file)
    format_POMDP_df <- function(x, file, type = c("T", "R", 
                                                  "O"), digits = 7) {
      type <- match.arg(type)
      if (type == "R") 
        cat(paste0(type, ": ", .to_0_idx(x[[1L]]), " : ", 
                   .to_0_idx(x[[2L]]), " : ", .to_0_idx(x[[3L]]), 
                   " : ", .to_0_idx(x[[4L]]), " ", sprintf(paste0("%.", 
                                                                  digits, "f"), x[[5L]]), collapse = "\n"), 
            file = file, append = TRUE)
      else cat(paste0(type, ": ", .to_0_idx(x[[1L]]), " : ", 
                      .to_0_idx(x[[2L]]), " : ", .to_0_idx(x[[3L]]), 
                      " ", sprintf(paste0("%.", digits, "f"), x[[4L]]), 
                      collapse = "\n"), file = file, append = TRUE)
      cat("\n", file = file)
    }
    format_POMDP_matrix <- function(x, file, type = c("T", 
                                                      "R", "O"), action, start = NULL, digits = 7) {
      type <- match.arg(type)
      if (inherits(x, "Matrix")) 
        return(format_POMDP_dgc(x, file, type, action, 
                                start, digits))
      action <- .to_0_idx(factor(action, levels = actions))
      if (type == "R") 
        code <- paste0(type, ": ", action, " : ", start, 
                       "\n")
      else code <- paste0(type, ": ", action, "\n")
      if (is.character(x) && length(x) == 1) 
        code <- paste0(code, x, "\n")
      else code <- paste0(code, .format_number_fixed(x, 
                                                     digits, type), "\n")
      cat(code, file = file, append = TRUE)
      cat("\n", file = file)
    }
    format_POMDP_dgc <- function(x, file, type = c("T", "R", 
                                                   "O"), action, start = NULL, digits = 7) {
      type <- match.arg(type)
      x <- as(x, "TsparseMatrix")
      if (nnzero(x) == 0L) 
        return()
      action <- .to_0_idx(factor(action, levels = actions))
      if (type == "R") 
        prefix <- paste(action, start, sep = " : ")
      else prefix <- action
      code <- paste0(type, ": ", prefix, " : ", x@i, " : ", 
                     x@j, " ", sprintf(paste0("%.", digits, "f"), 
                                       x@x))
      code <- paste0(code, collapse = "\n")
      cat(code, file = file, append = TRUE)
      cat("\n\n", file = file)
    }
    if (is.data.frame(transition_prob)) 
      format_POMDP_df(transition_prob, file, "T", digits)
    else {
      if (!identical(names(transition_prob), "*") && !setequal(names(transition_prob), 
                                                               actions)) 
        stop("the names of given transition probability matrices do not match the actions!")
      for (a in names(transition_prob)) format_POMDP_matrix(transition_prob[[a]], 
                                                            file, "T", a, digits = digits)
    }
    if (is.data.frame(observation_prob)) 
      format_POMDP_df(observation_prob, file, "O", digits)
    else {
      if (!identical(names(observation_prob), "*") && !setequal(names(observation_prob), 
                                                                actions)) 
        stop("the names of given observation probability matrices do not match the actions!")
      for (a in names(observation_prob)) format_POMDP_matrix(observation_prob[[a]], 
                                                             file, "O", a, digits = digits)
    }
    if (is.data.frame(reward)) {
      if (any((reward$value < -reward_max | reward$value > 
               reward_max) & reward$value != -Inf)) 
        stop("Reward values supported by the solver need to be in [", 
             -reward_max, ", ", +reward_max, "].")
      reward$value[reward$value == -Inf] <- reward_neg_Inf
      format_POMDP_df(reward, file, "R", digits)
    }
    else {
      if (!identical(names(reward), "*") && !setequal(names(reward), 
                                                      actions)) 
        stop("names of the rewards list do not match the actions!")
      for (a in names(reward)) {
        if (!identical(names(reward[[a]]), "*") && !setequal(names(reward[[a]]), 
                                                             states)) 
          stop("names of the second level of the rewards list do not match the states!")
        for (s in names(reward[[a]])) {
          reward_matrix <- reward[[a]][[s]]
          if (any((reward_matrix < -reward_max | reward_matrix > 
                   reward_max) & reward_matrix != -Inf)) 
            stop("Reward values supported by the solver need to be in [", 
                 -reward_max, ",", +reward_max, "]")
          reward_matrix[reward_matrix == -Inf] <- reward_neg_Inf
          format_POMDP_matrix(reward_matrix, file, "R", 
                              a, s, digits = digits)
        }
      }
    }
  })
}
solve_POMDP_kp<- function () 
{
  model<- sc
  horizon<- Inf
  discount = NULL
  initial_belief <- NULL
  terminal_values<- NULL
  method = "grid"
  digits<- 7
  parameter<- NULL
  timeout = Inf
  verbose = TRUE
  
  converged <- NA
  methods <- c("grid", "enum", "twopass", "witness", "incprune")
  method <- match.arg(method, methods)
  if (is.character(model)) 
    model <- read_POMDP(model)
  if (!inherits(model, "POMDP")) 
    stop("x needs to be a POMDP!")
  if (is.null(horizon)) 
    horizon <- model$horizon
  if (is.null(horizon)) 
    horizon <- Inf
  model$horizon <- horizon
  if (!is.null(initial_belief)) {
    model$start <- initial_belief
  }
  if (is.null(model$start)) 
    model$start <- "uniform"
  if (is_timedependent_POMDP(model)) 
    return(.solve_POMDP_time_dependent(model = model, horizon = horizon, 
                                       discount = discount, terminal_values = terminal_values, 
                                       method = method, parameter = parameter, verbose = verbose))
  if (horizon < 1) 
    horizon <- Inf
  if (is.null(terminal_values)) 
    terminal_values <- model$terminal_values
  if (!is.null(terminal_values) && length(terminal_values) == 
      1 && terminal_values == 0) 
    terminal_values <- NULL
  model$terminal_values <- terminal_values
  if (is.null(discount)) 
    discount <- model$discount
  if (is.null(discount)) {
    message("No discount rate specified. Using .9!")
    discount <- 0.9
  }
  model$discount <- discount
  file_prefix <- tempfile(pattern = "pomdp_")
  pomdp_filename <- paste0(file_prefix, ".POMDP")
  if (verbose) 
    cat("Writing the problem to", pomdp_filename, "\n")
  if (is.null(model$transition_prob) || is.null(model$observation_prob) || 
      is.null(model$reward)) {
    cat("Writing with writelines \n")
    writeLines(model$problem, con = pomdp_filename)
  } else {
    cat("Writing with writepomdp \n")
    write_POMDP(model, pomdp_filename, digits = digits)
  }
  if (!is.null(terminal_values)) {
    if (!is.matrix(terminal_values)) 
      terminal_values <- rbind(terminal_values)
    if (ncol(terminal_values) != length(model$states)) 
      stop("number of terminal values does not match the number of states.")
    colnames(terminal_values) <- as.character(model$states)
    terminal_values_filename <- .write_alpha_file(file_prefix, 
                                                  terminal_values)
  }
  if (!is.null(parameter$grid)) {
    if (method != "grid") 
      warning("Custom grids are ignored by all methods but 'grid'!")
    parameter$grid_filename <- .write_grid_file(file_prefix, 
                                                parameter$grid)
    parameter$grid <- NULL
    parameter$fg_type <- "file"
  }
  if (!is.null(parameter)) {
    paras <- as.vector(sapply(names(parameter), FUN = function(n) {
      c(paste0("-", n), 
        if (is.logical(parameter[[n]])) {
          if (parameter[[n]]) "true" else "false"
        } else parameter[[n]])
    }))
  } else {
    paras <- NULL
  }
  pomdp_args <- c("-pomdp", pomdp_filename, "-method", method)
  if (method == "grid") 
    pomdp_args <- append(pomdp_args, c("-fg_save", "true"))
  if (is.finite(horizon)) 
    pomdp_args <- append(pomdp_args, c("-horizon", horizon))
  if (timeout || is.finite(horizon)) 
    pomdp_args <- append(pomdp_args, c("-save_all", "true"))
  if (!is.null(discount)) 
    pomdp_args <- append(pomdp_args, c("-discount", discount))
  if (!is.null(terminal_values)) 
    pomdp_args <- append(pomdp_args, c("-terminal_values", 
                                       terminal_values_filename))
  if (!is.null(paras)) 
    pomdp_args <- append(pomdp_args, paras)
  if (verbose) 
    cat("Calling pomdp-solve with the following arguments:", 
        paste(pomdp_args, collapse = " "), "\nSolver output:", 
        sep = "\n")
  solver_output <- processx::run(pomdpSolve::find_pomdp_solve(), 
                                 args = pomdp_args, echo = verbose, timeout = timeout, 
                                 error_on_status = FALSE)
  if (solver_output$status != 0 && !verbose) {
    cat(paste(solver_output$stdout))
  }
  
  if (solver_output$timeout) {
    cat(" timeout reached!\n")
    if (is.finite(horizon))
      stop("Unfinished solutions cannot be used for finite horizon problems!")
    ep <- strsplit(solver_output$stdout, "\n")[[1]]
    ep <- ep[length(ep) - 1]
    ep <- as.integer(gsub("Epoch: (\\d+).*", "\\1", ep))
    if (is.na(ep))
      stop("Could not find a solved epoch. You may need to increase the timeout.")
    cat("Load last solved epoch: Epoch ", ep, "\n\n")
    converged <- FALSE
    alpha <- list(.get_alpha_file(file_prefix, model, ep))
    pg <- list(.get_pg_file(file_prefix, model, ep))
  } else if (!is.finite(horizon)) {
    converged <- TRUE
    alpha <- list(.get_alpha_file(file_prefix, model))
    pg <- list(.get_pg_file(file_prefix, model))
  } else {
    converged <- FALSE
    alpha <- list()
    pg <- list()
    for (i in 1:horizon) {
      r <- suppressWarnings(try({
        alpha[[i]] <- .get_alpha_file(file_prefix, model, 
                                      i)
        pg[[i]] <- .get_pg_file(file_prefix, model, i)
      }, silent = TRUE))
      if (inherits(r, "try-error")) {
        if (verbose) 
          cat("Convergence: Finite-horizon POMDP converged early at epoch:", 
              i - 1, "\n")
        converged <- TRUE
        pg <- tail(pg, n = 1L)
        alpha <- tail(alpha, n = 1L)
        break
      }
    }
    alpha <- rev(alpha)
    pg <- rev(pg)
  }
  belief <- .get_belief_file(file_prefix, model)
  model$solution <- structure(list(method = method, parameter = parameter, 
                                   converged = converged, total_expected_reward = NA, initial_belief = NA, 
                                   initial_pg_node = NA, belief_points_solver = belief, 
                                   pg = pg, alpha = alpha), class = "POMDP_solution")
  rew <- reward_node_action(model, belief = model$start)
  model$solution$initial_belief <- rew$belief
  model$solution$total_expected_reward <- rew$reward
  model$solution$initial_pg_node <- rew$pg_node
  model$solution$solver_output <- structure(solver_output, 
                                            class = "text")
  if (inherits(model, "MDP")) 
    model$solution$policy <- .MDP_policy_from_POMDP(model)
  model
}

POMDP_kp<- function (states, actions, observations, transition_prob, observation_prob, 
          reward, discount = 0.9, horizon = Inf, terminal_values = NULL, 
          start = "uniform", info = NULL, name = NA) 
{
  discount = 0.9
  horizon = Inf
  terminal_values = NULL
  info = NULL
  start = start
  name = "sc"
  transition_prob = transition_matrices
  observation_prob = observation_matrices
  observations = states
  reward = rewards3_fixed
  saved_x <- list(name = name, discount = discount, horizon = horizon, 
            states = states, actions = actions, observations = observations, 
            transition_prob = transition_prob, observation_prob = observation_prob, 
            reward = reward, start = start, terminal_values = terminal_values, 
            info = info)
  class(saved_x) <- c("POMDP", "list")
  x <- check_and_fix_MDP_kp(saved_x)
  x
}
sum1<- function (x, digits = getOption("digits")) 
{
  if (is.matrix(x)) 
    all(apply(x, MARGIN = 1, sum1))
  else zapsmall(sum(x), digits) == 1
}
.is_timedependent_field<- function (x, field) 
{
  field <- match.arg(field, c("transition_prob", "observation_prob", 
                              "reward"))
  m <- x[[field]]
  if (is.null(m)) 
    stop("Field ", field, " does not exist.")
  if (!is.list(m) || is.data.frame(m)) 
    return(FALSE)
  if (!is.list(m[[1]])) 
    return(FALSE)
  if (field == "reward" && !is.list(m[[1]][[1]])) 
    return(FALSE)
  if (length(m) != length(x$horizon)) 
    stop("Inconsistent POMDP specification. Field ", field, 
         " does not contain data for the appropriate number of episodes.")
  TRUE
}

check_and_fix_MDP_kp<- function (x) 
  
{
  x<-saved_x
  check_func <- function(x, func, name) {
    req_formals <- head(names(formals(func)), -1)
    if (!identical(names(formals(x)), req_formals)) 
      stop(name, " function needs formal arguments: ", 
           paste(sQuote(req_formals), collapse = ", "))
  }
  check_df <- function(x, field, func) {
    req_columns <- names(formals(func))
    if (is.null(colnames(field))) 
      colnames(field) <- req_columns
    if (!identical(colnames(field), req_columns)) 
      stop("The ", deparse(substitute(field)), " data.frame needs columns named: ", 
           paste(sQuote(req_columns), collapse = ", "))
    field[field == "*"] <- NA
    field <- type.convert(field, as.is = TRUE)
    for (i in grep("action", colnames(field))) {
      if (is.numeric(field[[i]])) 
        field[[i]] <- x$actions[field[[i]]]
      field[[i]] <- factor(field[[i]], levels = x$actions)
    }
    for (i in grep("state", colnames(field))) {
      if (is.numeric(field[[i]])) 
        field[[i]] <- x$states[field[[i]]]
      field[[i]] <- factor(field[[i]], levels = x$states)
    }
    for (i in grep("observation", colnames(field))) {
      if (is.numeric(field[[i]])) 
        field[[i]] <- x$observations[field[[i]]]
      field[[i]] <- factor(field[[i]], levels = x$observations)
    }
    field
  }
  if (is.numeric(x$states) && length(x$states) == 1L) 
    x$states <- paste0("s", seq_len(x$states))
  if (is.numeric(x$actions) && length(x$actions) == 1L) 
    x$actions <- paste0("a", seq_len(x$actions))
  if (inherits(x, "POMDP")) {
    if (is.numeric(x$observations) && length(x$observations) == 
        1L) 
      x$observations <- paste0("o", seq_len(x$observations))
  }
  x$discount <- as.numeric(x$discount)
  if (length(x$discount) != 1L || x$discount <= 0 || x$discount > 
      1) 
    stop("discount has to be a single value in the range (0,1].")
  if (is.null(x$horizon)) 
    x$horizon <- Inf
  x$horizon <- as.numeric(x$horizon)
  if (any(x$horizon != floor(x$horizon))) 
    stop("horizon needs to be an integer.")
  if (is.numeric(x$start) && length(x$start) == length(x$states)) {
    if (!sum1(x$start)) 
      stop("The start probability vector does not add up to 1.")
    if (is.null(names(x$start))) 
      names(x$start) <- x$states
    else x$start <- x$start[x$states]
  }
  if (any(is.na(x$start))) 
    stop("start containes undefined start states.")
  if (is.character(x$start)) {
    if (!(identical(x$start, "uniform") || all(x$start %in% 
                                               x$states))) 
      stop("when using characters for start, then it needs to be the keyword 'uniform' or a set of start states.")
  }
  if ((is.null(x$transition_prob) || (inherits(x, "POMDP") && 
                                      is.null(x$observation_prob)) || is.null(x$reward)) && 
      is.null(x$problem)) 
    stop("transition_prob, observation_prob or reward can only miss if the field problem is set!")
  if (!is.null(x$transition_prob) && !.is_timedependent_field(x, 
                                                              "transition_prob")) {
    if (is.function(x$transition_prob)) 
      check_func(x$transition_prob, T_, "transition_prob")
    else if (is.data.frame(x$transition_prob)) 
      x$transition_prob <- check_df(x, x$transition_prob, 
                                    T_)
    else {
      if (is.null(names(x$transition_prob))) 
        names(x$transition_prob) <- x$actions
      if (all(names(x$transition_prob) != x$actions)) 
        x$transition_prob <- x$transition_prob[x$actions]
      for (a in x$actions) {
        if (is.null(x$transition_prob[[a]])) 
          stop("transition_prob for action ", a, " is missing!")
        if (is.matrix(x$transition_prob[[a]])) {
          if (!identical(dim(x$transition_prob[[a]]), 
                         c(length(x$states), length(x$states)))) 
            stop("transition_prob matrix for action ", 
                 a, ": has not the right dimensions!")
          if (!sum1(x$transition_prob[[a]])) 
            stop("transition_prob matrix for action ", 
                 a, ": rows do not add up to 1!")
          if (is.null(dimnames(x$transition_prob[[a]]))) 
            dimnames(x$transition_prob[[a]]) <- list(x$states, 
                                                     x$states)
          else x$transition_prob[[a]][x$states, x$states]
        }
      }
    }
  }
  if (!is.null(x$transition_prob) && .is_timedependent_field(x, 
                                                             "transition_prob")) {
    for (e in seq_along(x$horizon)) {
      if (is.function(x$transition_prob[[e]])) 
        check_func(x$transition_prob[[e]], T_, "transition_prob")
      else if (is.data.frame(x$transition_prob[[e]])) 
        x$transition_prob[[e]] <- check_df(x, x$transition_prob[[e]], 
                                           T_)
      else {
        if (is.null(names(x$transition_prob[[e]]))) 
          names(x$transition_prob[[e]]) <- x$actions
        if (all(names(x$transition_prob[[e]]) != x$actions)) 
          x$transition_prob[[e]] <- x$transition_prob[[e]][x$actions]
        for (a in x$actions) {
          if (is.null(x$transition_prob[[e]][[a]])) 
            stop("transition_prob for action ", a, " is missing in epoch ", 
                 e, "!")
          if (is.matrix(x$transition_prob[[e]][[a]])) {
            if (!identical(dim(x$transition_prob[[e]][[a]]), 
                           c(length(x$states), length(x$states)))) 
              stop("transition_prob matrix for action ", 
                   a, " in epoch ", e, ": has not the right dimensions!")
            if (!sum1(x$transition_prob[[e]][[a]])) 
              stop("transition_prob matrix for action ", 
                   a, " in epoch ", e, ": rows do not add up to 1!")
            if (is.null(dimnames(x$transition_prob[[e]][[a]]))) 
              dimnames(x$transition_prob[[e]][[a]]) <- list(x$states, 
                                                            x$states)
            else x$transition_prob[[e]][[a]][x$states, 
                                             x$states]
          }
        }
      }
    }
  }
  if (!is.null(x$reward) && !.is_timedependent_field(x, "reward")) {
    if (!inherits(x, "POMDP")) {
      cat("not inherits pomdp\n")
      R_ <- function(action = NA, start.state = NA, end.state = NA, 
                     value) data.frame(action = action, start.state = start.state, 
                                       end.state = end.state, value = as.numeric(value), 
                                       stringsAsFactors = FALSE)
    }
    if (is.function(x$reward)) 
      cat("x$reward is func \n")
      check_func(x$reward, R_, "reward")
    if (is.data.frame(x$reward)) {
      cat("x$reward is df \n")
      x$reward <- check_df(x, x$reward, R_)
    }
  }
  else {
    if (is.null(names(x$reward))) 
      names(x$reward) <- x$actions
    if (all(names(x$reward) != x$actions)) 
      x$reward <- x$reward[x$actions]
    for (a in x$actions) {
      if (is.null(x$reward[[a]])) 
        stop("reward for action ", a, " is missing!")
      for (s in x$states) {
        if (is.null(x$reward[[a]][[s]])) 
          stop("reward for action ", a, " and state ", 
               s, " is missing!")
        if (is.matrix(x$reward[[a]][[s]])) {
          if (!identical(dim(x$reward[[a]][[s]]), c(length(x$states), 
                                                    length(x$observations)))) 
            stop("reward matrix for action ", a, " and start.state ", 
                 s, ": has not the right dimensions!")
          if (is.null(dimnames(x$reward[[a]][[s]]))) 
            dimnames(x$reward[[a]][[s]]) <- list(x$states, 
                                                 x$observations)
          else x$reward[[a]][[s]][x$states, x$observations]
        }
      }
    }
  }
  if (!is.null(x$reward) && .is_timedependent_field(x, "reward")) {
    for (e in seq_along(x$horizon)) {
      if (is.function(x$reward[[e]])) 
        check_func(x$reward[[e]], R_, "reward")
      if (is.data.frame(x$reward[[e]])) 
        x$reward[[e]] <- check_df(x, x$reward[[e]], R_)
      else {
        if (is.null(names(x$reward[[e]]))) 
          names(x$reward[[e]]) <- x$actions
        if (all(names(x$reward[[e]]) != x$actions)) 
          x$reward[[e]] <- x$reward[[e]][x$actions]
        for (a in x$actions) {
          if (is.null(x$reward[[e]][[a]])) 
            stop("reward for action ", a, " in episode ", 
                 e, " is missing!")
          for (s in x$states) {
            if (is.null(x$reward[[e]][[a]][[s]])) 
              stop("reward for action ", a, " and state ", 
                   s, " in episode ", e, " is missing!")
            if (is.matrix(x$reward[[e]][[a]][[s]])) {
              if (!identical(dim(x$reward[[e]][[a]][[s]]), 
                             c(length(x$states), length(x$observations)))) 
                stop("reward matrix for action ", a, 
                     " and start.state ", s, " in episode ", 
                     e, ": has not the right dimensions!")
              if (is.null(dimnames(x$reward[[e]][[a]][[s]]))) 
                dimnames(x$reward[[e]][[a]][[s]]) <- list(x$states, 
                                                          x$observations)
              else x$reward[[e]][[a]][[s]][x$states, 
                                           x$observations]
            }
          }
        }
      }
    }
  }
  if (inherits(x, "POMDP") && !is.null(x$terminal_values)) {
    if (length(x$terminal_values) != 1L && length(x$terminal_values) != 
        length(x$states) && (is.matrix(x$terminal_values) && 
                             ncol(x$terminal_values) != length(x$states))) 
      stop("Terminal values are not in the right format.")
  }
  if (inherits(x, "POMDP") && !is.null(x$solution)) {
    if (any(sapply(x$solution$alpha, ncol) != length(x$states))) 
      stop("Alpha vectors do not have the right dimension.")
    x$solution$pg <- lapply(x$solution$pg, FUN = function(y) {
      y$action <- factor(y$action, levels = x$actions)
      y
    })
  }
  if (inherits(x, "POMDP") && !is.null(x$observation_prob) && 
      !.is_timedependent_field(x, "observation_prob")) {
    if (is.function(x$observation_prob)) 
      check_func(x$observation_prob, O_, "observation_prob")
    else if (is.data.frame(x$observation_prob)) 
      x$observation_prob <- check_df(x, x$observation_prob, 
                                     O_)
    else {
      if (is.null(names(x$observation_prob))) 
        names(x$observation_prob) <- x$actions
      if (all(names(x$observation_prob) != x$actions)) 
        x$observation_prob <- x$observation_prob[x$actions]
      for (a in x$actions) {
        if (is.null(x$observation_prob[[a]])) 
          stop("observation_prob for action ", a, " is missing!")
        if (is.matrix(x$observation_prob[[a]])) {
          if (!identical(dim(x$observation_prob[[a]]), 
                         c(length(x$states), length(x$observations)))) 
            stop("observation_prob matrix for action ", 
                 a, ": has not the right dimensions!")
          if (!sum1(x$observation_prob[[a]])) 
            stop("observation_prob matrix for action ", 
                 a, ": rows do not add up to 1!")
          if (is.null(dimnames(x$observation_prob[[a]]))) 
            dimnames(x$observation_prob[[a]]) <- list(x$states, 
                                                      x$observations)
          else x$observation_prob[[a]][x$states, x$observations]
        }
      }
    }
  }
  if (inherits(x, "POMDP") && !is.null(x$observation_prob) && 
      .is_timedependent_field(x, "observation_prob")) {
    for (e in seq_along(x$horizon)) {
      if (is.function(x$observation_prob[[e]])) 
        check_func(x$observation_prob[[e]], O_, "observation_prob")
      if (is.data.frame(x$observation_prob[[e]])) 
        x$observation_prob[[e]] <- check_df(x, x$observation_prob[[e]], 
                                            O_)
      else {
        if (is.null(names(x$observation_prob[[e]]))) 
          names(x$observation_prob[[e]]) <- x$actions
        if (all(names(x$observation_prob[[e]]) != x$actions)) 
          x$observation_prob[[e]] <- x$observation_prob[[e]][x$actions]
        for (a in x$actions) {
          if (is.null(x$observation_prob[[e]][[a]])) 
            stop("observation_prob for action ", a, " is missing!")
          if (is.matrix(x$observation_prob[[e]][[a]])) {
            if (!identical(dim(x$observation_prob[[e]][[a]]), 
                           c(length(x$states), length(x$observations)))) 
              stop("observation_prob matrix for action ", 
                   a, ": has not the right dimensions!")
            if (!all(rowSums(x$observation_prob[[e]][[a]]) == 
                     1)) 
              stop("observation_prob matrix for action ", 
                   a, ": rows do not add up to 1!")
            if (is.null(dimnames(x$observation_prob[[e]][[a]]))) 
              dimnames(x$observation_prob[[e]][[a]]) <- list(x$states, 
                                                             x$observations)
            else x$observation_prob[[e]][[a]][x$states, 
                                              x$observations]
          }
        }
      }
    }
  }
  x
}
simulate_POMDP_kp<- function (model, n = 1000, belief = NULL, horizon = NULL, epsilon = NULL, 
          delta_horizon = 0.001, digits = 7L, return_beliefs = FALSE, 
          return_trajectories = FALSE, engine = "r", verbose = FALSE, 
          ...) 
{
  model = sc
  n = 1000
  belief = vec
  engine = "r"
  horizon = Inf
  digits = 7L
  epsilon = NULL
  return_trajectories = TRUE
  time_start <- proc.time()
  engine <- match.arg(tolower(engine), c("cpp", "r"))

  solved <- is_solved_POMDP(model)
  cat("Model solved:", solved)
  dt <- is_timedependent_POMDP(model)
  cat("Model timdedependent:", dt)
  if (is.null(belief)) 
    belief <- start_vector(model)
  if (!is.numeric(belief) || length(belief) != length(model$states) || 
      !sum1(belief)) 
    stop("Initial belief is misspecified!")
  n <- as.integer(n)
  digits <- as.integer(digits)
  if (is.null(horizon)) 
    horizon <- model$horizon
  if (is.null(horizon) || is.infinite(horizon)) {
    if (is.null(model$discount) || !(model$discount < 1)) 
      stop("Simulation needs a finite simulation horizon.")
    max_abs_R <- .max_abs_reward(model)
    horizon <- ceiling(log(delta_horizon/max_abs_R)/log(model$discount)) + 
      1
  }
  horizon <- as.integer(horizon)
  if (is.null(epsilon)) {
    if (!solved) 
      epsilon <- 1
    else epsilon <- 0
  }
  if (!solved && epsilon != 1) 
    stop("epsilon has to be 1 for unsolved models.")
  disc <- model$discount
  if (is.null(disc)) 
    disc <- 1

  states <- as.character(model$states)
  n_states <- length(states)
  states_absorbing <- which(absorbing_states(model))
  obs <- as.character(model$observations)
  n_obs <- length(obs)
  actions <- as.character(model$actions)
  current_episode <- 1L
  if (dt) {
    dt_horizon <- model$horizon
    dt_episodes <- cumsum(c(1, head(model$horizon, -1)))
    dt_trans_m <- lapply(1:length(dt_horizon), FUN = function(ep) transition_matrix(model, 
                                                                                    ep))
    dt_obs_m <- lapply(1:length(dt_horizon), FUN = function(ep) observation_matrix(model, 
                                                                                   ep))
    trans_m <- dt_trans_m[[current_episode]]
    obs_m <- dt_obs_m[[current_episode]]
  }
  else {
    trans_m <- transition_matrix(model, sparse = NULL)
    obs_m <- observation_matrix(model, sparse = NULL)
  }
  if (verbose) {
    cat("Simulating POMDP trajectories.\n")
    cat("- engine: r\n")
    cat("- horizon:", horizon, "\n")
    cat("- n:", n, "- parallel workers:", foreach::getDoParWorkers(), 
        "\n")
    cat("- epsilon:", epsilon, "\n")
    if (dt) 
      cat("- time-dependent:", length(dt_horizon), "episodes", 
          "\n")
    cat("- discount factor:", disc, "\n")
    cat("- starting belief:\n")
    print(head(belief, n = 10))
    if (length(belief) > 10) 
      cat("(Remaining belief components supressed)")
    cat("\n")
  }
  sim <- foreach(i = 1:n) %dopar% {
    alpha_vec_id <- NA_integer_
    s <- sample.int(length(states), 1L, prob = belief)
    b <- belief
    rew <- 0
    e <- 1L
    action_cnt <- rep(0L, length(actions))
    names(action_cnt) <- actions
    state_cnt <- rep(0L, length(states))
    names(state_cnt) <- states
    obs_cnt <- rep(0L, length(obs))
    names(obs_cnt) <- obs
    if (return_beliefs) 
      visited_belief_states <- matrix(NA, nrow = horizon, 
                                      ncol = n_states, dimnames = list(NULL, states))
    else visited_belief_states <- matrix(nrow = 0, ncol = 0)
    if (return_trajectories) 
      trajectory <- data.frame(episode = rep(NA_integer_, 
                                             horizon), time = rep(NA_integer_, horizon), simulation_state = NA_integer_, 
                               alpha_vector_id = NA_integer_, a = NA_integer_, 
                               r = NA_real_)
    else trajectory <- NULL
    for (j in 1:horizon) {
      if (dt) {
        if (length(current_episode <- which(j == dt_episodes)) == 
            1L) {
          if (verbose) 
            cat("- Switching to episode", current_episode, 
                "at epoch", j, "\n")
          obs_m <- dt_obs_m[[current_episode]]
          trans_m <- dt_trans_m[[current_episode]]
        }
      }
      if (runif(1) < epsilon) {
        a <- sample.int(length(actions), 1L)
      }
      else {
        if (!model$solution$converged) 
          e <- .get_pg_index(model, j)
        alpha_vec_id <- which.max(model$solution$alpha[[e]] %*% 
                                    b)
        a <- as.integer(model$solution$pg[[e]][["action"]])[alpha_vec_id]
      }
      s_prev <- s
      s <- sample.int(length(states), 1L, prob = trans_m[[a]][s, 
      ])
      o <- sample.int(length(obs), 1L, prob = obs_m[[a]][s, 
      ])
      action_cnt[a] <- action_cnt[a] + 1L
      state_cnt[s] <- state_cnt[s] + 1L
      obs_cnt[o] <- obs_cnt[o] + 1L
      r <- reward_matrix(model, a, s_prev, s, o, episode = current_episode)
      rew <- rew + r * disc^(j - 1L)
      b <- .update_belief(b, a, o, trans_m, obs_m, digits = digits)
      if (return_beliefs) 
        visited_belief_states[j, ] <- b
      if (return_trajectories) 
        trajectory[j, ] <- data.frame(episode = i, time = j - 
                                        1L, simulation_state = s_prev, alpha_vector_id = alpha_vec_id, 
                                      a = a, r = r)
      if (s %in% states_absorbing) {
        if (return_trajectories) 
          trajectory <- trajectory[1:j, , drop = FALSE]
        visited_belief_states <- visited_belief_states[1:j, 
                                                       , drop = FALSE]
        break
      }
    }
    if (j == sum(model$horizon) && !is.null(model$terminal_values)) {
      rew <- rew + model$terminal_values[s] * disc^j
    }
    rownames(visited_belief_states) <- NULL
    list(action_cnt = action_cnt, state_cnt = state_cnt, 
         obs_cnt = obs_cnt, reward = rew, belief_states = visited_belief_states, 
         trajectory = trajectory)
  }
  time_end <- proc.time()
  if (verbose) 
    print(time_end - time_start)
  rew <- Reduce(c, lapply(sim, "[[", "reward"))
  trajectories <- NULL
  if (return_trajectories) {
    trajectories <- Reduce(rbind, lapply(sim, "[[", "trajectory"))
    trajectories$simulation_state <- factor(trajectories$simulation_state, 
                                            levels = seq_along(states), labels = states)
    trajectories$a <- factor(trajectories$a, levels = seq_along(actions), 
                             labels = actions)
  }
  list(avg_reward = mean(rew, na.rm = TRUE), action_cnt = Reduce("+", 
                                                                 lapply(sim, "[[", "state_cnt")), state_cnt = Reduce("+", 
                                                                                                                     lapply(sim, "[[", "state_cnt")), obs_cnt = Reduce("+", 
                                                                                                                                                                       lapply(sim, "[[", "obs_cnt")), reward = rew, belief_states = Reduce(rbind, 
                                                                                                                                                                                                                                           lapply(sim, "[[", "belief_states")), trajectories = trajectories)
}