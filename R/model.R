#' check if system is ready for model run
check_model <- function(ip) {
  if (!is(ip, "isopath")) stop ("can only check model ready for an isopath", call. = FALSE)

  message(sprintf("Checking model for %d scenario(s)...", nrow(ip$parameters)))

  tryCatch({

    # do a sample calculation for each parameter set
    ip$parameters %>%
      group_by_(.dots = names(ip$parameters)) %>%
      do({
        # safety checks - will through errors if needed
        if (nrow(.) != 1) stop("there seem to be multiple identical run scenarios, please make sure each set of parameters is unique", call. = FALSE)
        dx_var <- get_ode_matrix(ip, eval = T, param = .)
        if (any(is.na(dx_var$`dx/dt`))) {
          dx_text <- get_ode_matrix(ip, eval = F)$`dx/dt`
          stop("some derivatives could not be computed, check the expressions: ",
               paste(dx_var$x[is.na(dx_var$`dx/dt`)], "=", dx_text[is.na(dx_var$`dx/dt`)]) %>%
                 paste(collapse = ", "), call. = FALSE)
        }
        tibble()
      })

  },

  # error handling
  error = function(e) {
    stop("encountered the following error during pre-check of the model run: <<< ", e$message, " >>>. If this is a missing object error, please make sure that all referenced columns are defined in the parameter set and that scenarios are unique. Currently available to the system are: ", names(ip$parameters) %>% paste(collapse = ", "), call. = FALSE)
  })

  return(invisible(ip))
}

#' compile the ODE expression for the isopath
get_ode_expression <- function(ip) {
  exp_text <- paste0("list(",
                     (ip %>% get_ode_matrix() %>% mutate(exp = paste(x, "=", `dx/dt`)))$exp %>%
                       paste(collapse = ", "), ")")
  interp(lazy(x), x = parse(text = exp_text, keep.source = F, n = NULL)[[1]])
}

#' get ode derivative
get_ode_function <- function(ip) {

  # lazy expression for optimized derivative calculation
  exp_lazy <- ip %>% get_ode_expression()

  # vector of variable component pools
  pools <- (ip %>% get_component_matrix()  %>% filter(variable == T))$component

  # derivatives function
  function(t, y, p) {

    # current ODE state
    ode_state <- c(as.list(y), p)

    # evaluate the expression for the variable dx/dt
    dx_var <- lazy_eval(exp_lazy, data = ode_state)

    # check if there are any missing values other computational problems
    if ( any(sapply(dx_var, length) == 0) || any(is.nan(unlist(dx_var))) )
      stop("some derivatives could not be computed: ",
           dx_var[sapply(dx_var, length) == 0 | is.nan(unlist(dx_var))] %>% names() %>% paste(collapse = ", "),
           " (state variables: ", paste(names(y), "=", y) %>% paste(collapse = ", "), ")",
           call. = FALSE)

    # make sure not running out of any component or running into infinity problems (with simple loop for speed)
    for (pool in pools){
      if (ode_state[[pool]] > 1e100)
        stop("overflowing pool: ", pool, " (size: ", ode_state[[pool]], " | dx/dt: ", signif(dx_var[[pool]],3), ")",
             ". Check your system to make sure that none of the component pools overflows (reaches a size of 10^100). ",
             "State variables: ", paste(names(y), "=", y) %>% paste(collapse = ", "),
             call. = FALSE)
      else if (ode_state[[pool]] + dx_var[[pool]] < 0)
        stop("depleted the following pool: ", pool, " (size: ", ode_state[[pool]], " | dx/dt: ", signif(dx_var[[pool]],3), ")",
             ". Please make sure that none of the component pools runs out completely. ",
             "State variables: ", paste(names(y), "=", y) %>% paste(collapse = ", "),
             call. = FALSE)
    }

    # all dx/dt, default value 0 for all non modified variables (=parameters)
    dx <- modifyList(setNames(as.list(rep(0, length(y))), names(y)), dx_var)

    # convert to correct list format
    return(list(unlist(dx)))
  }
}




#' Run the reaction model
#'
#' Run it for a certain number of time steps.
#'
#' @param time_steps the number of time steps, can be a number or expression (referring to a parameter in the isopath)
#' @param ... additional parameters passed on to the \link{ode} solver
#' @param make_state_var vector of parameters that should be included as state variables (so they can be changed as part of any special events passed to \code{...}, see \link{ode} for detail on the \code{events} parameter). All variable components' pools and isotopic compositions are always inclduded as state parameters and don't need to be added explicitly here.
#' @export
run_model <- function(ip, time_steps, ..., make_state_var = c()) {
  if (!is(ip, "isopath")) stop ("can only run model for an isopath", call. = FALSE)

  # model ready checks
  if (missing(time_steps))
    stop("time_steps is required when calling 'run_model' (either a number or an expression referencing a parameter)")
  check_model(ip)

  # info
  message(sprintf("Running model for %d scenario(s)...", nrow(ip$parameters)))

  # derivatives function
  func <- get_ode_function(ip)

  # steps evaluation
  steps_exp <- lazy(time_steps)

  # state variables
  state_vars <- c(get_ode_matrix(ip)$x, make_state_var) %>% unique()
  constants <- setdiff(names(ip$parameters), state_vars)

  # run each scenario by grouping by each row
  result <-
    ip$parameters %>%
    group_by_(.dots = names(ip$parameters)) %>%
    do({
      # time steps
      times <- seq(0, lazy_eval(steps_exp, data = .), by = 1)

      # attempt to solve ODE
      sln <- tibble()
      tryCatch({
        sln <- ode(y = unlist(.[state_vars]), times = times, func = func,
                   parms = as.list(.[constants]), ...)
      },
      error = function(e) {
        message("WARNING: encountered the following error while trying to solve one parameter set (skipping to the next set of parameters):\n<<< ", e$message, " >>>.\nParameters: ", paste(names(unlist(.)), "=", unlist(.)) %>% paste(collapse = ", "))
      })

      # output
      as.data.frame(sln)

    }) %>%
    ungroup()

  if (nrow(result) == 0) stop("None of the scenarios could be computed successfully", call. = FALSE)

  result %>% select_(.dots = c("time",c(state_vars, constants))) %>% return()
}

#' Run model to steady state
#'
#' Reports the time required to reach steady-state (in model time, with dt=1) and creates <variable>.t0 columns for each state variable. Uses \link{runsteady} internally for greatest flexibility although \link{stode} can work for some reaction systems and would be a little faster.
#'
#' @param ... additional parameters passed on to the steady state ODE (\link{runsteady}) solver, use \link{runsteady} parameter \code{verbose = TRUE} for additional output during computations, use parameters \code{rtol}, \code{atol} and \code{ctol} to adjust steady-state tolerances
#' @export
run_steady_state <- function(ip, ...) {
  if (!is(ip, "isopath")) stop ("can only run model for an isopath", call. = FALSE)

  # model ready checks
  check_model(ip)

  # info
  message(sprintf("Looking for steady-state for %d scenario(s)...", nrow(ip$parameters)))

  # ode function
  func <- get_ode_function(ip)

  # state variables
  state_vars <- get_ode_matrix(ip)$x
  state_vars_t0 <- paste0(state_vars, "_t0")
  constants <- setdiff(names(ip$parameters), state_vars)

  # run each scenario by grouping by each row
  result <-
    ip$parameters %>%
    group_by_(.dots = names(ip$parameters)) %>%
    do({
      sln <- tibble()
      # attempt to solve ODE
      tryCatch({
        timing <-
          system.time(
            out <- runsteady(y = unlist(.[state_vars]), func = func, times = c(0, Inf),
                         parms = as.list(.[constants]),  ...)
          )

        if (attributes(out)$steady) {
          sln <- bind_cols(
            tibble(time = attributes(out)$time),
            rename_(., .dots = setNames(state_vars, state_vars_t0))[state_vars_t0], # t0 values
            as_tibble(as.list(out$y)) # solutions
          )
        } else {
          message("WARNING: did not reach steady state (elapsed time: ", signif(timing[['elapsed']], 5), "s) for one parameter set (skipping to the next set of parameters). Run with verbose = TRUE for additional detail.\nParameters: ", paste(names(unlist(.)), "=", unlist(.)) %>% paste(collapse = ", "))
        }

      },
      error = function(e) {
        message("WARNING: encountered the following error while trying to find steady state for one reaction system (skipping to the next set of parameters):\n<<< ", e$message, " >>>.\nParameters: ",
                paste(names(unlist(.)), "=", unlist(.)) %>% paste(collapse = ", "))
      })

      # output
      sln

    }) %>%
    ungroup()

  if (nrow(result) == 0) stop("None of the scenarios could be run to steady-state.", call. = FALSE)

  result %>% select_(.dots = c("time", state_vars, state_vars_t0, constants)) %>% return()
}


