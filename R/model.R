#' check if system is ready for model run
check_model <- function(ip) {
  if (!is(ip, "isopath")) stop ("can only check model ready for an isopath", call. = FALSE)

  tryCatch({

    # parameters
    missing_params <- setdiff(ip$info$variables, names(ip$parameters))
    if (length(missing_params) > 0)
      stop("parameters required for minimal parameter set missing: ",
           missing_params %>% paste(collapse = ", "), call. = FALSE)

    # do a sample calculation for each parameter set
    ip$parameters %>%
      group_by_(.dots = names(ip$parameters)) %>%
      do({
        # safety checks - will through errors if needed
        if (nrow(.) != 1) stop("there seem to be multiple identical run scenarios, please make sure each set of parameters is unique", call. = FALSE)
        dx_var <- get_ode_matrix(ip, eval = T, param = .)$`dx/dt`
        if (any(is.na(dx_var))) {
          stop("some derivatives could not be computed: ", dx_var[is.na(dx_var)] %>% names() %>% paste(collapse = ", "), call. = FALSE)
        }
        data_frame()
      })

  },

  # error handling
  error = function(e) {
    stop("encountered the following error during pre-check of the model run: >>> ", e$message, " <<< Please make sure that all referenced columns are defined in the parameter set and that scenarios are unique. Currently available to the system are: ", names(ip$parameters) %>% paste(collapse = ", "), call. = FALSE)
  })

  return(invisible(ip))
}


#' Run the reaction model
#'
#' Run it for a certain number of time steps.
#'
#' This will purposefully include the entire parameter set in the state variables to allow for events to change
#' any parameter at specific time points or in response to functional values
#'
#' @param time_steps the number of time steps, can be a number or expression (referring to a parameter in the isopath)
#' @param ... additional parameters passed on to the \link{ode} solver
run_model <- function(ip, time_steps, ..., quiet = FALSE) {
  if (!is(ip, "isopath")) stop ("can only run model for an isopath", call. = FALSE)

  # model ready checks
  if (missing(time_steps))
    stop("time_steps is required when calling 'run_model' (either a number or an expression referencing a parameter)")
  check_model(ip)

  # steps evaluation
  steps_exp <- lazy(time_steps)

  # construct lazy expression for optimized derivative calculation
  exp_text <- paste0("list(",
                    (ip %>% get_ode_matrix() %>% mutate(exp = paste(x, "=", `dx/dt`)))$exp %>%
                      paste(collapse = ", "), ")")
  exp_lazy <- interp(lazy(x), x = parse(text = exp_text, keep.source = F, n = NULL)[[1]])

  # vector of variable component pools
  pools <- (ip %>% get_component_matrix()  %>% filter(variable == T))$component

  # derivatives function
  calc_derivs <- function(t, y, p) {

    # current ODE state
    ode_state <- as.list(y)

    # evaluate the expression for the variable dx/dt
    dx_var <- lazy_eval(exp_lazy, data = ode_state)

    # check if there are any missing values
    if ( any(sapply(dx_var, length) == 0) )
      stop("some derivatives could not be computed: ",
           dx_var[sapply(dx_var, length) == 0] %>% names() %>% paste(collapse = ", "), call. = FALSE)

    # make sure not running out of any component (with simple loop for speed)
    for (pool in pools){
      if (ode_state[[pool]] + dx_var[[pool]] < 0)
        stop("depleted the following pool: ", pool, " (size: ", ode_state[[pool]], " | dx/dt: ", signif(dx_var[[pool]],3), ")",
             ". Please make sure that none of the component pools runs out completely.", call. = FALSE)
    }

    # all dx/dt, default value 0 for all non modified variables (=parameters)
    dx <- modifyList(
      setNames(as.list(rep(0, length(ode_state))), names(ode_state)),
      dx_var)

    # convert to correct list format
    return(list(unlist(dx)))
  }

  # info
  if (!quiet) message(sprintf("Running model for %d scenario(s)...", nrow(ip$parameters)))

  # run each scenario by grouping by each row
  result <-
    ip$parameters %>%
    group_by_(.dots = names(ip$parameters)) %>%
    do({
      # time steps
      times <- seq(0, lazy_eval(steps_exp, data = .), by = 1)

      # attempt to solve ODE
      sln <- data_frame()
      tryCatch({
        sln <- ode(y = unlist(.), times = times, func = calc_derivs, ...)
      },
      error = function(e) {
        message("WARNING: encountered the following error while trying to solve one reaction system (skipping to the next set of parameters): ", e$message)
      })

      # output
      return(as.data.frame(sln))

    }) %>%
    ungroup()

  if (nrow(result) == 0) stop("None of the scenarios could be computed successfully", call. = FALSE)

  result %>% select_(.dots = c("time", names(ip$parameters))) %>% return()
}
