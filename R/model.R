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
        get_component_change_summary(ip, param = ., check_missing = TRUE)
        get_isotope_change_summary(ip, param = ., check_missing = TRUE)
      })

  },

  # error handling
  error = function(e) {
    stop("encountered the following error during pre-check of the model run: >>> ", e$message, " <<< Please make sure that all referenced columns are defined in the parameter set and that scenarios are unique. Currently available to the system are: ", names(ip$parameters) %>% paste(collapse = ", "))
  })

  return(invisible(ip))
}

#' Run the reaction model
#'
#' Run it for a certain number of time steps
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

  # derivatives function
  calc_derivs <- function(t, y, params) {

    # get ODE solutions
    ode_state <- (c(as.list(y), params))
    component_change <- ip %>% get_component_change_summary(param = ode_state, check_missing = FALSE)
    isotope_change <- ip %>% get_isotope_change_summary(param = ode_state, check_missing = FALSE)

    # make sure not running out of any component
    if (any(with(component_change, pool_size + `dx/dt` < 0))) {
      run_out <- component_change %>%
        filter(pool_size + `dx/dt` < 0) %>%
        mutate(label = paste0(component, " (pool: ", signif(pool_size, 3), " | dx/dt: ", signif(`dx/dt`, 3), ")"))
      stop("depleted the following pool(s): ", run_out$label %>% paste(collapse = ", "),
           ". Please make sure that none of the component pools runs out completely.", call. = FALSE)
    }

    # all values (combined changes in pool size and isotopic composition)
    dxdt <-
      bind_rows(
        component_change %>% rename(variable = component),
        isotope_change %>% mutate(variable = paste0(component, ".", isotope)) %>%
          select(variable, `dx/dt`)
      )

    # put the dx in the right order and fill in the constants (0s)
    dx <- c()
    for (name in names(y)) {
      i <- which(name == dxdt$variable)
      if (length(i) == 1)
        dx <- c(dx, dxdt$`dx/dt`[i])
      else
        dx <- c(dx, 0) # non-variable parameters
    }

    # convert to list
    return(list(dx))
  }

  # variables and constants in the system (i.e. non-variable parameters)
  # In order to allow events to affect parameters as well as constants, treating
  # them all as variable with dx/dts evaluating to 0 for all parameters.
  variables <- names(ip$parameters)
  constants <- c()

  # info
  if (!quiet) message(sprintf("Running model for %d scenarios...", nrow(ip$parameters)))

  # run each scenario by grouping by each row
  ip$parameters %>%
    group_by_(.dots = c(variables, constants)) %>%
    do({
      # time steps
      times <- seq(0, lazy_eval(steps_exp, data = .[constants]), by = 1)

      # attempt to solve ODE
      sln <- data_frame()
      tryCatch({
        sln <- ode(y = .[variables] %>% unlist(),
                   times = times, func = calc_derivs,
                   parms = .[constants] %>% as.list(), ...)
      },
      error = function(e) {
        message("WARNING: encountered the following error while trying to solve one reaction system (skipping to the next set of parameters): ", e$message)
      })

      # output
      sln %>% as.data.frame() %>% return()

    }) %>%
    ungroup() %>%
    select_(.dots = c("time", c(variables, constants)))
}
