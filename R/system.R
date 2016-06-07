#' Create a pathway for isotope fractionation
#'
#' @export
isopath <- function() {
  structure(list(
    isotopes = list(),
    components = list(),
    reactions = list(),
    parameters = data_frame(),
    info = list()
  ), class = "isopath")
}

#' get the component matrix for an isopath
#' @param na how to group compounds that have no isotopes specified
#' @family system information
#' @export
get_component_matrix <- function(ip, na = "unspecified") {
  if (!is(ip, "isopath")) stop ("can only get component matrix from an isopath")
  lapply(ip$components, function(i) {
    if ( length(i$isotopes) == 0 ) i$isotopes <- setNames(1, na)
    c(list(component = i$name, variable = i$variable), as.list(i$isotopes)) %>% as_data_frame()
  }) %>%
    bind_rows()
}

#' get the reaction matrix for an isopath
#'
#' @param evaluate if TRUE, tries to evaluate the expressions stored for the different fluxes (requires parameters provided)
#' @param parameters data to use for evaluating expressions
#' @family system information
#' @note TODO: write tests and replace get_reaction_matrix
#' @export
get_reaction_matrix2 <- function(ip, evaluate = FALSE, parameters = ip$parameters[1,]) {
  if (!is(ip, "isopath")) stop ("can only get reaction matrix from an isopath")

  lapply(ip$reactions, function(i) {
    c(list(reaction = i$name, abscissa = i$abscissa),
      as.list(i$components),
      list(flux =
             if (evaluate) {
               lazy_eval(i$flux, data = parameters)
             } else {
               deparse(i$flux$expr, width.cutoff = 500L) %>% paste0(collapse = "")
             }))%>%
      as_data_frame()
  }) %>%
    bind_rows()
}

#' get the combined reaction_component_matrix
#'
#' if \code{evaluate = TRUE}, excludes components that are not variable
#'
#' @inheritParams get_reaction_matrix
#' @family system information
#' @note TODO: write tests and replace get_reaction_component_matrix
#' @export
get_reaction_component_matrix2 <- function(ip, evaluate = FALSE, parameters = ip$parameters[1,]) {
  if (!is(ip, "isopath")) stop ("can only calculate reaction component matrix from an isopath")

  cps <- ip %>% get_component_matrix()
  rxn <- ip %>% get_reaction_matrix2(eval = evaluate, param = parameters)

  if (nrow(cps) == 0 || nrow (rxn) == 0) return(data_frame())

  # pool size for reaction component
  get_size <- function(component) {
    eq <- interp(lazy(pool_), pool_ = as.name(component))
    if (evaluate) lazy_eval(eq, parameters)
    else deparse(eq$expr, width.cutoff = 500L) %>% paste0(collapse = "")
  }

  # dx/dt for reaction component
  get_change <- function(rxn, comp_stoic, variable) {
    if (!variable) return(0)
    eq <- interp(lazy(stoi_ * flux_), stoi_ = comp_stoic, flux_ = ip$reactions[[rxn]]$flux)
    if (evaluate) lazy_eval(eq, parameters)
    else deparse(eq$expr, width.cutoff = 500L) %>% paste0(collapse = "")
  }

  # construct data frame
  df <- left_join(
    rxn %>% gather(component, comp_stoic, -reaction, -abscissa, -flux),
    cps %>% select(component, variable),
    by = "component") %>%
    # remove listings that don't have the component
    filter(!is.na(comp_stoic))

  # exclude invariable components if evaluating
  if (evaluate)
    df <- df %>% filter(variable)

  df %>%
    # get pool size and dx/dt
    mutate(
      pool_size = mapply(get_size, component),
      `dx/dt` = mapply(get_change, reaction, comp_stoic, variable)) %>%
    # calculate abscissa for each component
    mutate(abscissa = abscissa - ifelse(comp_stoic < 0, 1, 0 )) %>%
    # arrange by component, then abscissa
    arrange(component, abscissa, reaction) %>%
    # order columns
    select(component, abscissa, variable, reaction, comp_stoic, flux, pool_size, `dx/dt`)
}

#' get the combined reaction_isotope_matrix
#'
#' if \code{evaluate = TRUE}, excludes components that are not variable
#'
#' @inheritParams get_reaction_matrix
#' @family system information
#' @note TODO: write tests and replace all the get_flux... functions!!
#' @export
get_reaction_isotope_matrix <- function(ip, evaluate = FALSE, parameters = ip$parameters[1,]) {
  if (!is(ip, "isopath")) stop ("can only calculate reaction isotope matrix from an isopath")

  # component and reaction matrices
  cps <- ip %>% get_component_matrix()
  rxn <- ip %>% get_reaction_matrix2(eval = evaluate, param = parameters)

  if (nrow(cps) == 0 || nrow (rxn) == 0) return(data_frame())

  # pool size and isotopic composition
  get_size <- function(component) {
    eq <- interp(lazy(pool_), pool_ = as.name(component))
    if (evaluate) lazy_eval(eq, parameters)
    else deparse(eq$expr, width.cutoff = 500L) %>% paste0(collapse = "")
  }
  get_isotope <- function(component, isotope) {
    eq <- interp(lazy(pool_iso_), pool_iso_ = as.name(paste0(component, ".", isotope)))
    if (evaluate) lazy_eval(eq, parameters)
    else deparse(eq$expr, width.cutoff = 500L) %>% paste0(collapse = "")
  }

  # flux isotope
  get_flux_isotope_eq <- function(rxn, component, isotope) {
    iso <- ip$reactions[[rxn]]$isotopes[[isotope]]
    if (iso %>% is("lazy")) return(iso) # flux isotopes for all components
    else if (is.null(iso[[component]])) return(interp(lazy(x), x = as.name(paste0(component, ".", isotope)))) # no isotope effect
    else return(iso[[component]]) # component specific
  }
  get_flux_isotope <- function(rxn, component, isotope) {
    eq <- get_flux_isotope_eq(rxn, component, isotope)
    if (evaluate) lazy_eval(eq, parameters)
    else deparse(eq$expr, width.cutoff = 500L) %>% paste0(collapse = "")
  }

  # get dx/dt for isotope
  get_isotope_change <- function(rxn, component, isotope, comp_stoic, variable) {
    if (!variable) return(0)
    iso <- get_flux_isotope_eq(rxn, component, isotope)
    eq <- interp(lazy(stoi_ * flux_ / pool_ * (flux_iso_ - pool_iso_)),
                 stoi_ = comp_stoic, flux_ = ip$reactions[[rxn]]$flux,
                 pool_ = as.name(component),
                 flux_iso_ = iso,
                 pool_iso_ = as.name(paste0(component, ".", isotope)))
    if (evaluate) lazy_eval(eq, parameters)
    else deparse(eq$expr, width.cutoff = 500L) %>% paste0(collapse = "")
  }

  # construct data frame
  df <- left_join(
    rxn %>% gather(component, comp_stoic, -reaction, -abscissa, -flux),
    cps %>% gather(isotope, iso_stoic, -component, -variable),
    by = "component"
  ) %>%
    # remove listings that don't have component or isotope
    filter(!is.na(comp_stoic), !is.na(iso_stoic))

  # exclude invariable components if evaluating
  if (evaluate)
    df <- df %>% filter(variable)

  df %>%
    # flux and pool isotopic composition, plus dx/dt
    mutate(
      pool_size = mapply(get_size, component),
      pool_isotope = mapply(get_isotope, component, isotope),
      flux_isotope = mapply(get_flux_isotope, reaction, component, isotope),
      `dx/dt` = mapply(get_isotope_change, reaction, component, isotope, comp_stoic, variable)) %>%
    # calculate abscissa for each component
    mutate(abscissa = abscissa - ifelse(comp_stoic < 0, 1, 0 )) %>%
    # arrange
    arrange(isotope, abscissa, component, desc(comp_stoic)) %>%
    # order columns
    select(isotope, component, abscissa, variable, reaction, comp_stoic,
           flux, flux_isotope, pool_size, pool_isotope, `dx/dt`)
}

#' get the ODE
#'
#' Get the system of ordinary differential equations
#'
#' @inheritParams get_reaction_matrix
#' @family system information
#' @note TODO: write tests and replace all the get_flux... functions!!
#' @export
get_ode_matrix <- function(ip, evaluate = FALSE, parameters = ip$parameters[1,]) {
  if (!is(ip, "isopath")) stop ("can only calculate ode matrix from an isopath")

  # value
  get_value <- function(component) {
    eq <- interp(lazy(pool_), pool_ = as.name(component))
    if (evaluate) {
      value <- lazy_eval(eq, parameters)
      if (length(value) == 1) return(value)
      else if (length(value) == 0) return(NA)
      else stop("should not be possible to happen")
    } else deparse(eq$expr, width.cutoff = 500L) %>% paste0(collapse = "")
  }

  # get net isotope change
  get_net_change <- function(eq_text) {
    if (evaluate) {
      eq <- interp(lazy(x), x = parse(text = eq_text, keep.source = F, n = NULL)[[1]])
      value <- lazy_eval(eq, parameters)
      if (length(value) == 1) return(value)
      else if (length(value) == 0) return(NA)
      else stop("should not be possible to happen")
    } else return(eq_text)
  }

  bind_rows(
    ip %>% get_reaction_component_matrix2(eval = F) %>%
      select(x = component, variable, `dx/dt`),
    ip %>% get_reaction_isotope_matrix(eval = F) %>%
      mutate(x = paste0(component, ".", isotope)) %>%
      select(x, variable, `dx/dt`)
  ) %>%
    filter(variable == TRUE) %>%
    group_by(x) %>%
    summarize(eqn_net = paste(`dx/dt`, collapse = " + ")) %>%
    mutate(
      value = mapply(get_value, x) %>% unname(),
      `dx/dt` = mapply(get_net_change, eqn_net) %>% unname()
    ) %>%
    select(x, value, `dx/dt`)
}

#' @export
print.isopath <- function(x, ...) {
  # implement this eventually
  sprintf("Isopath [%d components, %d reactions, %d parameter sets]:\n",
          length(x$components), length(x$reactions), nrow(x$parameters)) %>% cat()
  cat("\nCOMPONENTS - ")
  print(x %>% get_component_matrix())
  cat("\nREACTIONS - ")
  print(x %>% get_reaction_matrix2() %>% select(-flux))
  cat("\nORDINARY DIFFERENTIAL EQUATIONS - ")
  odes <- x %>% get_ode_matrix(eval = F) %>% select(x, `dx/dt`)
  max_chars <- getOption("width") - max(nchar(odes$x))-15
  odes %>% mutate(
    `dx/dt` = substr(`dx/dt`, 1, max_chars) %>%
      paste0( ifelse(nchar(`dx/dt`) > max_chars, "...", ""))) %>%
    print()
  cat("\nPARAMETERS - ")
  print(x$parameters)
}

#' Merge multiple isopaths
#' @export
merge_isopaths <- function(...) {
  stop("not implemented yet")
}



