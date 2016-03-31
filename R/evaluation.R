#' retrieve flux matrix
#' @param evaluate if TRUE, tries to evaluate the expressions stored for the different fluxes (requires parameters provided)
#' @param parameters data to use for evaluating expressions
#' @export
get_flux_matrix <- function(ip, evaluate = FALSE, parameters = ip$parameters[1,]){
  if (!is(ip, "isopath")) stop ("can only get flux matrix from an isopath", call. = F)
  lapply(ip$reactions, function(i) {
    if (evaluate) lazy_eval(i$flux, parameters)
    else deparse(i$flux$expr)
  }) %>% as_data_frame() %>% gather(reaction, flux)
}

#' get component changes (dx/dt) for each reaction and component with the parameters provided
#' @param check_missing whether to check for any missing flux components
#' @export
get_component_change_matrix <- function(ip, parameters = ip$parameters[1,], check_missing = TRUE) {
  if (!is(ip, "isopath")) stop ("can only get component change for an isopath", call. = F)
  component_flux <-
    ip$info$variable_reaction_component_matrix %>%
    select(reaction, component, comp_stoic) %>%
    distinct() %>%
    # join in reaction
    left_join(
      ip %>% get_flux_matrix(eval = T, param = parameters),
      by = c("reaction")
    ) %>%
    # join in data (not actually needed for computation but good for safety checks)
    left_join(
      as_data_frame(parameters)[ip$info$variables] %>%
        gather(component, pool_size),
      by = c("component")
    ) %>%
    # calculate actual stoichiometry weighted dx/dt
    mutate(`dx/dt` = comp_stoic * flux)

  # safety checks
  if (check_missing) {
    missing_flux <- component_flux %>% filter(is.na(`dx/dt`))
    if (nrow(missing_flux) > 0) {
      stop("missing mass flux for the following reaction components: ",
           paste(missing_flux$reaction, " ", missing_flux$component) %>% paste(collapse = ", "), call. = F)
    }
  }

  return(component_flux)
}


#' get the summary of component changes (dx/dt) for the parameters provided
#' @export
get_component_change_summary <- function(ip, parameters = ip$parameters[1,], ...) {
  if (!is(ip, "isopath")) stop ("can only get component change for an isopath", call. = F)
  ip %>%
    get_component_change_matrix(param = parameters, ...) %>%
    group_by(component, pool_size) %>%
    summarize(`dx/dt` = sum(`dx/dt`)) %>%
    ungroup()
}

#' retrieve flux isotopes matrix
#' @param evaluate if TRUE, tries to evaluate the expressions stored for the different fluxes
#' @param parameters only needed if \code{evaluate=TRUE}, data frame with all the parameters needed to evaluate the flux isotope expressions
#' @export
get_flux_isotope_matrix <- function(ip, evaluate = FALSE, parameters = ip$parameters[1,]) {
  if (!is(ip, "isopath")) stop ("can only get flux isotope matrix from an isopath", call. = F)
  df <- list()

  for (rxn in ip$reactions) {
    if (length(rxn$isotopes) > 0) {

      df[[rxn$name]] <-
        mapply(function(iso_name, iso) {
          retval <-
            if (iso %>% is("lazy")) { # flux isotopes for all components
              data_frame(
                component =
                  filter(ip$info$reaction_component_matrix, reaction == rxn$name, isotope == iso_name)$component,
                flux_isotope =
                  if(evaluate) lazy_eval(iso, parameters) else deparse(iso$expr))
            } else
              lapply(iso, function(iso.component) # flux isotopes for individual components
                if (evaluate) lazy_eval(iso.component, parameters)
                else deparse(iso.component$expr)
              ) %>% as_data_frame() %>% gather(component, flux_isotope)
          return(mutate(retval, isotope = iso_name, reaction = rxn$name))
        },
        names(rxn$isotopes), rxn$isotopes, SIMPLIFY = F, USE.NAMES = F)

    }
  }

  # convert to data frame
  df <- df %>% unlist(recursive = F) %>% bind_rows()

  # arrange columns if any recovered
  if (nrow(df) > 0)
    df <- df %>% select(reaction, isotope, component, flux_isotope)

  return(df)
}

#' get isotope change matrix by reaction and component + isotope
#' @param check_missing whether to check for any missing flux components
#' @export
get_isotope_change_matrix <- function(ip, parameters = ip$parameters[1,], check_missing = TRUE) {
  if (!is(ip, "isopath")) stop ("can only get isotope change for an isopath", call. = F)

  # get just data values from the parameters
  data <- as_data_frame(parameters)[ip$info$variables]

  # combine information
  flux_isotopes <-
    ip$info$variable_reaction_component_matrix %>%
    select(-abscissa, -variable) %>%
    # join in the flux matrix to figure out flux values
    # NOTE: this + data join below is actually faster than joining in get_flux_change_summary! tested with microbenchmark
    left_join(
      ip %>% get_flux_matrix(eval = T, param = parameters),
      by = c("reaction")
    ) %>%
    # join in the flux isotope matrix to figure out flux isotope values
    left_join(
      ip %>% get_flux_isotope_matrix(eval = T, param = parameters),
      by = c("reaction", "component", "isotope")
    ) %>%
    # join in the actual data on existing pool sizes
    left_join(
      data %>% gather(component, pool_size),
      by = c("component")
    ) %>%
    # join in the actual data on existing isotope values
    mutate(comp.isotope = paste0(component, ".", isotope)) %>%
    left_join(
      data %>% gather(comp.isotope, pool_isotope),
      by = c("comp.isotope")
    ) %>%
    # calculate dx/dt
    # Note: iso_stoic goes both in numerator and denominator so cancels out!
    mutate(
      `dx/dt` = comp_stoic * flux / pool_size * (flux_isotope - pool_isotope)
    )

  # safety checks
  if (check_missing) {
    missing_flux <- flux_isotopes %>% filter(is.na(`dx/dt`))
    if (nrow(missing_flux) > 0) {
      stop("missing isotope flux for the following reaction + component + isotopes: ",
           paste(missing_flux$reaction, "+", missing_flux$component, "+",
                 missing_flux$isotope) %>% paste(collapse = ", "), call. = F)
    }
  }

  return(flux_isotopes %>% select(-iso_stoic, -comp.isotope))
}

#' get the isotope change summary
#' @export
get_isotope_change_summary <- function(ip, parameters = ip$parameters[1,], ...) {
  if (!is(ip, "isopath")) stop ("can only get isotope change for an isopath", call. = F)
  ip %>%
    get_isotope_change_matrix(param = parameters, ...) %>%
    group_by(isotope, component, pool_isotope) %>%
    summarize(`dx/dt` = sum(`dx/dt`)) %>%
    ungroup()
}
