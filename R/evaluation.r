#' retrieve flux matrix
#' @param evaluate if TRUE, tries to evaluate the expressions stored for the different fluxes
#' @param parameters
#' @export
get_flux_matrix <- function(ip, evaluate = FALSE, parameters = ip$parameters){
  if (!is(ip, "isopath")) stop ("can only get flux matrix from an isopath", call. = F)
  lapply(ip$reactions, function(i) {
    if (evaluate) lazy_eval(i$flux, parameters)
    else deparse(i$flux$expr)
  }) %>% as_data_frame() %>% gather(reaction, flux)
}

#' retrieve flux isotopes matrix
#' @param evaluate if TRUE, tries to evaluate the expressions stored for the different fluxes
#' @param parameters only needed if \code{evaluate=TRUE}, data frame with all the parameters needed to evaluate the flux isotope expressions
#' @export
get_flux_isotope_matrix <- function(ip, evaluate = FALSE, parameters = ip$parameters) {
  if (!is(ip, "isopath")) stop ("can only get flux isotope matrix from an isopath", call. = F)
  df <- list()
  for (rxn in ip$reactions) {
    if (length(rxn$isotopes) > 0) {
      df[[rxn$name]] <-
        mapply(function(rxn, iso_name, iso) {
          retval <-
            if (iso %>% is("lazy")) { # flux isotopes for all components
              data_frame(
                component = "<all>",
                flux = if(evaluate) lazy_eval(iso, parameters) else deparse(iso$expr))
            } else
              lapply(iso, function(iso.component) # flux isotopes for individual components
                if (evaluate) lazy_eval(iso.component, parameters)
                else deparse(iso.component$expr)) %>%
            as_data_frame() %>% gather(component, flux)
          return(mutate(retval, isotope = iso_name, reaction = rxn))
        },
        rxn$name, names(rxn$isotopes), rxn$isotopes, SIMPLIFY = F, USE.NAMES = F)
    }
  }

  # convert to data frame
  df <- df %>% unlist(recursive = F) %>% bind_rows()

  # sort if any recovered
  if (nrow(df) > 0)
    df <- df %>% select(reaction, isotope, component, flux)

  return(df)
}
