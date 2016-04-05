#' Create a reaction for a pathway
#' @param name of the component (use a vector to add multiple components with the same isotopes and \code{variable} attribute)
#' @param variable - whether this component has variable concentration and isotopic composition, default is that it can vary
#' @param ... the isotopes (refer to by name)
#' @note non-standard evaluation
#' @export
add_component <- function(ip, name, ..., variable = TRUE) {
  for (component in name) {
    ip <- add_component_(ip, component, variable = variable,
                   .dots = lazy_dots(...) %>% lapply(function(i) deparse(i$expr)))
  }
  return(invisible(ip))
}

#' add isotopes to a component
#' @export
#' @note standard evaluation
add_component_ <- function(ip, name, ..., .dots = list(), variable = FALSE) {
  if (!is(ip, "isopath")) stop ("component can only be added to an isopath", call. = FALSE)
  if (!grepl("^\\w+$", name)) stop("only alphanumeric component names allowed: ", name, call. = FALSE)
  ip$components[[name]] <- list(
    name = name,
    variable = variable,
    isotopes = sapply(c(list(...), .dots), parse_reaction_component)
  )
  missing <- setdiff(ip$components[[name]]$isotopes %>% names(), names(ip$isotopes))
  if (length(missing) > 0)
    stop("missing isotope definition(s), make sure to add this with add_isotope() first: ",
         missing %>% paste(collapse = ", "), call. = FALSE)

  # store info
  ip <- ip %>% store_info()

  return(invisible(ip))
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
