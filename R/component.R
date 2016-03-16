#' Create a reaction for a pathway
#' @param fixed - whether this component has fixed concentration and isotopic composition, default is that it can vary
#' @param ... the isotopes (refer to by name)
#' @note non-standard evaluation
#' @export
add_component <- function(ip, name, ..., fixed = FALSE) {
  add_component_(ip, name, fixed = fixed, .dots = lazy_dots(...) %>% lapply(function(i) deparse(i$expr)))
}

#' add isotopes to a component
#' @export
#' @note standard evaluation
add_component_ <- function(ip, name, ..., .dots = list(), fixed = FALSE) {
  if (!is(ip, "isopath")) stop ("component can only be added to an isopath", stop. = FALSE)
  if (!grepl("^\\w+$", name)) stop("only alphanumeric component names allowed: ", name, stop. = FALSE)
  ip$components[[name]] <- list(
    name = name,
    fixed = fixed,
    isotopes = sapply(c(list(...), .dots), parse_reaction_component)
  )
  missing <- setdiff(ip$components[[name]]$isotopes %>% names(), names(ip$isotopes))
  if (length(missing) > 0)
    stop("missing isotope definition(s), make sure to add this with add_isotope() first: ",
         missing %>% paste(collapse = ", "), stop. = FALSE)
  return(ip)
}

#' get the component matrix for an isopath
#' @param na how to group compounds that have no isotopes specified
#' @export
get_component_matrix <- function(ip, na = "unspecified") {
  if (!is(ip, "isopath")) stop ("can only get component matrix from an isopath")
  lapply(ip$components, function(i) {
    if ( length(i$isotopes) == 0 ) i$isotopes <- setNames(1, na)
    c(list(component = i$name), as.list(i$isotopes)) %>% as_data_frame()
  }) %>%
  bind_rows()
}
