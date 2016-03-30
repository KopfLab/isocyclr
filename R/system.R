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

#' Get the dynamic variables in the system
#' @export
get_variables <- function(ip){
  if (!is(ip, "isopath")) stop ("need an isopath to generate a list of dynamics variables")

  ip$components %>%
    lapply(function(i) {
      if (i$variable) # only variable components are included
        c(i$name, paste0(i$name, ".", names(i$isotopes)))
      else
        c()
    }) %>%
    unlist() %>% unname()
}

#' store the information about the system in the info list
#' (this is done for performance reasons to reduce time
#' intenstive calculations during ode integrations)
store_info <- function(ip) {
  stopifnot(is(ip, "isopath"))
  ip$info$reaction_component_matrix <- ip %>% get_reaction_component_matrix()
  ip$info$variable_reaction_component_matrix <- ip$info$reaction_component_matrix %>% filter(variable == T)
  ip$info$variables <- ip %>% get_variables()
  return(ip)
}

#' @export
print.isopath <- function(x, ...) {
  # implement this eventually
  sprintf("Isopath [%d components, %d reactions]:\n", length(x$components), length(x$reactions)) %>% cat()
  cat("\nCOMPONENTS - ")
  print(x %>% get_component_matrix())
  cat("\nREACTIONS - ")
  print(x %>% get_reaction_matrix())
}

#' Merge multiple isopaths
#' @export
merge_isopaths <- function(...) {
  stop("not implemented yet")
}



