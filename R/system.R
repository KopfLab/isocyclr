#' Create a pathway for isotope fractionation
#'
#' @export
isopath <- function() {
  structure(list(
    isotopes = list(),
    components = list(),
    reactions = list(),
    parameters = data_frame()
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



