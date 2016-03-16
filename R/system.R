#' Create a pathway for isotope fractionation
#'
#' @export
isopath <- function() {
  structure(list(
    isotopes = list(),
    components = list(),
    reactions = list(),
    parameters = list()
  ), class = "isopath")
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



