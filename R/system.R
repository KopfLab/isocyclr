#' Create a pathway for isotope fractionation
#'
#' @export
isopath <- function() {
  structure(list(
    isotopes = list(),
    components = list(),
    reactions = list()
  ), class = "isopath")
}

#' @export
print.isopath <- function(x, ...) {
  # implement this eventually
  print.default(x)
}

#' Merge multiple isopaths
#' @export
merge_isopaths <- function(...) {
  stop("not implemented yet")
}



