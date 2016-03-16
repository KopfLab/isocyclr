#' Add an isotope
#' @export
add_isotope <- function(ip, name) {
  if (!is(ip, "isopath")) stop ("isotope can only be added to an isopath")
  if (!grepl("^\\w+$", name)) stop("only alphanumeric isotope names allowed: ", name)
  ip$isotopes[[name]] <- list()

  return(ip)
}
