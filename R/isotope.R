#' Add an isotope
#' @export
add_isotope <- function(ip, name) {
  if (!is(ip, "isopath")) stop ("isotope can only be added to an isopath")
  ip$isotopes[[name]] <- list()

  return(ip)
}
