#' Add an isotope
#' @param name the name of the isotope to add (alpha numeric names only), can be a vector to add multiple isotopes at once
#' @export
add_isotope <- function(ip, name) {
  if (!is(ip, "isopath")) stop ("isotope can only be added to an isopath")
  if (any(!grepl("^\\w+$", name))) stop("only alphanumeric isotope names allowed: ",
                                        name %>% paste(collapse = ", "))
  for (isotope in name) {
    ip$isotopes[[isotope]] <- list()
  }

  # store info
  ip <- ip %>% store_info()

  return(ip)
}
