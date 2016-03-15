#' Create a reaction for a pathway
#' @param ... the isotopes (refer to by name)
#' @export
add_component <- function(ip, name, ...) {
  if (!is(ip, "isopath")) stop ("component can only add be added to an isopath")
  ip$components[[name]] <- parse_component_isotopes(...)
  return(ip)
}

#' parse component isotope
parse_component_isotopes <- function(...) {
  lazyeval::lazy_dots(...) %>%
    sapply(function(i) {
      iso <- deparse(i$expr)
      m <- regexec("^(\\d*)( \\* )?(.+)$", iso)
      parts <- regmatches(iso, m)[[1]]
      if ( length(parts) == 0)
        stop("Can't parse component isotope: ", iso, call. = F)
      list(
        if (parts[2] == "") 1 else as.numeric(parts[2])) %>%
        setNames(parts[4])
    })
}
