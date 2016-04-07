#' add parameter data frame to the isopath
#' @param ... parameters, either a whole set of parameters as a single data frame, or a list of individual values, will try to merge with existing parameters and if neceesary overwrite existing parameters with the same name
#' @export
set_parameters <- function(ip, ...) {
  if (!is(ip, "isopath")) stop ("parameters can only be set for an isopath")

  ldots <- lazy_dots(...)
  if (length(list(...)) == 1 && is.data.frame(..1)) {
    # overwrite if single data frame passed in
    ip$parameters <- ..1
  } else if (nrow(ip$parameters) == 0) {
    # overwrite if nothing set yet
    ip$parameters <- data_frame(...)
  } else {
    # mutate otherwise
    tryCatch({
      ip$parameters <- mutate_(ip$parameters, .dots = ldots)
    },
    error = function(e) {
      stop("something went wrong trying to merge the new parameters with the existing ones: '", e$message, "'.", call. = FALSE)
    })
  }

  return(invisible(ip))
}


#' @rdname set_parameters
#' @export
get_parameters <- function(ip) {
  if (!is(ip, "isopath")) stop ("can only get parameters from an isopath")
  ip$parameters
}
