#' add parameter data frame to the isopath
#' @param ... parameters, either a whole set of parameters as a single data frame, or a list of individual values, will try to merge with existing parameters and if neceesary overwrite existing parameters with the same name
#' @export
set_parameters <- function(ip, ...) {
  if (!is(ip, "isopath")) stop ("parameters can only be set for an isopath")

  args <- list(...)
  if (length(args) == 1 && is.data.frame(..1)) params <- ..1
  else params <- args %>% as_data_frame()

  if (nrow(ip$parameters) > 0) {
    # attempt to merge
    tryCatch(
      ip$parameters <- bind_cols(ip$parameters[!names(ip$parameters) %in% names(params)], params),
      error = function(e) {
        stop("something went wrong trying to merge the new parameters with the existing ones: '", e$message, "'.")
      })
  } else
    ip$parameters <- params

  return(ip)
}

#' get the parameter data frame
#' @export
get_parameters <- function(ip) {
  if (!is(ip, "isopath")) stop ("can only get parameters from an isopath")
  ip$parameters
}
