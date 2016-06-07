#' add parameter data frame to the isopath
#'
#' to expand on an existing parameter set with multiple variations, use \code{\link{expand_parameters}}
#'
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

#' expand existing parameter set
#'
#' this expands each existing parameter set with a combinatoric
#' expansion of the provided parameters. provided parameters can
#' be numeric vectors or expressions that can reference existing parameters
#' @param ... named vectors or expressions to expand existing parameter sets
#' @export
expand_parameters <- function(ip, ...) {

  if (!is(ip, "isopath")) stop ("parameters can only be expanded for an isopath")
  if (nrow(ip$parameters) == 0) stop("no parameters set yet for this isopath")

  ldots <- lazy_dots(...)
  all_params <- names(ip$parameters)

  ip$parameters <-
    ip$parameters %>%
    # group by everything to get unique sets
    group_by_ (.dots = all_params) %>%
    # expand for each row
    do({
      if (nrow(.) != 1) stop("existing parameter sets are not unique and cannot be expanded, please make sure each set of parameters is unique", call. = FALSE)
      # evaluate the expansion parameters in the data frame
      # in case it contains derived parameters
      lapply(ldots, function(i) lazy_eval(i, data = .)) %>%
        expand.grid(stringsAsFactors = FALSE)

    }) %>%
    ungroup()

  return(invisible(ip))
}

#' mutate parameters
#'
#' modify existing parameters with mutate statements
#' @param ... mutate expressions
#' @export
mutate_parameters <- function(ip, ...) {
  if (!is(ip, "isopath")) stop ("parameters can only be mutate for an isopath")
  if (nrow(ip$parameters) == 0) stop("no parameters set yet for this isopath")

  ldots <- lazy_dots(...)
  ip$parameters <-
    ip$parameters %>%
    mutate_(.dots = ldots)

  return(invisible(ip))
}

#' @rdname set_parameters
#' @export
get_parameters <- function(ip) {
  if (!is(ip, "isopath")) stop ("can only get parameters from an isopath")
  ip$parameters
}
