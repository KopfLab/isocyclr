#' add parameter data frame to the isopath
#'
#' to expand on an existing parameter set with multiple variations, use \code{\link{expand_parameters}}
#'
#' @param ... parameters, either a whole set of parameters as a single data frame, or a list of individual values, will try to merge with existing parameters and if neceesary overwrite existing parameters with the same name
#' @family parameters
#' @export
set_parameters <- function(ip, ...) {
  if (!is(ip, "isopath")) stop ("parameters can only be set for an isopath")

  ldots <- lazy_dots(...)
  if (length(ldots) == 1 && is.data.frame(df <- lazy_eval(ldots[[1]]))) {
    # overwrite if single data frame passed in
    ip$parameters <- df
  } else if (nrow(ip$parameters) == 0) {
    # overwrite if nothing set yet
    ip$parameters <- tibble(...)
  } else {
    # mutate otherwise
    tryCatch({
      ip$parameters <- mutate(ip$parameters, ...)
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
#' @family parameters
#' @export
expand_parameters <- function(ip, ...) {

  if (!is(ip, "isopath")) stop ("parameters can only be expanded for an isopath")
  if (nrow(ip$parameters) == 0) stop("no parameters set yet for this isopath")
  ip$parameters <- ip$parameters %>% expand_data_frame(...)
  return(invisible(ip))
}

#' expand data frame
#' utility function to expand data frames with non-standard evaluation
#' used by \code{\link{expand_parameters}}
#' @param ... use for non standard evaluation
#' @param .dots use for standard evaluation
#' @export
expand_data_frame <- function(df, ..., .dots = list()) {

  # dots
  ldots <- c(lazy_dots(...), .dots)

  # unique row ids
  df <- df %>% mutate(..rowid = 1:nrow(df))

  df %>%
    # group by everything to get unique sets
    group_by(!!!purrr::map(names(df), rlang::sym)) %>%
    # expand for each row
    do({
      # evaluate the expansion in the data frame in case it contains derived fields
      lapply(ldots, function(i) lazy_eval(i, data = .)) %>%
        expand.grid(stringsAsFactors = FALSE)

    }) %>%
    ungroup() %>%
    select(-..rowid)
}

#' mutate parameters
#'
#' modify existing parameters with mutate statements
#' @param ... mutate expressions
#' @family parameters
#' @export
mutate_parameters <- function(ip, ...) {
  if (!is(ip, "isopath")) stop ("parameters can only be mutate for an isopath")
  if (nrow(ip$parameters) == 0) stop("no parameters set yet for this isopath")

  ip$parameters <-
    ip$parameters %>%
    mutate(...)

  return(invisible(ip))
}

#' @rdname set_parameters
#' @export
get_parameters <- function(ip) {
  if (!is(ip, "isopath")) stop ("can only get parameters from an isopath")
  ip$parameters
}
