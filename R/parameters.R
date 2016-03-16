#' add parameter data frame to the isopath
#' @export
set_parameters <- function(ip, params) {
  if (!is(ip, "isopath")) stop ("parameters can only be set for an isopath")
  if (!is.data.frame(params))
    stop("parameters must be provided as a data frame. to see what columns are required minimally, call 'get_parameters_template()' for your isopath.") # note: also all the things called out in the function calls matter!

  missing <- setdiff(names(params), names(ip %>% get_parameters_template()))
  if (length(missing) > 0)
    stop("parameters required for minimal parameter set missing: ",
         missing %>% paste(collapse = ", "), call. = FALSE)

  ip$parameters <- params
  return(ip)
}


#' get parameter template
#' @export
get_parameters_template <- function(ip) {
  if (!is(ip, "isopath")) stop ("need an isopath to generate a parameter template")

  df <- list()
  df[ip$components %>%
      sapply(function(i) {
        c(i$name, paste0(i$name, ".", names(i$isotopes)))
      }) %>% unlist()] <- 0
  df %>% as_data_frame()
}
