#' add parameter data frame to the isopath
#' @export
set_parameters <- function(ip, params) {
  if (!is(ip, "isopath")) stop ("parameters can only be set for an isopath")
  if (!is.data.frame(params))
    stop("parameters must be provided as a data frame. The minimal set of parameters comprise the starting values for all the variables in the system, call 'get_variables(isopath)' to see a list of minimally required columns. Additionally, all parameters referenced in calculations can be included.") # note: also all the things called out in the function calls matter!

  missing <- setdiff(ip %>% get_variables(), names(params))
  if (length(missing) > 0)
    stop("parameters required for minimal parameter set missing: ",
         missing %>% paste(collapse = ", "), call. = FALSE)

  ip$parameters <- params
  return(ip)
}
