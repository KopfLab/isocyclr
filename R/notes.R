# reversible reaction: A + B == C + D
#
# irreversible reaction: A + B -> C + D
#
#
#

# construction notes
#' Modify an isopath by adding components or reactions or another isopath
#'
# "+.isocyclr" <- function(e1, e2) {
#   # Get the name of what was passed in as e2, and pass along so that it
#   # can be displayed in error messages
#   e2name <- deparse(substitute(e2))
#
#   if      (is.theme(e1))  add_theme(e1, e2, e2name)
#   else if (is.ggplot(e1)) add_ggplot(e1, e2, e2name)
# }
