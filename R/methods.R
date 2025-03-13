
#'@export
# The print method (shows just the best design)
# print.design_list <- function(object, ...) {
print.design_list <- function(x, ...) {
  print(x$BestDesign, ...)
}

#'@export
# The summary method (to show all the designs)
# summary.design_list <- function(object, ...) {
summary.design_list <- function(object, ...) {
  print(object$AllDesigns, ...)
}