#' Launch the app
#'
#' Shiny Web App using Rfastp.
#' @export
#' @rdname launch_app
launch_app <- function() {
  renv::restore("renv.lock")
  shiny::runApp("R")
}
