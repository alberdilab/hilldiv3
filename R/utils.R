# Internal utilities.

# Marker for engine/feature functions that are scaffolded but not yet
# implemented. Keeps the package loadable and the API discoverable while the
# corresponding maths/feature is being filled in.
.NotYetImplemented <- function() {
  call <- sys.call(-1)
  fn <- if (!is.null(call)) deparse(call[[1]]) else "this function"
  cli::cli_abort(c(
    "{.fn {fn}} is not implemented yet in hilldiv3.",
    "i" = "This is a scaffolded stub; see the TODO in the source."
  ))
}
