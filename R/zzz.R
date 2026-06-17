# Register S3 methods for generics that live in Suggested packages, so those
# packages remain optional. ggplot2::autoplot() is only needed when a user
# actually calls autoplot() on a result.
.onLoad <- function(libname, pkgname) {
  for (cls in c("hill_diversity", "hill_profile", "hill_evenness",
                "hill_partition", "hill_hierarchy", "hill_dissimilarity",
                "hill_similarity")) {
    s3_register("ggplot2::autoplot", cls)
  }
}

# Standalone S3-method registration for generics in Suggested packages, vendored
# from the rlang/vctrs "import-standalone-s3-register" snippet so we do not
# depend on a particular rlang version exporting it.
s3_register <- function(generic, class, method = NULL) {
  stopifnot(is.character(generic), length(generic) == 1)
  stopifnot(is.character(class), length(class) == 1)

  pieces <- strsplit(generic, "::")[[1]]
  stopifnot(length(pieces) == 2)
  package <- pieces[[1]]
  generic <- pieces[[2]]

  caller <- parent.frame()

  get_method_env <- function() {
    top <- topenv(caller)
    if (isNamespace(top)) asNamespace(environmentName(top)) else caller
  }
  get_method <- function(method) {
    if (is.null(method)) {
      get(paste0(generic, ".", class), envir = get_method_env())
    } else {
      method
    }
  }

  register <- function(...) {
    envir <- asNamespace(package)
    method_fn <- get_method(method)
    registerS3method(generic, class, method_fn, envir = envir)
  }

  # Avoid registration failures during loading (pkgload or regular).
  if (isNamespaceLoaded(package)) {
    register()
  }

  setHook(packageEvent(package, "onLoad"), function(...) register())
}
