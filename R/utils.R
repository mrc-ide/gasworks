
`%||%` <- function(a, b) { # nolint
  if (is.null(a)) b else a
}

named_list <- function(names, x = NULL) {
  ret <- vector("list", length(names))
  names(ret) <- names
  if (!is.null(x)) ret[] <- x
  ret
}
