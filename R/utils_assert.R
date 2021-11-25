assert_is <- function(x, what, name = deparse(substitute(x))) {
  if (!inherits(x, what)) {
    stop(sprintf("'%s' must be a %s", name,
                 paste(what, collapse = " / ")), call. = FALSE)
  }
  invisible(x)
}

assert_named <- function(x, unique = FALSE, name = deparse(substitute(x))) {
  if (is.null(names(x))) {
    stop(sprintf("'%s' must be named", name), call. = FALSE)
  }
  if (unique && any(duplicated(names(x)))) {
    stop(sprintf("'%s' must have unique names", name), call. = FALSE)
  }
}

assert_character <- function(x, name = deparse(substitute(x))) {
  if (!(is.character(x))) {
    stop(sprintf("'%s' must be a character", name), call. = FALSE)
  }
  invisible(x)
}


assert_integer <- function(x, len = length(x), name = deparse(substitute(x)),
                           what = "integer") {
  if (!(is.integer(x))) {
    eps <- sqrt(.Machine$double.eps)
    usable_as_integer <- is.numeric(x) && (max(abs(round(x) - x)) < eps)
    if (!usable_as_integer) {
      stop(sprintf("'%s' must be an %s", name, what), call. = FALSE)
    }
    x <- as.integer(round(x))
  }
  assert_length(x, len, name)
  invisible(x)
}


assert_length <- function(x, len = length(x), name = deparse(substitute(x))) {
  if (length(x) != len) {
    stop(sprintf("'%s' must be of length %s", name, len), call. = FALSE)
  }
  invisible(x)
}


assert_strictly_increasing <- function(x, len = length(x),
                                       name = deparse(substitute(x))) {
  force(name)
  if (any(diff(x) <= 0)) {
    stop(sprintf("'%s' must be strictly increasing", name), call. = FALSE)
  }
  assert_length(x, len, name)
  invisible(x)
}

assert_unit_interval <- function(x, len = length(x),
                                 name = deparse(substitute(x))) {
  force(name)
  if (any(x < 0) | any(x > 1)) {
    stop(sprintf("'%s' must be between 0 and 1", name), call. = FALSE)
  }
  assert_length(x, len, name)
  invisible(x)
}

assert_positive <- function(x, len = length(x),
                                    name = deparse(substitute(x))) {
  force(name)
  if (any(x <= 0)) {
    stop(sprintf("'%s' must be greater than 0", name), call. = FALSE)
  }
  assert_length(x, len, name)
  invisible(x)
}

assert_nonnegative <- function(x, len = length(x),
                            name = deparse(substitute(x))) {
  force(name)
  if (any(x < 0L)) {
    stop(sprintf("'%s' must be greater than or equal to 0", name),
         call. = FALSE)
  }
  assert_length(x, len, name)
  invisible(x)
}

assert_positive_integer <- function(x, len = length(x),
                                    name = deparse(substitute(x))) {
  force(name)
  x <- assert_integer(x, len, name)
  x <- assert_positive(x, len, name)
  invisible(x)
}

assert_nonnegative_integer <- function(x, len = length(x),
                                       name = deparse(substitute(x))) {
  force(name)
  x <- assert_integer(x, len, name)
  x <- assert_nonnegative(x, len, name)
  invisible(x)
}
