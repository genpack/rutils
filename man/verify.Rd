% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rutils.R
\name{verify}
\alias{verify}
\title{Verifies that a given isvariable is valid.}
\usage{
verify(
  var,
  allowed = NULL,
  domain = NULL,
  lengths = NULL,
  dims = NULL,
  null_allowed = T,
  names_domain = NULL,
  rownames_domain = NULL,
  names_include = NULL,
  names_identical = NULL,
  rownames_include = NULL,
  rownames_identical = NULL,
  fix = F,
  err_src = paste(deparse(sys.calls()[[max(sys.nframe() - 1, 1)]]), collapse = ""),
  err_msg = "Error in Operation!",
  default = NULL,
  varname = deparse(substitute(var))
)
}
\arguments{
\item{var}{Any variable}

\item{allowed}{Vector of characters containing valid classes. For example: c('numeric', 'integer')}

\item{domain}{Valid domain for the given variable. If \code{var} is numeric or integer,
\code{domain} must be a numeric vector of two elements specifying the lower and upper bounds.
If \code{var} is character or factor, \code{domain} is a list of valid values.
If \code{domain} is \code{NULL} (default), domain compliance will not be chacked.}

\item{lengths}{Valid lengths for the given variable. Only used for vectors
or classes for which generic function length() is defined.}

\item{dims}{Valid dimensions for the given variable. Only used for data.frames and matrices.
or classes for which generic function dim() is defined.}

\item{null_allowed}{A logical True if NULL values can be passed to argument \code{var}}

\item{names_domain}{Vector of characters specifying valid domain for the names of the given variable.
Only used for classes for which generic function names() is defined, like vectors, data.frames and matrices.}

\item{rownames_domain}{Vector of characters specifying valid domain for the rownames of the given variable.
Only used for classes for which generic function rownames() is defined, like data.frames and matrices.}

\item{names_include}{Vector of characters specifying names that the given variable must include.
Only used for classes for which generic function names() is defined, like vectors, data.frames and matrices.}

\item{names_identical}{Vector of characters. Names of the given variable is checked to be identical to this vector.
Only used for classes for which generic function names() is defined, like vectors, data.frames and matrices.}

\item{rownames_include}{Vector of characters specifying rownames that the given variable must include.
Only used for classes for which generic function rownames() is defined, like data.frames and matrices.}

\item{rownames_identical}{Vector of characters. Row names of the given variable is checked to be identical to this vector.
Only used for classes for which generic function rownames() is defined, like data.frames and matrices.}

\item{err_src}{A character string specifying the source generating the error if verification fails.}

\item{err_msg}{A character string specifying the error message if input if the input variable is the result of a failed operation.}
}
\description{
Verifies that a given isvariable is valid.
}
