# Header
# Filename:     gener.R
# Description:  Contains some general functions useful for R programming
# Author:       Nicolas Berta 
# Email :       nicolas.berta@gmail.com
# Start Date:   21 October 2013
# Last change:  05 August 2019
# Version:      3.0.9

# Version   Date               Action
# -----------------------------------
# 1.2.3     18 July 2016       Function no.missing() added
# 1.2.4     18 July 2016       Function prettyDate() added
# 1.2.5     19 July 2016       Function yesterday() added
# 1.3.0     26 August 2016     Argument 'centralize' added to function extract.seasonality()
# 1.3.1     07 September 2016  Function verify() modified: Checks class inheritance rather than being exactly the same class specified by argument 'allowed'
# 1.3.2     08 September 2016  Action operators: %-%, %U% and %^% created for difference, union and intersect of sets
# 1.3.3     15 September 2016  Action operator: %<% created for subset
# 1.3.4     29 September 2016  Functions numerics() and nominals() added
# 1.3.5     29 September 2016  Library 'visgen.R' transferred to package viser
# 1.3.6     29 September 2016  Documentation for functions high.pass.mean() and high.pass.moving.mean() added
# 1.3.7     06 October 2016    Documentation added for functions nominals() and numerics()
# 1.3.8     14 October 2016    Some time formats added to valid.time.formats
# 1.4.0     27 October 2016    function fortday modified: given time argument is converted to timeDate
# 1.4.1     27 October 2016    function is.clean() added
# 1.4.2     24 November 2016   function datetimes() added
# 1.4.3     25 November 2016   function verify() modified: Logical argument 'null_allowed' added, defaulted by TRUE, set to FALSE if you want the verification to fail if argument 'var' passed to the function is NULL.
# 1.5.0     01 December 2016   File io.R (ver 1.2.2) added to the package
# 1.5.1     08 December 2016   Function sample.time() added
# 1.5.2     20 December 2016   Function numerics() modified: class 'difftime' added as a numeric type
# 1.5.3     20 December 2016   Function delib() added. delib() is used to de-library or detach a loaded package.
# 1.5.4     03 February 2017   Comperative operator %==% added which checks equality for sets
# 1.5.6     08 February 2017   Function smartMap() added
# 1.5.7     13 March 2017      Function vec2Col() added
# 1.5.8     13 March 2017      Small modifications to function verify()
# 1.5.9     16 March 2017      Function high.pass.prod() addded.
# 1.6.0     16 March 2017      Six functions added: mean.narm(), sd.narm(), sum.narm(), prod.narm(), rowSums.narm() and colSums.narm() added.
# 1.6.1     26 March 2017      Functions list.extend() and list.extract() added
# 1.6.2     26 March 2017      Functions coerce() added
# 1.6.3     27 March 2017      Function list.edit() added
# 1.6.4     27 March 2017      Action merge %<==>% added
# 1.6.5     30 March 2017      Functions extend.char() and repeat.char() modified: Now support for vector arguments
# 1.6.6     31 March 2017      Function chooseif() renamed to chif() and a small update: returns NULL if first argument is NULL or empty
# 1.6.7     31 March 2017      Function is.empty() modified: code improved
# 1.6.8     31 March 2017      Function colorise() added.
# 1.7.0     14 April 2017      Function color.mean() added.
# 1.7.1     14 April 2017      Function verify() modified: stops silent
# 1.7.2     19 April 2017      Function rename.items() added
# 1.7.3     22 April 2017      Function is.empty() modified: Returns TRUE if a list of NULLs or NAs is given
# 1.7.4     10 May 2017        Function column2Rownames() added. Works better than column_to_rownames from tibble package (makes unique, overwitres previous rownames, returns a data.frame)
# 1.7.5     21 May 2017        Function coerce() modified: retains names in the output
# 1.7.6     26 May 2017        Function zero.omit() added
# 1.7.7     26 May 2017        Function sandwich() added
# 1.7.8     26 May 2017        Function unname() added
# 1.7.9     26 May 2017        Function na2zero() added
# 1.7.10    26 May 2017        Function column2Rownames() modified: Default value for argument 'column'
# 1.8.0     01 June 2017       Function list.extract.field() added. If L is a named list of lists, instead of x = c('name', 'address'); L[[x]]$jojo, use: L %>% list.extract.field(x, field_name = 'jojo')
# 1.8.1     06 June 2017       Fixed the bug for factors in function vect.extend()
# 1.8.2     09 June 2017       Function vect.map() modified: Fixed the bug when input is a vector with equal values
# 1.8.3     12 July 2017       Function support() added
# 1.8.4     20 July 2017       Bug in function support() fixed!
# 1.8.5     03 August 2017     Small bug in function appendCol() fixed!
# 1.8.6     03 August 2017     'difftime' added to valid.numeric.classes
# 1.8.7     22 August 2017     Function toVectorList() added
# 1.8.9     28 August 2017     nibetree.R updated to version 0.0.2
# 2.2.2     21 September 2017  clust.R added to the package
# 2.2.3     27 October 2017    Function compareTables() added to gener.R
# 2.2.4     30 October 2017    Function compareTables() modified: small bugs fixed regarding notrig flag
# 2.2.5     31 October 2017    Function check() modified: You can pass any argument to this function
# 2.2.6     01 November 2017   Function swap() modified: Works for data.frame and matrix as well as list
# 2.2.7     08 November 2017   Function as.time() modified: Calls function coerce() instead of as()
# 2.2.8     08 November 2017   Function as.time() modified. calls char2Time() when input is character or factor
# 2.2.9     08 November 2017   Function setTZ() added and exported. Changes the time value by changing the time zone
# 2.3.0     10 November 2017   Function list.flatten() added and exported.
# 2.3.1     13 November 2017   Function time2Char() modifed: Returns time with desired format.
# 2.3.2     22 November 2017   Function col2Hex() added
# 2.3.4     28 November 2017   Functions diffDate() and diffTime() added
# 2.3.5     30 November 2017   Small bug in function diffDate() rectified: Did not work for empty input arguments
# 2.3.6     30 November 2017   Small bug in function as.time() rectified: returned error when given arguments is NA
# 2.3.7     11 December 2017   Small modification in function diffDate(): extends both d1 and d2 to max(length(d1), length(d2))
# 2.3.8     05 January 2018    Bug fixed in function diffTime(): t2 and t1 swapped
# 2.4.0     05 January 2018    Functions addDate() and addTime() added.
# 2.4.1     08 January 2018    Functions addDate() and addTime() modified: A smarter method used
# 2.4.2     16 January 2018    Functions addDate() and addTime() modified
# 2.4.3     06 Febryary 2018   package timeDate is not imported and is no more a dependency
# 2.4.4     12 Febryary 2018   Functions assert() and verify() modified: Detects varname and err_src arguments by itself if these arguments are null
# 2.4.5     21 Febryary 2018   Functions assert() and verify() modified: concatenates err_src to a single character if deparse returns multiple lables
# 2.4.6     21 Febryary 2018   Function nameColumns() added
# 2.4.7     27 Febryary 2018   Function column2Rownames() modified: Argument remove added. Option whether or not remove the original column
# 2.4.8     09 March 2018      Function is.empty() modified: Bug rectified: When a character verctor of NA was passed, would return NA
# 2.4.9     03 May 2018        Function nameColumns() modified: Argument classes can be NULL
# 2.5.0     03 May 2018        Documentation added for function nameColumns()
# 2.5.1     15 May 2018        io.R changed to version 1.3.2
# 2.5.2     21 May 2018        Function verify modified: Argument 'fix' added to fix the input variable to meet required conditions. Currently Only works for arguments domain, names_identical and rownames_identical
# 2.5.3     21 May 2018        Action function '%<==>%' modified: works with left side object to be NULL
# 2.5.4     25 May 2018        Function equal() modified: Does not call all.equal() any more. returns TRUE if the absolute difference is higher than given tolerance
# 2.5.6     01 June 2018       Functions diffTime() and diffDate() modified. Bugs rectified!
# 2.5.7     05 June 2018       Functions diffTime() modified. Bug rectified!
# 2.5.9     13 June 2018       Functions getCallingFunctionName() and contrastColor() added and exported.
# 2.6.0     13 June 2018       Action function merge '%<==>%' modified.
# 2.6.1     13 June 2018       Function verify() modified: calls getCallingFunctionName() to get a default value for argument err_src
# 2.6.2     15 June 2018       Function contrastColors() added. Finds contrasting color (black or white) for a vector of colors: todo: need to write a better function
# 2.6.3     20 June 2018       More formats added to global variable valid.time.formats, formats with '/' seperators eliminated.
# 2.6.4     20 June 2018       Function char2Time() modified: Converts all '/' to '-' before searching for correct format.
# 2.6.5     20 June 2018       Function as.time() modified: does not return error if all the timeset can not be converted
# 2.6.6     22 June 2018       Function findChainParts() added.
# 2.6.7     22 June 2018       Function verify() modified: Finds varname and err_src even if called via magrittr pipe
# 2.6.8     22 June 2018       Function verify() modified: error messages changed. varname and err_src are built only if verifucation failed.
# 2.6.9     25 June 2018       Function vect.extend() modified: returns input vector with no change if N = length(v)
# 2.7.0     25 June 2018       Global variable valid.date.formats added and exported.
# 2.7.1     25 June 2018       Function char2Date() modified: converts '/' to '-' before reading formats
# 2.7.2     26 June 2018       Functions numerics() and nominals() modified: A faster and more professional code. integer is now only included in numerics
# 2.7.3     02 July 2018       Functions nameColumns() modified: now works well for tibbles too.
# 2.7.5     03 July 2018       Functions na2zero() and cumulative() modified: does not change integer columns to double
# 2.7.7     05 July 2018       linalg.R updated to ver 1.2.3
# 2.7.8     24 July 2018       Functions list.edit() and list.add() modified: works with NULL as the first input argument
# 2.7.9     10 October 2018    Function appendCol() modified. Small change to rectify a bug
# 2.8.0     24 February 2019   Functions assert() and verify() modified: default value for argument err_src updated.
# 2.8.1     14 March 2019      Function list.default() added.
# 2.8.2     26 March 2019      Function charFilter() added. Test it and see how it works.
# 2.8.3     11 April 2019      Function partition() added. 
# 3.0.0     26 June 2019       io.r changes to version 2.0.0 
# 3.0.1     17 July 2019       io.r changes to version 2.0.1 
# 3.0.9     05 August 2019     Functions pdf(), cdf(), cdf.inv() and gen.random() plus 4 type conversion functions added.


# --------------------------------------------



#' @import magrittr
#' @include linalg.R
#' @include io.R
#' @include nibetree.R
#' @include clust.R

#' @export
support = function(...){
  packages = as.character(c(...))
  # verify(assert, 'character')
  if (length(packages) < 1){return(NULL)}
  if (length(packages) > 1){lapply(packages, support)}
  else {
    if (!require(packages, character.only = T)){
      res = try(install.packages(packages))
      if(inherits(res, 'try-error')){
        stop(paste("\n","\n", "Package", packages , "is not available and cannot be installed from cran! Please install it manually!", "\n", "\n"))
      }
      library(packages, character.only = T)
    }
  }
}

support('magrittr')

#' @export
wdlabel  = c(Mon = 'Monday', Tue = 'Tuesday', Wed = 'Wednesday', Thu = 'Thursday', Fri = 'Friday', Sat = 'Saturday', Sun = 'Sunday')

#' @export
fdlabel  = c(Mon.1 = 'Monday.1'  , Tue.1 = 'Tuesday.1', Wed.1 = 'Wednesday.1', Thu.1 = 'Thursday.1', Fri.1 = 'Friday.1',
             Sat.2 = 'Saturday.2', Sun.2 = 'Sunday.2' , Mon.2 = 'Monday.2'   , Tue.2 = 'Tuesday.2' , Wed.2 = 'Wednesday.2',
             Thu.2 = 'Thursday.2', Fri.2 = 'Friday.2' , Sat.1 = 'Saturday.1' , Sun.1 = 'Sunday.1')

#' @export
mntlabel = c(Jan = 'January', Feb = 'February' , Mar = 'March',   Apr = 'April',    May = 'May', Jun = 'June', Jul = 'July',
             Aug = 'August' , Sep = 'September', Oct = 'October', Nov = 'November', Dec = 'December')

#' @export
valid.colors = c('green', 'purple', 'aqua', 'blue', 'red', 'yellow', 'magenta', 'cyan', 'black', 'grey', 'orange')

#' @export
valid.time.classes = c("character", "factor", "Date", "POSIXct", "POSIXlt", "timeDate")

#' @export
valid.numeric.classes = c('numeric', 'integer', 'difftime')

#' @export
valid.nominal.classes = c('character', 'factor', 'logical')

#' @export
valid.ordinal.classes = c('integer')

#' @export
valid.classes = c(valid.nominal.classes, valid.time.classes, valid.numeric.classes) %>% unique

#' @export
valid.date.formats = c(
  "%Y-%m-%d", "%d-%m-%Y", "%d-%m-%y", "%B %d, %Y", "%B %d, %y", "%d%b%Y", "%d%b%y", "%A, %d %B %Y", "%A, %d %B %y", "%d%m%Y", "%d%m%y", "%m-%d-%Y", "%m-%d-%y", "%d-%b-%Y", "%d-%b-%y")

#' @export
valid.time.formats = c(
  "%Y-%m-%d %I:%M:%S %p", "%d-%m-%y %I:%M:%S %p", "%d-%m-%Y %I:%M:%S %p", "%d-%m-%y %I:%M %p", "%d-%m-%Y %I:%M %p", "%d-%b-%y %I:%M:%S %p", "%d-%b-%Y %I:%M:%S %p",
  "%Y-%m-%d %H:%M:%S", "%d-%m-%y %H:%M:%S", "%d-%m-%Y %H:%M:%S", "%d-%m-%y %H:%M", "%d-%m-%Y %H:%M", "%d-%b-%y %H:%M:%S", "%d-%b-%Y %H:%M:%S",
  valid.date.formats)


valid.time.formats.old = c("%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M", "%d/%m/%Y", "%B %d, %Y", "%d%b%y", "%A, %d %B %Y", "%d%m%y", "%m/%d/%Y", "%d/%b/%Y %H:%M:%S", "%d/%b/%Y", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d", '%Y/%m/%d %H:%M:%S', '%Y/%m/%d')

#' @export
valid.time.periods = c('sec' = 1, 'min' = 60, 'hour' = 3600, 'day' = 86400, 'week' = 7*86400, 'fortnight' = 14*86400, 'month' = 86400*365/12, 'quarter' = 86400*365/4, 'year' = 86400*365)

#' Use this function to check if a variable is empty or not
#'
#' @param x A variable. Can be a vector,matrix,list or data.frame
#' @return A boolean: TRUE if argument \code{x} is empty, FALSE otherwise
#' @examples
#' a = matrix()
#' is.empty(a)
#' [1] TRUE
#' a = c(2, 4)
#' is.empty(a)
#' [1] FALSE
#'
#' @export
is.empty = function(x){
  if (is.null(x)) {return(T)}
  # x %<>% na.omit causes great problems
  if      (inherits(x, c('matrix', 'data.frame', 'tibble', 'data.table'))){flag = (dim(x)[1] == 0) | (dim(x)[2] == 0) | (sum(!is.na(x)) == 0)}
  else if (inherits(x, 'list')){
    #y = x %>% unlist
    #if(inherits(y, 'list')){flag = y %>% list.clean %>% is.empty} else (flag = length(x) == 0)
    flag = (length(x) == 0)
  }
  else if (inherits(x, 'character')){flag = (sum(x != '', na.rm = T) == 0) | (sum(!is.na(x)) == 0)}
  else if (length(x) == 0){flag = T}
  else if (inherits(x, 'function')){flag = F} else {flag = sum(!is.na(x)) == 0}
  return(flag)
}

#' Rounds \code{x} to the nearest multiple of \code{N}
#'
#' @param x A numeric or integer
#' @param N A numeric or integer
#' @param adjust A character string:
#' 'closest': Returns the closest multiple of \code{N} \cr
#' 'top'   ': Returns the closest multiple of \code{N} that is greater than \code{x} \cr
#' 'bottom' : Returns the closest multiple of \code{N} that is smaller than \code{x} \cr
#' @return A numeric
#' @examples
#' roundto.multiple(345.45, 1)
#' [1] 345
#' roundto.multiple(345.45, 7)
#' [1] 343
#' roundto.multiple(345.45, 7.2)/7.2
#' [1] 48
#'
#' @export
roundto.multiple = function(x, N, adjust = 'closest'){
  # Check arguments:
  adjust <- tolower(adjust)
  verify(adjust, allowed = 'character', domain = c('closest', 'top', 'bottom'))
  
  k = x %/% N
  switch(adjust,
         'closest' = {flg    = (x - k*N > 0.5*N)
         inc(k[flg])},
         'top'     = {k = k + 1})
  
  return (k*N)
}

inc = function(v, increment = 1){
  v + increment
}

#' Terminates the execution with an error message if the given conditional statement is \code{FALSE}
#'
#' @param flag A boolean: Condition to be checked
#' @param err_msg A character string: Message to be displayed
#' @return NULL
#' @examples
#' x = 1
#' assert(x > 1, "Argument x can not be greater than 1")
#'
#' Error in assert(x > 1, "Argument x can not be greater than 1") :
#'   Argument x can not be greater than 1
#'
#' @export
assert <- function(flag, err_msg = 'Assertion Error !', err_src = paste(deparse(sys.calls()[[max(sys.nframe()-1,1)]]), collapse = "")){

  #
  if (length(flag) > 1){flag = (sum(!flag) == 0)}
  if (!flag) {
    if(is.null(err_src)){err_src = paste(deparse(sys.calls()[[sys.nframe()-1]]), collapse = "")}
    if (typeof(err_src) == "language"){err_src = as.character(as.expression(err_src))}
    if (!is.null(err_src)){err_msg %<>% make.err.msg(err_src)}
    stop(err_msg, call. = F)}
}

#' @export
warnif <- function(flag, wrn_msg, wrn_src = NULL){
  if (flag){cat(paste('Warning:', wrn_src, '\n \n', wrn_msg, '\n \n'))}
}

#' Makes a pretty error message to be displayed on the console
#'
#' @param err_msg A character string: Message to be displayed
#' @param err_src A character string: Source from which the error message is issued.
#' @return A character string containing the modified error message
#' @examples
#' stop(make.err.msg(err_msg = "Some error message", err_src = match.call()[[1]]))
#'
#' Error:
#'
#' Error from: match.call:
#' ------------------------
#'
#' Some error message
#'
#' @export
make.err.msg <- function(err_msg = '', err_src = 'gener::make.err.msg'){
  if (typeof(err_src) == "language"){err_src = as.character(as.expression(err_src))}
  paste0('\n from ', err_src, ':', '\n', '\n',
         #repeat.char('-',9 + nchar(err_src)),'\n', '\n',
         err_msg, '\n', '\n')
}

# Returns the name of the calling function. Useful when you want to know which function has called the function you are in
#' @export
getCallingFunctionName <- function(){
  calling_fcn <- deparse(sys.call(-1))
  stringr::str_replace_all(calling_fcn, pattern = "([a-z0-9_]*)(.*)",
                           replacement = "\\1")
}


#' @export
contrastColors = function(colors)
{
  colors %>% sapply(contrastColor)
}


contrastColor = function(color)
{
  rgb_colors <- ((grDevices::col2rgb(color) %>% as.numeric())/255)^2.2
  luminance  <- (0.2126 * rgb_colors[1]) + (0.7152 * rgb_colors[2]) + (0.0722 * rgb_colors[3])
  saturation <- (max(rgb_colors) - min(rgb_colors) + 1e-05)/(max(rgb_colors) + 1e-05)
  
  if (saturation < 0.35) {
    if (luminance > 0.3) {
      contrasting_color <- "#000000"
    }
    else {contrasting_color <- "#FFFFFF"}
  }
  else {
    contrasting_color <- "#FFFFFF"
  }
  contrasting_color
}

# The pipe function creates an environment that keeps track of the chain parts.
# I tried walking up the current execution environments looking for this variable and then use the lhs info stored there to find the symbol at the start of the pipe.
# This isn't well tested.
findChainParts <- function(){
  i <- 1
  while(!("chain_parts" %in% ls(envir=parent.frame(i))) && i < sys.nframe()) {
    i <- i+1
  }
  parent.frame(i)
}


#' Verifies that a given isvariable is valid.
#'
#' @param var Any variable
#' @param allowed Vector of characters containing valid classes. For example: c('numeric', 'integer')
#' @param domain Valid domain for the given variable. If \code{var} is numeric or integer,
#'        \code{domain} must be a numeric vector of two elements specifying the lower and upper bounds.
#'        If \code{var} is character or factor, \code{domain} is a list of valid values.
#'        If \code{domain} is \code{NULL} (default), domain compliance will not be chacked.
#' @param lengths Valid lengths for the given variable. Only used for vectors
#'          or classes for which generic function length() is defined.
#' @param dims Valid dimensions for the given variable. Only used for data.frames and matrices.
#'       or classes for which generic function dim() is defined.
#' @param null_allowed A logical True if NULL values can be passed to argument \code{var}
#' @param names_domain  Vector of characters specifying valid domain for the names of the given variable.
#'               Only used for classes for which generic function names() is defined, like vectors, data.frames and matrices.
#' @param rownames_domain Vector of characters specifying valid domain for the rownames of the given variable.
#'               Only used for classes for which generic function rownames() is defined, like data.frames and matrices.
#' @param names_include  Vector of characters specifying names that the given variable must include.
#'               Only used for classes for which generic function names() is defined, like vectors, data.frames and matrices.
#' @param names_identical  Vector of characters. Names of the given variable is checked to be identical to this vector.
#'               Only used for classes for which generic function names() is defined, like vectors, data.frames and matrices.
#' @param rownames_include  Vector of characters specifying rownames that the given variable must include.
#'               Only used for classes for which generic function rownames() is defined, like data.frames and matrices.
#' @param rownames_identical  Vector of characters. Row names of the given variable is checked to be identical to this vector.
#'               Only used for classes for which generic function rownames() is defined, like data.frames and matrices.
#' @param err_src A character string specifying the source generating the error if verification fails.
#' @param err_msg A character string specifying the error message if input if the input variable is the result of a failed operation.
#' @export
verify <- function(var, allowed = NULL, domain = NULL, lengths = NULL, dims = NULL, null_allowed = T,
                   names_domain = NULL, rownames_domain = NULL, names_include = NULL, names_identical = NULL,
                   rownames_include = NULL, rownames_identical = NULL, fix = F,
                   err_src = paste(deparse(sys.calls()[[max(sys.nframe()-1,1)]]), collapse = ""), err_msg = 'Error in Operation!', default = NULL, varname = deparse(substitute(var))){

  if (is.null(var)){if(null_allowed){return(default)} else {stop("NULL value is not allowed for argument '" %++% varname  %++% "'", call. = F)}}
  
  clsv = class(var)
  fail = (clsv == 'try-error')
  
  if (sum(fail, na.rm = T) != 0){stop(make.err.msg(err_msg, err_src = err_src), call. = F)}
  
  if (!is.null(allowed)){fail = !inherits(var, allowed)}
  
  if (sum(fail, na.rm = T) != 0){
    if (is.null(allowed)){s = ""} else {s = "Argument " %++% varname %++% " must be class " %++% paste(allowed, collapse = " or ") %++% "."}
    err.class.fail  = make.err.msg(paste("Class verification failed!", s, " was class", paste(class(var), collapse = " and ")), err_src = err_src)
    stop(err.class.fail, call. = F)
  }
  
  if (!is.null(lengths)){
    fail = !(length(var) %in% lengths)
    if (sum(fail, na.rm = T) != 0){stop(make.err.msg(paste0("Length verification failed!", " Argument " , varname , " must be of length " , paste(lengths, collapse = " or ") , " was of length " , length(var)), err_src = err_src), call. = F)}
  }
  if (!is.null(names_include)){
    fail  = !(names_include  %in% names(var))
    if (sum(fail, na.rm = T) != 0){
      erstr = paste0("Field inclusion verification failed! Argument ", varname , " must have these fields but it does not: ", paste(names_include[fail], collapse = ", "))
      stop(make.err.msg(erstr, err_src = err_src), call. = F)
    }
  }
  if (!is.null(names_domain)){
    fail = !(names(var) %in% names_domain)
    if (sum(fail, na.rm = T) != 0){
      erstr = paste0("Field domain verification failed! Argument ", varname , " must not have these fields but it does: ", paste(names(var)[fail], collapse = ", "))
      stop(make.err.msg(erstr, err_src = err_src), call. = F)
    }
  }
  if (!is.null(rownames_domain)){
    fail = !(rownames(var) %in% rownames_domain)
    if (sum(fail, na.rm = T) != 0){
      erstr = paste0("Row-name domain verification failed! Argument ", varname, " must not have these row-names but it does: ", paste(rownames(var)[fail], collapse = " , "))
      stop(make.err.msg(erstr, err_src = err_src), call. = F)
      # todo: fix should be done for other arguments as well
    }
  }
  if (!is.null(rownames_include)){
    fail = !rownames_include  %in%  (rownames(var))
    if (sum(fail, na.rm = T) != 0){
      erstr = paste0("Row-name inclusion verification failed! Argument ", varname, " must have these row-names but it does not: ", paste(rownames_include[fail], collapse = " , "))
      stop(make.err.msg(erstr, err_src = err_src), call. = F)
    }
  }
  if (!is.null(names_identical)){
    if(fix){
      fail = F
      if(var %>% dim %>% length == 2){
        var = var[, names %^% names_identical]
      } else if (inherits(var, 'list')){
        var %<>% list.extract(names %^% names_identical)
      } else {
        var = var[names %^% names_identical]
      }} else {
        fail = sum(!identical(sort(names(var)) , sort(as.character(names_identical))), na.rm = T) > 0
      }
    
    if (fail){stop(make.err.msg(paste0("Field verification failed!", "Argument ", varname , " must only have these fields ", paste(names_identical, collapse = ", "), " , but contained these fields: ", paste(names(var), collapse =  ", ")), err_src = err_src))}
  }
  if (!is.null(rownames_identical)){
    if(fix){fail = F; var = var[rownames_identical %^% rownames(var),]} else {
      fail = !identical(rownames(var), rownames_identical)
    }
    if (sum(fail, na.rm = T) != 0){stop(make.err.msg(paste0("Row-name verification failed!"), err_src = err_src), call. = F)}
  }
  
  if (!is.null(dims)){
    fail = sum(!(dim(var) %in% dims), na.rm = T) > 0
    if (fail){stop(make.err.msg(paste0("Dimension verification failed!", "Argument ", varname ," must be dim ", paste(dims, collapse = " or "), " was dim c(", paste(dim(var), collapse = ','), ")"), err_src = err_src), call. = F)}
  }
  
  if (!is.null(domain)){
    # assert(class(domain) %in% clsv, "Argument 'domain' should have the same class as argument 'var'!", err_src = "gener::verify")
    if (inherits(domain, c('character', 'factor', 'logical'))){
      if(fix){fail = F; var = var[which(var %in% domain)]} else {fail = sum(!(var %in% domain), na.rm = T) > 0}
    }
    else if (inherits(domain, c('numeric', 'integer', valid.time.classes))){
      if(fix){fail = F; var = ifelse(var < domain[1], domain[1], ifelse(var > domain[2], domain[2], var))}
      else {
        fail = sum((var < domain[1]) | (var > domain[2]), na.rm = T) > 0
      }
    }
    else {
      err.domain.class = make.err.msg(paste("Argument 'domain' is of class",class(domain)[1]," which is not supported!"), err_src = err_src)
      stop(err.domain.class, call. = F)
    }
    
    if (fail) {
      err.domain.fail  = make.err.msg(paste("domain verification failed! \n", "Argument ", varname, " must be in domain (", paste(domain, collapse = ' , '), ") \n",
                                            " Item(s): ", paste(var[fail], collapse = ' , ') ," were not in domain!") , err_src = err_src)
      stop(err.domain.fail, call. = F)
    }
  }
  return(var)
}


verify.new <- function(var, allowed = NULL, domain = NULL, lengths = NULL, dims = NULL,
                       names_domain = NULL, rownames_domain = NULL, names_include = NULL, names_identical = NULL,
                       rownames_include = NULL, rownames_identical = NULL, fix = F, null_allowed = T,
                       err_src = NULL, err_msg = 'Error in Operation!', default = NULL, varname = NULL){
  
  errcode = 0
  
  clsv = class(var)
  
  # Check to see if var is try-error
  if(inherits(var, 'try-error')){
    errcode   = 1
    errprefix = "Operation "
    errsuffix = "failed with error message " %++% as.character(clsv)
  }
  
  # Check to see if var is NULL
  if (!errcode & is.null(var)){
    if(null_allowed){return(default)}
    else {
      errcode   = 2
      errprefix = "NULL value is not allowed for variable "
      errsuffix = "!"
    }
  }
  
  # Check to see if var is NULL
  if(!errcode & inherits(allowed, 'character')){
    if(!inherits(var, allowed)){
      if(fix){
        var = try(var %>% coerce(allowed[1]), silent = T)
        if(!inherits(var, allowed)){
          errcode   = 4
          errprefix = "Argument "
          errsuffix = " must inherit class " %++%  paste(allowed, collapse = " or ") %++%  " but inherits " %++% paste(clsv, collapse = " --> ") %++% "! Furthermore, the value cannot be coerced to class " %++% allowed[1] %++% "."
        }
      } else {
        errcode   = 3
        errprefix = "Argument "
        errsuffix = " must inherit class " %++%  paste(allowed, collapse = " or ") %++%  " but inherits " %++% paste(clsv, collapse = " --> ") %++% "!"
      }
    }
  }
  
  # Check to see if var is in domain:
  if(!errcode & !is.null(domain)){
    # assert(class(domain) %in% clsv, "Argument 'domain' should have the same class as argument 'var'!", err_src = "gener::verify")
    if (inherits(domain, c('character', 'factor', 'logical', 'integer'))){
      if(sum(!(var %in% domain), na.rm = T) > 0){
        if(fix){var = var[which(var %in% domain)]}
        else{
          var_domain = var %-% domain
          errcode    = 5
          errprefix  = "Argument "
          errsuffix  = paste0(" must be in domain (",
                              paste(domain, collapse = " , "),
                              "). Element",
                              chif(length(var_domain) > 1, "s ", " "),
                              paste(var_domain, collapse = " and "),
                              chif(length(var_domain) > 1, " are ", " is "),
                              "not in the domain!")
        }
      }
    }
    
    else if (inherits(domain, c('numeric', valid.time.classes))){
      if(sum((var < domain[1]) | (var > domain[2]), na.rm = T) > 0){
        if(fix){var = ifelse(var < domain[1], domain[1], ifelse(var > domain[2], domain[2], var))}
        else {errcode = 5}
      }
    }
    
    else {stop(make.err.msg(paste("Argument 'domain' is of class", class(domain)[1], " which is not supported!"), err_src = 'gener::verify'))}
  }
  
  # Check to see if var lengths matches requirements:
  if(!errcode & inherits(length, c('numeric', 'integer'))){
    lengths = as.integer(lengths)
    if(!(length(var) %in% lengths)){
      if(fix){
        var = try(var %>% extend(lengths[1]), silent = T)
        if(inherits(var, 'try-error')){
          errcode   = 8
          errprefix = 'Length of Variable '
          errsuffix = ' must be ' %++% smartConcat(lengths, " or ") %++% ". Extension attempt failed with error " %++% as.character
        }
      } else {
        errcode   = 7
        errprefix = 'Length of Variable '
        errsuffix = ' must be ' %++% smartConcat(lengths, " or ") %++% ". Given argument has length " %++% length(var) %++% " which is not allowed!"
      }
    }
  }
  
  # Check to see if var dimension matches requirements:
  if (!errcode & inherits(dims, 'list')){
    dimvar = dim(var)
    pass   = (length(dims) == 0)
    for (e in dims){pass = pass | e == dimvar}
    if(!pass){
      if(fix){
        var = var[dims[[1]][1], dims[[1]][2]]
      } else {
        errcode   = 9
        errprefix = "Dimension verification failed for argument "
        errsuffix = "!"
      }
    }
  }
  
  # Check to see if names(var) is in names_domain:
  if(!errcode & inherits(names_domain, 'character')){
    namesvar = names(var)
    if(!is.null(namesvar)){
      if(!(namesvar %in% names_domain)){
        if(fix){var = var[names_domain]} else {
          errcode   = 11
          errprefix = "Argument "
          errsuffix = paste0(
            " must not have any of these names but it does: ",
            paste(namesvar %-% names_domain , collapse = ", "),
            "!")
        }
      }
    }
  }
  
  # Check to see if rownames(var) is in rownames_domain:
  if(!errcode & inherits(rownames_domain, 'character')){
    rownamesvar = rownames(var)
    if(!is.null(rownamesvar)){
      if(!(rownamesvar %in% rownames_domain)){
        if(fix){var = var[rownames_domain, ]} else {
          errcode   = 13
          errprefix = "Argument "
          errsuffix = paste0(
            " must not have any of these rownames but it does: ",
            paste(rownamesvar %-% rownames_domain , collapse = ", "),"!")
        }
      }
    }
  }
  
  # Check to see if names(var) include all elements in names_include:
  if(!errcode & inherits(names_include, 'character')){
    namesvar = names(var)
    if(is.null(namesvar)){errcode = 16} else {
      if(sum(!(names_include  %in% names(var))) > 0){
        errcode   = 15
        errprefix = "Argument "
        errsuffix = paste0(
          " must contain these names but it does not: ",
          paste(names_include %-% namesvar, collapse = ", "), "!")
      }
    }
  }
  
  # Check to see if names(var) is identical to names_identical:
  if(!errcode & inherits(names_identical, 'character')){
    namesvar = names(var)
    if(!identical(namesvar, names_identical)){
      errcode   = 17
      errprefix = 'Argument '
      errsuffix = paste0(
        " names must be identical to ",
        paste(names_identical, collapse = ", "), "!")
    }
  }
  
  # Check to see if rownames(var) include all elements in rownames_include:
  if(!errcode & inherits(rownames_include, 'character')){
    rownamesvar = rownames(var)
    if(is.null(rownamesvar)){errcode = 20} else {
      if(sum(!(rownames_include  %in% rownamesvar)) > 0){
        errcode   = 19
        errprefix = "Argument "
        errsuffix = paste0(
          " must contain these rownames but it does not: ",
          paste(rownames_include %-% rownamesvar, collapse = ", "), "!")
      }
    }
  }
  
  # Check to see if rownames(var) is identical to rownames_identical:
  if(!errcode & inherits(rownames_identical, 'character')){
    rownamesvar = rownames(var)
    if(!identical(rownamesvar, rownames_identical)){
      errcode   = 21
      errprefix = 'Rownames of argument '
      errsuffix = paste0(
        "must be identical to ",
        smartconcat(rownames_identical, ", "), "!")
    }
  }
  
  
  # Build proper error message:
  if(errcode > 0){
    # build varname
    if(is.null(varname)){varname = deparse(substitute(var))}
    if(varname == "."){ee <- findChainParts(); varname = deparse(ee$lhs)}
    # build source
    if(inherits(err_src, 'character')){err_src %<>% paste(collapse = ' --> ')} else {
      syscalls = sys.calls()
      err_src  = syscalls[syscalls %>% as.character %>% grep(pattern = varname) %>% min]
    }
    stop(paste0("Error from: ", err_src, "\n", errprefix, varname, errsuffix, "\n"))
  }
  
  return(var)
}

# if(2 > 1) {stop(lfeed("No", "This is not a good way", "be careful!"))}
lfeed <- function(...){
  paste('\n', ..., sep = '\n')
}

#' Replicates the given charachter to make a longer character string
#'
#' @param char a character: Character to be replicated
#' @param N An integer: number of times \code{char} should be repeated
#' @return A character string
#' @examples
#' repeat.char('X', 3)
#' [1] "XXX"
#'
#' @export
repeat.char <- function(char, N){
  if(length(N) > length(char)){char %<>% vect.extend(length(N))} else {N %<>% vect.extend(length(char))}
  N[N < 0] <- 0
  out = character()
  for(i in seq(char)){out = c(out, paste(replicate(N[i], char[i]), collapse = ""))}
  return(out)
}

#' Maps a given vector into a given range
#'
#' @param x A numeric vector: Vector to be mapped
#' @param minimum A numeric scalar: Specifies the lower bound of the output
#' @param maximum A numeric scalar: Specifies the higher bound of the output
#' @return A numeric vector
#' @examples
#' vect.map(c(1, 2, 3), minimum = 100, maximum = 300)
#' [1] 100 200 300
#'
#' @export
vect.map = function(x, minimum = 0, maximum = 1){
  #
  xl = min(x, na.rm = T)
  xh = max(x, na.rm = T)
  w  = xh - xl
  if (equal(w,0)){return(rep(0.5*(maximum + minimum), length(x)))}
  return(minimum + (maximum - minimum)*(x - xl)/w)
}

#' Replaces missing values of a vector by a given value
#'
#' @param v vector of any type: Vector containing missing values (\code{NA})
#' @param x Scalar value same type as \code{v}: The value by which the missing values are replaced (default is \code{0})
#' @return Vector same type as \code{v}: in which missing values are replaced by \code{x}
#' @examples
#' vect.clean(c(1, NA, 3), 2)
#' [1] 1 2 3
#'
#' @export
vect.clean <- function(v, x = 0){
  cv = v
  for (i in 1:length(v)){
    if (is.na(v[i])){
      cv[i] = x
    }
  }
  return(cv)
}

#' @export
vect.unclean <- function(v, x = 0){
  cv = v
  for (i in 1:length(v)){
    if (v[i] == x){
      cv[i] = NA
    }
  }
  return(cv)
}


#' @export
mat.clean <- function(M, x = 0){
  apply(M, 2, vect.clean, x = x)
}

#' @export
mat.unclean <- function(M, x = 0){
  apply(M, 2, vect.unclean, x = x)
}


#' @export
cumulative = function(v, by.row = T){
  verify(v, c('matrix', 'data.frame', 'numeric', 'integer', 'logical'), varname = 'v', null_allowed = F)
  if (inherits(v, c('matrix', 'data.frame'))){
    N   = ifelse(by.row, nrow(v), ncol(v))
    s   = ifelse(by.row, rep(as.integer(0), ncol(v)), rep(as.integer(0), nrow(v)))
    cv  = v
    if (by.row){
      for (i in sequence(N)){
        s      = s + v[i,]
        cv[i,] = s
      }
    } else {
      for (i in sequence(N)){
        s      = s + v[,i]
        cv[,i] = s
      }
    }
  } else if (inherits(v, c('numeric', 'integer', 'logical'))){
    N = length(v)
    s  = as.integer(0)
    cv = v
    for (i in sequence(N)){
      s     = s + v[i]
      cv[i] = s
    }
  }
  return(cv)
}

#' @export
tab.cumulative  = function(v){
  ct = c()
  t  = table(v)
  n  = length(t)
  s  = 0
  for (i in 1:n){
    s  = s + t[i]
    ct = c(ct, s)
  }
  names(ct) = names(t)
  return(ct)
}


#' @export
tab.top = function(t, threshold =0.5){
  tp = sort(t/sum(t), decreasing = TRUE)
  return(names(tp)[tp > threshold])
}

#' @export
tab.top.cumulative = function(t, threshold =0.5){
  if (sum(t) == 0){return (c())}
  if (is.null(names(t))){names(t) = as.character(t)}
  tp = sort(t/sum(t), decreasing = TRUE)
  i  = 1
  s  = tp[i]
  while (s < threshold){
    i = i + 1
    s = s + tp[i]
  }
  return(names(tp)[1:i])
}

#' @export
equal = function(x1, x2, tolerance = .Machine$double.eps){
  absdiff = abs(x1 - x2)
  return(absdiff < tolerance)
}

#' @export
vect.dim.equal <- function(v1,v2){
  #returns TRUE if the two given vectors have the same length
  n1 = length(v1)
  n2 = length(v2)
  return ((n1==n2))
}

#' @export
mat.dim.equal <- function(m1,m2){
  #returns TRUE if the two given matrices have the same dimensions
  d1 = dim(m1)
  d2 = dim(m2)
  return ((d1[1]==d2[1]) & (d1[2]==d2[2]))
}

#' Returns mean of a vector where values lower than a given threshold are filtered
#'
#' @param v numeric or integer: Given vector for which the filtered mean is computed
#' @param threshold numeric (scalar): value specifying the filtering threshold (default is \code{0})
#' @param na.rm logical: Should missing values be filtered? (default is TRUE)
#' @return numeric (scalar): Filtered mean computed
#' @examples
#' high.pass.mean(c(1, NA, 3, 2, -1, -2), threshold = 1)
#' [1] 2.5
#'
#' @export
high.pass.mean <- function(v, threshold = 0, na.rm = T){
  x = mean(v[v > threshold], na.rm = na.rm)
  if (is.na(x)){x = 0}
  return(x)
}

#' @export
mean_narm <- function(v){mean(v, na.rm = T)}

#' @export
sd_narm <- function(v){sd(v, na.rm = T)}

#' @export
sum_narm <- function(v){sum(v, na.rm = T)}

#' @export
prod_narm <- function(v){prod(v, na.rm = T)}

#' @export
rowSums_narm <- function(v){rowSums(v, na.rm = T)}

#' @export
colSums_narm <- function(v){colSums(v, na.rm = T)}



#' @export
high.pass.prod <- function(v, threshold = 0, na.rm = T){
  x = prod(v[v > threshold], na.rm = na.rm)
  if (is.na(x)){x = 0}
  return(x)
}


#' Returns moving average(mean) of a vector where values lower than a given threshold are filtered
#'
#' @param v numeric or integer: Given vector for which the filtered moving mean is computed
#' @param threshold numeric (scalar): value specifying the filtering threshold (default is \code{0})
#' @param weight integer (scalar): Inclusion window size for computing the moving average
#' @param na.rm logical: Should missing values be filtered? (default is TRUE)
#' @return numeric (scalar): Filtered mean computed
#' @examples
#' high.pass.moving.mean(1:100, weight = 10, threshold = 95)
#' [1] 98
#' # (96 + 97 + 98 + 99 + 100)/5 = 98
#'
#' @export
high.pass.moving.mean <- function(v, threshold = 0, weight = 14, na.rm = T){
  if (na.rm){v = v[!is.na(v)]}
  x = mean(tail(v[v > threshold], n = weight))
  if (is.na(x)){x = 0}
  return(x)
}

#' @export
high.pass.moving.sd <- function(v, threshold = 0.001, weight = 14, na.rm = T){
  if (na.rm){v = v[!is.na(v)]}
  x = sd(tail(v[v > threshold], n = weight))
  if (is.na(x)){x = 0}
  return(x)
}

#' @export
high.pass.sd <- function(v, threshold = 0.001, na.rm = T){
  x = sd(v[v > threshold], na.rm = na.rm)
  if (is.na(x)){x = 0}
  return(x)
}

low.pass.mean <- function(v, threshold, na.rm = T){
  x = mean(v[v < threshold], na.rm = na.rm)
  if (is.na(x)){x = 0}
  return(x)
}

low.pass.sd <- function(v, threshold, na.rm = T){
  x = sd(v[v < threshold], na.rm = na.rm)
  if (is.na(x)){x = 0}
  return(x)
}

#' @export
floorto.precision <- function(a, precision = 0.01){
  return(floor(a/precision)*precision)
}

#' @export
clear.workspace <- function(){
  rm(list = ls())
}

## Returns a matrix of zeros
zeros <- function (m,n){
  return(matrix(rep(0, m*n), m , n))
}

#' @export
mod = function(m, d){
  return(m - (m %/% d)*d)
}

#' @export
mod.seq = function(start, end, num){
  s = start
  v = c()
  for (i in 1:num){
    v = c(v, s)
    s = s + 1
    if (s > end){
      s = 0
    }
  }
  return (v)
}
# todo: convert it to an action

#' @export
fday <- function(t, base = '2000-01-01', zone = "GMT"){
  t = as.POSIXlt(t, tz = zone)
  b = as.POSIXlt(base, tz = zone)
  ndif = as.integer(difftime(t, b, units = 'day'))
  return(t$wday + 7*((ndif %/% 7) %% 2))
}

#' Returns the fortnight day of the given time
#' @export
fortday <- function(time, base = '2000-01-01'){
  support('timeDate')
  base = timeDate(base)
  time = as.timeDate(time)
  ndif = as.integer(time - base)
  x    = (ndif %/% 7) %% 2
  return(fdlabel[dayOfWeek(time) %++% '.'  %++% (x + 1)])
}

#' Extends the given vector \code{v} to length \code{N}
#'
#' @param v Vector of any type
#' @param N An integer
#' @return Vector same type as \code{v}
#' @examples
#' a = c(1, 2, 3)
#' vect.extend(a, 5)
#' [1] 1 2 3 1 2
#'
#' @export
vect.extend <- function(v, N){
  if(is.null(v)){return(v)}
  if(N == length(v)){return(v)}
  # in R, concatenating two factors makes an integer!
  if (inherits(v, 'factor')){was.factor = T; v %<>% as.character} else {was.factor = F}
  m = length(v)
  if (m >= N){
    return(v[sequence(N)])
  }
  ve = v
  while (N > m + length(ve)){ve = c(ve, v)}
  out = c(ve, v[sequence(N - length(ve))])
  if(was.factor){return (out %>% as.factor)} else {return(out)}
}


#' @export
list.extract   = function(obj, ...){
  N = c(...)
  l = list()
  if(inherits(N, 'character')){N = which(names(obj) %in% N)}
  for(i in N){
    l = c(l,list(obj[[i]]))
  }
  names(l) <- names(obj)[N]
  return(l)
}


#' @export
list.extract.field   = function(obj, field_name, ...){
  out = c()
  a   = obj %>% list.extract(...)
  for(i in a){out = c(out, i[[field_name]])}
  return(out)
}

#' @export
list.remove = function(obj, ...){
  fields = c(...)
  for (fig in fields){obj[[fig]] <- NULL}
  obj
}

#' @export
list.extend <- function(v, N){
  if(is.null(v)){return(v)}
  m = length(v)
  if (m >= N){
    return(v  %>% list.extract(sequence(N)))
  }
  ve = v
  while (N > m + length(ve)){ve = c(ve, v)}
  return (c(ve, v %>% list.extract(sequence(N - length(ve)))))
}


#' Returns the most frequent element of given vector \code{v}
#' @param v Vector
#' @return The most frequent element of given vector \code{v}
#' @export
most.common <- function(v){
  tbl = table(v)
  nms = names(tbl)
  mc  = nms[order(tbl, decreasing = T)[1]]
  if (is.null(mc)){mc = NA}
  return(mc)
}


#' @export
date.adjust = function(x, t, zone = "GMT"){
  t = as.POSIXlt(t, tz = zone)
  if (length(x) == 1){return(rep(x, length(t)))} else
    if (length(x) == 7){
      aux = t$wday
      aux[aux == 0] = 7
      return(x[aux])
    } else
      if (length(x) == 14){
        aux = fday(t)
        aux[aux == 0] = 14
        return(x[aux])
      } else
        if (length(x) == length(t)){return(x)} else {stop(make.err.msg("Arguments 'x' and 't' must have the same length", err_src = match.call()[[1]]))}
}

#' Returns if a given value meets a given feasibility condition. The feasibility condition is specified by
#' a feasibility range and a set of non-feasible values.
#' @param v A numeric or a vector of numerics or integers
#' @param range Vector of two numerics containing the lower and upper bounds of the feasibility range.
#'        If \code{range} is NULL (default), range condition checking is skipped.
#' @param exclude Vector of numerics containing non-feasible values
#'        If \code{exclude} is NULL (default), non-feasibility condition checking is skipped.
#' @return Logical TRUE if the given value is not missing, within the given range and not in \code{exclude},
#'                 FALSE otherwise.
#' @export
feasible = function(v, range = NULL, excludes = NULL){
  index = !is.na(v)
  if (!is.null(range)){index = index & (v <= range[2]) & (v >= range[1])}
  if (!is.null(excludes)){index = index & !(v %in% excludes)}
  return(index)
}

#' @export
first.feasible <- function(...){
  w = which(feasible(...))
  if (length(w) == 0){return(NA)} else {return(w[1])}
}

#' @export
last.feasible <- function(...){
  w = which(feasible(...))
  if (length(w) == 0){return(NA)} else {return(w[length(w)])}
}

#' @export
data.frame.na <- function(nrow, ncol, row_names = c(), col_names = c()){
  M = matrix(NA, nrow, ncol)
  D = as.data.frame(M)
  if (!is.empty(col_names)){
    assert(length(col_names) == ncol, "Argument 'col_names' does not match argument 'ncol'")
    colnames(D) <- col_names
  }
  if (!is.empty(row_names)){
    assert(length(row_names) == nrow, "Argument 'row_names' does not match argument 'nrow'")
    rownames(D) <- row_names
  }
  return(D)
}

#' @export
make.unique = function(v, char_duplicated = '.', char_non_duplicated = '', add_prefix = F){
  verify(v, allowed = 'character')
  dup.ind = duplicated(v)
  while (sum(dup.ind) > 0){
    if (add_prefix){
      v[dup.ind]  = paste0(char_duplicated, v[dup.ind], char_duplicated)
      v[!dup.ind] = paste0(char_non_duplicated, v[!dup.ind], char_non_duplicated)
    } else {
      v[dup.ind]  = paste0(v[dup.ind], char_duplicated)
      v[!dup.ind] = paste0(v[!dup.ind], char_non_duplicated)
    }
    dup.ind = duplicated(v)
  }
  return(v)
}

add.equivalent <- function(X, Y){
  # Inputs (X, Y): data.frame with two columns:
  #              Column 1 is any type
  #              Column 2 is numeric contains the values
  # X and Y should not necessarily have the same number of rows
  # adds values from Y to their equivalent values in X.
  # (equivalency is determined by matching values in column 1)
  Ny = dim(Y)[1]
  for (k in sequence(Ny)){
    if (Y[k,2] != 0){
      j      = which(X[,1] == Y[k, 1])
      X[j,2] = X[j,2] + Y[k,2]
    }
  }
  return (X)
}


#' @export
check = function(...){
  l = list(...)
}

# Given a time vector and a seasonal effect, this function
# applies the seasonal effect on each time interval accordingly"
#' @export
distribute.seasonality = function(time, season.values, seasonality = 'dow'){
  data.seas = rep(0, length(time))
  switch(seasonality,
         'dow' = {
           dow.lbls = dayOfWeek(time)
           for (x in unique(dow.lbls)){data.seas[x == dow.lbls] = season.values[x ,1]}},
         'dof' = {
           dof.lbls = fortday(time)
           for (x in unique(dof.lbls)){data.seas[x == dof.lbls] = season.values[x ,1]}},
         'moy' = {
           moy.lbls = mntlabel[months(time)]
           for (x in unique(moy.lbls)){data.seas[x == moy.lbls] = season.values[x ,1]}},
         'doy' = {
           tt = as.POSIXlt(time)
           # doy.lbls = tt$mon*31 + tt$mday
           doy.lbls = paste(tt$mday, tt$mon)
           for (x in unique(doy.lbls)){data.seas[x == doy.lbls] = season.values[x ,1]}}
  )
  data.seas[is.na(data.seas)] <- 0
  return(data.seas)
}

#' @export
add = function(x,y){
  if (is.na(x)){x = 0}
  if (is.na(y)){y = 0}
  return(x+y)
}

#' @export
repeat.row <- function(v, n){
  # Returns a matrix in which vector v repeats as multiple rows
  R = c()
  for (i in sequence(n)){
    R = rbind(R, v)
  }
  return(R)
}

#' @export
repeat.col <- function(v, n){
  # Returns a matrix in which vector v repeats as multiple columns
  R = c()
  for (i in sequence(n)){
    R = cbind(R, v)
  }
  return(R)
}

#' @export
'%++%' = function(obj1, obj2){
  paste0(obj1, obj2)
}

#' @export
'%-%' = function(obj1, obj2){
  setdiff(obj1, obj2)
}

# subset
#' @export
'%<%' = function(obj1, obj2){
  sum(!(obj1 %in% obj2)) == 0
}

#' @export
'%==%' = function(obj1, obj2){
  (obj1 %<% obj2) & (obj2 %<% obj1)
}

#' @export
'%U%' = function(obj1, obj2){
  if(inherits(obj1, 'factor')){obj1 = as.character(obj1)}
  if(inherits(obj2, 'factor')){obj2 = as.character(obj2)}
  unique(c(obj1, obj2))
}

#' @export
'%^%' = function(obj1, obj2){
  intersect(obj1, obj2)
}

#' @export
as.time = function(timeset, target_class = "POSIXct", formats = valid.time.formats, ...){
  if (inherits(timeset, target_class)) {return(timeset)}
  
  if (inherits(timeset, c('factor', 'character'))){
    tt = timeset %>% as.character %>% char2Time(formats = formats) %>% coerce(target_class, ...)
    if (inherits(tt,target_class) & (sum(is.na(tt)) == 0)){return(tt)}
  }
  
  if (target_class == 'character'){tt = timeset %>% time2Char}
  else if (target_class == 'factor'){tt = timeset %>% time2Char %>% as.factor}
  else if (target_class == "timeDate"){
    support('timeDate')
    tt <- try(as.timeDate(timeset, ...), silent = T)} else {tt <- try(timeset %>% coerce(target_class, ...), silent = T)}
  
  if (inherits(tt,target_class)){return(tt)}
  
  assert(F, "Given timeset can not be converted to a timeDate object. Check the class of timeset and the format!", match.call()[[1]])
}

as.time.new = function(timeset, target_class = "POSIXct", formats = valid.time.formats, ...){
  if (inherits(timeset, target_class)) {return(timeset)}
  
  if (inherits(timeset, c('factor', 'character'))){
    tt = timeset %>% as.character %>% char2Time(formats = formats) %>% coerce(target_class, ...)
    if (inherits(tt,target_class)){
      # if (inherits(tt,target_class) & (sum(is.na(tt)) == 0)){return(tt)}
      warnif(sum(is.na(tt)) > 0, 'Could not convert all given timeset!')
      return(tt)}
  }
  
  if (target_class == 'character'){tt = timeset %>% time2Char}
  else if (target_class == 'factor'){tt = timeset %>% time2Char %>% as.factor}
  else if (target_class == "timeDate"){
    support('timeDate')
    tt <- try(as.timeDate(timeset, ...), silent = T)} else {tt <- try(timeset %>% coerce(target_class, ...), silent = T)}
  
  if (inherits(tt,target_class)){return(tt)}
  
  assert(F, "Given timeset can not be converted to a timeDate object. Check the class of timeset and the format!", match.call()[[1]])
}

#' @export
extract.seasonality = function(time, data, seasonality = 'dow', centralize = F){
  # Verifications
  assert(seasonality %in% c('dow', 'moy', 'doy'), err_msg = "Unknown value for argument 'seasonality'. Must be in c('dow', 'moy', 'doy')", match.call()[[1]])
  assert(length(time) == length(data), "Arguments 'time' and 'data' must have the same lengths", match.call()[[1]])
  assert(require(timeDate), "Package 'timeDate' is not installed!")
  
  time = try(as.timeDate(time), silent = T)
  verify(time, err_msg = "Argument 'time' can not be converted to 'timeDate' object")
  
  switch(seasonality,
         'dow' = {
           S   <- aggregate(data, by = list(dayOfWeek(time)), FUN = mean)
           if (centralize){S$x   = S$x - mean(S$x)}
           rownames(S) = S[, 1]
           S           = S[,-1, drop = F]
         },
         'moy' = {
           mlb   = mntlabel[months(time)]
           S     = aggregate(data, by = list(mlb), FUN = mean)
           if (centralize){S$x   = S$x - mean(S$x)}
           rownames(S) = S[, 1]
           S           = S[,-1, drop = F]
         },
         'doy' = {
           tt    = as.POSIXlt(time)
           dylb  = paste(tt$mday, tt$mon)
           # dylb  = tt$mon*31 + tt$mday
           S     = aggregate(data, by = list(dylb), FUN = mean)
           if (centralize){S$x   = S$x - mean(S$x)}
           rownames(S) = S[, 1]
           S   = S[,-1, drop = F]
         })
  return(S)
}

#' @export
char2Date = function(datechar, formats = valid.date.formats){
  nf    = sequence(length(datechar))
  dates = as.Date(rep(NA, length(nf)))
  datechar %<>% gsub(pattern = '/', replacement = '-')
  i = 1
  while ((length(nf) > 0) & (i <= length(formats))){
    dates[nf] = as.Date(datechar[nf], format = formats[i])
    yearcount = dates[nf] %>% format('%Y') %>% as.integer
    dates[nf][which(yearcount < 100)] <- NA
    nf        = which(is.na(dates))
    i         = i + 1
  }
  return(dates)
}

char2Time.old = function(timechar, formats = valid.time.formats){
  nf    = sequence(length(timechar))
  times = as.POSIXlt(rep(NA, length(nf)))
  i     = 1
  while ((length(nf) > 0) & (i <= length(formats))){
    times[nf] = strptime(timechar[nf], format = formats[i])
    nf        = which(is.na(times))
    i         = i + 1
  }
  return(times)
}

#' @export
char2Time = function(timechar, formats = valid.time.formats){
  nf    = sequence(length(timechar))
  times = as.POSIXlt(rep(NA, length(nf)))
  timechar %<>% gsub(pattern = '/', replacement = '-')
  i     = 1
  while ((length(nf) > 0) & (i <= length(formats))){
    times[nf] = strptime(timechar[nf], format = formats[i])
    yearcount = times[nf] %>% format('%Y') %>% as.integer
    times[nf][which(yearcount < 100)] <- NA
    nf        = which(is.na(times))
    i         = i + 1
  }
  if(length(nf) > 0){
    if(require('anytime'))
      times[nf] = timechar[nf] %>% anytime::anytime()
  }
  
  nf = which(is.na(times))
  warnif(length(nf) > 0, 'Could not convert all given strings!')
  return(times)
}

# http://statistics.berkeley.edu/computing/r-dates-times
#' @export
time2Char = function(time, format = '%Y-%m-%d %H:%M:%S', omit_zero_time = F, make_unique = F){
  tc     <- as.character.POSIXt(time, format = format)
  if(omit_zero_time){
    change <- substr(tc, 12, 19) == "00:00:00"
    tc[change] <- substr(tc[change], 1, 10)
  }
  if (make_unique){tc = make.unique(tc)}
  return(tc)
}


# as.character.POSIXt('2017/12/01 13:28:32' %>% as.POSIXlt, format = '%d/%m/%Y %H:%M')


#' @export
# todo: Should change function name to time2DateChar see where this function is used!
time2Date = function(time, make_unique = F){
  ds <- as.character(as.Date(time))
  if (make_unique){ds = make.unique(ds)}
  return(ds)
}

# Example:
# script.list.fields(list.str = "default.dygraphs.tsline.settings", sublist.str = 'RangeSelector')
#' @export
script.list.fields = function(list.str, sublist.str){
  nms = eval(parse(text = paste0("names(" , list.str , "[['", sublist.str, "']])")))
  lst = nms[length(nms)]
  scr = ""
  for (p in nms){
    scr = paste0(scr , p, " = ", list.str , "$" , sublist.str, "$", p)
    if (p != lst){scr = scr %++% ', '}
  }
  return(scr)
}

#' @export
list2Json = function(L = list(), fields = names(L), fields_remove = c(), quotations = F){
  verify(L, "list")
  verify(fields, 'character', domain = names(L))
  verify(fields_remove, 'character', domain = fields)
  fields = setdiff(fields, fields_remove)
  
  scr  = "{"
  N    = length(fields)
  for (i in sequence(N)){
    x = fields[i]
    u = L[[x]]
    if (!is.null(u)){
      if (quotations){scr = paste0(scr,"'", x, "': ")} else {scr = paste0(scr, x, ": ")}
      
      if      (class(u) == 'list')     {scr = scr %++% list2Json(u, quotations = quotations)}
      else if (class(u) == 'character'){scr = paste0(scr, "'", u, "'")}
      else if (class(u) == 'logical')  {scr = scr %++% ifelse(u, 'true', 'false')}
      else                             {scr = scr %++% as.character(u)}
      if (i != N){
        scr = scr %++% ", "
      }
    }
  }
  scr = scr %++% "}"
  return(scr)
}

#' @export
list2Script = function(l, fields = NULL, fields_remove = NULL, arguments = NULL){
  if (is.null(fields)){fields = names(l)}
  if (  is.null(arguments)  |  !identical(names(arguments), fields)){
    assert(sum(!names(arguments) %in% fields) == 0, "fields must include names of argument 'arguments'!")
    args = fields
    names(args) <- fields
    if (! is.null(arguments)){
      for (j in names(arguments)){
        args[[j]] <- arguments[[j]]
      }
    }
    arguments <- args
  }
  
  fields = setdiff(intersect(fields, names(l)), fields_remove)
  
  scr = ""
  for (i in fields){
    if (inherits(l[[i]], "character")){
      if (length(l[[i]]) > 1){
        scr <- paste0(scr, arguments[i], " = c(", paste("'" %++% l[[i]] %++% "'", collapse = ', '), "), ")
      } else if (length(l[[i]]) == 1){scr <- paste0(scr, arguments[i], " = '", paste(l[[i]], collapse = ', '), "', ")}
    }
    else if (inherits(l[[i]], "list")){scr = paste0(scr, arguments[i], " = ", 'list(', list2Script(l[[i]]), '), ')}
    else if (class(l[[i]]) %in% c('interger', 'numeric', 'logical')){
      if (length(l[[i]]) > 1){
        scr <- paste0(scr, arguments[i], " = c(", paste(l[[i]], collapse = ', '), "), ")
      } else if (length(l[[i]]) == 1){scr <- paste0(scr, arguments[i], " = ", paste(l[[i]], collapse = ', '), ", ")}
    }
  }
  scr = substr(scr, 1, nchar(scr) - 2)
  return(scr)
}

#' @export
today      = function(){as.Date(timeDate())}

#' @export
yesterday  = function(){as.Date(timeDate()) - 1}


#' @export
today.char = function(){time.to.char(today())}

#' returns TRUE if none of the elements of given argument v is missing
#' @export
no.missing = function(v){
  return(sum(is.na(v)) == 0)
}

#' @export
prettyDate = function(x){
  if (!inherits(x, 'Date')){x = as.Date(x)}
  # paste(weekdays(x), mday(x), months(x), year(x))
  paste(weekdays(x), as.POSIXlt(x)$mday, months(x), as.numeric(format(x,'%Y')))
}

#' Returns numeric column labels of a data frame
#' @field df data.frame: The table for which numeric columns are required
#' @return vector of characters containing numeric column labels of the given table
#' @export
numerics = function(df){
  classes = sapply(df, class)
  classes[classes %in% c('numeric', 'integer', 'double', 'difftime')] %>% names
}

#' Returns nominal(categorical) column labels of a data frame
#' @field df data.frame: The table for which nominal columns are required
#' @return vector of characters containing nominal column labels of the given table
#' @export
nominals = function(df){
  classes = sapply(df, class)
  classes[classes %in% c('factor', 'character', 'logical', 'ordered', 'integer')] %>% names
}

#' Returns labels of the columns of the given data frame containing date-time values
#' @field df data.frame: The table for which date-time columns are required
#' @return vector of characters containing date-time column labels of the given data.frame
#' @export
datetimes = function(df){
  figures = names(df)
  tims = c()
  vtc  = valid.time.classes %-% c('character', 'factor')
  for (k in figures){
    if (inherits(df[, k], vtc)){tims = c(tims, k)}
  }
  return(tims)
}

#' @export
is.clean = function(v){
  sum(is.na(v)) == 0
}

#' @export
sample.time <- function(from, to, size) {
  dayseq <- seq.Date(as.Date(from),as.Date(to),by="day")
  dayselect <- sample(dayseq,size,replace=TRUE)
  hourselect <- sample(1:24,size,replace=TRUE)
  minselect <- sample(0:59,size,replace=TRUE)
  as.POSIXlt(paste(dayselect, hourselect,":",minselect,sep="") )
}



#' @export
delib <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

# If vector x contain outliers (extremely high or low values), then linear mapping will contain those outliers as well.
# smartMap() maps given vector 'x' to range (0,1) but also gives a more uniform distribution by limiting the effect of outliers
# The range [0,1] is divided into segments with lengths weighted as counts of values of x after being binned in 'n' bins
# Each segment of range (0,1) is for one bin
# The values in each bin are then linearly mapped to the range of their segment for which the width depends on the count of values that fall into that bin
#' @export
smartMap = function(x, n = 10){
  y = x - min(x)
  y = as.integer(vect.map(y, 1, n))
  t = tabulate(y, nbins = n)
  s = c(0, cumulative(t/sum(t)))
  z = rep(NA, length(x))
  for (i in 1:n){
    a = which(y == i)
    if (length(a) > 1){z[a] = vect.map(x[a], s[i], s[i + 1])}
    else if (length(a) == 1){z[a] = s[i]}
  }
  assert(sum(is.na(z)) == 0)
  return(z)
}

#' @export
colorise = function(x, palette){
  nuniques <- length(unique(x))
  palcols <- (grDevices::colorRampPalette(palette))(nuniques)
  if (!is.numeric(x) | nuniques < 10) {
    y <- as.numeric(as.factor(x))
    xcols <- palcols[y]
  }
  else {
    ecum <- ecdf(x)
    xcols <- palcols[ceiling(nuniques * ecum(x))]
  }
  xcols
}

# Converts a numeric vector into colors where color density reflects value
#' @export
vec2Col = function(v, base = c(red = 1, green = 1, blue = 1)){
  verify(v, valid.numeric.classes, varname = 'v')
  v = v[!is.na(v)]
  
  if (inherits(base, 'character')){
    verify(base, 'character', domain = colours(), lengths = 1, varname = 'base')
    base = vect.normalize(as.numeric(col2rgb(base)))
    names(base) = c('red', 'green', 'blue')
  }
  verify(base, 'numeric', domain = c(0, 1), names_identical = c('red', 'green', 'blue'), lengths = 3, varname = 'base')
  
  vn = vect.normalize(v)
  res = rgb(vn*base['red'], vn*base['green'], vn*base['blue'])
  names(res) <- names(v)
  return(res)
}

#' @export
color.mean = function(...){
  mn = (try(c(...) %>% col2rgb, silent = T) %>%
          verify(err_msg = 'Invalid colors given! Cannot convert to rgb.', err_src = 'gener::color.mean') %>%
          rowMeans)/255
  rgb(red = mn['red'], green = mn['green'], blue = mn['blue'])
}

#' @export
appendCol = function(data, columns, column_names = NULL){
  "
  Appends Vector or data.frame or matrix to the data
  Arguments:
  - columns Vector or data.frame. length or dim must match
  - column_names [Optional] Vector of character strings containing the column names to be appended
  "
  # todo: add verifications
  old.ns = names(data)
  if(columns %>% inherits(c('data.frame', 'matrix'))){
    
    N   = nrow(data)
    M   = ncol(columns)
    
    assert(nrow(columns) == N, 'arguments columns and data must have the same number of rows!')
    column_names %<>% verify('character', lengths = M, default = names(columns), varname = 'column_names')
    if(is.empty(column_names) & M > 0){column_names = 'X' %++% sequence(M)}
    names(columns) <- column_names
    
    for (fig in column_names){
      data[, fig] = columns[, fig]
    }
  } else {
    # Old Code is now used only for vectors:
    if(ncol(data) == 0){data = data.frame(columns, stringsAsFactors = inherits(columns, 'factor'))}
    else               {data %<>% cbind(columns %>% vect.extend(nrow(data)), stringsAsFactors = inherits(columns, 'factor'))}
    
    if (!is.null(column_names)){names(data) = c(old.ns, column_names)}
  }
  
  return(data)
  
}

#' @export
chif = function(a, x, y){
  if (is.empty(a)){return(NULL)}
  if (a) return(x) else {return(y)}
}

#' @export
extend.char = function (s, n, fillChar = ' ', left = T){
  M = nchar(s)
  return(ifelse(M >= n,
                chif(left, substr(s,1,n), substr(s, M - n + 1,M)),
                chif(left, paste0(s, repeat.char(fillChar, n - nchar(s))), paste0(repeat.char(fillChar, n - M), s))))
}

#' @export
swap = function(L, x = NULL, y = NULL, ..., vf = T){
  if(is.null(x) & is.null(y)){return(L)}
  if(vf){
    verify(L, c('list', 'data.frame', 'matrix'), varname = 'L', null_allowed = F)
    verify(x, 'character', varname = 'x', null_allowed = F)
    verify(y, 'character', varname = 'y', null_allowed = F)
  }
  if(inherits(L, 'list')){
    aux    <- L[[x]]
    L[[x]] <- L[[y]]
    L[[y]] <- aux
  } else if (inherits(L, c('data.frame', 'matrix'))){
    wx = which(colnames(L) == x)
    wy = which(colnames(L) == y)
    assert((length(wx) == 1) & (length(wx) == 1), 'Given column labels x and y should refer to a single column')
    aux     <- L[, x]
    L[, x]  <- L[, y]
    L[, y] <- aux
    colnames(L)[wx] <- y
    colnames(L)[wy] <- x
  } else {stop('I should not be here!')}
  return(L %>% swap(...))
}

#' @export
coerce = function(v, class_name, ...){
  # todo: complete verification: add domain, ...
  verify(class_name, 'character')
  if(is.empty(v) & (class_name %in% valid.classes)){v %<>% as.logical()}
  
  if(!inherits(v, class_name)){
    errmsg = paste("Given variable", deparse(substitute(v)), "cannot be coerced to class " %++% class_name)
    scr = paste0('as.', class_name, '(v, ...)')
    nms = names(v)
    if(inherits(v, 'POSIXct') & class_name == 'Date'){v %<>% lubridate::force_tz('GMT')}
    out = try(eval(parse(text = scr)), silent = T) %>% verify(err_msg = errmsg, err_src = 'gener::coerce')
    assert(inherits(out, class_name), err_msg = errmsg, err_src = 'gener::coerce')
    if(!is.null(nms)){names(out) <- nms}
    return(out)
  } else {return(v)}
}


#' @export
list.add = function(l = NULL, ...){
  if(is.empty(l)){l = list()}
  return(c(l, list(...)))
}

# list add
ladd = function(L, ...){
  # todo: write for classes other than list
  LL  = list(...)
  nms = names(LL)
  for (i in seq(LL)){
    item = LL[[i]]
    name = nms[i]
    if(!is.empty(item)){
      if(is.null(name)){
        L = c(L, item)
      } else {
        L[[name]] = item
      }
    }
  }
  return(L)
}

#' @export
list.clean = function(L){
  i = 1
  for (item in L){
    if(inherits(item, 'list')){item %<>% list.clean}
    if (item %>% is.empty){
      L[[i]] <- NULL
    } else {
      L[[i]] = item
      i = i + 1
    }
  }
  return(L)
}

# Merging Action:
#' @export
'%<==>%' = function(obj1, obj2){
  # Verifications:
  msg = 'The two objects must be of the same class!'
  src = match.call()[[1]]
  if(is.null(obj1)){if(inherits(obj2, 'list')){obj1 = list()}} # todo: what about other classes?!
  if (inherits(obj1, 'list') | typeof(obj1) == 'list'){
    assert(inherits(obj2, 'list') | typeof(obj2) == 'list', msg, src)
    
    ns   = names(obj2)
    if(is.null(ns)){return(obj1 %>% c(obj2))}
    w  = (ns %in% names(obj1)) & (ns != '')
    obj1 %<>% c(obj2 %>% list.extract(which(!w)))
    w    %<>% which
    for (i in w){
      if(inherits(obj1[[ns[i]]], 'list')){
        if(inherits(obj2[[i]], 'list')){
          obj1[[ns[i]]] = obj1[[ns[i]]] %<==>% obj2[[i]]
        } else {
          obj1[[ns[i]]] =  c(obj1[[ns[i]]], obj2[[i]])
        }
      }
      else {obj1[[ns[i]]] = obj2[[i]]}
    }
  }
  # todo: write similarly for named vectors & data.frames
  return(obj1)
}

#' @export
list.edit = function(l = NULL, ...){
  if(is.empty(l)){l = list()}
  return(l %<==>% list(...))
}

#' @export
aggr = function(df, id_cols = NULL, value_cols = NULL, func = 'mean', vf = T){
  # Verifications:
  if(vf){
    verify(df, c('data.frame', 'tibble', 'matrix'))
    id_cols    %<>% verify('character', domain = names(df), default = nominals(df), varname = 'id_cols')
    value_cols %<>% verify('character', domain = names(df), default = numerics(df), varname = 'value_cols')
  }
  
  grpscr = paste0('dplyr::group_by(', paste(id_cols, collapse = ', '), ')')
  varscr = paste0(value_cols, ' = ', func, '(', value_cols, ')')
  sumscr = paste0('dplyr::summarize(', paste(varscr, collapse = ', '), ')')
  
  eval(parse(text = paste('df', grpscr, sumscr, sep = ' %>% '))) %>% as.data.frame
}

reorder.tbc = function(v, new_order){
  if(!inherits(v, 'factor')){v %<>% as.factor}
  verify(new_order, 'character', domain = levels(v), varname = 'new_order')
  a = seq(new_order)
  names(a) <- new_order
  names(a)[a[levels(v)[v %>% as.numeric]]] %>% factor(levels = names(a))
} # todo: check with function reorder() in package stats

rename.items = function(lst, ...){
  if(is.empty(lst)){return(lst)}
  to   = c(...) %>% verify('character')
  from = to %>% names
  
  names(to) <- NULL
  
  for(i in seq(from)){
    w = which(names(lst) == from[i])
    names(lst)[w] <- to[i]
  }
  return(lst)
}

#' @export
column2Rownames = function(df, colname = NULL, remove = T){
  df %<>% as.data.frame
  nms = names(df)
  colname %<>% verify('character', domain = nms, default = names(df)[1], lengths = 1, varname = 'colname')
  ww = which(nms == colname)
  df = df[!duplicated(df[, colname]),, drop = F] # can add option to summarize(aggregate) later
  rownames(df) <- df[, colname]
  if(remove){return(df[, - ww, drop = F])} else {return(df)}
}

#' @export
rownames2Column = function(df, colname = 'rownames', colpos = 1){
  nms = colnames(df) %>% verify(default = 'V' %++% as.character(1:ncol(df)))
  verify(colname, 'character', lengths = 1, varname = 'colname')
  df  = df[, nms != colname, drop = F]
  m   = ncol(df)
  assert(colpos <= m + 1, "Maximum value for colpos is " %++% (m + 1))
  df %<>% cbind(rownames(df))
  colnames(df) <- c(nms, colname)
  rownames(df) <- NULL
  if (colpos == m + 1){return(df)}
  if (colpos == 1){return(df[, c(m + 1, 1:m)])}
  else {return(df[c(1:(colpos-1), m + 1, colpose:m)])}
}

#' @export
zero.omit = function(M, colname = NULL){
  if(inherits(M, c('matrix', 'data.frame'))){
    if(is.null(colname)){colname = M %>% ncol %>% sequence}
    for (col in colname){
      M = M[M[,col] != 0,, drop = F]
    }
    return(M)
  }
  
  keep = M != 0
  nms = names(M)[keep]
  M = M[keep]
  names(M) <- nms
  return(M)
}

swapNamesAndContents = function(v){
  nms = names(v)
  names(nms) <- v %>% as.character
  return(nms)
}

#' @export
sandwich = function(v, ss = "'"){
  ss %++% v %++% ss
}


#' @export
unname = function(v){
  names(v) <- NULL
  v
}


#' @export
na2zero = function(v){
  v[is.na(v)] <- as.integer(0)
  return(v)
}


#' @export
trim = function(v, lb = 0, ub = Inf){
  verify(v, 'numeric')
  
  v[(v >= lb) & (v <= ub)]
}


#' @export
removeRownames = function(v){
  rownames(v) <- NULL
  return(v)
}


factor2Character = function(df){
  i <- sapply(df, is.factor)
  df[i] <- lapply(df[i], as.character)
  return(df)
}


offset2GMT = function(tt){
  tt %>% as.character %>% strptime(format = '%Y-%m-%d %H:%M:%S', tz = 'GMT') %>% as.POSIXct
}


#' Converts a data.frame or matrix into a list of named vectors where each column is a named vector
#' @export
toVectorList = function(df){
  lst = list()
  rwnms = rownames(df)
  for(i in df %>% ncol %>% sequence){
    lst[[i]] = df[,i]
    if(!is.null(rwnms)){names(lst[[i]]) <- rwnms}
  }
  names(lst) <- names(df)
  return(lst)
}


# Compares two data.frames and returns a list of data.frames indicating which rows and columns should change from what to what
#' @export
compareTables = function(t1, t2){
  if(identical(t1,t2)){return(list())}
  # if number of rows or columns are different, trigger redraw
  nn = nrow(t1)
  mm = ncol(t1)
  notrig = (nrow(t2) == nn) & (ncol(t2) == mm)
  
  cold = colnames(t1)
  if(!is.null(cold)){
    cnew   = colnames(t2)
    notrig = notrig & (cnew %==% cold)
    if(notrig){if(!identical(cnew, cold)){t2 = t2[, cold]}}
  }
  
  i = 0
  while (notrig & (i < mm)){
    i = i + 1
    if(!inherits(t2[,i], t1[,i] %>% class)){
      t2i = try(coerce(t2[,i], class(t1[,i])), silent = T)
      notrig = notrig & inherits(t2i, t1[, i] %>% class)
      if(notrig){t2[,i] = t2i}
    }
    if(inherits(t1[,i], 'factor')){notrig = notrig & levels(t1[,i]) %==% levels(t2[,i])}
    notrig = notrig & sum(is.na(t1[,i]) >= is.na(t2[,i])) == nn
  }
  
  if(notrig){
    changes = list()
    
    rnew = rownames(t2)
    rold = rownames(t1)
    
    if(!is.null(rnew) & !is.null(rold)){
      diff = which(rnew != rold)
      if(length(diff) > 0){changes %<>% list.add(data.frame(row = diff, col = 0, from = rold[diff], to = rnew[diff], stringsAsFactors = F))}
    }
    
    for(i in sequence(ncol(t1))){
      if(inherits(t1[,i], 'numeric')){diff = which(abs(t2[,i] - t1[,i]) > 0.001)}
      else {diff = which(t2[,i] != t1[,i])}
      diff = try(c(diff, which(is.na(t1[,i]) & !is.na(t2[,i]))), silent = T) %>% verify(default = diff)
      if(length(diff) > 0){changes %<>% list.add(data.frame(row = diff, col = i, from = t1[diff,i], to = t2[diff,i], stringsAsFactors = F))}
    }
    
    return(changes)
  }
  
  return(NULL)
}

#' @export
setTZ = function(time, tz){
  time %>% as.character %>% as.POSIXct(tz = tz)
  # look at lubridate::force_tz() and lubridate::with_tz()
}


#' @export
list.flatten = function(lst){
  nms = list()
  for (i in lst){
    if(inherits(i, 'list')){nms = c(nms, i %>% list.flatten)} else {nms = c(nms, i)}
  }
  return(nms)
}

# Converts color name to equivalent hex code:
#' @export
col2Hex = function(color){
  col = try(col2rgb(color)) %>% verify
  rgb(red = col['red',], green = col['green',], blue = col['blue',], maxColorValue = 255)
}


is.holiday = function(d, holidays = c()){
  return((format(d, "%u") %in% c("6", "7")) | (d %in% holidays))
}

# Gives difference between two given dates excluding weekends and list of given holidays
#' @export
diffDate <- function(d1, d2, holidays = c()) {
  mm = max(length(d1), length(d2))
  d1 = try(as.Date(d1), silent = T) %>% verify('Date', err_msg = "Could not convert argument 'd1' to date", err_src = 'gener::diffDate') %>% vect.extend(mm)
  d2 = try(as.Date(d2), silent = T) %>% verify('Date', err_msg = "Could not convert argument 'd2' to date", err_src = 'gener::diffDate') %>% vect.extend(mm)
  if(is.empty(d1) | is.empty(d2)){return(numeric(length(d1)))}
  if(!is.empty(holidays)){holidays = try(as.Date(holidays), silent = T) %>% verify('Date', err_msg = "Could not convert argument 'holidays' to date", err_src = 'gener::diffDate')}
  
  w1 = which(d1 %>% is.holiday)
  while(length(w1) > 0){
    d1[w1] <- d1[w1] + 1
    w1 = which(d1 %>% is.holiday)
  }
  
  w2 = which(d2 %>% is.holiday)
  while(length(w2) > 0){
    d2[w2] <- d2[w2] + 1
    w2 = which(d2 %>% is.holiday)
  }
  
  dl <- min(d1, d2)
  # while (dl %>% is.holiday){dl = dl + 1}
  #
  # d2[d2 < dl] <- dl
  # d1[d1 < dl] <- dl
  #
  dh <- max(d1, d2)
  # while (dh %>% is.holiday){dh = dh + 1}
  
  d <- seq(dl, dh, 1)
  d <- d[!(d %>% is.holiday)]
  
  func = function(x){max(which(x >= d))}
  
  n1 = d1 %>% sapply(func)
  n2 = d2 %>% sapply(func)
  
  return(n1 - n2)
}

# Gives difference between two given times excluding weekends and list of given holidays
#' @export
diffTime = function(t1, t2, units = 'hours', holidays = c()){
  units %>% verify('character', domain = c('hours', 'mins', 'secs', 'days'), default = 'hours', varname = 'units', err_src = 'gener::diffTime')
  unitsPerDay = c(hours = 24, mins = 24*60, secs = 24*3600, days = 1)
  
  t1 %<>% as.time(target_class = 'POSIXct') %>% setTZ('GMT')
  t2 %<>% as.time(target_class = 'POSIXct') %>% setTZ('GMT')
  dd1 = as.Date(t1)
  dd2 = as.Date(t2)
  tt1 = ifelse(dd1 %>% is.holiday, 0, 3600*lubridate::hour(t1) +  60*lubridate::minute(t1) + lubridate::second(t1))
  tt2 = ifelse(dd2 %>% is.holiday, 0, 3600*lubridate::hour(t2) +  60*lubridate::minute(t2) + lubridate::second(t2))
  diffDate(dd1, dd2, holidays = holidays)*unitsPerDay[units] +  (tt1 - tt2)*unitsPerDay[units]/(24*3600)
}


#' @export
addDate <- function(ds, nd, holidays = c()) {
  mm = max(length(ds), length(nd))
  d1 = try(as.Date(ds), silent = T) %>% verify('Date', err_msg = "Could not convert argument 'ds' to date", err_src = 'gener::addDate') %>% vect.extend(mm)
  nd = try(as.integer(nd), silent = T) %>% verify('integer', err_msg = "Could not convert argument 'nd' to integer", err_src = 'gener::addDate') %>% vect.extend(mm) %>% na2zero
  
  
  if(is.empty(d1) | is.empty(nd)){return(d1)}
  if(!is.empty(holidays)){holidays = try(as.Date(holidays), silent = T) %>% verify('Date', err_msg = "Could not convert argument 'holidays' to date", err_src = 'gener::addDate')}
  
  out = d1 + 2*nd - diffDate(d1 + nd, d1, holidays = holidays)
  
  for (i in sequence(length(out))){
    while ((format(out[i], "%u") %in% c("6", "7")) | (out[i] %in% holidays)){out[i] = out[i] + 1}
  }
  return(out)
}


#' @export
addTime <- function(ts, nt, units = 'hours', holidays = c()){
  units %>% verify('character', domain = c('hours', 'mins', 'secs', 'days'), default = 'hours', varname = 'units', err_src = 'gener::diffTime')
  mm = max(length(ts), length(nt))
  ts = try(as.time(ts), silent = T) %>% verify('POSIXct', err_msg = "Could not convert argument 'ts' to time", err_src = 'gener::addTime') %>% vect.extend(mm)
  nt = try(as.numeric(nt), silent = T) %>% verify('numeric', err_msg = "Could not convert argument 'nt' to numeric", err_src = 'gener::addTime') %>% vect.extend(mm) %>% na2zero
  
  if(is.empty(ts) | is.empty(nt)){return(ts)}
  
  unitsPerDay = c(hours = 24, mins = 24*60, secs = 86400, days = 1)
  
  t1 = ts %>% as.time(target_class = 'POSIXct')
  t2 = t1 + nt*86400/unitsPerDay[units]
  
  d1 = t1 %>% format('%Y-%m-%d') %>% as.Date
  d2 = t2 %>% format('%Y-%m-%d') %>% as.Date
  
  dd = as.numeric(d2 - d1)
  d2 = d1 %>% addDate(as.numeric(d2 - d1), holidays = holidays)
  dd = as.numeric(d2 - d1) - dd
  
  out = t2 + 86400*dd
  
  return(out)
}


#' This function, changes the column names of the given table (data.frame) and checks for the column classes
#'
#' @param dataset A table(data.frame) on which the column names is going to change
#' @param columns A named list of characters indicating which names should change to what.
#' List values should be a subset of original column names and list names are names to which the original column names should change
#' @param classes a named list of characters inducating the required classes of each column.
#' List names should be the new column names and values are classes required for each column
#' @param convert A boolean indicating if columns should be coerced to the required class. Function returns error is coercing is not possible.
#' Defaulted is TRUE. If FALSE, the function verifies the column classes and returns error if don't match.
#' @return A table(data.frame) same class as the input table \code{dataset}
#' @examples
#' Data = data.frame(x = c(1,2,3), y = c('10', '20', '30'), z = c('One', 'two', 'three'), stringsAsFactors = F)
#' Data %>% nameColumns(columns = list(XX = 'x', YY = 'y', ZZ = 'z'), classes = list(XX = 'integer', YY = 'numeric'))
#'
#' @export
nameColumns = function(dataset, columns, classes = NULL, convert = T){
  verify(dataset, 'data.frame')
  verify(classes, 'list')
  nms = names(columns)
  assert(length(nms) == length(columns), "Argument 'column' must be a named list! All elements must be named.")
  
  for (i in nms){
    e = columns[[i]] %>% verify(c('character', 'factor', 'numeric', 'integer'), lengths = 1)
    if(!is.empty(e)){
      if((i %in% names(dataset)) & (e != i)){dataset[, i] <- NULL}
      if (inherits(e, c('character', 'factor'))){
        w = which(names(dataset) == columns[[i]])
        assert(length(w) > 0, 'Column ' %++% e %++% ' does not exist in the dataset!')
        assert(length(w) == 1, 'Column ' %++% e %++% ' appears more than once in the dataset!')
        names(dataset)[w] <- i
      } else if(inherits(e, c('numeric', 'integer'))){
        assert(ncol(dataset) >= e, 'Given dataset does not have ' %++% e %++% ' columns!')
        names(dataset)[e] <- i
      }
    }
  }
  
  dataset = dataset[, !duplicated(names(dataset)), drop = F]
  # Test classes:
  if(!is.null(classes)){
    classes %<>% verify('list') %>% list.extract(names(classes) %^% colnames(dataset))
    for(i in names(classes)){
      e = classes[[i]]
      correct = dataset %>% pull(i) %>% inherits(e)
      if(!correct){
        ermsg = 'Column ' %++% i %++% ' is expected to be of class ' %++% paste(e, collapse = ' or ') %++% '. Was class ' %++% paste(class(dataset[, i]), collapse = ' --> ') %++% '!'
        if(convert){
          dataset[, i] <- try(dataset %>% pull(i) %>% coerce(e[1]), silent = T) %>% verify(e, err_msg = ermsg %++% ' Coerction attempt failed.')
        } else {
          assert(FALSE, ermsg)
        }
      }
    }
  }
  
  return(dataset)
}

# Converts a data.frame to list. Respects rownames and colnames
#' @export
df2List = function(df){
  lll  = list()
  rows = rownames(df)
  # if(is.null(rows)){rows = df %>% nrow %>% sequence}
  for(i in rows){lll[[i]] <- df[i,] %>% as.list}
  return(lll)
}

#' @export
applyFunctionList = function(x, L){y = x; for (fun in L){y = fun(y)}; return(y)}


# secondsFromDayStart = function(time){
#   if(require('lubridate')){
#     return(3600*hour(time) + 60*minute(time) + second(time))
#   } else {}
# }

#' @export
list.default = function(l, ...){
  lst = list(...)
  
  for(i in names(lst)){
    if(is.null(l[[i]])){
      l[[i]] = lst[[i]]
    }
  }
  
  return(l)
}

# todo: add reverse as argument
#' @export
charFilter = function(str, ..., and = T, match_case = F){
  pat = c(...) %>% verify('character')
  if(!match_case){
    strl = str %>% tolower
    patl = pat %>% tolower
  } else {
    strl = strl
    patl = pat
  }
  w = chif(and, strl %>% length %>% sequence, c())
  for(p in patl){
    w = chif(and, intersect(w, grep(p, strl)), union(w, grep(p, strl)))
  }
  str[w]
}

#' @export
partition = function(tbl, ratio = 0.7){
  N   = nrow(tbl)
  ind = N %>% sequence %>% sample(size = floor(ratio*N), replace = F)
  list(part1 = tbl[ind, ], part2 = tbl[- ind, ])
}



#' @export
fmap = c(normal = 'norm', gaussian = 'norm', 
         beta = 'beta', binomial = 'binom', 
         'chi-squared' = 'chisq', chi = 'chisq', chi_squared = 'chisq', chisquared = 'chisq', 
         exponential = 'exp', 
         geometric = 'geom', hypergeomrtric = 'hyper', hyper_geometric = 'hyper', 'hyper-geometric' = 'hyper',
         logistic = 'logis',
         lognormal = 'lnorm', 'log-normal' = 'lnorm', log_normal = 'lnorm',
         poisson = 'pois', 
         student_t = 't', 'student_t' = 't',
         studentized_range = 'tukey', 'studentized-range' = 'tukey',
         uniform = 'unif', 
         wilcoxon_rank_sum = 'wilcox', wilcoxon = 'wilcox', wilcoxon_signed_rank = 'signrank',
         wilcoxonranksum = 'wilcox', ranksom = 'wilcox', rank_som = 'wilcox',
         wilcoxonsignedrank = 'signrank',signed_rank = 'signrank', signedrank = 'signrank')

#' @export
pdf = function(family, ...){
  family  %<>% tolower 
  if(family %in% names(fmap)){family = fmap[family]}
  parse(text = paste0('d', family, '(...)')) %>% eval
}

#' @export
cdf = function(family, ...){
  family  %<>% tolower 
  if(family %in% names(fmap)){family = fmap[family]}
  parse(text = paste0('p', family, '(...)')) %>% eval
}

#' @export
cdf.inv = function(family, ...){
  family  %<>% tolower 
  if(family %in% names(fmap)){family = fmap[family]}
  parse(text = paste0('q', family, '(...)')) %>% eval
}

#' @export
gen.random = function(family, ...){
  family  %<>% tolower 
  if(family %in% names(fmap)){family = fmap[family]}
  parse(text = paste0('r', family, '(...)')) %>% eval
}

#' @export  
string2factor = function(df){
  cls = colnames(df) %>% sapply(function(i) class(df[,i])[1])
  for(i in which(cls == 'character')) df[,i] %<>% as.factor 
  return(df)
}

#' @export
factor2integer = function(df){
  cls = colnames(df) %>% sapply(function(i) class(df[,i])[1])
  for(i in which(cls == 'factor')) df[,i] %<>% as.integer 
  return(df)
}

#' @export
numeric2time = function(df){
  cls = colnames(df) %>% sapply(function(i) class(df[,i])[1])
  for(i in which(cls %in% c('integer', 'numeric'))) df[,i] %<>% as.POSIXct(origin = '1970-01-01')
  return(df)
}

#' @export
integer2numeric = function(df){
  cls = colnames(df) %>% sapply(function(i) class(df[,i])[1])
  for(i in which(cls %in% c('integer'))) df[,i] %<>% as.numeric
  return(df)
}

#' @export
logical2integer = function(df){
  cls = colnames(df) %>% sapply(function(i) class(df[,i])[1])
  for(i in which(cls == 'logical')) df[,i] %<>% as.integer 
  return(df)
}
