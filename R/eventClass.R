#' Eventlist
#' 
#' An eventlist is a data.frame with the necessary parameters to define an event as columns and specific events as rows.
#' Event time and value can be passed as parameters, which can also be estimated.
#' 
#' The function `addEvent` is pipe-friendly
#'  
#' @param var Character, the state to which the event is applied
#' @param time Character or Numeric, the time at which the event happens
#' @param value Character or Numeric, the value of the event
#' @param root Character or NA, condition to trigger the event instead of time
#' @param method Character, options are "replace", "add" or "multiply"
#'
#' @return data.frame with class eventlist
#' @export
#'
#' @examples
#' \dontrun{
#'   eventlist(var = "A", time = "5", value = 1, method = "add")
#'   
#'   events <- addEvent(NULL, var = "A", time = "5", value = 1, method = "add")
#'   events <- addEvent(events, var = "A", time = "10", value = 1, method = "add")
#'   
#'   # With symbols
#'   events <- eventlist()
#'   # Set "A" to "value_switch" at time "time_Switch"
#'   events <- addEvent(events, var = "A", time = "time_switch",
#'                      value = "value_switch", method = "replace")
#'   # Set "B" to 2 when "A" reaches "A_target". The time-parameter for
#'   # internal use will be "time_root".
#'   events <- addEvent(events, var = "B", time = "time_root", value = 2,
#'                      root = "A - A_target", method = "replace")
#' }
eventlist <- function(var = NULL, time = NULL, value = NULL, root = NULL, method = NULL) {

  # If any non-root field is supplied, default root to NA so all columns share
  # the same length. Mirrors as.eventlist.list() which fills root with NA when
  # absent. Without this, eventlist(var = "A", time = 5, value = "x", method = "add")
  # crashes in data.frame() with a 0-vs-1 length conflict.
  if (is.null(root) &&
      (!is.null(var) || !is.null(time) || !is.null(value) || !is.null(method))) {
    root <- NA
  }

  out <- data.frame(var = var,
                    time = time,
                    value = value,
                    root = root,
                    method = method,
                    stringsAsFactors = FALSE)
  
  class(out) <- c("eventlist", "data.frame")
  return(out)
  
}

#' Coerce to eventlist
#'
#' @param ... not used
#' @export
as.eventlist <- function(x, ...) {
  UseMethod("as.eventlist", x)
}


#' @export
#' @rdname as.eventlist
#' @param x list, data.frame
as.eventlist.list <- function(x, ...) {
  
  # Check names
  required <- c("var", "time", "value", "method")
  if (!all(required %in% names(x)))
    stop("x needs to provide var, time, value, and method.")
  
  # Check for optional root entry
  if (!"root" %in% names(x))
    x[["root"]] <- rep(NA, length(x[["var"]]))
  
  # Convert to list of characters and numeric
  x <- lapply(x[c(required, "root")], function(element) {
    if (!is.numeric(element)) as.character(element) else element
  })
  
  # Return eventlist object
  do.call(eventlist, x)
  
  
}

#' @export
#' @rdname as.eventlist
as.eventlist.data.frame <- function(x, ...) {

  as.eventlist.list(as.list(x))
    
}

#' @rdname eventlist
#' @param event object of class `eventlist`
#' @param ... not used
#' @export
addEvent <- function(event, var, time = 0, value = 0, root = NA, method = "replace", ...) {
  
  UseMethod("addEvent", event)
  
}

#' @export
addEvent.eventlist <- function(event, var, time = 0, value = 0, root = NA, method = "replace", ...) {
  
  event.new <- data.frame(var = var, time = time, value = value, root = root, method = method, stringsAsFactors = FALSE)
  as.eventlist(rbind(event, event.new))
  
}

#' @export
addEvent.NULL <- function(event, var, time = 0, value = 0, root = NA, method = "replace", ...) {
  
  event.new <- data.frame(var = var, time = time, value = value, root = root, method = method, stringsAsFactors = FALSE)
  as.eventlist(rbind(event, event.new))
  
}
