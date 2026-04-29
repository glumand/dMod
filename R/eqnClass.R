
## Class "eqnlist" and its constructor ------------------------------------------


#' Coerce to an equation list
#' @description Translates a reaction network, e.g. defined by a data.frame, into an equation list object.
#' @param ... additional arguments to be passed to or from methods.
#' @details If `data` is a `data.frame`, it must contain columns "Description" (character),
#' "Rate" (character), and one column per ODE state with the state names.
#' The state columns correspond to the stoichiometric matrix.
#' @return Object of class [eqnlist]
#' @rdname eqnlist
#' @export
as.eqnlist <- function(data, volumes, ...) {
  UseMethod("as.eqnlist", data)
}

#' @export
#' @param data data.frame with columns Description, Rate, and one colum for each state
#' reflecting the stoichiometric matrix
#' @rdname eqnlist
as.eqnlist.data.frame <- function(data, volumes = NULL, compartments = NULL, compartmentOf = NULL,
                                   reactionCompartment = NULL, ...) {
  description <- as.character(data$Description)
  rates <- as.character(data$Rate)
  states <- setdiff(colnames(data), c("Description", "Rate"))
  smatrix <- as.matrix(data[, states]); colnames(smatrix) <- states

  if (is.null(volumes))             volumes             <- attr(data, "volumes")
  if (is.null(compartments))        compartments        <- attr(data, "compartments")
  if (is.null(compartmentOf))       compartmentOf       <- attr(data, "compartmentOf")
  if (is.null(reactionCompartment)) reactionCompartment <- attr(data, "reactionCompartment")

  eqnlist(smatrix, states, rates, volumes, description,
          compartments = compartments, compartmentOf = compartmentOf,
          reactionCompartment = reactionCompartment)

}


#' @export
#' @rdname eqnlist
#' @param x object of class `eqnlist`
is.eqnlist <- function(x) {

  expected_names <- c("smatrix", "states", "rates", "volumes", "description",
                      "compartments", "compartmentOf", "reactionCompartment")

  #Empty list
  if (is.null(x$smatrix)) {
    if (length(x$states) == 0 &&
        length(x$rates) == 0 &&
        is.null(x$volumes) &&
        length(x$description) == 0 &&
        is.null(x$compartments) &&
        is.null(x$compartmentOf) &&
        is.null(x$reactionCompartment)
    ) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    #Non-empty list
    rc_ok <- is.null(x$reactionCompartment) ||
             (length(x$reactionCompartment) == length(x$rates) &&
              all(is.na(x$reactionCompartment) | x$reactionCompartment %in% names(x$compartments)))
    if (inherits(x, "eqnlist") &&
        all(names(x) == expected_names) &&
        all(names(x$smatrix) == names(x$states)) &&
        dim(x$smatrix)[1] == length(x$rates) &&
        dim(x$smatrix)[2] == length(x$states) &&
        is.matrix(x$smatrix) &&
        !is.null(x$compartments) && !is.null(x$compartmentOf) &&
        all(x$compartmentOf %in% names(x$compartments)) &&
        all(x$states %in% names(x$compartmentOf)) &&
        rc_ok
    ) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}


## Class "eqnlist" and its methods ------------------------------------------

#' Determine conserved quantites by finding the kernel of the stoichiometric
#' matrix
#'
#' @param S Stoichiometric matrix
#' @param weight One of `"none"` (default) or `"volume"`. When `"volume"`, the
#'   columns of `S` are multiplied by their compartment volume before the kernel
#'   is computed, so the returned quantities are conserved in *amount* rather
#'   than concentration. Requires `volumes` to be supplied as numeric.
#' @param volumes Optional named numeric vector of volume values keyed by state,
#'   aligned with `colnames(S)`. Only consulted when `weight = "volume"`.
#' @return Data frame with conserved quantities carrying an attribute with the
#'   number of conserved quantities.
#' @author Malenke Mader, \email{Malenka.Mader@@fdm.uni-freiburg.de}
#'
#' @example inst/examples/equations.R
#' @export
conservedQuantities <- function(S, weight = c("none", "volume"), volumes = NULL) {
  weight <- match.arg(weight)
  # Get kernel of S
  S[is.na(S)] <- 0
  if (weight == "volume") {
    if (is.null(volumes))
      stop("`weight = \"volume\"` requires a named numeric `volumes` argument.")
    if (is.null(colnames(S)))
      stop("`S` must have column names when `weight = \"volume\"`.")
    missing_vol <- setdiff(colnames(S), names(volumes))
    if (length(missing_vol) > 0L)
      stop("`volumes` missing entries for: ", paste(missing_vol, collapse = ", "))
    v_num <- suppressWarnings(as.numeric(volumes[colnames(S)]))
    if (anyNA(v_num))
      stop("`weight = \"volume\"` requires all volumes to be numeric; got symbolic expression(s).")
    S <- sweep(S, 2, v_num, "*")
  }
  v <- nullZ(S)
  n_cq <-  ncol(v)
  
  # Iterate over conserved quantities, removes 0s, etc.
  if (n_cq > 0) {
    if (is.null(colnames(S))) stop("Columns of stoichiometric matrix not named.") else variables <- colnames(S)
    cq <- matrix(nrow = ncol(v), ncol = 1)
    for (iCol in 1:ncol(v)) {
      is.zero <- v[, iCol] == 0
      cq[iCol, 1] <- sub("+-", "-", paste0(v[!is.zero, iCol], "*", variables[!is.zero], collapse = "+"), fixed = TRUE)
    }
    
    colnames(cq) <- paste0("Conserved quantities: ", n_cq)
    cq <- as.data.frame(cq)
    attr(x = cq, which = "n") <- n_cq
    
  } else {
    cq <- c()
  }
  
  return(cq)
}



#' Generate a table of reactions (data.frame) from an equation list
#' 
#' @param eqnlist object of class [eqnlist]
#' @return `data.frame` with educts, products, rate and description. The first
#' column is a check if the reactions comply with reaction kinetics.
#' 
#' @example inst/examples/equations.R
#' @export
getReactions <- function(eqnlist) {
  
  # Extract information from eqnlist
  S <- eqnlist$smatrix
  rates <- eqnlist$rates
  description <- eqnlist$description
  variables <- eqnlist$states
  
  # Determine lhs and rhs of reactions
  if(is.null(S)) return()
  
  reactions <- apply(S, 1, function(v) {
    
    numbers <- v[which(!is.na(v))]
    educts <- -numbers[numbers < 0]
    #educts[which(educts == 1)] <- " "
    products <- numbers[numbers > 0]
    #products[which(products == 1)] <- 
    educts <- paste(paste(educts, names(educts), sep = "*"), collapse=" + ")
    products <- paste(paste(products, names(products), sep = "*"), collapse=" + ")
    educts <- gsub("1*", "", educts, fixed = TRUE)
    products <- gsub("1*", "", products, fixed = TRUE)
    
    reaction <- paste(educts, "->", products)
    return(c(educts, products))
    
  })
  educts <- reactions[1,]
  products <- reactions[2,]
  
  # Check for consistency
  exclMarks.logical <- unlist(lapply(1:length(rates), function(i) {
    
    myrate <- rates[i]
    parsedRate <- getParseData(parse(text=myrate, keep.source = TRUE))
    symbols <- parsedRate$text[parsedRate$token=="SYMBOL"]
    
    educts <- variables[which(S[i,]<0)]
    
    !all(unlist(lapply(educts, function(e) any(e==symbols))))
    
  }))
  exclMarks <- rep(" ", ncol(reactions))
  exclMarks[exclMarks.logical] <- "!"
  

  # Generate data.frame  
  out <- data.frame(exclMarks, educts, "->", products, rates, description, stringsAsFactors = FALSE)
  colnames(out) <- c("Check", "Educt",  "->",  "Product", "Rate", "Description")
  rownames(out) <- 1:nrow(out)
  
  return(out)
  
}


#' Add reaction to reaction table
#'
#' @param eqnlist equation list, see [eqnlist]
#' @param from character with the left hand side of the reaction, e.g. "2*A + B"
#' @param to character with the right hand side of the reaction, e.g. "C + 2*D"
#' @param rate character. The rate associated with the reaction. The name is employed as a description
#' of the reaction.
#' @param description Optional description instead of `names(rate)`.
#' @param compartment Character, compartment ID to which any *new* states introduced
#' by this reaction are assigned. Defaults to `"defaultComp"`. If the compartment does
#' not yet exist on `eqnlist`, it is created with volume `"1"`.
#' @param rateCompartment Optional compartment ID naming the frame in which `rate`
#' is a concentration-rate. Required when educts span multiple compartments
#' (e.g. membrane binding `L_ext + R_cyt -> Complex`); leave as `NA` (the
#' default) to let [getFluxes()] infer the frame from the educts.
#' @return An object of class [eqnlist].
#' @examples
#' f <- eqnlist()
#' f <- addReaction(f, "2*A+B", "C + 2*D", "k1*B*A^2")
#' f <- addReaction(f, "C + A", "B + A", "k2*C*A")
#'
#'
#' @example inst/examples/equations.R
#' @export
#' @rdname addReaction
addReaction <- function(eqnlist, from, to, rate, description = names(rate),
                         compartment = "defaultComp", rateCompartment = NA_character_) {


  if (missing(eqnlist)) eqnlist <- eqnlist()

  volumes <- eqnlist$volumes
  compartments_in <- eqnlist$compartments
  compartmentOf_in <- eqnlist$compartmentOf
  reactionCompartment_in <- eqnlist$reactionCompartment

  # Analyze the reaction character expressions
  educts <- getSymbols(from)
  eductCoef <- 0
  if(length(educts) > 0) eductCoef <- sapply(educts, function(e) sum(getCoefficients(from, e)))
  products <- getSymbols(to)
  productCoef <- 0
  if(length(products) > 0) productCoef <- sapply(products, function(p) sum(getCoefficients(to, p)))


  # States introduced by this reaction
  states <- unique(c(educts, products))

  # Description
  if(is.null(description)) description <- ""

  # Stoichiometric matrix
  smatrix <- matrix(NA, nrow = 1, ncol=length(states)); colnames(smatrix) <- states
  if(length(educts)>0) smatrix[,educts] <- -eductCoef
  if(length(products)>0) {
    filled <- !is.na(smatrix[,products])
    smatrix[,products[filled]] <- smatrix[,products[filled]] + productCoef[filled]
    smatrix[,products[!filled]] <- productCoef[!filled]
  }


  smatrix[smatrix == "0"] <- NA


  # data.frame
  mydata <- cbind(data.frame(Description = description, Rate = as.character(rate)), as.data.frame(smatrix))
  row.names(mydata) <- NULL


  if(!is.null(eqnlist)) {
    mydata0 <- as.data.frame(eqnlist)
    mydata <- combine(mydata0, mydata)
  }

  # Extend compartment assignment for brand-new states with the `compartment` arg.
  new_states <- setdiff(states, names(compartmentOf_in))
  compartments_out <- compartments_in
  compartmentOf_out <- compartmentOf_in
  if (length(new_states) > 0L) {
    if (is.null(compartments_out)) compartments_out <- list()
    if (is.null(compartmentOf_out)) compartmentOf_out <- character(0)
    if (!compartment %in% names(compartments_out)) {
      compartments_out[[compartment]] <- list(volume = "1", rule = NULL)
    }
    compartmentOf_out <- c(compartmentOf_out,
                           setNames(rep(compartment, length(new_states)), new_states))
  }

  # Extend reactionCompartment with the value for this new reaction. When the
  # input list has no annotations (NULL), pad with NA for the existing rates
  # so the final vector lines up with the combined data.frame rows.
  existing_n <- length(eqnlist$rates)
  if (is.null(reactionCompartment_in)) reactionCompartment_in <- rep(NA_character_, existing_n)
  reactionCompartment_out <- c(reactionCompartment_in, as.character(rateCompartment))
  if (all(is.na(reactionCompartment_out))) reactionCompartment_out <- NULL

  as.eqnlist(mydata, volumes = volumes,
             compartments = compartments_out, compartmentOf = compartmentOf_out,
             reactionCompartment = reactionCompartment_out)

}


#' Generate list of fluxes from equation list
#' 
#' @param eqnlist object of class [eqnlist].
#' @param type "conc." or "amount" for fluxes in units of concentrations or
#' number of molecules. 
#' @return list of named characters, the in- and out-fluxes for each state.
#' @example inst/examples/equations.R
#' @export
getFluxes <- function(eqnlist, type = c("conc", "amount")) {

  type <- match.arg(type)

  description <- eqnlist$description
  rate <- eqnlist$rates
  variables <- eqnlist$states
  SMatrix <- eqnlist$smatrix
  compartments <- eqnlist$compartments
  compartmentOf <- eqnlist$compartmentOf
  reactionCompartment <- eqnlist$reactionCompartment

  if (is.null(SMatrix)) return()

  # Defensive fallback: an eqnlist constructed outside our constructor may have
  # NULL compartment info. Treat every state as living in an implicit "defaultComp"
  # compartment with volume "1" — matches legacy behavior.
  if (is.null(compartments) || is.null(compartmentOf)) {
    compOf <- setNames(rep("defaultComp", length(variables)), variables)
    compartments <- list(defaultComp = list(volume = "1", rule = NULL))
  } else {
    compOf <- compartmentOf[variables]
  }
  volumes <- vapply(compOf, function(cid) compartments[[cid]]$volume, character(1))
  names(volumes) <- variables

  # Resolve per-reaction reference compartment V_ref (concentration-rate frame).
  # Priority: (1) user-supplied reactionCompartment[i] if non-NA, (2) unique
  # educt compartment, (3) unique product compartment for pure synthesis.
  # When educts span multiple compartments and no annotation is given, we
  # error with a clear message pointing the user at `reactionCompartment`.
  nR <- nrow(SMatrix)
  vref_cid <- character(nR)
  for (i in seq_len(nR)) {
    if (!is.null(reactionCompartment) && !is.na(reactionCompartment[i])) {
      vref_cid[i] <- reactionCompartment[i]
      next
    }
    row_i <- SMatrix[i, ]
    educt_idx <- which(!is.na(row_i) & row_i < 0)
    product_idx <- which(!is.na(row_i) & row_i > 0)
    cand <- if (length(educt_idx) > 0) unique(compOf[educt_idx])
            else if (length(product_idx) > 0) unique(compOf[product_idx])
            else character(0)
    if (length(cand) == 1L) {
      vref_cid[i] <- cand
    } else if (length(cand) > 1L) {
      stop(sprintf(
        "Reaction %d (\"%s\") spans compartments (%s). Pass `reactionCompartment` to name the frame in which the rate is a concentration-rate.",
        i, description[i], paste(cand, collapse = ", ")))
    } else {
      stop(sprintf("Reaction %d (\"%s\") has no species; cannot determine reference compartment.",
                   i, description[i]))
    }
  }
  vref_vol <- vapply(vref_cid, function(cid) compartments[[cid]]$volume, character(1))

  # generate equation expressions
  terme <- lapply(1:length(variables), function(j) {
    v <- SMatrix[,j]
    nonZeros <- which(!is.na(v))
    var.description <- description[nonZeros]
    positives <- which(v > 0)
    destin_cid <- compOf[[j]]
    destin_vol <- volumes[[j]]

    # Uniform flux formula: flux_X = stoich_X * rate * V_ref / V_X for every
    # species in every reaction. For single-educt-compartment reactions this
    # is equivalent to the legacy asymmetric formula.
    switch(type,
           conc = {
             volumes.ratios <- paste0("*(", vref_vol, "/", destin_vol, ")")
             volumes.ratios[vref_cid == destin_cid] <- ""
           },
           amount = {
             volumes.ratios <- paste0("*(", vref_vol, ")")
           }
    )

    numberchar <- as.character(v)
    if (nonZeros[1] %in% positives) {
      numberchar[positives] <- paste(c("", rep("+", length(positives)-1)), numberchar[positives], sep = "")
    } else {
      numberchar[positives] <- paste("+", numberchar[positives], sep = "")
    }
    var.flux <- paste0(numberchar[nonZeros], "*(", rate[nonZeros], ")", volumes.ratios[nonZeros])
    names(var.flux) <- var.description

    # Dilution term: if state j's compartment has a non-null volume rule,
    # d[X]/dt picks up -[X]*(dV/dt)/V (SBML concentration-correction). Always
    # zero in the constant-volume case because `rule` is NULL there.
    r <- compartments[[destin_cid]]$rule
    if (!is.null(r) && nzchar(r)) {
      dilution <- paste0("-(", variables[j], ")*(", r, ")/(", destin_vol, ")")
      names(dilution) <- paste0("dilution_", destin_cid)
      var.flux <- c(var.flux, dilution)
    }

    return(var.flux)
  })

  fluxes <- terme
  names(fluxes) <- variables

  return(fluxes)


}



#' Symbolic time derivative of equation vector given an equation list
#' 
#' The time evolution of the internal states is defined in the equation list.
#' Time derivatives of observation functions are expressed in terms of the
#' rates of the internal states.
#' 
#' @param observable named character vector or object of type [eqnvec]
#' @param eqnlist equation list
#' @details Observables are translated into an ODE
#' @return An object of class [eqnvec]
#' @example inst/examples/equations.R
#' @export
dot <- function(observable, eqnlist) {
  
 
  # Analyze the observable character expression
  symbols <- getSymbols(observable)
  states <- intersect(symbols, eqnlist$states)
  derivatives <- lapply(observable, function(obs) {
    out <- lapply(as.list(states), function(x) paste(deparse(D(parse(text=obs), x), width.cutoff = 500),collapse=""))
    names(out) <- states
    return(out)
  })
  
  # Generate equations from eqnist
  f <- as.eqnvec(eqnlist)
  
  newodes <- sapply(derivatives, function(der) {
    
    prodSymb(matrix(der, nrow = 1), matrix(f[names(der)], ncol = 1))
    
#     
#     out <- sapply(names(der), function(n) {
#       d <- der[n]
#       
#       if (d != "0") {
#         prodSymb(matrix(d, nrow = 1), matrix(f[names(d)], ncol = 1))
#       } else  {
#         return("0")
#       }
#         
#       
#       
#       #paste( paste("(", d, ")", sep="") , paste("(", f[names(d)], ")",sep=""), sep="*") else return("0")
#     })
#     out <- paste(out, collapse = "+")
#     
#     return(out)
    
  })
  
  as.eqnvec(newodes)
}



#' Coerce equation list into a data frame
#' 
#' @param x object of class [eqnlist]
#' @param ... other arguments
#' @return a `data.frame` with columns "Description" (character), 
#' "Rate" (character), and one column per ODE state with the state names. 
#' The state columns correspond to the stoichiometric matrix.
#' @export
as.data.frame.eqnlist <- function(x, ...) {

  eqnlist <- x

  if(is.null(eqnlist$smatrix)) return()

  data <- data.frame(Description = eqnlist$description,
                     Rate = eqnlist$rate,
                     eqnlist$smatrix,
                     stringsAsFactors = FALSE)

  attr(data, "volumes") <- eqnlist$volumes
  attr(data, "compartments") <- eqnlist$compartments
  attr(data, "compartmentOf") <- eqnlist$compartmentOf
  attr(data, "reactionCompartment") <- eqnlist$reactionCompartment

  return(data)
}

#' Write equation list into a csv file
#' 
#' @param eqnlist object of class [eqnlist]
#' @param ... Arguments going to [write.table][utils::write.table]
#' 
#' @export
#' @importFrom utils file.edit getParseData install.packages installed.packages read.csv str tail write.csv
write.eqnlist <- function(eqnlist, ...) {
  
  
  arglist <- list(...)
  argnames <- names(arglist)
  if (!"row.names" %in% argnames) arglist$row.names <- FALSE
  if (!"na" %in% argnames) arglist$na <- ""
  
  arglist$x <- as.data.frame(eqnlist)
  
  do.call(write.csv, arglist)
  
}


#' subset of an equation list
#' 
#' @param x the equation list
#' @param ... logical expression for subsetting
#' @details The argument `...` can contain "Educt", "Product", "Rate" and "Description".
#' The "%in%" operator is modified to allow searches in Educt and Product (see examples).
#' 
#' @return An object of class [eqnlist]
#' @examples
#' reactions <- data.frame(Description = c("Activation", "Deactivation"), 
#'                         Rate = c("act*A", "deact*pA"), A=c(-1,1), pA=c(1, -1) )
#' f <- as.eqnlist(reactions)
#' subset(f, "A" %in% Educt)
#' subset(f, "pA" %in% Product)
#' subset(f, grepl("act", Rate))
#' @export subset.eqnlist
#' @export
subset.eqnlist <- function(x, ...) {
  
  eqnlist <- x
  
  # Do selection on data.frame
  data <- getReactions(eqnlist)
  if(is.null(data)) return()
  
  data.list <- list(Educt = lapply(data$Educt, getSymbols), 
                    Product = lapply(data$Product, getSymbols),
                    Rate = data$Rate,
                    Description = data$Description,
                    Check = data$Check)
  
  "%in%" <- function(x, table) sapply(table, function(mytable) any(x == mytable))
  select <- which(eval(substitute(...), data.list))
  if (length(select) == 0) return(NULL)
  
  # Translate subsetting on eqnlist entries
  # smatrix
  smatrix <- submatrix(eqnlist$smatrix, rows = select)
  empty <- sapply(1:ncol(smatrix), function(i) all(is.na(smatrix[, i])))
  smatrix <- submatrix(smatrix, cols = !empty)
  
  # states and rates
  states <- colnames(smatrix)
  rates <- eqnlist$rates[select]

  # volumes (derived view; filter to surviving states)
  volumes <- eqnlist$volumes
  if(!is.null(volumes)) volumes <- volumes[intersect(names(volumes),  states)]

  # compartments/compartmentOf: restrict to surviving states and drop unreferenced compartments.
  # `%in%` is locally shadowed above; use base::`%in%` explicitly.
  compartmentOf <- eqnlist$compartmentOf
  compartments <- eqnlist$compartments
  if (!is.null(compartmentOf)) {
    compartmentOf <- compartmentOf[intersect(names(compartmentOf), states)]
    if (!is.null(compartments)) {
      used_cids <- unique(compartmentOf)
      compartments <- compartments[base::`%in%`(names(compartments), used_cids)]
    }
  }

  # description
  description <- eqnlist$description[select]

  reactionCompartment <- if (!is.null(eqnlist$reactionCompartment)) eqnlist$reactionCompartment[select] else NULL
  if (!is.null(reactionCompartment) && all(is.na(reactionCompartment))) reactionCompartment <- NULL

  eqnlist(smatrix, states, rates, volumes, description,
          compartments = compartments, compartmentOf = compartmentOf,
          reactionCompartment = reactionCompartment)


}


#' Print or pander equation list
#' 
#' @param x object of class [eqnlist]
#' @param pander logical, use pander for output (used with R markdown)
#' @param ... additional arguments
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' @author Daniel Kaschek, \email{daniel.kaschek@@physik.uni-freiburg.de}
#' 
#' @export
print.eqnlist <- function(x, pander = FALSE, ...) {

  eqnlist <- x

  # Entities to print and pander
  cq <- conservedQuantities(eqnlist$smatrix)
  r <- getReactions(eqnlist)

  # Compartment block: only show when there is a meaningful assignment to
  # surface (i.e. more than one compartment, or a compartment with non-unit
  # volume, or any compartment with a rule).
  comp_lines <- .format_compartments(eqnlist$compartments, eqnlist$compartmentOf)

  # Print or pander?
  if (!pander) {
    print(cq)
    if (length(comp_lines) > 0L) {
      cat("\n")
      cat(comp_lines, sep = "\n")
    }
    cat("\n")
    print(r)
  } else {
    pander::panderOptions("table.alignment.default", "left")
    pander::panderOptions("table.split.table", Inf)
    pander::panderOptions("table.split.cells", Inf)
    exclude <- "Check"
    r <- r[, setdiff(colnames(r), exclude)]
    r$Rate <- paste0(format.eqnvec(as.character(r$Rate)))
    pander::pander(r)
  }
}


# Internal: render a compact compartment summary for print.eqnlist.
# Returns character(0) when the model has exactly one compartment whose volume
# is "1" and no rule (the implicit-default case for models that never used
# the compartment feature), so legacy output stays unchanged.
.format_compartments <- function(compartments, compartmentOf) {
  if (is.null(compartments) || is.null(compartmentOf)) return(character(0))
  if (length(compartments) == 1L) {
    only <- compartments[[1L]]
    if (identical(only$volume, "1") && is.null(only$rule)) return(character(0))
  }
  header <- "Compartments:"
  comp_entries <- vapply(names(compartments), function(cid) {
    entry <- compartments[[cid]]
    rule_txt <- if (!is.null(entry$rule) && nzchar(entry$rule)) paste0(", rule=", entry$rule) else ""
    sprintf("  %s (V=%s%s)", cid, entry$volume, rule_txt)
  }, character(1))
  assign_lines <- sprintf("  %s: %s",
                          names(split(names(compartmentOf), compartmentOf)),
                          vapply(split(names(compartmentOf), compartmentOf),
                                 function(sts) paste(sts, collapse = ", "), character(1)))
  c(header, comp_entries, "States by compartment:", assign_lines)
}



## Class "eqnvec" and its constructors --------------------------------------------



#' Coerce to an equation vector
#' 
#' @param x object of class `character` or `eqnlist`
#' @param ... arguments going to the corresponding methods
#' @details If `x` is of class `eqnlist`, [getFluxes] is called and coerced
#' into a vector of equations.
#' @return object of class [eqnvec].
#' @export
as.eqnvec <- function(x, ...) {
  UseMethod("as.eqnvec", x)
}

#' Generate equation vector object
#'
#' @param names character, the left-hand sides of the equation
#' @rdname as.eqnvec
#' @export
as.eqnvec.character <- function(x = NULL, names = NULL, ...) {
  
  equations <- x
  
  if (is.null(equations)) return(NULL)
  
  if (is.null(names)) names <- names(equations)
  if (is.null(names)) stop("equations need names")
  if (length(names) != length(equations)) stop("Length of names and equations do not coincide")
  try.parse <- try(parse(text = equations), silent = TRUE)
  if (inherits(try.parse, "try-error")) stop("equations cannot be parsed: ", try.parse)
  
  out <- structure(equations, names = names)
  class(out) <- c("eqnvec", "character")
  
  return(out)
  
}



#' Transform equation list into vector of equations
#' 
#' @description An equation list stores an ODE in a list format. The function
#' translates this list into the right-hand sides of the ODE.
#' @rdname as.eqnvec
#' @export
as.eqnvec.eqnlist <- function(x, ...) {
  
  eqnlist <- x
  
  terme <- getFluxes(eqnlist, ...)
  if(is.null(terme)) return()
  terme <- lapply(terme, function(t) paste(t, collapse=" "))
  
  
  terme <- do.call(c, terme)
  
  as.eqnvec(terme, names(terme))
  
}

#' @export
c.eqnlist <- function(...) {

  inputs <- list(...)
  inputs <- inputs[!vapply(inputs, function(x) is.null(x) || is.null(x$smatrix), logical(1))]
  if (length(inputs) == 0L) return(eqnlist())

  # Merge stoichiometry / rates / description via the data.frame path
  out <- Reduce(combine, lapply(inputs, as.data.frame))

  # Merge compartments with conflict detection
  all_compartments <- list()
  for (el in inputs) {
    if (is.null(el$compartments)) next
    for (cid in names(el$compartments)) {
      new_entry <- el$compartments[[cid]]
      if (cid %in% names(all_compartments)) {
        old_entry <- all_compartments[[cid]]
        if (!identical(old_entry$volume, new_entry$volume)) {
          stop(sprintf("Compartment conflict: '%s' has volume '%s' in one eqnlist and '%s' in another.",
                       cid, old_entry$volume, new_entry$volume))
        }
        if (!identical(old_entry$rule, new_entry$rule)) {
          stop(sprintf("Compartment conflict: '%s' has different `rule` in the input eqnlists.", cid))
        }
      } else {
        all_compartments[[cid]] <- new_entry
      }
    }
  }

  all_compartmentOf <- character(0)
  for (el in inputs) {
    if (is.null(el$compartmentOf)) next
    for (st in names(el$compartmentOf)) {
      cid <- unname(el$compartmentOf[[st]])
      if (st %in% names(all_compartmentOf)) {
        if (!identical(unname(all_compartmentOf[[st]]), cid)) {
          stop(sprintf("State '%s' assigned to different compartments across input eqnlists.", st))
        }
      } else {
        all_compartmentOf[st] <- cid
      }
    }
  }

  if (length(all_compartments) == 0L) all_compartments <- NULL
  if (length(all_compartmentOf) == 0L) all_compartmentOf <- NULL

  # Concatenate reactionCompartment annotations. If any input has them, we need
  # to produce a vector of length nrow(combined). Missing entries become NA.
  any_rc <- any(vapply(inputs, function(el) !is.null(el$reactionCompartment), logical(1)))
  if (any_rc) {
    all_rc <- unlist(lapply(inputs, function(el) {
      if (is.null(el$reactionCompartment)) rep(NA_character_, length(el$rates))
      else el$reactionCompartment
    }))
    if (all(is.na(all_rc))) all_rc <- NULL
  } else {
    all_rc <- NULL
  }

  as.eqnlist(out, compartments = all_compartments, compartmentOf = all_compartmentOf,
             reactionCompartment = all_rc)

}


#' @export
#' @param x obect of any class
#' @rdname eqnvec
is.eqnvec <- function(x) {
  if (inherits(x, "eqnvec") &&
      length(x) == length(names(x))
  )
    return(TRUE)
  
  else
    return(FALSE)
}


## Class "eqnvec" and its methods --------------------------------------------





#' Encode equation vector in format with sufficient spaces
#' 
#' @param x object of class [eqnvec]. Alternatively, a named parsable character vector.
#' @param ... additional arguments
#' @return named character
#' @export format.eqnvec
#' @export
format.eqnvec <- function(x, ...) {
  
  eqnvec <- x
  
  eqns <- sapply(eqnvec, function(eqn) {
    parser.out <- getParseData(parse(text = eqn, keep.source = TRUE))
    parser.out <- subset(parser.out, terminal == TRUE)
    # parser.out$text[parser.out$text == "*"] <- "*" (avoid non-ASCII characters for CRAN)
    out <- paste(parser.out$text, collapse = "")
    return(out)
  })
  
  patterns <- c("+", "-", "*", "/")
  for (p in patterns) eqns <- gsub(p, paste0(" ", p, " "), eqns, fixed = TRUE)
  
  return(eqns)
    
  
}

#' Print equation vector
#' 
#' @param x object of class [eqnvec].
#' @param width numeric, width of the print-out
#' @param pander logical, use pander for output (used with R markdown)
#' @param ... not used right now
#' 
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' 
#' @import stringr
#' @export
print.eqnvec <- function(x, width = 140, pander = FALSE, ...) {
  
  eqnvec <- x

  # Stuff to print
  m_odr <- "Idx"
  m_rel <- " <- "
  m_sep <- " "
  m_species <- names(eqnvec)
  
  # Width of stuff to print
  m_odrWidth <- max(3, nchar(m_odr))
  m_speciesWidth <- max(nchar(m_species), nchar("outer"))
  m_lineWidth <- max(width, m_speciesWidth + 10)
  m_relWidth <- nchar(m_rel)
  m_sepWidth <- nchar(m_sep)
  
  # Compound widths
  m_frontWidth <- m_odrWidth + m_speciesWidth + m_relWidth + m_sepWidth
  m_eqnWidth <- m_lineWidth - m_frontWidth
  
  # Order of states for alphabetical for print out
  m_eqnOrder <- order(m_species)
  
  # Iterate over species
  m_msgEqn <- do.call(c, mapply(function(eqn, spec, odr) {
    return(paste0(
      str_pad(string = odr, side = "left", width = m_odrWidth),
      m_sep,
      str_pad(string = spec, side = "left", width = m_speciesWidth),
      m_rel,
      str_wrap(string = gsub(x = eqn, pattern = " ", replacement = "", fixed = TRUE),
               width = m_eqnWidth, exdent = m_frontWidth)
    ))
  }, eqn = eqnvec[m_eqnOrder], spec = m_species[m_eqnOrder], odr = m_eqnOrder, SIMPLIFY = FALSE))
  
  # Print to command line or to pander
  if (!pander) {
    cat(paste0(str_pad(string = m_odr, side = "left", width = m_odrWidth),
               m_sep,
               str_pad(string = "Inner", side = "left", width = m_speciesWidth),
               m_rel,
               "Outer\n"))
    cat(m_msgEqn, sep = "\n")
  } else {
    pander::panderOptions("table.alignment.default", "left")
    pander::panderOptions("table.split.table", Inf)
    pander::panderOptions("table.split.cells", Inf)
    out <- as.data.frame(unclass(eqnvec), stringsAsFactors = FALSE)
    colnames(out) <- "" #  as.character(substitute(eqnvec))
    out[, 1] <- format.eqnvec(out[, 1])
    pander::pander(out)
    
  }
  

}



#' Summary of an equation vector
#' 
#' @param object of class [eqnvec].
#' @param ... additional arguments
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' 
#' @export
summary.eqnvec <- function(object, ...) {
  
  eqnvec <- object
  
  m_msg <- mapply(function(name, eqn) {
    m_symb <- paste0(getSymbols(eqn), sep = ", ", collapse = "")
    m_msg <- paste0(name, " = f( ", m_symb, ")")
    }, name = names(eqnvec), eqn = eqnvec)
  cat(m_msg, sep = "\n")
}


#' @export
c.eqnvec <- function(...) {
 
  out <- lapply(list(...), unclass)
  out <- do.call(c, out)
  if (any(duplicated(names(out)))) {
    stop("Names must be unique")
  }
  
  as.eqnvec(out)
}

#' @export
"[.eqnvec" <- function(x, ...) {
  out <- unclass(x)[...]
  class(out) <- c("eqnvec", "character")
  return(out)
}



#' Identify linear variables in an equation vector using sympy
#'
#' @param eqnvec An object of class `eqnvec`, representing a set of equations.
#' @details This function calls Python's `sympy` library via `reticulate` to symbolically analyze equations and determine if variables appear linearly in all equations.
#'
#' @return A character vector of variables that occur linearly in all equations.
#'
#' @examples
#' eqnvec <- as.eqnvec(
#'   c("-k1*A", "k1*A - k2*B", "-k3*B*C/(Km+C) + k4*pC", "k3*B*C/(Km+C) - k4*pC"),
#'   names = c("A", "B", "C", "pC")
#' )
#' getLinVars(eqnvec)
#'
#' @importFrom reticulate import py_run_string
#' @export
getLinVars <- function(eqnvec) {
  if (!inherits(eqnvec, "eqnvec")) {
    stop("Input 'eqnvec' must be of class 'eqnvec'.")
  }
  
  sympy <- reticulate::import("sympy")
  sympy_zero <- sympy$Integer(0)
  
  variables <- names(eqnvec)
  sympy_vars <- lapply(variables, sympy$Symbol)
  sympy_eqns <- lapply(as.character(eqnvec), sympy$simplify)
  
  is_linear_in_eq <- function(eqn, var) {
    first_derivative <- sympy$diff(eqn, var)
    second_derivative <- sympy$diff(first_derivative, var)
    is_second_derivative_zero <- sympy$simplify(second_derivative) == sympy_zero
    is_first_derivative_nonzero <- sympy$simplify(first_derivative) != sympy_zero
    is_second_derivative_zero && is_first_derivative_nonzero
  }
  
  linear_vars <- sapply(seq_along(sympy_vars), function(i) {
    var <- sympy_vars[[i]]
    all(sapply(sympy_eqns, function(eqn) is_linear_in_eq(eqn, var)))
  })
  
  variables[linear_vars]
}

#' Log-transform variables in an equation vector using SymPy
#'
#' @param eqnvec An object of class `eqnvec`, representing a set of equations.
#' @param whichVar A character vector specifying the variables to be log-transformed.
#' @details The function applies a logarithmic transformation to the specified variables in the equation vector.
#' For a variable `var`, the transformation is `log(var)` and the derivative `d/dt log(var)` is replaced by 
#' `(1/exp(log(var))) * d/dt var`, substituting `var` with `exp(log(var))` in the original equation.
#' The original variable equations are replaced or removed, and the transformed equations are added with updated names prefixed by `log_`.
#' The equations are simplified using SymPy to ensure mathematical correctness and compactness.
#'
#' @return An updated `eqnvec` object with transformed equations.
#'
#' @examples
#' eqnvec <- as.eqnvec(
#'   c("-k1*A + k2*B", "k1*A - k2*B"),
#'   names = c("A", "B")
#' )
#' log_transformed_eqns <- x2logx(eqnvec, c("A", "B"))
#'
#' @importFrom reticulate import
#' @export
x2logx <- function(eqnvec, whichVar) {
  if (!inherits(eqnvec, "eqnvec")) {
    stop("Input 'eqnvec' must be of class 'eqnvec'.")
  }
  
  if (!all(whichVar %in% names(eqnvec))) {
    stop("All variables in whichVar must be present in eqnvec names.")
  }

  # Import SymPy
  sympy <- reticulate::import("sympy")
  
  # Create symbolic variables for all variables in eqnvec
  variables <- names(eqnvec)
  sympy_vars <- lapply(variables, sympy$Symbol)
  names(sympy_vars) <- variables
  
  # Create log-transformed symbolic variables
  log_vars <- paste0("log_", whichVar)
  sympy_log_vars <- lapply(log_vars, sympy$Symbol)
  names(sympy_log_vars) <- log_vars
  
  # Convert equations to SymPy expressions
  sympy_eqns <- lapply(as.character(eqnvec), sympy$simplify)
  names(sympy_eqns) <- variables
  
  # Create substitution list for exp(log(x)) -> x
  subs_list <- lapply(whichVar, function(var) {
    list(
      sympy_vars[[var]], 
      sympy$exp(sympy_log_vars[[paste0("log_", var)]])
    )
  })
  
  # Transform equations
  new_eqns <- list()
  
  # Process non-transformed variables
  for (var in setdiff(variables, whichVar)) {
    expr <- sympy_eqns[[var]]
    # Substitute exp(log(x)) for each transformed variable
    for (sub in subs_list) {
      expr <- sympy$simplify(expr$subs(sub[[1]], sub[[2]]))
    }
    new_eqns[[var]] <- as.character(expr)
  }
  
  # Process transformed variables
  for (var in whichVar) {
    log_var <- paste0("log_", var)
    expr <- sympy_eqns[[var]]
    
    # Apply chain rule: d/dt log(x) = (1/x) * dx/dt
    # First substitute all other transformations
    for (sub in subs_list) {
      expr <- sympy$simplify(expr$subs(sub[[1]], sub[[2]]))
    }
    
    # Then apply the chain rule
    expr <- sympy$simplify(expr / sympy$exp(sympy_log_vars[[log_var]]))
    new_eqns[[log_var]] <- as.character(expr)
  }
  
  # Create new eqnvec with transformed equations
  result_names <- c(setdiff(variables, whichVar), paste0("log_", whichVar))
  result_eqns <- unlist(new_eqns[result_names])
  
  as.eqnvec(result_eqns, result_names)
}



