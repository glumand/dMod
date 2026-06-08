#' Compare two objects and return differences
#' 
#' Works eigher on a list or on two arguments. In case of a list,
#' comparison is done with respect to a reference entry. Besides the
#' objects themselves also some of their attributes are compared,
#' i.e. "equations", "parameters" and "events" and "forcings".
#' 
#' @param vec1 object of class [eqnvec], `character` or
#' `data.frame`. Alternatively, a list of such objects.
#' @param vec2 same as vec1. Not used if vec1 is a list.
#' @param reference numeric of length one, the reference entry.
#' @param ... arguments going to the corresponding methods
#' @return `data.frame` or list of data.frames with the differences. 
#' 
#' @export
#' @examples
#' ## Compare equation vectors
#' eq1 <- eqnvec(a = "-k1*a + k2*b", b = "k2*a - k2*b")
#' eq2 <- eqnvec(a = "-k1*a", b = "k2*a - k2*b", c = "k2*b")
#' compare(eq1, eq2)
#' 
#' ## Compare character vectors
#' c1 <- c("a", "b")
#' c2 <- c("b", "c")
#' compare(c1, c2)
#' 
#' ## Compare data.frames
#' d1 <- data.frame(var = "a", time = 1, value = 1:3, method = "replace")
#' d2 <- data.frame(var = "a", time = 1, value = 2:4, method = "replace")
#' compare(d1, d2)
#' 
#' ## Compare structures like prediction functions
#' fn1 <- function(x) x^2
#' attr(fn1, "equations") <- eq1
#' attr(fn1, "parameters") <- c1
#' attr(fn1, "events") <- d1
#' 
#' fn2 <- function(x) x^3
#' attr(fn2, "equations") <- eq2
#' attr(fn2, "parameters") <- c2
#' attr(fn2, "events") <- d2
#' 
#' mylist <- list(f1 = fn1, f2 = fn2)
#' compare(mylist)

compare <- function(vec1, ...) {
  UseMethod("compare", vec1)
}

#' @export
#' @rdname compare
compare.list <- function(vec1, vec2 = NULL, reference = 1, ...) {
  
  index <- (1:length(vec1))[-reference]
  diffable.attributes <- c("equations", "parameters", "forcings", "events")
  
  
  out.total <- lapply(index, function(i) {
    
    # Compare objects if possible
    vec1.inner <- vec1[[reference]]
    vec2.inner <- vec1[[i]]
    out1 <- NULL
    if(any(class(vec1.inner) %in% c("eqnvec", "data.frame"))) {
      out1 <- list(compare(vec1.inner, vec2.inner))
      names(out1) <- "object"
    }
      
    # Compare comparable attributes of the object if available
    out2 <- NULL
    attributes1 <- attributes(vec1.inner)[diffable.attributes]
    attributes2 <- attributes(vec2.inner)[diffable.attributes]
    slots <- names(attributes1)[!is.na(names(attributes1))]
    out2 <- lapply(slots, function(n) {
      compare(attributes1[[n]], attributes2[[n]])
    })
    names(out2) <- slots
    
    c(out1, out2)
    
  })
  names(out.total) <- names(vec1)[index]
  
  ## Do resorting of the list
  innernames <- names(out.total[[1]])
  out.total <- lapply(innernames, function(n) {
    out <- lapply(out.total, function(out) out[[n]])
    out[!sapply(out, is.null)]
  })
  names(out.total) <- innernames
  
  
  return(out.total)
  
  
  
}

#' @export
#' @rdname compare
compare.character <- function(vec1, vec2 = NULL, ...) {
  missing <- setdiff(vec1, vec2)
  additional <- setdiff(vec2, vec1)
  
  out <- do.call(rbind, 
          list(different = NULL, 
               missing = data.frame(name = missing), 
               additional = data.frame(name = additional)
          )
  )
  
  if(nrow(out) == 0) out <- NULL
  return(out)
  
  
  
}

#' @export
#' @rdname compare
compare.eqnvec <- function(vec1, vec2 = NULL, ...) {

  names1 <- names(vec1)
  names2 <- names(vec2)
  
  missing <- setdiff(names1, names2)
  additional <- setdiff(names2, names1)
  joint <- intersect(names1, names2)
  
  # Compare joint equations
  v1 <- format(vec1)
  v2 <- format(vec2)
  not.coincide <- which(as.character(v1[joint]) != as.character(v2[joint]))
  
  different <- data.frame(name = names(v2[not.coincide]), equation = as.character(v2[not.coincide]))
  missing <- data.frame(name = names(v2[missing]), equation = as.character(v2[missing]))
  additional <- data.frame(name = names(v2[additional]), equation = as.character(v2[additional]))
  
  out <- do.call(rbind, list(different = different, missing = missing, additional = additional))
  if(nrow(out) == 0) out <- NULL
  return(out)
  
  
}

#' @export
#' @rdname compare
compare.data.frame <- function(vec1, vec2 = NULL, ...) {
  
  additional <- !duplicated(rbind(vec1, vec2))[-(1:nrow(vec1))]
  missing <- !duplicated(rbind(vec2, vec1))[-(1:nrow(vec2))]
  
  out <- do.call(rbind, list(different = character(0), missing = vec1[missing, ], additional = vec2[additional, ]))
  if(nrow(out) == 0) out <- NULL
  return(out)
  
  
}


#' Combine several data.frames by rowbind
#' 
#' @param ... data.frames or matrices with not necessarily overlapping colnames
#' @details This function is useful when separating models into independent csv model files,
#' e.g.~a receptor model and several downstream pathways. Then, the models can be recombined 
#' into one model by `combine()`.
#' 
#' @return A `data.frame`
#' @export
#' @examples
#' data1 <- data.frame(Description = "reaction 1", Rate = "k1*A", A = -1, B = 1)
#' data2 <- data.frame(Description = "reaction 2", Rate = "k2*B", B = -1, C = 1)
#' combine(data1, data2)
#' @export
combine <- function(...) {
  
  # List of input data.frames
  mylist <- list(...)
  # Remove empty slots
  is.empty <- sapply(mylist, is.null)
  mylist <- mylist[!is.empty]
  
  mynames <- unique(unlist(lapply(mylist, function(S) colnames(S))))
  
  mylist <- lapply(mylist, function(l) {
    
    if(is.data.frame(l)) {
      i <- sapply(l, is.factor)
      l[i] <- lapply(l[i], as.character)
      present.list <- as.list(l)
      missing.names <- setdiff(mynames, names(present.list))
      missing.list <- structure(as.list(rep(NA, length(missing.names))), names = missing.names)
      combined.data <- do.call(function(...) cbind.data.frame(..., stringsAsFactors = FALSE), c(present.list, missing.list))
      rownames(combined.data) <- rownames(l)
    }
    if(is.matrix(l)) {
      present.matrix <- as.matrix(l)
      missing.names <- setdiff(mynames, colnames(present.matrix))
      missing.matrix <- matrix(0, nrow = nrow(present.matrix), ncol = length(missing.names), 
                             dimnames = list(NULL, missing.names))
      combined.data <- submatrix(cbind(present.matrix, missing.matrix), cols = mynames)
      rownames(combined.data) <- rownames(l)
    }
    
    return(combined.data)
  })
  
  out <- do.call(rbind, mylist)
  
  return(out)
  
  
}

#' Submatrix of a matrix returning ALWAYS a matrix
#' 
#' @param M matrix
#' @param rows Index vector
#' @param cols Index vector
#' @return The matrix `M[rows, cols]`, keeping/adjusting attributes like ncol nrow and dimnames.
#' @export
submatrix <- function(M, rows = 1:nrow(M), cols = 1:ncol(M)) {
  
 M[rows, cols, drop = FALSE] 
  
  # myrows <- (structure(1:nrow(M), names = rownames(M)))[rows]
  # mycols <- (structure(1:ncol(M), names = colnames(M)))[cols]
  # 
  # if(any(is.na(myrows)) | any(is.na(mycols))) stop("subscript out of bounds")
  # 
  # matrix(M[myrows, mycols], 
  #        nrow = length(myrows), ncol = length(mycols), 
  #        dimnames = list(rownames(M)[myrows], colnames(M)[mycols]))

}


#' Embed two matrices into one blockdiagonal matrix
#' 
#' @param M matrix of type character
#' @param N matrix of type character
#' @return Matrix of type character containing M and N as upper left and lower right block
#' @examples
#' M <- matrix(1:9, 3, 3, dimnames = list(letters[1:3], letters[1:3]))
#' N <- matrix(1:4, 2, 2, dimnames = list(LETTERS[1:2], LETTERS[1:2]))
#' blockdiagSymb(M, N)
#' @export
blockdiagSymb <- function(M, N) {
  
  red <- sapply(list(M, N), is.null)
  if(all(red)) {
    return()
  } else if(red[1]) {
    return(N)
  } else if(red[2]) {
    return(M)
  }
  
  A <- matrix(0, ncol=dim(N)[2], nrow=dim(M)[1])
  B <- matrix(0, ncol=dim(M)[2], nrow=dim(N)[1])
  result <- rbind(cbind(M, A), cbind(B, N))
  colnames(result) <- c(colnames(M), colnames(N))
  rownames(result) <- c(rownames(M), rownames(N))
  
  return(result)
  
}



#' Translate wide output format (e.g., from ODE solver) into long format
#'
#' Converts simulation output in wide format into a tidy long format suitable for
#' plotting or further analysis (e.g., with \pkg{ggplot2}). The function assumes
#' that the first column of \code{out} represents a time-like variable and the
#' remaining columns contain values.
#'
#' @param out A \code{data.frame}, \code{matrix}, or a \code{list} of matrices in wide format.
#' @param keep Integer vector specifying the column indices to keep (default is \code{1}).
#' @param na.rm Logical. If \code{TRUE}, missing values are removed in the long-format output.
#'
#' @details
#' If \code{out} is a list, the list names are added as an additional column named
#' \code{"condition"}. This is particularly useful for plotting results from multiple
#' simulation conditions with \pkg{ggplot2}.
#'
#' @return A \code{data.frame} in long format with the following columns:
#' \itemize{
#'   \item \code{"time"} — values from \code{out[, 1]}.
#'   \item \code{"name"} — column names from \code{out[, -1]}.
#'   \item \code{"value"} — corresponding numeric values.
#'   \item \code{"condition"} — if \code{out} was a list, contains the list names.
#' }
#'
#' @export
wide2long <- function(out, keep = 1, na.rm = FALSE) {
  
  UseMethod("wide2long", out)
  
  
}

#' @rdname wide2long
#' @export
wide2long.data.frame <- function(out, keep = 1, na.rm = FALSE) {
  
  wide2long.matrix(out, keep = keep, na.rm = na.rm)
  
}

#' @rdname wide2long
#' @export
wide2long.matrix <- function(out, keep = 1, na.rm = FALSE) {
  
  timenames <- colnames(out)[keep]
  allnames <- colnames(out)[-keep]
  if (any(duplicated(allnames))) warning("Found duplicated colnames in out. Duplicates were removed.")
  times <- out[,keep]
  ntimes <- nrow(out)
  values <- unlist(out[,allnames])
  outlong <- data.frame(times, 
                        name = factor(rep(allnames, each = ntimes), levels = allnames), 
                        value = as.numeric(values))
  colnames(outlong)[1:length(keep)] <- timenames
  
  if (na.rm) outlong <- outlong[!is.na(outlong$value),]
  
  return(outlong)
  
}

#' @rdname wide2long
#' @export
wide2long.list <- function(out, keep = 1, na.rm = FALSE) {
  
  conditions <- names(out)
  
  outlong <- do.call(rbind, lapply(1:max(c(length(conditions), 1)), function(cond) {
    
    cbind(wide2long.matrix(out[[cond]]), condition = conditions[cond])
    
  }))
  
  
  
  return(outlong)
  
}


#' Translate long to wide format (inverse of wide2long.matrix) 
#' 
#' @param out data.frame in long format 
#' @return data.frame in wide format 
#' @export
long2wide <- function(out) {
  
  timename <- colnames(out)[1]
  times <- unique(out[,1])
  allnames <- unique(as.character(out[,2]))
  M <- matrix(out[,3], nrow=length(times), ncol=length(allnames))
  M <- cbind(times, M)
  colnames(M) <- c(timename, allnames)
  
  return(M)
  
}


#' Bind named list of data.frames into one data.frame
#' 
#' @param mylist A named list of data.frame. The data.frames are expected to have the same structure.
#' @details Each data.frame ist augented by a "condition" column containing the name attributed of
#' the list entry. Subsequently, the augmented data.frames are bound together by `rbind`.
#' @return data.frame with the originial columns augmented by a "condition" column.
#' @export
lbind <- function(mylist) {
  
  conditions <- names(mylist)
  #numconditions <- suppressWarnings(as.numeric(conditions))
  #
  # if(!any(is.na(numconditions))) 
  #   numconditions <- as.numeric(numconditions) 
  # else 
  numconditions <- conditions

  
  outlong <- do.call(rbind, lapply(1:length(conditions), function(cond) {
    
    myout <- mylist[[cond]]
    if (nrow(myout) > 0)
      myout[["condition"]] <- numconditions[cond]
    else
      myout[["condition"]] <- character(0)
    
    return(myout)
    
  }))
  
  return(outlong)
  
}

#' Alternative version of expand.grid
#' @param seq1 Vector, numeric or character
#' @param seq2 Vector, numeric or character
#' @return Matrix ob combinations of elemens of `seq1` and `seq2`
expand.grid.alt <- function(seq1, seq2) {
  cbind(Var1=rep.int(seq1, length(seq2)), Var2=rep(seq2, each=length(seq1)))
}


## Windows-only: temp Makevars (existing user Makevars + `lines`) for R_MAKEVARS_USER.
.compileMakevarsUser <- function(lines) {
  f <- Sys.getenv("R_MAKEVARS_USER", unset = NA)
  if (is.na(f) || !file.exists(f)) {
    cand <- path.expand(c("~/.R/Makevars.ucrt", "~/.R/Makevars.win64",
                          "~/.R/Makevars.win", "~/.R/Makevars"))
    cand <- cand[file.exists(cand)]
    f <- if (length(cand)) cand[1] else NA
  }
  prev <- if (!is.na(f) && file.exists(f)) readLines(f, warn = FALSE) else character()
  mv <- tempfile(fileext = ".mk")
  writeLines(c(prev, lines), mv)
  mv
}


#' Compile model-related C/C++ code
#'
#' @description
#' Compiles model objects ([parfn], [obsfn], [prdfn]) related C/C++ files into shared libraries via `R CMD SHLIB`.
#'
#' @details
#' Per-file compile and link flags are taken from the `"compileInfo"`
#' attribute that [odemodel()], [Xs()], [Xf()], [Y()] and [P()] attach to
#' their return values. Each entry carries the source file together with the
#' `compileArgs` and `linkArgs` reported by the backend that produced it
#' (`cOde::funC`, `CppODE::CppODE`, `CppODE::CVODE`, ...), so solver-specific
#' libraries reach only the files that need them. Objects without
#' `compileInfo` fall back to modelname-based file discovery in the current
#' working directory.
#'
#' @param ... One or more model objects.
#' @param output Optional name for a combined shared library. When set, all
#'   files are linked into one object and the union of their `linkArgs` is
#'   applied.
#' @param args Additional compiler/linker flags applied to every file.
#' @param cores Parallel compilation jobs (Unix only, requires `cores > 1`).
#' @param verbose If `TRUE`, print compiler commands.
#'
#' @return Invisibly `TRUE` on success.
#' @export
compile <- function(..., output = NULL, args = NULL, cores = 1, verbose = FALSE) {

  ## save & restore env
  old <- Sys.getenv(c("PKG_CFLAGS","PKG_CXXFLAGS","PKG_CPPFLAGS","PKG_LIBS"), unset = NA)
  on.exit({
    for (n in names(old))
      if (is.na(old[n])) Sys.unsetenv(n) else Sys.setenv(structure(old[n], names = n))
  }, add = TRUE)

  objs <- list(...)
  if (!length(objs)) stop("No objects")
  obj.names <- as.character(substitute(list(...)))[-1]
  Rbin  <- shQuote(file.path(R.home("bin"), "R"))
  so    <- .Platform$dynlib.ext
  cfg   <- function(x) trimws(system(paste(Rbin, "CMD config", x), intern = TRUE))
  strip <- function(x) trimws(gsub("(^| )-std=[^ ]+", "", x))

  ## classify objects
  is_dmod <- vapply(objs, inherits, logical(1), c("obsfn","parfn","prdfn"))
  is_cpp  <- vapply(objs, function(o) !is.null(attr(o, "srcfile")), logical(1))

  ## Collect per-file build info.
  ## Primary source is `attr(o, "compileInfo")` carrying
  ## (srcfile, compileArgs, linkArgs) as reported by cOde/CppODE/CVODE.
  ## Falls back to modelname-based file discovery for objects that lack
  ## compileInfo, and to the bare `srcfile` attribute for raw CppODE objects.
  info_from_compileInfo <- unlist(
    lapply(objs, function(o) attr(o, "compileInfo")),
    recursive = FALSE
  )

  info_fallback <- list()
  for (i in seq_along(objs)) {
    o <- objs[[i]]
    if (!is.null(attr(o, "compileInfo"))) next
    if (is_dmod[i]) {
      b <- outer(modelname(o), c("","_deriv","_s","_s2","_sdcv","_dfdx","_dfdp"), paste0)
      cand <- c(paste0(b, ".c"), paste0(b, ".cpp"))
      src <- cand[file.exists(cand)]
      for (s in src)
        info_fallback[[length(info_fallback) + 1]] <-
          list(srcfile = normalizePath(s, winslash = "/", mustWork = TRUE),
               compileArgs = "", linkArgs = "")
    } else if (is_cpp[i]) {
      s <- attr(o, "srcfile")
      if (length(s) && nzchar(s) && file.exists(s))
        info_fallback[[length(info_fallback) + 1]] <- list(
          srcfile     = normalizePath(s, winslash = "/", mustWork = TRUE),
          compileArgs = attr(o, "compileArgs") %||% "",
          linkArgs    = attr(o, "linkArgs")    %||% "")
    }
  }

  info <- c(info_from_compileInfo, info_fallback)

  ## Expand entries with multiple srcfiles (e.g. cOde spills _deriv.c
  ## alongside the main .c) into one entry per file.
  info <- unlist(lapply(info, function(e) {
    if (!length(e$srcfile)) return(list())
    if (length(e$srcfile) == 1L) return(list(e))
    lapply(e$srcfile, function(s) list(srcfile = s, compileArgs = e$compileArgs, linkArgs = e$linkArgs))
  }), recursive = FALSE)

  info <- Filter(function(e) length(e$srcfile) == 1L && nzchar(e$srcfile) && file.exists(e$srcfile), info)
  if (!length(info)) stop("No source files found")

  ## Deduplicate by srcfile, keeping the first (non-empty) flags we saw.
  ord <- order(vapply(info, function(e) e$srcfile, character(1)))
  info <- info[ord]
  keep <- !duplicated(vapply(info, function(e) e$srcfile, character(1)))
  info <- info[keep]

  files      <- vapply(info, function(e) e$srcfile, character(1))
  roots      <- sub("\\.[^.]+$", "", basename(files))
  roots_full <- sub("\\.[^.]+$", "", files)

  ## compiler flags
  if (.Platform$OS.type == "windows") cores <- 1
  pic  <- if (.Platform$OS.type == "windows") "" else "-fPIC"
  base <- paste("-O2 -DNDEBUG -w", pic)

  ## KLU detection
  uses_klu <- any(vapply(objs, function(o) isTRUE(attr(o, "sparse")), logical(1)))
  klu_flag <- ""; klu_lib <- ""
  if (uses_klu) {
    sd <- if (nzchar(.Platform$r_arch)) file.path("lib", .Platform$r_arch) else "lib"
    lp <- system.file(sd, "libcppode_ss.a", package = "CppODE")
    if (file.exists(lp)) { klu_flag <- "-DKLU"; klu_lib <- shQuote(lp) }
  }

  ## shared pieces (compiler/linker) that apply to every file
  cxx_base <- paste(base, klu_flag)
  extra_args <- paste(c(args), collapse = " ")
  if (nzchar(extra_args)) {
    base     <- paste(base,     extra_args)
    cxx_base <- paste(cxx_base, extra_args)
  }
  ## BLAS/LAPACK: on Windows `R CMD config BLAS_LIBS` returns a value with
  ## unexpanded `$(R_HOME)`/`$(R_ARCH)` references. Those go into PKG_LIBS as
  ## an env var, and make should re-expand them, but in practice the
  ## expansion is unreliable inside SHLIB-generated link commands -- the
  ## final g++ invocation comes out without any BLAS libs. We sidestep that
  ## by building an absolute -L path here and skipping `R CMD config`.
  if (.Platform$OS.type == "windows") {
    r_bin   <- file.path(R.home("bin"), .Platform$r_arch)
    blaslapack <- paste0("-L", shQuote(r_bin), " -lRlapack -lRblas")
  } else {
    blaslapack <- paste(cfg("LAPACK_LIBS"), cfg("BLAS_LIBS"))
  }
  base_libs <- paste(klu_lib, blaslapack)
  cppflags  <- paste0("-I", system.file("include", package = "CppODE"))

  ## Compiler invocation bits cached up front so parallel forks don't each
  ## re-spawn R-CMD-config. Used by compile_one_obj() for the direct
  ## $CC/$CXX -c path.
  cc_bin      <- cfg("CC")
  cxx_bin     <- cfg("CXX")
  cflags_R    <- cfg("CFLAGS")
  cxxflags_R  <- cfg("CXXFLAGS")
  cpicflags   <- cfg("CPICFLAGS")
  cxxpicflags <- cfg("CXXPICFLAGS")
  r_inc       <- paste0("-I", shQuote(R.home("include")))

  ## toolchain report (use a representative entry for display)
  if (any(grepl("\\.c$",   files))) cat(sprintf("using C compiler:   %s [%s]\n", strip(cc_bin),  trimws(base)))
  if (any(grepl("\\.cpp$", files))) cat(sprintf("using C++ compiler: %s [%s]\n", strip(cxx_bin), trimws(cxx_base)))

  ## unload stale DLLs
  loaded <- getLoadedDLLs()
  for (i in seq_along(roots))
    if (roots[i] %in% names(loaded)) try(dyn.unload(loaded[[roots[i]]][["path"]]), silent = TRUE)
  if (!is.null(output)) try(dyn.unload(paste0(output, so)), silent = TRUE)

  ## Compile one file with its own compile/link flags applied via PKG_*.
  ## Each invocation sets the env just before shelling out to R CMD SHLIB,
  ## so per-file linkArgs (e.g. Sundials libs for CVODE) reach only the
  ## files that need them. Works inside mclapply because each fork has its
  ## own env.
  compile_one <- function(entry) {
    extra_c <- entry$compileArgs %||% ""
    pkg_c  <- trimws(paste(base,     extra_c))
    pkg_cx <- trimws(paste(cxx_base, extra_c))
    pkg_l  <- trimws(paste(base_libs, entry$linkArgs %||% ""))
    Sys.setenv(
      PKG_CFLAGS   = pkg_c,
      PKG_CXXFLAGS = pkg_cx,
      PKG_CPPFLAGS = cppflags,
      PKG_LIBS     = pkg_l
    )
    if (.Platform$OS.type == "windows") {
      mv <- .compileMakevarsUser(c(
        paste("PKG_CFLAGS =",   pkg_c),
        paste("PKG_CXXFLAGS =", pkg_cx),
        paste("PKG_CPPFLAGS =", cppflags),
        paste("PKG_LIBS =",     pkg_l)
      ))
      old_mu <- Sys.getenv("R_MAKEVARS_USER", unset = NA)
      Sys.setenv(R_MAKEVARS_USER = mv)
      on.exit({
        if (is.na(old_mu)) Sys.unsetenv("R_MAKEVARS_USER")
        else Sys.setenv(R_MAKEVARS_USER = old_mu)
        unlink(mv)
      }, add = TRUE)
    }
    cmd <- paste(Rbin, "CMD SHLIB", shQuote(entry$srcfile))
    if (verbose) cat(cmd, "\n")
    if (system(cmd, ignore.stdout = !verbose, ignore.stderr = !verbose) != 0)
      stop("Compilation failed: ", entry$srcfile)
  }

  ## Compile a single source to a .o object only, via a direct $CC/$CXX -c
  ## invocation. Used by the combined-output path so the per-file compile
  ## phase can run in parallel; the subsequent R CMD SHLIB link then sees
  ## the .o files are up-to-date and skips recompilation.
  compile_one_obj <- function(entry) {
    src     <- entry$srcfile
    extra_c <- entry$compileArgs %||% ""
    is_cpp  <- grepl("\\.cpp$", src, ignore.case = TRUE)
    obj     <- sub("\\.[^.]+$", ".o", src)

    if (is_cpp) {
      cmd <- paste(
        cxx_bin, r_inc, cppflags,
        trimws(paste(cxx_base, extra_c)),
        cxxpicflags, cxxflags_R,
        "-c", shQuote(src),
        "-o", shQuote(obj)
      )
    } else {
      cmd <- paste(
        cc_bin, r_inc, cppflags,
        trimws(paste(base, extra_c)),
        cpicflags, cflags_R,
        "-c", shQuote(src),
        "-o", shQuote(obj)
      )
    }

    if (verbose) cat(cmd, "\n")
    if (system(cmd, ignore.stdout = !verbose, ignore.stderr = !verbose) != 0)
      stop("Compilation failed: ", src)
    obj
  }

  if (is.null(output)) {
    if (.Platform$OS.type == "unix" && cores > 1)
      parallel::mclapply(info, compile_one, mc.cores = cores)
    else for (e in info) compile_one(e)
    for (r in roots_full) dyn.load(paste0(r, so))
  } else {
    ## Combined output: per-file compile to .o (parallel on Unix when cores>1,
    ## serial otherwise — including on Windows), then a single R CMD SHLIB
    ## link over the original sources. Because every .o is freshly written
    ## above, make sees them as up-to-date and only runs the link recipe;
    ## passing the source list lets SHLIB pick the C++ linker when any
    ## source is .cpp, which a .o-only invocation would miss. The pre-compile
    ## also has to run on Windows: the single-call SHLIB (compile + link in
    ## one go) was occasionally producing .dll files that LoadLibrary
    ## couldn't resolve when the source pulled in BLAS via the symbolic-
    ## mode chain wrapper -- splitting compile and link sidesteps that.
    if (.Platform$OS.type == "unix" && cores > 1)
      parallel::mclapply(info, compile_one_obj, mc.cores = cores)
    else
      for (e in info) compile_one_obj(e)

    ## Link step: union of every entry's linkArgs (dedup) so Sundials-dependent
    ## files still pull their libs.
    all_link <- unique(unlist(lapply(info, function(e) strsplit(trimws(e$linkArgs %||% ""), "\\s+")[[1]])))
    all_link <- all_link[nzchar(all_link)]
    all_compile <- unique(unlist(lapply(info, function(e) strsplit(trimws(e$compileArgs %||% ""), "\\s+")[[1]])))
    all_compile <- all_compile[nzchar(all_compile)]

    pkg_cflags   <- trimws(paste(base,     paste(all_compile, collapse = " ")))
    pkg_cxxflags <- trimws(paste(cxx_base, paste(all_compile, collapse = " ")))
    pkg_libs     <- trimws(paste(base_libs, paste(all_link, collapse = " ")))
    Sys.setenv(
      PKG_CFLAGS   = pkg_cflags,
      PKG_CXXFLAGS = pkg_cxxflags,
      PKG_CPPFLAGS = cppflags,
      PKG_LIBS     = pkg_libs
    )

    ## Belt-and-suspenders for Windows: env-imported PKG_LIBS has been
    ## observed to vanish from SHLIB's generated link command on some R/rtools
    ## combinations, leaving the .dll unlinked against BLAS/LAPACK. Drop a
    ## per-link Makevars(.win) alongside the source files so make picks it up
    ## even if the environment doesn't make it through. We clean it up after
    ## the link so the directory state stays hermetic.
    mv_dir  <- dirname(files[1])
    mv_name <- if (.Platform$OS.type == "windows") "Makevars.win" else "Makevars"
    mv_path <- file.path(mv_dir, mv_name)
    mv_pre  <- if (file.exists(mv_path)) readLines(mv_path, warn = FALSE) else NULL
    writeLines(c(
      paste("PKG_CFLAGS =",   pkg_cflags),
      paste("PKG_CXXFLAGS =", pkg_cxxflags),
      paste("PKG_CPPFLAGS =", cppflags),
      paste("PKG_LIBS =",     pkg_libs)
    ), mv_path)
    on.exit({
      if (is.null(mv_pre)) try(unlink(mv_path), silent = TRUE)
      else                 try(writeLines(mv_pre, mv_path), silent = TRUE)
    }, add = TRUE)

    ## Windows fallback for BLAS/LAPACK: inject PKG_* via R_MAKEVARS_USER.
    if (.Platform$OS.type == "windows") {
      mv <- .compileMakevarsUser(c(
        paste("PKG_CFLAGS =",   pkg_cflags),
        paste("PKG_CXXFLAGS =", pkg_cxxflags),
        paste("PKG_CPPFLAGS =", cppflags),
        paste("PKG_LIBS =",     pkg_libs)
      ))
      old_mu <- Sys.getenv("R_MAKEVARS_USER", unset = NA)
      Sys.setenv(R_MAKEVARS_USER = mv)
      on.exit({
        if (is.na(old_mu)) Sys.unsetenv("R_MAKEVARS_USER")
        else Sys.setenv(R_MAKEVARS_USER = old_mu)
        unlink(mv)
      }, add = TRUE)
    }

    output <- sub(paste0("\\", so, "$"), "", output)
    out <- file.path(dirname(files[1]), paste0(output, so))
    try(dyn.unload(out), silent = TRUE)
    if (file.exists(out)) unlink(out)
    cmd <- paste(Rbin, "CMD SHLIB", paste(shQuote(files), collapse = " "), "-o", shQuote(out))
    if (verbose) cat(cmd, "\n")
    ## Capture SHLIB output and strip its compiler banner: on this path the
    ## .o files are already fresh from compile_one_obj, so make only runs
    ## the link recipe and the "using C/C++ compiler:" line would be
    ## misleading. Stderr is folded into stdout via shell redirection so
    ## error messages still surface when verbose = TRUE.
    out_lines <- suppressWarnings(
      system(paste(cmd, "2>&1"), intern = TRUE)
    )
    status <- attr(out_lines, "status")
    if (verbose) {
      out_lines <- out_lines[!grepl("^using (C|C\\+\\+) compiler:", out_lines)]
      writeLines(out_lines)
    }
    if (!is.null(status) && status != 0L)
      stop("Compilation failed:\n", paste(out_lines, collapse = "\n"))
    if (!file.exists(out))
      stop("R CMD SHLIB returned exit 0 but did not produce ", out, ":\n",
           paste(out_lines, collapse = "\n"))
    dyn.load(out)
    for (i in which(is_dmod))
      eval.parent(parse(text = paste0("modelname(", obj.names[i], ") <- '", output, "'")))
  }

  invisible(TRUE)
}


#' Determine loaded DLLs available in working directory
#' 
#' @return Character vector with the names of the loaded DLLs available in the working directory
#' @export
getLocalDLLs <- function() {
  
  all.dlls <- getLoadedDLLs()
  is.local <- sapply(all.dlls, function(x) grepl(getwd(), unclass(x)$path, fixed = TRUE))
  names(is.local)[is.local]
  
}


#' Load shared object for a dMod object
#' 
#' Usually when restarting the R session, although all objects are saved in
#' the workspace, the dynamic libraries are not linked any more. `loadDLL`
#' is a wrapper for `dyn.load` that uses the "modelname" attribute of
#' dMod objects like prediction functions, observation functions, etc. to
#' load the corresponding shared object.
#' 
#' @param ... objects of class prdfn, obsfn, parfn, objfn, ...
#' 
#' @export
loadDLL <- function(...) {
  
  .so <- .Platform$dynlib.ext
  models <- modelname(...)
  files <- paste0(outer(models, c("", "_s", "_s2", "_sdcv", "_deriv", "_dfdx", "_dfdp"), paste0), .so)
  files <- files[file.exists(files)]
  
  for (f in files) {
    try(dyn.unload(f), silent = TRUE)
    dyn.load(f)
  }
  message("The following local files were dynamically loaded: ", paste(files, collapse = ", "))
}


sanitizeCores <- function(cores)  {
  
  max.cores <- parallel::detectCores()
  min(max.cores, cores)
 #  
 # if (Sys.info()[['sysname']] == "Windows") cores <- 1
 # return(cores)
  
}

sanitizeConditions <- function(conditions) {
  
  new <- str_replace_all(conditions, "[[:punct:]]", "_")
  new <- str_replace_all(new, "\\s+", "_")
  return(new)
  
}

sanitizePars <- function(pars = NULL, fixed = NULL) {
  
  # Convert fixed to named numeric
  if (!is.null(fixed)) fixed <- structure(as.numeric(fixed), names = names(fixed))
  
  # Convert pars to named numeric
  if (!is.null(pars)) {
    pars <- structure(as.numeric(pars), names = names(pars))
    # remove fixed from pars
    pars <- pars[setdiff(names(pars), names(fixed))]
  }
    
  
  return(list(pars = pars, fixed = fixed))
  
}


sanitizeData <- function(x, required = c("name", "time", "value"), imputed = c(sigma = NA, lloq = -Inf)) {
  
  all.names <- names(x)
  
  missing.required <- setdiff(required, all.names)
  missing.imputed <- setdiff(names(imputed), all.names)
  
  if (length(missing.required) > 0)
      stop("These mandatory columns are missing: ", paste(missing.required, collapse = ", "))
  
  if (length(missing.imputed) > 0) {
    
    for (n in missing.imputed) x[[n]] <- imputed[n]
    
  }
  
  list(data = x, columns = c(required, names(imputed)))
  
}


#' Print list of dMod objects in .GlobalEnv
#' 
#' @description Lists the objects for a set of classes.
#'   
#' @param classlist List of object classes to print.
#' @param envir Alternative environment to search for objects.
#' @examples 
#' \dontrun{
#' lsdMod()
#' lsdMod(classlist = "prdfn", envir = environment(obj)) 
#' }
#' 
#' @export
lsdMod <- function(classlist = c("odemodel", "parfn", "prdfn", "obsfn", "objfn", "datalist"), envir = .GlobalEnv){
  glist <- as.list(envir)
  out <- list()
  for (a in classlist) {
    flist <- which(sapply(glist, function(f) any(class(f) == a)))
    out[[a]] <- names(glist[flist])
    #cat(a,": ")
    #cat(paste(out[[a]], collapse = ", "),"\n")
  }
  
  unlist(out)
  
  
  
}





#' Select attributes.
#' 
#' @description Select or discard attributes from an object.
#'   
#' @param x The object to work on
#' @param atr An optional list of attributes which are either kept or removed. 
#'   This parameter defaults to dim, dimnames, names,  col.names, and row.names.
#' @param keep For keep = TRUE, atr is a positive list on attributes which are 
#'   kept, for keep = FALSE, \option{atr} are removed.
#'   
#' @return x with selected attributes.
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' @author Mirjam Fehling-Kaschek, \email{mirjam.fehling@@physik.uni-freiburg.de}
#'   
#' @export
attrs <- function(x, atr = NULL, keep = TRUE) {

  if (is.null(atr)) {
    atr <- c("class", "dim", "dimnames", "names", "col.names", "row.names")
  }
  
  xattr <- names(attributes(x))
  if (keep == TRUE) {
    attributes(x)[!xattr %in% atr] <- NULL
  } else {
    attributes(x)[xattr %in% atr] <- NULL
  }
  
  return(x)
}




#' Print object and its "default" attributes only.
#' 
#' @param x Object to be printed
#' @param list_attributes Prints the names of all attribute of x, defaults to 
#'   TRUE
#'   
#' @details Before the \option{x} is printed by print.default, all its arguments
#'   not in the default list of [attrs()] are removed.
#'   
#' @author Wolfgang Mader, \email{Wolfgang.Mader@@fdm.uni-freiburg.de}
#' @author Mirjam Fehling-Kaschek, 
#'   \email{mirjam.fehling@@physik.uni-freiburg.de}
#'   
#' @export
print0 <- function(x, list_attributes = TRUE ) {
  if (list_attributes == TRUE) {
    cat("List of all attributes: ", names(attributes(x)), "\n")
  }

  print.default(attrs(x))
}


# Cross-platform parallel-apply. Unix forks via doParallel; Windows uses a
# PSOCK cluster with explicit library-path + variable export. Driven through
# foreach::%dopar% so the caller body is the same on both platforms.
.parallelLapply <- function(X, FUN, cores = 1L,
                            extraExports = character(0),
                            envir = parent.frame()) {

  cores <- max(1L, as.integer(cores))
  if (cores == 1L)
    return(lapply(X, FUN))

  is_windows <- Sys.info()[["sysname"]] == "Windows"

  if (is_windows) {
    cluster <- parallel::makeCluster(cores)
    on.exit({
      parallel::stopCluster(cluster)
      doParallel::stopImplicitCluster()
    }, add = TRUE)
    doParallel::registerDoParallel(cl = cluster)
    parallel::clusterCall(cl = cluster, function(x) .libPaths(x), .libPaths())
    parallel::clusterEvalQ(cluster, suppressMessages(library(dMod)))
    if (length(extraExports) > 0L)
      parallel::clusterExport(cluster, envir = envir,
                              varlist = extraExports)
  } else {
    doParallel::registerDoParallel(cores = cores)
    on.exit(doParallel::stopImplicitCluster(), add = TRUE)
  }

  i <- NULL  # silence R CMD check NSE warning
  loaded_packages <- .packages()
  out <- foreach::foreach(i = seq_along(X),
                          .packages = loaded_packages) %dopar% {
    FUN(X[[i]])
  }
  out
}


`%dopar%` <- foreach::`%dopar%`


