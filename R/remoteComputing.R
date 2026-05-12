#' Detect number of free cores
#'
#' @description Estimates free cores from the 1-min load average.
#' Supports Linux, macOS, and remote machines via SSH.
#' On Windows, returns 1 with a warning (no load average available).
#' Result is floored at 1 so it can be fed straight into
#' `mclapply(mc.cores = ...)` without crashing under heavy load.
#' @param machine character vector of SSH hosts, e.g. "user@@localhost".
#' NULL (default) for the local machine.
#' @return numeric vector of free cores (>= 1) with attributes "ncores" and "used".
#' @export
detectFreeCores <- function(machine = NULL) {
  
  .getLoadAndCores <- function(prefix = NULL) {
    cmd <- function(x) {
      if (!is.null(prefix)) x <- paste("ssh", prefix, x)
      system(x, intern = TRUE)
    }
    
    # Detect OS: use uname for remote, Sys.info() for local
    os <- if (!is.null(prefix)) cmd("uname") else Sys.info()[["sysname"]]
    
    if (grepl("Windows", os, ignore.case = TRUE)) {
      # No load average on Windows — return 1 free core as safe default
      warning("detectFreeCores: load average not available on Windows, returning 1")
      nCores <- parallel::detectCores()
      return(list(free = 1, nCores = nCores, occupied = NA_real_))
    }
    
    if (grepl("Darwin", os, ignore.case = TRUE)) {
      # macOS: sysctl provides load avg as "{ x.xx x.xx x.xx }"
      occupied <- as.numeric(strsplit(cmd("sysctl -n vm.loadavg"), " ")[[1]][2])
      nCores <- as.numeric(cmd("sysctl -n hw.ncpu"))
    } else {
      # Linux: read 1-min load average from /proc/loadavg
      occupied <- as.numeric(strsplit(cmd("cat /proc/loadavg"), " ", fixed = TRUE)[[1]][1])
      nCores <- as.numeric(cmd("nproc --all"))
    }
    
    # Floor at 1: callers feed `free` straight into mclapply(mc.cores = ...)
    # which rejects 0. On a CI runner under heavy compile load the 1-min load
    # average can exceed nCores, which previously gave free = 0 and crashed
    # downstream. Reporting 1 instead of 0 here means "serialise, don't die".
    list(free = max(1L, round(nCores - occupied)), nCores = nCores, occupied = occupied)
  }
  
  if (!is.null(machine)) {
    output <- lapply(machine, .getLoadAndCores)
    freeCores <- vapply(output, `[[`, numeric(1), "free")
    attr(freeCores, "ncores") <- vapply(output, `[[`, numeric(1), "nCores")
    attr(freeCores, "used") <- vapply(output, `[[`, numeric(1), "occupied")
  } else {
    res <- .getLoadAndCores()
    freeCores <- res$free
    attr(freeCores, "ncores") <- res$nCores
    attr(freeCores, "used") <- res$occupied
  }
  
  freeCores
}


#' Run an R expression in the background (only on UNIX)
#' 
#' @description Generate an R code of the expression that is copied via `scp`
#' to any machine (ssh-key needed). Then collect the results.
#' @details `runbg()` generates a workspace from the `input` argument
#' and copies the workspace to the remote machines via `scp`. This will only
#' work if *an ssh-key had been generated and added to the authorized keys
#' on the remote machine*. The code snippet, i.e. the `...` argument, can
#' include several intermediate results but only the last call which is not
#' redirected into a variable is returned via the variable `.runbgOutput`,
#' see example below.
#'
#' Depending on the `compile` and `link` arguments, build-related files are
#' handled as follows:
#' \itemize{
#'   \item `compile = TRUE`: C/C++ source files are transferred and compiled
#'         remotely via `R CMD SHLIB`.
#'   \item `link = TRUE`: Pre-compiled object files (`.o`) are transferred and
#'         linked remotely.
#'   \item Both `FALSE` (default): Existing shared objects (`.so`) are copied
#'         directly.
#' }
#' When `compile` or `link` is `TRUE`, objects of class `obsfn`, `parfn`, or
#' `prdfn` in the workspace automatically get their `modelname` updated to
#' point to the newly built shared object.
#' @param ... Some R code.
#' @param machine Character vector, e.g. `"localhost"` or `"knecht1.fdm.uni-freiburg.de"`
#' or `c("localhost", "localhost")`.
#' @param filename Character, defining the filename of the temporary file. A random
#' file name is chosen if `NULL`.
#' @param input Character vector, the objects in the workspace that are stored
#' into an R data file and copied to the remote machine.
#' @param compile Logical. If `TRUE`, C/C++ source files (`.c`, `.cpp`) are
#' transferred to the remote machine and fully recompiled into a shared object
#' (`.so`). If set to `TRUE`, this overrides `link = TRUE`.
#' @param link Logical. If `TRUE`, only existing object files (`.o`) are
#' transferred to the remote machine and linked into a shared object (`.so`),
#' skipping compilation. If no `.o` files are found, an error is raised.
#' This option is ignored if `compile = TRUE`.
#' @param wait Logical. Wait until executed. If `TRUE`, the code checks if the
#' result file is already present in which case it is loaded. If not present,
#' `runbg()` starts, produces the result and loads it as `.runbgOutput` directly
#' into the workspace. If `wait = FALSE`, `runbg()` starts in the background
#' and the result is only loaded into the workspace when the `get()` function
#' is called, see Value section.
#' @param recover Logical. This option is useful to recover the functions
#' `check()`, `get()`, `purge()` and `terminate()`, e.g. when a session has
#' crashed. Then, the functions are recreated without restarting the job. They
#' can then be used to get the results of a job without having to do it manually.
#' Requires the correct filename, so if the previous `runbg()` was run with
#' `filename = NULL`, you have to specify the filename manually.
#' @param walltime Optional character. Maximum runtime in the format `"HH:MM:SS"`.
#' If exceeded, the job will be terminated.
#' @return List of functions `check()`, `get()`, `purge()` and `terminate()`. 
#' `check()` checks if the result is ready.
#' `get()` copies the result file
#' to the working directory and loads it into the workspace as an object called `.runbgOutput`. 
#' This object is a list named according to the machines that contains the results returned by each
#' machine.
#' `purge()` deletes the temporary folder
#' from the working directory and the remote machines.
#' `terminate()` kills all running processes associated with this job on the remote machines.
#' @export
#' @examples
#' \dontrun{
#' out_job1 <- runbg({
#'          M <- matrix(rnorm(1e2), 10, 10)
#'          solve(M)
#'          }, machine = c("localhost", "localhost"), filename = "job1")
#' out_job1$check()          
#' out_job1$get()
#' result <- .runbgOutput
#' print(result)
#' out_job1$purge()
#' }
#' \dontrun{
#' # Recover a runbg job with the option "recover"
#' out_job1 <- runbg({
#'          M <- matrix(rnorm(1e2), 10, 10)
#'          solve(M)
#'          }, machine = c("localhost", "localhost"), filename = "job1")
#' Sys.sleep(1)
#' remove(out_job1)
#' try(out_job1$check())
#' out_job1 <- runbg({
#'   "This code is not run"
#' }, machine = c("localhost", "localhost"), filename = "job1", recover = TRUE)
#' out_job1$get()
#' result <- .runbgOutput
#' print(result)
#' out_job1$purge()
#' }
runbg <- function(..., machine = "localhost", filename = NULL, input = ls(.GlobalEnv), compile = FALSE, link = FALSE, wait = FALSE, recover = FALSE, walltime = NULL) {
  
  expr <- as.expression(substitute(...))
  nmachines <- length(machine)
  
  # compile takes precedence over link
  if (compile) link <- FALSE
  
  # Generate a random filename if none is provided
  if (is.null(filename))
    filename <- paste0("tmp_", paste(sample(c(0:9, letters), 5, replace = TRUE), collapse = ""))
  
  filename0 <- filename
  filename <- paste(filename, 1:nmachines, sep = "_")
  
  # Initialize output list
  out <- structure(vector("list", 4), names = c("check", "get", "purge", "terminate"))
  
  # Check whether results are ready on all machines
  out[[1]] <- function() {
    
    check.out <- sapply(1:nmachines, function(m) length(suppressWarnings(
      system(paste0("ssh ", machine[m], " ls ", filename[m], "_folder/ | grep -x ", filename[m], "_result.RData"), 
             intern = TRUE))))
    
    if (all(check.out > 0)) {
      cat("Result is ready!\n")
      return(TRUE)
    } else if (any(check.out > 0)) {
      cat("Result from machines", paste(which(check.out > 0), collapse = ", "), "are ready.")
      return(FALSE)
    } else {
      cat("Not ready!\n") 
      return(FALSE)
    }
    
  }
  
  # Fetch result files from remote machines and load into workspace
  out[[2]] <- function() {
    
    result <- structure(vector(mode = "list", length = nmachines), names = machine)
    for (m in 1:nmachines) {
      .runbgOutput <- NULL
      system(paste0("scp ", machine[m], ":", filename[m], "_folder/", filename[m], "_result.RData ./"), ignore.stdout = TRUE, ignore.stderr = TRUE)
      check <- try(load(file = paste0(filename[m], "_result.RData")), silent = TRUE) 
      if (!inherits(check, "try-error")) result[[m]] <- .runbgOutput
    }
    
    .GlobalEnv$.runbgOutput <- result
    
  }
  
  # Remove temporary folders and files on remote machines and locally
  out[[3]] <- function() {
    
    for (m in 1:nmachines) {
      folder_exists <- suppressWarnings(
        system(paste0("ssh ", machine[m], " '[ -d ", filename[m], "_folder ] && echo 1 || echo 0'"), 
               intern = TRUE)
      )
      if (folder_exists == "1") {
        system(paste0("ssh ", machine[m], " rm -r ", filename[m], "_folder"))
      }
      
      rout_exists <- suppressWarnings(
        system(paste0("ssh ", machine[m], " '[ -f ", filename[m], ".Rout ] && echo 1 || echo 0'"), 
               intern = TRUE)
      )
      if (rout_exists == "1") {
        system(paste0("ssh ", machine[m], " rm ", filename[m], ".Rout"))
      }
    }
    
    local_files <- list.files(pattern = paste0(filename0, ".*"))
    if (length(local_files) > 0) {
      system(paste0("rm ", filename0, "*"))
    }
  }
  
  # Kill all processes associated with this job on the remote machines
  out[[4]] <- function() {
    for (m in 1:nmachines) {
      pids <- suppressWarnings(
        system(paste0("ssh ", machine[m], 
                      " 'ps aux | grep \"", filename[m], "\" | grep -v grep'"), 
               intern = TRUE)
      )
      
      if (length(pids) > 0) {
        running_pids <- sapply(strsplit(pids, "\\s+"), function(x) {
          if (grepl("R", x[8])) x[2] else NULL
        })
        running_pids <- running_pids[!sapply(running_pids, is.null)]
        
        if (length(running_pids) > 0) {
          system(paste0("ssh ", machine[m], " 'kill ", paste(running_pids, collapse = " "), "'"))
          cat("Terminated", length(running_pids), "running processes on", machine[m], "\n")
        } else {
          cat("No running processes found on", machine[m], "\n")
        }
      }
    }
  }
  
  # Recover control functions without re-submitting the job
  if (recover) return(out)
  
  # If result files already exist locally and wait = TRUE, load them directly
  resultfile <- paste(filename, "result.RData", sep = "_")
  if (all(file.exists(resultfile)) & wait) {
    
    result <- structure(vector(mode = "list", length = nmachines), names = machine)
    for (m in 1:nmachines) {
      load(file = resultfile[m])
      result[[m]] <- .runbgOutput
    }
    .GlobalEnv$.runbgOutput <- result
    return(out)
  }
  
  # Save current workspace to be transferred to remote machines
  save(list = input, file = paste0(filename0, ".RData"), envir = .GlobalEnv)
  
  # Collect currently loaded packages to replicate the library state remotely
  pack <- sapply(strsplit(search(), "package:", fixed = TRUE), function(v) v[2])
  pack <- pack[!is.na(pack)]
  pack <- paste(paste0("try(library(", pack, "))"), collapse = "\n")
  
  output <- ".runbgOutput"
  
  # Compiler flags mirroring compile() in tools.R
  # Written into a shell script to avoid quoting issues with nested SSH commands
  if (compile || link) {
    
    if (compile) {
      sourcefiles <- paste(
        c(list.files(pattern = glob2rx("*.c")), list.files(pattern = glob2rx("*.cpp"))),
        collapse = " "
      )
    } else {
      object_files <- Sys.glob("*.o")
      if (length(object_files) == 0)
        stop("No .o files found for linking! You must compile first.")
      sourcefiles <- paste(object_files, collapse = " ")
    }
    
  }
  
  if (compile || link) {
    # R code to load the newly built shared object and update modelnames of known function objects
    objfns <- 'obj.fns <- ls()[sapply(ls(), function(nm) inherits(get(nm, envir=.GlobalEnv), c("obsfn", "parfn", "prdfn")))]'
    setmn <- sprintf('for (o in obj.fns) eval(parse(text=paste0("modelname(", o, ") <- \'%s\'")))', paste0(filename0, "_shared_object"))
    load_so <- paste0("dyn.load('", filename0, "_shared_object.so')")
  } else {
    objfns <- NULL
    setmn <- NULL
    load_so <- NULL
  }
  
  # Assemble the R script to be executed on each remote machine
  program <- lapply(1:nmachines, function(m) paste(
    c(
      pack,
      paste0("setwd('~/", filename[m], "_folder')"),
      "rm(list = ls())",
      if (!is.null(walltime)) paste0("Sys.setenv(R_TIMEOUT='", walltime, "')"),
      paste0("load('", filename0, ".RData')"),
      objfns,
      setmn,
      if (!is.null(load_so)) {
        load_so
      } else {
        c("files <- list.files(pattern = '\\\\.so$')", "for (f in files) dyn.load(f)")
      },
      paste0(".node <- ", m),
      if (!is.null(walltime)) {
        paste0(".runbgOutput <- try(tools::pskill(Sys.getpid(), tools::SIGALRM, ", walltime, "); ", as.character(expr), ")")
      } else {
        paste0(".runbgOutput <- try(", as.character(expr), ")")
      },
      paste0("save(", output ,", file = '", filename[m], "_result.RData')")
    ),
    collapse = "\n"
  ))
  
  for (m in 1:nmachines) {
    
    # Write the R script for this machine
    cat(program[[m]], file = paste0(filename[m], ".R"))
    
    # Create remote working directory
    system(paste0("ssh ", machine[m], " mkdir -p ", filename[m], "_folder/"), 
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    # Clear any leftover files from previous runs
    system(paste0("ssh ", machine[m], " rm -r ", filename[m], "_folder/*"), 
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    # Transfer workspace
    system(paste0("scp ", getwd(), "/", filename0, ".RData* ", machine[m], ":", filename[m], "_folder/"))
    
    # Transfer R script
    system(paste0("scp ", getwd(), "/", filename[m], ".R* ", machine[m], ":", filename[m], "_folder/"))
    
    # Transfer C/C++ source and object files (always; harmless if none exist)
    system(paste0("scp ", getwd(), "/*.c ", machine[m], ":", filename[m], "_folder/"),
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    system(paste0("scp ", getwd(), "/*.cpp ", machine[m], ":", filename[m], "_folder/"),
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    system(paste0("scp ", getwd(), "/*.o ", machine[m], ":", filename[m], "_folder/"),
           ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    # Transfer shared objects only when no remote build is requested
    if (!compile && !link) {
      system(paste0("scp ", getwd(), "/*.so ", machine[m], ":", filename[m], "_folder/"),
             ignore.stdout = TRUE, ignore.stderr = TRUE)
    }
    
    # When recompiling, remove stale .o files so R CMD SHLIB starts clean
    if (compile) {
      system(paste0("ssh ", machine[m], " 'rm -f ", filename[m], "_folder/*.o'"),
             ignore.stdout = TRUE, ignore.stderr = TRUE)
    }
    
    # Write and transfer compile script per machine, then build the remote command
    compile_cmd <- ""
    if (compile || link) {
      compile_script_content <- paste(
        "#!/bin/bash",
        "export PKG_CFLAGS='-O2 -DNDEBUG -w -fPIC'",
        "export PKG_CXXFLAGS='-O2 -DNDEBUG -w -fPIC'",
        paste0("export PKG_CPPFLAGS=-I$(Rscript -e ", "'cat(system.file(\"include\", package=\"CppODE\"))'", ")"),
        paste0("cd ", filename[m], "_folder"),
        paste0("R CMD SHLIB ", sourcefiles, " -o ", filename0, "_shared_object.so"),
        sep = "\n"
      )
      compile_script_file <- paste0(filename[m], "_compile.sh")
      cat(compile_script_content, file = compile_script_file)
      system(paste0("scp ", getwd(), "/", compile_script_file, " ", machine[m], ":", filename[m], "_folder/"))
      compile_cmd <- paste0("bash ", filename[m], "_folder/", compile_script_file)
    }
    
    # Execute: compile/link if requested, then run the R script
    # OMP_NUM_THREADS=1 and MKL_NUM_THREADS=1 ensure single-threaded execution per job
    system(paste0(
      "ssh ", machine[m], 
      " 'export OMP_NUM_THREADS=1 && export MKL_NUM_THREADS=1",
      if (nzchar(compile_cmd)) paste0(" && ", compile_cmd),
      " && R CMD BATCH --vanilla ", filename[m], "_folder/", filename[m], ".R'"
    ), intern = FALSE, wait = wait)
  }
  
  if (wait) {
    out$get()
    out$purge()
  } else {
    return(out)
  }
  
}

