#' Run any R function on a remote HPC system with SLURM
#'
#' @description
#' Generates R and bash scripts, transfers them to a remote HPC system via SSH,
#' and executes the given R code in parallel using the SLURM batch manager.
#' The function handles workspace export, job submission, and result retrieval.
#'
#' @details
#' `distributed_computing()` generates R and bash scripts designed to run
#' on an HPC system managed by SLURM. The current R workspace together with the
#' scripts are exported and transferred to the remote system via SSH.  
#' If ssh-key authentication is not possible, the SSH password can be provided and
#' is used by `sshpass` (which must be installed locally).
#'
#' The code to be executed remotely is passed to the `...` argument; its final
#' output is stored in `cluster_result`, which can be loaded in the local
#' workspace via the `get()` function.
#'
#' It is possible to run multiple repetitions of the same program (via `no_rep`)
#' or to pass a list of parameter arrays through `var_values`. Parameters that
#' vary between runs must be named `var_i`, where *i* matches the index
#' of the corresponding array in `var_values`.
#'
#' @param ... R code to be remotely executed. Parameters to be changed for each run
#' must be named `var_i` (see "Details").
#' @param jobname Unique name (character) for the run. Existing runs with the same
#' name will be overwritten. Must not contain the string "Minus".
#' @param partition SLURM partition name to use. Default is `"single"`.
#' @param cores Number of cores per node. Values above 16 may limit available nodes.
#' @param nodes Number of nodes per task. Default is 1; typically should not be changed.
#' @param mem_per_core Memory per CPU core in GB. Default is 2 GB.
#' @param walltime Maximum runtime in format `"hh:mm:ss"`. Default is 1 hour.
#' @param ssh_passwd Password string for SSH authentication via `sshpass`.
#' Optional, and only used if key-based authentication is unavailable.
#' @param machine SSH address of the remote HPC system, e.g. `"user@@cluster"`.
#' @param var_values List of parameter arrays. Each array corresponds to one variable
#' `var_i`. The length of each array determines the number of SLURM array jobs.
#' Mutually exclusive with `no_rep`.
#' @param no_rep Integer number of repetitions (mutually exclusive with `var_values`).
#' @param recover Logical; if `TRUE`, no computation is performed. Instead,
#' the returned list of functions `check()`, `get()`, and `purge()`
#' can be used to interact with previously submitted jobs.
#' @param purge_local Logical; if `TRUE`, the `purge()` function also
#' deletes local result files.
#' @param compile Logical; if `TRUE`, all C/C++ source files (`*.c`, `*.cpp`)
#' are transferred to the cluster and fully recompiled into shared objects (`.so`).
#' If set to `TRUE`, this overrides `link = TRUE`.
#' @param link Logical; if `TRUE`, only existing object files (`*.o`) are
#' transferred to the cluster and linked into shared objects (`.so`),
#' skipping compilation. If no `.o` files are found, an error is raised.
#' This option is ignored if `compile = TRUE`.
#' @param custom_folders Named vector with exactly three elements: `"compiled"`,
#' `"output"`, and `"tmp"`. Each value is a relative path specifying where
#' compiled files, temporary data, and output results should be stored.
#' If `NULL`, all operations occur in the current working directory.
#' @param resetSeeds Logical; if `TRUE` (default), removes `.Random.seed`
#' from the transferred workspace to ensure each node has independent random seeds.
#' @param returnAll Logical; if `TRUE` (default), retrieves all remote files.
#' If `FALSE`, only result files (`*result.RData`) are fetched.
#'
#' @return
#' A list containing three functions:
#' \itemize{
#'   \item `check()` – Checks whether all remote results are complete.
#'   \item `get()` – Downloads results and loads them into
#'         `cluster_result` in the local workspace.
#'   \item `purge()` – Deletes temporary remote files; optionally removes local ones.
#' }
#'
#' @examples
#' \dontrun{
#' out_distributed_computing <- distributed_computing(
#' {
#'   mstrust(
#'     objfun=objective_function,
#'     center=outer_pars,
#'     studyname = "study",
#'     rinit = 1,
#'     rmax = 10,
#'     fits = 48,
#'     cores = 16,
#'     iterlim = 700,
#'     sd = 4
#'   )
#' },
#' jobname = "my_name",
#' partition = "single",
#' cores = 16,
#' nodes = 1,
#' mem_per_core = 2,
#' walltime = "02:00:00",
#' ssh_passwd = "password",
#' machine = "cluster",
#' var_values = NULL,
#' no_rep = 20,
#' recover = F,
#' compile = F,
#' link = F
#' )
#' out_distributed_computing$check()
#' out_distributed_computing$get()
#' out_distributed_computing$purge()
#' result <- cluster_result
#' print(result)
#' 
#' 
#' # calculate profiles
#' var_list <- profile_pars_per_node(best_fit, 4)
#' profile_jobname <- paste0(fit_filename,"_profiles_opt")
#' method <- "optimize"
#' profiles_distributed_computing <- distributed_computing(
#'   {
#'     profile(
#'       obj = obj,
#'       pars =  best_fit,
#'       whichPar = (as.numeric(var_1):as.numeric(var_2)),
#'       limits = c(-5, 5),
#'       cores = 16,
#'       method = method,
#'       stepControl = list(
#'         stepsize = 1e-6,
#'         min = 1e-4, 
#'         max = Inf, 
#'         atol = 1e-2,
#'         rtol = 1e-2, 
#'         limit = 100
#'       ),
#'       optControl = list(iterlim = 20)
#'     )
#'   },
#'   jobname = profile_jobname,
#'   partition = "single",
#'   cores = 16,
#'   nodes = 1,
#'   walltime = "02:00:00",
#'   ssh_passwd = "password",
#'   machine = "cluster",
#'   var_values = var_list,
#'   no_rep = NULL,
#'   recover = F,
#'   compile = F,
#'   link = F
#' )
#' profiles_distributed_computing$check()
#' profiles_distributed_computing$get()
#' profiles_distributed_computing$purge()
#' profiles  <- NULL
#' for (i in cluster_result) {
#'   profiles <- rbind(profiles, i)
#' }
#' }
#' 
#' @export
distributed_computing <- function(
    ...,
    jobname,
    partition = "single",
    cores = 16,
    nodes = 1,
    mem_per_core = 2,
    walltime = "01:00:00",
    ssh_passwd = NULL,
    machine = "cluster",
    var_values = NULL,
    no_rep = NULL,
    recover = TRUE,
    purge_local = FALSE,
    compile = FALSE,
    link = FALSE,
    custom_folders = NULL,
    resetSeeds = TRUE,
    returnAll = TRUE
){
  original_wd <- getwd()
  if (is.null(custom_folders)) {
    output_folder_abs <- "./"
  } else if(!is.null(custom_folders) & !all(length(custom_folders) == 3 & sort(names(custom_folders)) == c("compiled", "output", "tmp"))) {
    warning("'custom_folders' must be named vector with exact three elements:\n
            'compiled', 'output', 'tmp', containing relative paths to the resp folders\n
            input is wrong, ignored.\n")
  } else {
    compiled_folder <- custom_folders["compiled"]
    output_folder <- custom_folders["output"]
    tmp_folder <- custom_folders["tmp"]
    
    system(paste0("cp ", compiled_folder, "* ", tmp_folder))
    
    setwd(output_folder)
    output_folder_abs <- getwd()
    
    setwd(original_wd)
    setwd(tmp_folder)
  }
  
  # Safety rule: never compile and link at the same time
  if (compile) link <- FALSE
  
  
  on.exit(setwd(original_wd))
  
  # - definitions - #
  
  # relative path to the working directory, will now allways be used
  
  wd_path <- paste0("./",jobname, "_folder/")
  data_path <- paste0(getwd(),"/",jobname, "_folder/")
  
  # number of repetitions
  if(!is.null(no_rep) & is.null(var_values)) {
    num_nodes <- no_rep - 1
  } else if(is.null(no_rep) & !is.null(var_values)) {
    num_nodes <- length(var_values[[1]]) - 1
  } else {
    stop("I dont know what you want how often done. Please set either 'no_rep' or pass 'var_values' (_not_ both!)")
  }
  
  # define the ssh command depending on 'sshpass' being used
  if(is.null(ssh_passwd)){
    ssh_command <- "ssh "
    scp_command <- "scp "
  } else {
    ssh_command <- paste0("sshpass -p ", ssh_passwd, " ssh ")
    scp_command <- paste0("sshpass -p ", ssh_passwd, " scp ")
  }
  
  # - output functions - #
  # Structure of the output 
  out <- structure(vector("list", 3), names = c("check", "get", "purge"))
  
  # check function
  out[[1]] <- function() {
    
    result_length <- length(
      suppressWarnings(
        system(
          paste0(ssh_command, machine, " 'ls ", jobname, "_folder/ | egrep *result.RData'"),
          intern = TRUE)
      )
    )
    
    if (result_length == num_nodes +1) {
      cat("Result is ready!\n")
      return(TRUE)
    }
    else if (result_length  < num_nodes +1) {
      cat("Result from", result_length, "out of", (num_nodes +1), "nodes are ready.")
      return(FALSE)
    }
    setwd(original_wd)
    
  }
  
  # get function
  out[[2]] <- function () {
    # copy all files back
    if (returnAll == T) {
      system(
        paste0(
          "mkdir -p ", output_folder_abs, "/", jobname,"_folder/results/; ",
          ssh_command, "-n ", machine, # go to remote
          " 'tar -C ", jobname, "_folder", " -czf - ./'", # compress all files on remote)
          " | ", # pipe to local
          "",
          "tar -C ", output_folder_abs, "/", jobname,"_folder/results/ -xzf -"
        )
      )
    } else {
      # copy only result files back
      system(
        paste0(
          "mkdir -p ", output_folder_abs, "/", jobname,"_folder/results/; ",
          ssh_command, "-n ", machine, " '",
          "find ", jobname, "_folder -type f -name \"*result.RData\" -exec tar -czf - {} +'",
          " | tar --strip-components=1 -xz -C ", shQuote(paste0(output_folder_abs, "/", jobname,"_folder/results/"))
        )
      )
    }
    
    
    # get list of all currently available output files
    # setwd(paste0(jobname,"_folder/results"))
    result_list <- structure(vector(mode = "list", length = num_nodes+1))
    result_files <- list.files(path=paste0(output_folder_abs,"/",jobname,"_folder/results/"),pattern = glob2rx("*result.RData"))
    # setwd("../../")
    
    # result_files <- Sys.glob(file.path(paste0(wd_path, "/results/*RData")))
    
    for (i in seq(1, length(result_files))) {
      cluster_result <- NULL
      check <- try(load(file = paste0(output_folder_abs,"/",jobname,"_folder/results/",result_files[i])), silent = TRUE) 
      if (!inherits("try-error", check)) result_list[[i]] <- cluster_result
    }
    
    if (length(result_files) != num_nodes +1) {
      cat("\n\tNot all results ready\n")
    }
    setwd(original_wd)
    # results_cluster <- Sys.glob(paste0(wd_path, "*.RData")) %>% map_dfr(load)
    .GlobalEnv$cluster_result <- result_list
  }
  
  
  
  # purge function
  out[[3]] <- function (purge_local = FALSE) {
    # remove files remote
    system(
      paste0(
        ssh_command, machine, " rm -rf ", jobname, "_folder"
      )
    )
    # also remove local files if want so
    if (purge_local) {
      system(
        paste0("rm -rf ", output_folder_abs, "/", jobname,"_folder")
      )
    }
    setwd(original_wd)
  }
  
  # if recover == T, stop here
  if(recover) {
    setwd(original_wd)
    return(out)
  } 
  
  
  
  
  
  
  
  
  
  # - calculation - #
  # create wd for this run
  system(
    paste0("rm -rf ",jobname,"_folder")
  )
  system(
    paste0("mkdir ",jobname,"_folder/")
  )
  
  
  
  # export current workspace
  save.image(file = paste0(wd_path,jobname, "_workspace.RData")) 
  
  
  # WRITE R
  
  
  # generate list of currently loaded packages
  package_list <- sapply(strsplit(search(), "package:", fixed = TRUE), function(v) v[2])
  package_list <- package_list[!is.na(package_list)]
  package_list <- paste(paste0("try(library(", package_list, "))"), collapse = "\n")
  if (compile | link) {
    objfns <- 'obj.fns <- ls()[sapply(ls(), function(nm) inherits(get(nm, envir=.GlobalEnv), c("obsfn", "parfn", "prdfn")))]'
    setmn <- sprintf('for (o in obj.fns) eval(parse(text=paste0("modelname(", o, ") <- \'%s\'")))\n', paste0(jobname, "_shared_object"))
    load_so <- paste0("dyn.load('",jobname,"_shared_object.so')")
  } else {
    setmodelname <- ""
    setmn <-""
    load_so <- ""
  }
  
  
  
  # generate parameter lists
  if (!is.null(var_values)) {
    var_list <- paste(
      lapply(
        seq(1,length(var_values)),
        function(i) {
          if( class(var_values[[i]]) == "character") {
            paste0("var_values_", i, "=c('", paste(var_values[[i]], collapse="','"),"')")
          } else {
            paste0("var_values_", i, "=c(", paste(var_values[[i]], collapse=","),")")
          }
          
          
        }
      ),
      collapse = "\n"
    )
    # cat(variable_list)
    
    # List of all names of parameters that will be changes between runs
    var_names <- paste(lapply(seq(1,length(var_values)), function(i) paste0("var_",i)))
    
    # Variables per run
    var_per_run <- paste(
      lapply(
        seq(1, length(var_values)),
        function(i) {
          paste0("var_", i, "=var_values_",i,"[(as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + 1)]")
        }
      ),
      collapse = "\n"
    )
  } else {
    var_list <- ""
    var_per_run <- ""
  }
  
  # cat(var_per_run)
  
  
  # define fixed pars
  fixedpars <- paste(
    "node_ID = Sys.getenv('SLURM_ARRAY_TASK_ID')",
    "job_ID = Sys.getenv('SLURM_JOB_ID')",
    paste0("jobname = ", "'",jobname,"'"),
    sep = "\n"
  )
  
  
  # WRITE R
  expr <- as.expression(substitute(...))
  # expr <- as.expression(substitute(called_function))
  cat(
    paste(
      "#!/usr/bin/env Rscript",
      "",
      "# Load packages",
      package_list,
      "try(library(tidyverse))",
      "",
      "# Load environment",
      paste0("load('",jobname,"_workspace.RData')"),
      "",
      if (resetSeeds == TRUE & exists(".Random.seed")) {
        paste0("# remove random seeds\nrm(.Random.seed)\nset.seed(as.numeric(Sys.getenv('SLURM_JOB_ID')))")
        
      },
      "",
      "# load shared object if precompiled",
      setmodelname,
      setmn,
      load_so,
      "",
      "# List of variablevalues",
      var_list,
      "",
      "# Define variable values per run",
      var_per_run,
      "",
      "# Fixed parameters",
      fixedpars,
      "",
      "",
      "",
      "# Paste function call",
      paste0("cluster_result <- try(", as.character(expr),")"),
      sep = "\n",
      "save(cluster_result, file = paste0(jobname,'_', node_ID, '_result.RData'))" #'_',job_ID, 
    ),
    file = paste0(wd_path, jobname,".R")
  )
  
  
  
  
  # WRITE BASH
  cat(
    paste(
      "#!/bin/bash",
      "",
      "# Job name",
      paste0("#SBATCH --job-name=",jobname),
      "# Define format of output, deactivated",
      paste0("#SBATCH --output=",jobname,"_%j-%a.out"),
      "# Define format of errorfile, deactivated",
      paste0("#SBATCH --error=",jobname,"_%j-%a.err"),
      "# Define partition",
      paste0("#SBATCH --partition=", partition),
      "# Define number of nodes per task",
      paste0("#SBATCH --nodes=", nodes),
      "# Define number of cores per node",
      paste0("#SBATCH --ntasks-per-node=",cores),
      "# Define walltime",
      paste0("#SBATCH --time=",walltime),
      "# Define of repetition",
      paste0("#SBATCH -a 0-", num_nodes),
      "# memory per CPU core",
      paste0("#SBATCH --mem-per-cpu=", mem_per_core, "gb"),
      "",
      "",
      "# Load compiler modules",
      "module load compiler/gnu/13.3",
      "# Load R modules",
      "module load math/R",
      # paste0("export OPENBLAS_NUM_THREADS=",cores),
      paste0("export OMP_NUM_THREADS=","1"), # paste0("export OMP_NUM_THREADS=",cores),
      paste0("export MKL_NUM_THREADS=", "1"), # paste0("export MKL_NUM_THREADS=",cores),
      "",
      "# Run R script",
      paste0("Rscript ", jobname, ".R"),
      sep = "\n" 
    ),
    file = paste0(wd_path,jobname,".sh")
  )
  
  
  if (compile) {
    # --- FULL RECOMPILATION (.cpp/.c -> .so) ---
    compile_files <- Sys.glob(paste0("*.cpp"))
    compile_files <- append(compile_files, Sys.glob(paste0("*.c")))
    compile_files <- paste(compile_files, collapse = " ")
    
    tar_locale <- paste0("tar -jcf - ", compile_files, " ", wd_path, "*")
    tar_remote <- paste0("tar -C ./ -jxf - ; mv -t ./", jobname, "_folder ", compile_files, "; ")
    
    sourcefiles <- paste(
      c(list.files(pattern = glob2rx("*.c")), list.files(pattern = glob2rx("*.cpp"))),
      collapse = " "
    )
    
    compile_remote <- paste0(
      "module load math/R; R CMD SHLIB ", sourcefiles, " -o ", jobname, "_shared_object.so; "
    )
    
  } else if (link) {
    # --- LINK ONLY (.o -> .so) ---
    object_files <- Sys.glob("*.o")
    if (length(object_files) == 0)
      stop("No .o files found for linking! You must compile first.")
    
    # Remove any old .so files before linking
    # unlink(list.files(pattern = "(\\.so)$"))
    
    compile_files <- paste(object_files, collapse = " ")
    tar_locale <- paste0("tar -jcf - ", compile_files, " ", wd_path, "*")
    tar_remote <- paste0("tar -C ./ -jxf - ; mv -t ./", jobname, "_folder ", compile_files, "; ")
    
    compile_remote <- paste0(
      "module load math/R; R CMD SHLIB ", paste(object_files, collapse = " "),
      " -o ", jobname, "_shared_object.so; "
    )
    
  } else {
    # --- NO BUILD ACTION (.so/.o already available) ---
    compile_files <- Sys.glob(paste0("*.so"))
    compile_files <- append(compile_files, Sys.glob(paste0("*.o")))
    compile_files <- paste(compile_files, collapse = " ")
    
    tar_locale <- paste0("tar -jcf - ", compile_files, " ", wd_path, "*")
    tar_remote <- paste0("tar -C ./ -jxf - ; mv -t ./", jobname, "_folder ", compile_files, "; ")
    compile_remote <- ""
  }
  
  ##
  # transfer and run files
  system(
    paste0(
      tar_locale, # compress all files in the local working dir
      " | ", ssh_command, machine, # pipe to ssh session on remote
      " 'if [ -d ", jobname, "_folder ]; then rm -Rf ", jobname,"_folder; fi ;", # remove folder if it exists
      " mkdir -p ", jobname,"_folder; ", # create new wd on remote
      tar_remote, # uncompress files in wd on remote, if necessary move files
      "cd ", jobname, "_folder; ", # change in said wd
      compile_remote, # compile files if said so, if not nothing happen
      "sbatch ", jobname, ".sh'" # start bash script
    )
  )
  
  setwd(original_wd)
  return(out)
}



#' Generate parameter list for distributed profile calculation
#' 
#' @description Generates list of `WhichPar` entries to facillitate distribute
#' profile calculation.
#' @details Lists to split the parameters for which the profiles are calculated
#' on the different nodes.
#' 
#' @param parameters list of parameters 
#' @param fits_per_node numerical, number of parameters that will be send to each node.
#' @param side determine if both sides are calculated (default) or if the profiles are split in 'left' and 'right' for calculation
#' 
#' @return List with two arrays: `from` contains the number of the starting
#' parameter, while `to` stores the respective upper end of the parameter list
#' per node.
#' @examples
#' \dontrun{
#' parameter_list <- setNames(1:10, letters[1:10])
#' var_list <- profile_pars_per_node(parameter_list, 4)
#' }
#' 
#' @export
profile_pars_per_node <- function(parameters, fits_per_node, side = c("both", "split")[1]) {
  # sanitize side input: must be either "left", "right" or "both"
  if (!(side %in% c("both", "split"))) {
    stop("'side' must be either 'both' or 'split'")
  }
  
  # get the number of parameters
  n_pars <- length(parameters)
  
  # Get number of fits per node
  fits_per_node <- fits_per_node
  
  # determine the number of nodes necessary
  no_nodes <- 1:ceiling(n_pars/fits_per_node)
  
  # generate the lists which parameters are send to which node
  pars_from <- fits_per_node
  pars_to_vec <- fits_per_node
  while (pars_from < (n_pars)) {
    pars_from <- pars_from + fits_per_node
    pars_to_vec <- c(pars_to_vec, pars_from)
  }
  pars_to_vec[length(pars_to_vec)] <- n_pars
  
  pars_from_vec <- c(1, pars_to_vec+1)
  pars_from_vec <- head(pars_from_vec, -1)
  
  # adjust for sides 
  if (side == "both") {
    side_vec <-  rep("both", length(no_nodes))
  } else {
    # split pars_to_vec and pars_from_vec by repeating each element twice
    pars_to_vec <- rep(pars_to_vec, each = 2)
    pars_from_vec <- rep(pars_from_vec, each = 2)
    side_vec <- rep(c("left", "right"), length(no_nodes))
  }
  
  out <- list(from=pars_from_vec, to=pars_to_vec, side = side_vec)
  
  return(out)
}



## Use Julia to calculate steady states -----------------------------------------

#' Install the julia setup  
#' 
#' @description Installs Julia and the necessary julia packages
#' 
#' 
#' @param installJulia boolean, default `false`. If set to true, juliaup and via this then Julia is installed. 
#' @param installJuliaPackages boolean, default `true`. If set to true, the necessary packages are installed.
#' @param JuliaProjectPath string, default `NULL`. Allows for installing the required packages to a separate project environment instead of the global package environment. Also need to specify same JuliaProjectPath to steadyStateToolJulia().
#' 
#' @return nothing
#' 
#' @export
installJuliaForSteadyStates <- function(installJulia = FALSE, installJuliaPackages = TRUE, JuliaProjectPath = NULL) {
  
  tryCatch(
    {
      system("git clone git@github.com:SeverinBang/JuliaSteadyStates.git ~/.JuliaSteadyStates/")
    },
    finally = {
      cat("github.com:SeverinBang/JuliaSteadyStates.git could not be cloned (again), check if ~/.JuliaSteadyStates already exists, if not write Severin your github username to be added to the repository")
    }
  )
  
  # install Julia
  if (installJulia) {
    system("sh -i ~/.JuliaSteadyStates/installJuliaUp.sh -y")
    system("juliaup add release")
    try(system("source ~/.bashrc"))
  }

  # install packages
  if (installJuliaPackages) {
    # Sensible to use JuliaProjectPath = "~/.JuliaSteadyStates/"
    if(!is.null(JuliaProjectPath)){
      ProjectPathSetter = paste0("--project='", JuliaProjectPath, "' ")
    } else {
      ProjectPathSetter = ""
    }
    system(paste0("julia ", ProjectPathSetter, "-e 'using Pkg; Pkg.add([\"CSV\", \"DataFrames\", \"Symbolics\", \"SymbolicUtils\", \"Catalyst\", \"Graphs\"])'"))
  }
  
}


#' Calculate the steady states of a given model  
#' 
#' @description Uses julia to calculate the steady state transformations
#' 
#' @param el the equation list
#' @param forcings vector of strings, default `c("","")`. The names of the forcings which will be set to zero before solving for the steady state equations.
#' @param neglect vector of strings, default `c("","")`. The names of the variables which will be neglected as fluxParameters and therefore will not be solved for.
#' @param verboseLevel integer, default `1`. The level of verbosity of the output, right now only 1 (all) and 0 (no) is implemented.
#' @param testSteadyState boolean, default `true`
#' @param JuliaPath string, default `NULL`. If specified, uses julia executable from given path.
#' @param FileExportPath string, default `getwd()`. Directory to which .csv files for transfer of equations between R and Julia are saved.
#' @param JuliaProjectPath string, default `NULL`. If specified, uses julia local project environment from given path, otherwise uses global package environment for code execution.
#' 
#' @return named vector with the steady state transformations. The names are the inner, the values are the outer parameters
#' 
#' @export
steadyStateToolJulia <- function(
    el,
    forcings = NULL,
    neglect = NULL,
    verboseLevel = 1,
    testSteadyState = TRUE,
    JuliaPath = NULL,
    FileExportPath = NULL,
    JuliaProjectPath = NULL
) {
  # prepare things:
  if (is.null(FileExportPath)) {
    myWD <- getwd()
  } else {
    myWD <- FileExportPath
  }
  dModEqnFileName = "EquationsForSteadyStates"
  
  dMod::write.eqnlist(el, file = file.path(myWD, paste0(dModEqnFileName, ".csv")))
  
  inputPath <- file.path(myWD, paste0(dModEqnFileName, ".csv"))
  fileName <- "SteadyStatesFromJulia"
  
  if (is.null(forcings)) {
    forcings <- c("","")
  }
  
  if (is.null(neglect)) {
    neglect <- c("","")
  }
  
  # load julia
  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    warning("The 'JuliaCall' package must be installed.")
    return(NULL)
  }
  
  if (!is.null(JuliaProjectPath)) {
    Sys.setenv(JULIA_PROJECT = JuliaProjectPath)
  }
  
  if (!is.null(JuliaPath)) {
    JuliaCall::julia_setup(JULIA_HOME = JuliaPath)
  } else if (dir.exists(file.path(Sys.getenv("HOME"),".juliaup/bin"))) {
    JuliaCall::julia_setup(JULIA_HOME = file.path(Sys.getenv("HOME"),".juliaup/bin"))
  } else if (file.exists("/usr/bin/julia")) {
    JuliaCall::julia_setup(JULIA_HOME = "/usr/bin/")
  } else {
    stop("No Julia installation found, please use juliaup to install julia.")
    return(NULL)
  }
 
  
  
  # call the julia steady state tool:
  JuliaCall::julia_source(file.path(Sys.getenv("HOME"),".JuliaSteadyStates/ODESteadyStateTrafo_function.jl"))
  
  JuliaCall::julia_call("determineSteadyStateTrafos", inputPath, forcings, neglect, myWD, fileName, verboseLevel = JuliaCall::julia_eval(paste0("Int(", verboseLevel, ")")), testSteadyState = ifelse(testSteadyState, 1, 0)) #JuliaCall::julia_eval(paste0(ifelse(testSteadyState, "true", "false")))
  
  # load the results
  steadyStatesFile = read.csv(paste0(myWD,"/",fileName, ".csv" ), dec = ".", sep = ",")
  
  steadyStates = data.table(keys = steadyStatesFile$Keys, values = steadyStatesFile$Values)
  steadyStates = steadyStates[!(keys %in% forcings)]
  
  
  sstates <- steadyStates$values
  
  sstates <- str_replace_all(sstates, "_initBySteadyStateTool", "")
  
  names(sstates) <- steadyStates$keys
  
  
  return(sstates)
}

## apply transformation to parameter sets ---------------------------------------

#' transform parametersets from the profiles for path plotting  
#' 
#' @description while using non-trivial steady states, parameters can couple trough the steady states. This function applies the transformation to the parameter sets from the profiles, to account for that.
#' 
#' 
#' @param profs parframe with the profiles, as returned from the \link{dMod::profile} function
#' @param trafo parameter transformation for the steady states, as returned by `P(steadystateTrafo)`. Currently no ther formulation is supported.
#' @param rescale character, default `"lin"` (no rescaling). The rescaling of the transformed parameters to the model scale, can be `"lin"`, `"log"`, `"log10"` or `"log2"`.
#' 
#' @return `parframe` of the input `profs` with the added columns of `trafo` applied to the parameters.
#' 
#' @export
addTrafoForPaths <- function(
    profs,
    trafo,
    rescale = c("lin", "log", "log10", "log2")[1]
) {
  
  # check if 'rescale' is a valid input
  if (!(rescale %in% c("lin", "log", "log10", "log2"))) {
    stop("'rescale' must be one of 'lin', 'log', 'log10' or 'log2'")
  }
  
  # build data.frame from trafo applied row wise to the entries (i.e. each parameterset)
  tDF <- do.call(
    rbind,
    lapply(
      seq_len(nrow(profs)),
      function (i) {
        parset <- do.call(c,profs[i,9:ncol(profs)])
        tParset <- trafo(parset)
        flattened <- flatten(tParset[1])
        if (is.list(flattened)) {
          namedTParset <- do.call(c, flatten(tParset[1]))
        } else {
          namedTParset <- flattened
        }
        # perform rescaling if necessary
        if (rescale == "log") {
          return(log(namedTParset))
        } else if (rescale == "log10") {
          return(log10(namedTParset))
        } else if (rescale == "log2") {
          return(log2(namedTParset))
        }
        setNames(as.data.frame(t(namedTParset)), paste0(names(namedTParset),"_trafo"))
      }
    )
  )
  
  # cast original profs to data.frame for joining  
  profsDF <- as.data.frame(profs)
  
  
  # add the new columns of the transformed profs parameters
  profsCombined <- cbind(profsDF, tDF)
  
  # make a parframe out of the combined data.frame and return it
  profsCombinedPF <- parframe(profsCombined)
  attr(profsCombinedPF, "metanames") <- attr(profs, "metanames")
  attr(profsCombinedPF, "obj.attributes") <- attr(profs, "obj.attributes")
  return(profsCombinedPF)
}
