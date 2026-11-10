# ====================================================================== #
# run_coevolution.R
#
# Coevolutionary GMF analysis of the empirical mutualistic networks.
#
# The script simulates coevolution in each empirical network with fewer than
# max_species species. At each time step, species traits reach their
# coevolutionary steady state under the exact coevolutionary model, and the
# temporal average of exact log fitness gives the numerical GMF. The
# analytical approximation is computed for the same networks, so the two can
# be compared.
#
# The scenario with independent fluctuations is computed by
# run_random_scenario.R, once per number of partners.
#
# Inputs
#   src/empirical_networks/*.csv
#
# Outputs
#   results/coevolution_gmf.csv     one row per species
#   results/network_subset.csv      one row per analysed network
#
# Each network is seeded with seed + net_id, where net_id is its row index in
# results/network_subset.csv. main_fig3.R simulates again the coevolution of
# its two focal species from these seeds, so the row order of that table must
# stay fixed.
# ====================================================================== #


#---------------------------------------------#
#---------Loading packages and functions------#
#---------------------------------------------#

library(here)
library(parallel)

source(here("src", "aux_functions.R"))
source(here("src", "main_functions.R"))


#---------------------------------------------#
#------------Analysis parameters--------------#
#---------------------------------------------#

# Network set. The empirical files are plant-by-animal incidence matrices.
# Networks with max_species or more species after cleaning are skipped.
network_dir <- here("src", "empirical_networks")
max_species <- 300

# Model parameters (Methods). alpha controls how strongly the benefits of
# mutualisms decay with trait mismatches, and m is the proportion of the
# selective pressures that comes from mutualisms.
alpha <- 0.2
m <- 0.5

# Coevolutionary scenario. Species traits start at mean_theta. At each time
# step, species traits reach their coevolutionary steady state after at most
# n_burnin iterations of the coevolutionary model, or earlier when the mean
# absolute change in species traits falls below eps. All simulation settings
# are written into the output table.
nsteps <- 5000
n_burnin <- 300
eps <- 1e-3
phi <- 0.5
mean_theta <- 0
sd_theta <- 2.88

# Reproducibility and parallelization.
seed <- 42
n_cores <- detectCores(logical=FALSE)

out_results <- here("results", "coevolution_gmf.csv")
out_network_subset <- here("results", "network_subset.csv")


#---------------------------------------------#
#---------Selecting empirical networks--------#
#---------------------------------------------#

network_files <- list.files(network_dir, pattern="\\.csv$")

if (length(network_files) == 0) {
  stop("No CSV files found in ", network_dir)
}

# Sorting fixes the order of the analysis set, which fixes net_id and hence
# the random seed of every network.
network_files <- sort(network_files)

# Paths are kept relative to the project root so that the stored table can be
# read on any computer without pointing outside the repository.
network_paths <- file.path("src", "empirical_networks", network_files)

subset_list <- list()
included_paths <- character(0)

for (path in network_paths) {

  # The network is read here only to decide whether it belongs to the analysis
  # set. Each worker reads it again before running the model, which keeps every
  # parallel job self-contained.
  net <- read_network(file=path, max_species=max_species)

  if (!net$included) next

  included_paths <- c(included_paths, path)

  subset_list[[length(subset_list) + 1]] <- data.frame(
    network=net$network,
    network_file=net$network_file,
    n_plants=net$n_plants,
    n_animals=net$n_animals,
    n_species=net$n_species,
    n_interactions=sum(net$incidence),
    connectance=net$connectance,
    leading_eigenvalue=net$leading_eigenvalue,
    stringsAsFactors=FALSE
  )
}

if (length(included_paths) == 0) {
  stop("No empirical networks passed the max_species filter.")
}

network_subset <- do.call(rbind, subset_list)
rownames(network_subset) <- NULL

dir.create(here("results"), showWarnings=FALSE, recursive=TRUE)

utils::write.csv(network_subset, file=out_network_subset, row.names=FALSE)

cat("\nSelected empirical network set.\n")
cat("Networks included:", nrow(network_subset), "\n")
cat("Maximum species threshold:", max_species, "\n")
cat("Network subset table:", out_network_subset, "\n\n")


#---------------------------------------------#
#---------Parallel coevolutionary run---------#
#---------------------------------------------#

# Each worker receives one network index. That index also sets the random seed,
# which gives every network a distinct and reproducible sequence of
# environmental optima.
network_ids <- seq_along(included_paths)

n_cores <- min(n_cores, length(network_ids))

cat("Running coevolutionary GMF analysis with", n_cores, "cores.\n\n")

cl <- makeCluster(n_cores)

# Workers do not inherit the sourced functions, so each worker sources both
# files again. The project root is resolved once here with here() and
# exported, because workers build every path from project_root.
project_root <- here()

clusterExport(cl, varlist=c("project_root"), envir=environment())

invisible(clusterEvalQ(cl, {
  source(file.path(project_root, "src", "aux_functions.R"))
  source(file.path(project_root, "src", "main_functions.R"))
  NULL
}))

clusterExport(cl, varlist=c(
  "included_paths", "project_root",
  "alpha", "m", "max_species", "seed",
  "nsteps", "n_burnin", "eps", "phi",
  "mean_theta", "sd_theta"
), envir=environment())

results_list <- tryCatch({

  parLapply(cl, network_ids, function(net_id) {

    set.seed(seed + net_id)

    net <- read_network(file=file.path(project_root, included_paths[net_id]),
                        max_species=max_species)

    if (!net$included) return(data.frame())

    # The reader stores the path it was given, which is absolute inside a
    # worker. The relative form is restored so the table stays portable.
    net$network_file <- included_paths[net_id]

    out <- simulate_coevolution(
      net=net,
      alpha=alpha,
      m=m,
      nsteps=nsteps,
      n_burnin=n_burnin,
      eps=eps,
      phi=phi,
      mean_theta=mean_theta,
      sd_theta=sd_theta
    )

    return(out)
  })

}, finally={
  stopCluster(cl)
})


#---------------------------------------------#
#-------------Saving output-------------------#
#---------------------------------------------#

keep <- vapply(results_list, nrow, numeric(1)) > 0
results_list <- results_list[keep]

if (length(results_list) == 0) {
  stop("No model results were returned.")
}

results <- do.call(rbind, results_list)
rownames(results) <- NULL

utils::write.csv(results, file=out_results, row.names=FALSE)


#---------------------------------------------#
#-------------Console diagnostics-------------#
#---------------------------------------------#

cat("\nCoevolutionary GMF analysis finished.\n")
cat("Rows:", nrow(results), "\n")
cat("Output:", out_results, "\n")
cat("alpha:", alpha, "\n")
cat("m:", m, "\n")
cat("nsteps:", nsteps, "\n")
cat("n_burnin:", n_burnin, "\n")
cat("eps:", eps, "\n")
cat("phi:", phi, "\n")
cat("sd_theta:", sd_theta, "\n")
