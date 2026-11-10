# ====================================================================== #
# run_random_scenario.R
#
# Scenario with independent fluctuations for the empirical mutualistic
# networks.
#
# In this scenario, the trait of the focal species remains fixed across
# generations, while partner traits and the environmental optimum fluctuate
# independently. Log fitness then depends on the network only through the
# number of partners of the focal species. Therefore, the analysis runs once
# for each number of partners observed across the empirical networks, rather
# than once per species.
#
# Two approximations are reported next to the exact simulation. The first
# evaluates the small-alpha approximation of log fitness on the partner
# traits and environmental optima drawn in the simulation. The second is the
# closed-form expectation of this approximation (equation 3 of the main
# text), which depends only on the number of partners and on the
# distributions of traits and environmental optima.
#
# Inputs
#   src/empirical_networks/*.csv     read only to collect numbers of partners
#
# Outputs
#   results/random_gmf_by_degree.csv   one row per number of partners
# ====================================================================== #


#---------------------------------------------#
#---------Loading packages and functions------#
#---------------------------------------------#

library(here)

source(here("src", "aux_functions.R"))
source(here("src", "main_functions.R"))


#---------------------------------------------#
#------------Analysis parameters--------------#
#---------------------------------------------#

network_dir <- here("src", "empirical_networks")
max_species <- 300

alpha <- 0.2
m <- 0.5

# Parameters of the scenario with independent fluctuations, shared with
# main_fig2.R.
#
# The trait of the focal species remains fixed one unit away from the mean of
# partner traits and environmental optima. Consequently, the population pays
# a persistent cost for its distance to the average position of the moving
# peak. Partner traits also vary more than the environmental optimum, which
# makes the portfolio effect of generalists visible across numbers of
# partners.
mean_zi <- 1
sd_zi <- 0
mean_zj <- 0
sd_zj <- 1.5
mean_theta <- 0
sd_theta <- 1

# Number of generations used to compute the temporal average of log fitness
# for each number of partners.
nsteps <- 10000

seed <- 42

out_results <- here("results", "random_gmf_by_degree.csv")


#---------------------------------------------#
#-----Collecting the numbers of partners------#
#---------------------------------------------#

# The networks are read only to learn which numbers of partners occur in the
# empirical networks, and how many species have each of them. Nothing else
# about network structure enters this scenario.

network_files <- sort(list.files(network_dir, pattern="\\.csv$"))

if (length(network_files) == 0) {
  stop("No CSV files found in ", network_dir)
}

observed_degrees <- lapply(network_files, function(file) {

  net <- read_network(file=file.path("src", "empirical_networks", file),
                      max_species=max_species)

  if (!net$included) return(NULL)

  return(net$species$degree)
})

degree_counts <- table(unlist(observed_degrees))
observed_degree <- sort(as.numeric(names(degree_counts)))

cat("\nDistinct degrees observed:", length(observed_degree), "\n")
cat("Degree range:", min(observed_degree), "-", max(observed_degree), "\n")
cat("Species represented:", sum(degree_counts), "\n\n")


#---------------------------------------------#
#-----Simulating each number of partners------#
#---------------------------------------------#

set.seed(seed)

degree_results <- vector("list", length(observed_degree))

for (k in seq_along(observed_degree)) {

  degree_k <- observed_degree[k]

  # Exact simulation, together with the small-alpha approximation evaluated
  # on the same partner traits and environmental optima.
  simulated <- simulate_random_gmf(
    degree=degree_k,
    alpha=alpha,
    m=m,
    nsteps=nsteps,
    mean_zi=mean_zi,
    sd_zi=sd_zi,
    mean_zj=mean_zj,
    sd_zj=sd_zj,
    mean_theta=mean_theta,
    sd_theta=sd_theta
  )

  # Closed-form expectation for the same number of partners.
  expected <- compute_random_gmf_terms(
    degree=degree_k,
    alpha=alpha,
    m=m,
    mean_zi=mean_zi,
    sd_zi=sd_zi,
    mean_zj=mean_zj,
    sd_zj=sd_zj,
    mean_theta=mean_theta,
    sd_theta=sd_theta
  )

  degree_results[[k]] <- data.frame(
    degree=degree_k,
    n_species=as.integer(degree_counts[as.character(degree_k)]),
    alpha=alpha,
    m=m,
    ln_gmf_exact=simulated$ln_gmf_exact,
    gmf_exact=simulated$gmf_exact,
    ln_gmf_approx_sim=simulated$ln_gmf_approx_sim,
    gmf_approx_sim=simulated$gmf_approx_sim,
    peak_height_sim=simulated$peak_height_sim,
    peak_distance_sim=simulated$peak_distance_sim,
    var_selective_target_sim=simulated$var_selective_target_sim,
    generalist_benefit_sim=simulated$generalist_benefit_sim,
    ln_gmf_approx=expected$ln_gmf_approx,
    gmf_approx=expected$gmf_approx,
    peak_height=expected$peak_height,
    peak_distance=expected$peak_distance,
    var_selective_target=expected$var_selective_target,
    generalist_benefit=expected$generalist_benefit,
    nsteps=nsteps,
    mean_zi=mean_zi,
    sd_zi=sd_zi,
    mean_zj=mean_zj,
    sd_zj=sd_zj,
    mean_theta=mean_theta,
    sd_theta=sd_theta,
    stringsAsFactors=FALSE
  )
}

results <- do.call(rbind, degree_results)
rownames(results) <- NULL

dir.create(here("results"), showWarnings=FALSE, recursive=TRUE)

utils::write.csv(results, file=out_results, row.names=FALSE)


#---------------------------------------------#
#-------------Console diagnostics-------------#
#---------------------------------------------#

cat("Random-scenario analysis finished.\n")
cat("Rows:", nrow(results), "\n")
cat("Output:", out_results, "\n")
cat("alpha:", alpha, "\n")
cat("m:", m, "\n")
cat("nsteps:", nsteps, "\n")
