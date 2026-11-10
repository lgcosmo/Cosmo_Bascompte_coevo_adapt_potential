# ====================================================================== #
# main_functions.R
#
# Functions that run the numerical simulations of the exact fitness function
# and coevolutionary model, which test the analytical approximations of the
# geometric mean fitness (GMF). They use the functions in aux_functions.R,
# which must be sourced first.
#
# The simulations follow the two scenarios of the main text:
#
#   1. Independent fluctuations ("random" in the function names). Partner
#      traits and environmental optima fluctuate independently across
#      generations. This scenario depends on the network only through the
#      number of partners of the focal species;
#
#   2. Coevolution. Environmental optima fluctuate across time steps, and
#      species traits reach their coevolutionary steady state within each
#      time step. This scenario depends on the whole network.
# ====================================================================== #


#---------------------------------------------#
#-------------Exact log fitness---------------#
#---------------------------------------------#

# Compute the exact log fitness of one species (log of equation 1 of the
# main text), given its trait z_i, the traits of its partners z_j and its
# environmental optimum theta_i:
#
#   ln w_i =
#     (m / alpha) ln[ sum_j (a_ij / k_i) exp{-alpha (z_j - z_i)^2} ]
#     - (1 - m)(theta_i - z_i)^2.

compute_ln_w <- function(z_i, z_j, theta_i, alpha, m) {

  if (length(z_j) == 0) return(NA_real_)

  sq_dist <- (z_j - z_i)^2

  # Compute log(mean(exp(-alpha * d^2))) with a stable log-sum-exp form.
  # This is algebraically identical to the expression in the equation, but
  # avoids numerical underflow when alpha or trait differences are large.
  log_kernel <- -alpha * sq_dist
  max_kernel <- max(log_kernel)
  log_S <- max_kernel + log(mean(exp(log_kernel - max_kernel)))

  ln_w <- (m / alpha) * log_S -
    (1 - m) * (theta_i - z_i)^2

  return(ln_w)
}


# Compute the exact log fitness of every species in a network, given the
# vectors of species traits z and environmental optima theta.

compute_ln_w_all <- function(A, z, theta, alpha, m) {

  N <- nrow(A)
  ln_w <- rep(NA_real_, N)

  for (i in seq_len(N)) {

    partners <- which(A[i, ] == 1)

    ln_w[i] <- compute_ln_w(z_i=z[i],
                            z_j=z[partners],
                            theta_i=theta[i],
                            alpha=alpha,
                            m=m)
  }

  return(ln_w)
}


#---------------------------------------------#
#-----Simulating independent fluctuations-----#
#---------------------------------------------#

# Simulate the scenario with independent fluctuations for a species with a
# given number of partners (degree). Partner traits and the environmental
# optimum are drawn at each of nsteps generations, and the trait of the
# focal species remains fixed when sd_zi = 0. Returns the exact log GMF and
# the temporal averages of the components of equation 2 of the main text,
# computed from the same draws.

simulate_random_gmf <- function(degree, alpha, m, nsteps,
                                mean_zi=0, sd_zi=0,
                                mean_zj=0, sd_zj=1,
                                mean_theta=0, sd_theta=1) {

  if (degree <= 0) {
    stop("degree must be positive for the random GMF simulation.")
  }

  ln_w <- rep(NA_real_, nsteps)
  ln_w_approx <- rep(NA_real_, nsteps)
  height <- rep(NA_real_, nsteps)
  distance <- rep(NA_real_, nsteps)
  target_var <- rep(NA_real_, nsteps)
  benefit <- rep(NA_real_, nsteps)

  for (t in seq_len(nsteps)) {

    z_i <- rnorm(1, mean_zi, sd_zi)
    z_j <- rnorm(degree, mean_zj, sd_zj)
    theta_i <- rnorm(1, mean_theta, sd_theta)

    # Exact log fitness from equation 1 of the main text.
    ln_w[t] <- compute_ln_w(z_i=z_i,
                            z_j=z_j,
                            theta_i=theta_i,
                            alpha=alpha,
                            m=m)

    # V_D is the variance in trait matching among partners at this
    # generation, the term that fuels the portfolio effect of generalists.
    D <- (z_j - z_i)^2
    D_bar <- mean(D)
    V_D <- mean(D^2) - D_bar^2

    benefit[t] <- (m * alpha / 2) * V_D

    # Equation 2 of the main text decomposes log fitness into the height of
    # the adaptive peak minus the squared distance to the peak. The adaptive
    # peak is the weighted average of partner traits and the environmental
    # optimum, and V_i is the variance of these selective targets.
    z_opt <- m * mean(z_j) + (1 - m) * theta_i

    target_var[t] <- m * mean(z_j^2) +
      (1 - m) * theta_i^2 -
      z_opt^2

    distance[t] <- (z_i - z_opt)^2
    height[t] <- -target_var[t] + benefit[t]

    ln_w_approx[t] <- height[t] - distance[t]
  }

  ln_gmf_exact <- mean(ln_w, na.rm=TRUE)
  ln_gmf_approx_sim <- mean(ln_w_approx, na.rm=TRUE)

  out <- data.frame(
    ln_gmf_exact=ln_gmf_exact,
    gmf_exact=exp(ln_gmf_exact),
    ln_gmf_approx_sim=ln_gmf_approx_sim,
    gmf_approx_sim=exp(ln_gmf_approx_sim),
    peak_height_sim=mean(height, na.rm=TRUE),
    peak_distance_sim=mean(distance, na.rm=TRUE),
    var_selective_target_sim=mean(target_var, na.rm=TRUE),
    generalist_benefit_sim=mean(benefit, na.rm=TRUE),
    stringsAsFactors=FALSE
  )

  return(out)
}


#---------------------------------------------#
#-------Coevolutionary steady-state model-----#
#---------------------------------------------#

# Iterate the coevolutionary model of the Methods until species traits reach
# their steady state under the environmental optima theta of one time step:
#
#   z_i(t+1) = z_i(t) + phi [
#       m sum_j q_ij(t)(z_j(t) - z_i(t))
#       + (1-m)(theta_i - z_i(t))
#   ]
#
# The contributions of partners, q_ij, follow the trait-matching rule and
# are updated at every iteration. Iterations stop after n_burnin steps, or
# when the mean absolute change in species traits falls below eps.

compute_z_eq <- function(A, z, theta, phi, m, alpha,
                         n_burnin, eps=1e-3) {

  stopifnot(all(rowSums(A) > 0))

  for (b in seq_len(n_burnin)) {

    # z_dif stores z_j - z_i for interacting species. Multiplying by A make sure that
    # non-interacting species do not contribute directly to coevolution (only indirectly through time).
    z_dif <- t(A * z) - A * z

    # Direct evolutionary effects, Q_ij = m q_ij, where the contributions of
    # partners, q_ij, follow the trait-matching rule and sum to one across
    # the partners of each species.
    Q_raw <- A * exp(-alpha * z_dif^2)
    Q_coev <- m * (Q_raw / rowSums(Q_raw))

    mut_sel <- rowSums(Q_coev * z_dif)
    env_sel <- (1 - m) * (theta - z)

    z_old <- z
    z <- z + phi * (mut_sel + env_sel)

    mean_change <- mean(abs(z - z_old))

    if (mean_change < eps) break
  }

  return(z)
}


#---------------------------------------------#
#-------Coevolutionary network simulation-----#
#---------------------------------------------#

# Compute the GMF of all species coevolving in one empirical network. At
# each of nsteps time steps, new environmental optima are drawn and species
# traits reach their coevolutionary steady state. The temporal average of
# exact log fitness gives the numerical GMF, which is returned together with
# the analytical approximation (equation 5 of the main text), with one row
# per species.

simulate_coevolution <- function(net, alpha, m,
                                 nsteps=5000,
                                 n_burnin=300,
                                 eps=1e-3,
                                 phi=0.5,
                                 mean_theta=0,
                                 sd_theta=2.88) {

  A <- net$A
  N <- nrow(A)

  # Species traits start at the mean of the environmental optima. At each
  # later time step, traits start from the steady state of the previous one.
  z <- rep(mean_theta, N)
  ln_w_mat <- matrix(NA_real_, nrow=nsteps, ncol=N)

  for (t in seq_len(nsteps)) {

    theta <- rnorm(N, mean_theta, sd_theta)

    z <- compute_z_eq(A=A,
                      z=z,
                      theta=theta,
                      phi=phi,
                      m=m,
                      alpha=alpha,
                      n_burnin=n_burnin,
                      eps=eps)

    ln_w_mat[t, ] <- compute_ln_w_all(A=A,
                                      z=z,
                                      theta=theta,
                                      alpha=alpha,
                                      m=m)
  }

  ln_gmf_exact <- colMeans(ln_w_mat, na.rm=TRUE)

  T_mat <- compute_T_matrix(A, m)

  approx <- compute_coevo_gmf_terms(A=A,
                                    T_mat=T_mat,
                                    sd_theta=sd_theta,
                                    alpha=alpha,
                                    m=m)

  out <- data.frame(
    network=net$network,
    network_file=net$network_file,
    type=net$species$type,
    species=net$species$species,
    degree=net$species$degree,
    n_plants=net$n_plants,
    n_animals=net$n_animals,
    n_species=net$n_species,
    connectance=net$connectance,
    leading_eigenvalue=net$leading_eigenvalue,
    alpha=alpha,
    m=m,
    ln_gmf_exact=ln_gmf_exact,
    gmf_exact=exp(ln_gmf_exact),
    ln_gmf_approx=approx$ln_gmf_approx,
    gmf_approx=approx$gmf_approx,
    peak_height=approx$peak_height,
    peak_distance=approx$peak_distance,
    var_selective_target=approx$var_selective_target,
    generalist_benefit=approx$generalist_benefit,
    nu_bar=approx$nu_bar,
    rho_sq_bar=approx$rho_sq_bar,
    r2_partner=approx$r2_partner,
    sigma_ii=approx$sigma_ii,
    sigma_jj=approx$sigma_jj,
    benefit_kernel=approx$benefit_kernel,
    indirect_network=approx$indirect_network,
    indirect_species=approx$indirect_species,
    indirect_partners=approx$indirect_partners,
    R_i=approx$R_i,
    R_bar=approx$R_bar,
    nsteps=nsteps,
    n_burnin=n_burnin,
    eps=eps,
    phi=phi,
    sd_theta=sd_theta,
    stringsAsFactors=FALSE
  )

  rownames(out) <- NULL

  return(out)
}
