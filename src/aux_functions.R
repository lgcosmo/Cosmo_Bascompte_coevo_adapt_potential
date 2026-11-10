# ====================================================================== #
# aux_functions.R
#
# Auxiliary functions for the analyses of geometric mean fitness (GMF) in
# empirical mutualistic networks. These functions read the empirical
# networks, compute the matrix of total evolutionary effects, and compute
# the analytical GMF of the two scenarios of the main text.
#
# The empirical networks are plant-by-animal incidence matrices stored as
# CSV files, with plants in rows, animals in columns and plant names in the
# first column. Any positive entry is treated as an interaction.
# ====================================================================== #


#---------------------------------------------#
#---------Small numerical utilities-----------#
#---------------------------------------------#

# Return the name of a network file without its extension, which is used as
# the network identifier.

get_network_name <- function(file) {

  network <- tools::file_path_sans_ext(basename(file))

  return(network)
}


# Compute the leading eigenvalue of a square adjacency matrix A. The number
# of pathways of length l connecting any two species grows as lambda^l, so
# the Methods use the leading eigenvalue to measure the number of pathways
# connecting species in each network. The real part is taken to discard the
# small imaginary components that eigen() can produce.

compute_leading_eigenvalue <- function(A) {

  if (nrow(A) == 0) return(NA_real_)

  eig_values <- eigen(A, only.values=TRUE)$values
  leading_eigenvalue <- max(Re(eig_values))

  return(leading_eigenvalue)
}


# Average a species-level vector x over the partners of each species,
# sum_j (a_ij / k_i) x_j, where A is the square binary adjacency matrix.

partner_mean <- function(A, x) {

  out <- as.vector(A %*% x) / rowSums(A)

  return(out)
}


#---------------------------------------------#
#-------Network reading and conversion--------#
#---------------------------------------------#

# Convert a plant-by-animal incidence matrix B into the square adjacency
# matrix A used by the model equations, with plants first and animals
# second:
#
#       plants    animals
#   [     0          B    ] plants
#   [   B^T          0    ] animals

incidence_to_adjacency <- function(incidence) {

  n_plants  <- nrow(incidence)
  n_animals <- ncol(incidence)
  n_species <- n_plants + n_animals

  # Start from an empty square matrix and fill the two off-diagonal
  # bipartite blocks. The diagonal blocks remain zero because species from
  # the same guild do not interact in the empirical incidence matrix.
  A <- matrix(0, nrow=n_species, ncol=n_species)

  plant_id  <- seq_len(n_plants)
  animal_id <- n_plants + seq_len(n_animals)

  A[plant_id, animal_id] <- incidence
  A[animal_id, plant_id] <- t(incidence)

  sp_names <- c(rownames(incidence), colnames(incidence))
  rownames(A) <- sp_names
  colnames(A) <- sp_names

  return(A)
}


# Read and clean one empirical network. The CSV files have plant names in
# the first column and animal names in the column headers. The matrix is
# binarized and species without interactions are dropped. Networks with
# max_species or more species are returned with included=FALSE and are not
# simulated. Returns a list with the incidence matrix, the square adjacency
# matrix, the species table and the descriptors of the network.

read_network <- function(file, max_species=300) {

  if (!file.exists(file)) {
    stop("Network file not found: ", file)
  }

  network <- get_network_name(file)
  network_file <- gsub("\\\\", "/", file)

  # check.names=FALSE preserves species names exactly as they appear in the
  # CSV headers. This is important because animal names are stored as column
  # names and may contain spaces or punctuation.
  raw_net <- utils::read.csv(file, check.names=FALSE,
                             stringsAsFactors=FALSE)

  if (ncol(raw_net) < 2) {
    stop("Network file must contain plant names plus at least one animal column: ",
         file)
  }

  plant_names  <- trimws(as.character(raw_net[[1]]))
  animal_names <- trimws(colnames(raw_net)[-1])

  # If a file contains duplicated or missing names, make.unique() gives each
  # species a stable identifier while preserving the original name whenever
  # possible. Empty names are replaced before make.unique() so they do not
  # collapse into indistinguishable blanks.
  plant_names[plant_names == "" | is.na(plant_names)] <- "unnamed_plant"
  animal_names[animal_names == "" | is.na(animal_names)] <- "unnamed_animal"

  plant_names  <- make.unique(plant_names)
  animal_names <- make.unique(animal_names)

  # Convert the interaction block to numeric values. Any missing or
  # non-numeric entry becomes zero. Then binarize: any positive value is an
  # interaction and all other values are no interaction.
  interaction_block <- as.matrix(raw_net[, -1, drop=FALSE])
  interaction_values <- suppressWarnings(as.numeric(interaction_block))
  interaction_values[is.na(interaction_values)] <- 0

  incidence <- matrix(interaction_values, nrow=nrow(raw_net),
                      ncol=ncol(raw_net) - 1)
  incidence <- ifelse(incidence > 0, 1, 0)

  rownames(incidence) <- plant_names
  colnames(incidence) <- animal_names

  # Empty species have no partners and cannot be evaluated by the fitness
  # equations, because the mutualistic term is undefined when k_i = 0.
  # Dropping them also makes the species count reflect the connected part of
  # the empirical incidence matrix.
  keep_plants  <- rowSums(incidence) > 0
  keep_animals <- colSums(incidence) > 0

  incidence <- incidence[keep_plants, keep_animals, drop=FALSE]

  n_plants  <- nrow(incidence)
  n_animals <- ncol(incidence)
  n_species <- n_plants + n_animals

  included <- n_plants > 0 && n_animals > 0 && n_species < max_species

  # Networks that are empty after cleaning, or too large for the requested
  # analysis, keep their metadata but leave the adjacency-related objects
  # empty. This lets wrapper functions report which networks were skipped.
  if (included) {

    A <- incidence_to_adjacency(incidence)

    plant_degree  <- rowSums(incidence)
    animal_degree <- colSums(incidence)
    degree <- c(plant_degree, animal_degree)

    species <- data.frame(
      type=c(rep("plant", n_plants), rep("animal", n_animals)),
      species=c(rownames(incidence), colnames(incidence)),
      degree=as.numeric(degree),
      stringsAsFactors=FALSE
    )

    # Connectance is computed on the original bipartite incidence matrix,
    # because possible interactions are plant-animal pairs, not all pairs in
    # the square adjacency matrix.
    connectance <- sum(incidence) / (n_plants * n_animals)

    leading_eigenvalue <- compute_leading_eigenvalue(A)

  } else {

    A <- NULL
    species <- data.frame()
    connectance <- NA_real_
    leading_eigenvalue <- NA_real_
  }

  return(list(
    network=network,
    network_file=network_file,
    included=included,
    incidence=incidence,
    A=A,
    species=species,
    n_plants=n_plants,
    n_animals=n_animals,
    n_species=n_species,
    connectance=connectance,
    leading_eigenvalue=leading_eigenvalue
  ))
}


#---------------------------------------------#
#------Coevolutionary analytical objects------#
#---------------------------------------------#

# Compute the matrix of total evolutionary effects, T, from the adjacency
# matrix A and the proportion of the selective pressures that comes from
# mutualisms, m (equation 4 of the main text):
#
#   z* = (I - Q)^(-1) Psi theta = T theta
#
# In the small-alpha approximation each partner contributes equally to
# selection, so Q_ij = m a_ij / k_i, and Psi has 1-m on the diagonal.

compute_T_matrix <- function(A, m) {

  N <- nrow(A)

  stopifnot(all(rowSums(A) > 0))

  # Direct evolutionary effects among mutualistic partners, Q_ij = m a_ij / k_i.
  Q <- m * (A / rowSums(A))

  Psi <- diag(1 - m, nrow=N, ncol=N)
  I   <- diag(1, nrow=N, ncol=N)

  T_mat <- solve(I - Q, Psi)

  rownames(T_mat) <- rownames(A)
  colnames(T_mat) <- colnames(A)

  return(T_mat)
}


# Compute the temporal covariance matrix of coevolved traits. Coevolved
# traits are linear combinations of the environmental optima, z* = T theta.
# If the optima are independent across species with standard deviation
# sd_theta, then:
#
#   Sigma = sd_theta^2 T T^T.

compute_trait_covariance <- function(T_mat, sd_theta) {

  Sigma <- sd_theta^2 * T_mat %*% t(T_mat)

  rownames(Sigma) <- rownames(T_mat)
  colnames(Sigma) <- rownames(T_mat)

  return(Sigma)
}


#---------------------------------------------#
#----T-matrix pathways and variance terms-----#
#---------------------------------------------#

# Compute the contribution of indirect evolutionary effects from the
# T-matrix. The off-diagonal entries of row i are the effects that species i
# receives from the other species. The entries of its mutualistic partners
# (a_ij = 1) are direct effects, and the remaining ones are indirect. The
# contribution of indirect effects of the network is F in the Methods:
#
#   F = sum_i sum_{j != i, a_ij = 0} T_ij / sum_i sum_{j != i} T_ij.
#
# Returns the contribution for each species (indirect_species), its average
# among the partners of each species (indirect_partners), and F
# (indirect_network).

compute_indirect_effects <- function(A, T_mat) {

  # The off-diagonal entries of row i carry every effect that species i
  # receives from the rest of the network. Zeroing the diagonal removes the
  # contribution of the focal species to itself.
  T_off <- unname(T_mat)
  diag(T_off) <- 0

  # Zeroing the entries of direct partners leaves only the effects that reach
  # the focal species through indirect pathways.
  T_indirect <- T_off
  T_indirect[A == 1] <- 0

  total_effects <- rowSums(T_off)
  indirect_effects <- rowSums(T_indirect)

  indirect_species <- indirect_effects / total_effects

  # Average contribution of indirect effects among the direct partners of
  # each focal species.
  indirect_partners <- partner_mean(A, indirect_species)

  # Network-level contribution of indirect effects, F in the Methods.
  indirect_network <- sum(indirect_effects) / sum(total_effects)

  out <- data.frame(
    indirect_species=indirect_species,
    indirect_partners=indirect_partners,
    indirect_network=indirect_network,
    stringsAsFactors=FALSE
  )

  return(out)
}


# Compute the Simpson diversity of the coevolutionary pathways reaching each
# species, S_i in the main text (stored as R_i), and its average across the
# partners of each species, S_bar (stored as R_bar):
#
#   S_i = 1 - sum_k T_ik^2.
#
# These terms give the conflict of evolutionary effects of equation 5 of the
# main text, E[V_i] = sd_theta^2 (S_i - m S_bar).

compute_R_terms <- function(A, T_mat) {

  R_i <- 1 - rowSums(T_mat^2)

  out <- data.frame(
    R_i=R_i,
    R_bar=partner_mean(A, R_i),
    stringsAsFactors=FALSE
  )

  return(out)
}


#---------------------------------------------#
#---Approximation: independent fluctuations---#
#---------------------------------------------#

# Compute the analytical GMF of a species with a given number of partners
# (degree) in the scenario with independent fluctuations. For small values
# of alpha, the log GMF is:
#
#   ln G(w_i) =
#     -m sd_zj^2
#     -(1-m) sd_theta^2
#     -sd_zi^2
#     -m(1-m)(mean_zj - mean_theta)^2
#     -(mean_zi - zbar_opt)^2
#     +(m alpha / 2) sigma_D^2 (1 - 1/k_i)
#
# With sd_zi = 0 this is equation 3 of the main text. With sd_zi > 0 the
# trait of the focal species also fluctuates independently. Returns the log
# GMF and its components: the height of the peak, the distance to the peak,
# the variance of selective targets and the portfolio effect of generalists.

compute_random_gmf_terms <- function(degree, alpha, m,
                                     mean_zi, sd_zi,
                                     mean_zj, sd_zj,
                                     mean_theta, sd_theta) {

  if (degree <= 0) {
    stop("degree must be positive for the random GMF approximation.")
  }

  zbar_opt <- m * mean_zj + (1 - m) * mean_theta

  var_selective_target <- m * sd_zj^2 -
    (m^2 * sd_zj^2 / degree) +
    m * (1 - m) * sd_theta^2 +
    m * (1 - m) * (mean_zj - mean_theta)^2

  # Variance of the squared trait differences, sigma_D^2 = 2 sigma_Dz^2 in
  # the SI. The focal-trait terms enter by conditioning on z_i, which is
  # shared by all partners within a generation, and then averaging over its
  # fluctuations.
  sigma_D2 <- 2 * sd_zj^4 +
    4 * sd_zj^2 * (sd_zi^2 + (mean_zj - mean_zi)^2)

  generalist_benefit <- (m * alpha / 2) * sigma_D2 * (1 - 1 / degree)

  peak_height <- -var_selective_target + generalist_benefit

  peak_distance <- sd_zi^2 +
    (m^2 * sd_zj^2 / degree) +
    (1 - m)^2 * sd_theta^2 +
    (mean_zi - zbar_opt)^2

  ln_gmf_approx <- peak_height - peak_distance

  out <- data.frame(
    ln_gmf_approx=ln_gmf_approx,
    gmf_approx=exp(ln_gmf_approx),
    peak_height=peak_height,
    peak_distance=peak_distance,
    var_selective_target=var_selective_target,
    generalist_benefit=generalist_benefit,
    stringsAsFactors=FALSE
  )

  return(out)
}


#---------------------------------------------#
#--Analytical approximation: coevolutionary---#
#---------------------------------------------#

# Compute the analytical GMF of coevolving species (equation 5 of the main
# text). Coevolving species sit at their adaptive peaks, so the log GMF is
# the temporal average of the height of the peak:
#
#   ln G(w_i) = -E[V_i] + m alpha (1 - 1/k_i)(nu_bar^2 - rho_sq_bar)
#             = -E[V_i] + m alpha (1 - 1/k_i) nu_bar^2 (1 - r_i^2).
#
# nu_bar is the average temporal variance of the trait differences between
# species i and its partners, and rho_sq_bar is the average squared
# covariance among these trait differences. Their ratio, r_i^2 (r2_partner),
# is the correlation among partners, that is, the temporal correlation
# between the squared trait differences of species i with its partners,
# averaged across all pairs of partners. All terms follow from the
# covariance matrix of coevolved traits, Sigma = sd_theta^2 T T^T
# (Supplementary Information, section S6). Returns one row per species with
# the log GMF, its components, the contribution of indirect effects and the
# diversity of coevolutionary pathways (R_i and R_bar).

compute_coevo_gmf_terms <- function(A, T_mat, sd_theta, alpha, m) {

  N <- nrow(A)

  stopifnot(all(rowSums(A) > 0))

  Sigma <- compute_trait_covariance(T_mat, sd_theta)
  Sigma_diag <- unname(diag(Sigma))

  # sigma_ii is the temporal variance of the focal trait induced by the T
  # matrix. sigma_jj is the m-weighted average partner trait variance in the
  # E[V_i] expression. Both depend only on the diagonal of Sigma, so the
  # function evaluates them for every species at once.
  sigma_ii <- Sigma_diag
  sigma_jj <- m * partner_mean(A, Sigma_diag)

  var_selective_target <- sigma_jj + (1 - m) * sd_theta^2 - sigma_ii

  # The remaining components compare each focal species with its own set of
  # partners, so they are computed one species at a time.
  nu_bar <- rep(NA_real_, N)
  rho_sq_bar <- rep(NA_real_, N)
  r2_partner <- rep(NA_real_, N)
  benefit_kernel <- rep(NA_real_, N)
  generalist_benefit <- rep(NA_real_, N)

  for (i in seq_len(N)) {

    partners <- which(A[i, ] == 1)
    ki <- length(partners)

    # nu_j is the temporal variance of the trait difference to partner j,
    # Var(z_j - z_i) = Sigma_jj + Sigma_ii - 2 Sigma_ij. The mean-field
    # approximation of the SI replaces each nu_j by their average, nu_bar.
    nu_j <- Sigma_diag[partners] + Sigma_diag[i] - 2 * Sigma[i, partners]
    nu_bar[i] <- mean(nu_j)

    # A species with a single partner has no variance in trait matching, so
    # the portfolio effect vanishes with the factor (1 - 1/k_i). The partner
    # correlation needs pairs of partners and stays NA.
    benefit_kernel[i] <- 0
    generalist_benefit[i] <- 0

    if (ki > 1) {

      # rho_vals holds the covariance between the trait differences to each
      # pair of distinct partners, rho_jj' = Cov(z_j - z_i, z_j' - z_i). The
      # covariance is symmetric, so averaging over the unordered pairs from
      # combn() equals the average over ordered pairs used in the SI.
      pairs <- utils::combn(partners, 2)

      rho_vals <- Sigma[cbind(pairs[1, ], pairs[2, ])] -
        Sigma[i, pairs[1, ]] -
        Sigma[i, pairs[2, ]] +
        Sigma_diag[i]

      # When samples covary, the expected sample variance subtracts the
      # average squared covariance (SI section S6). Dividing this average by
      # nu_bar^2 gives the correlation among partners, r_i^2 in the main
      # text, which is the temporal correlation between the squared trait
      # differences of the focal species with its partners, averaged across
      # all pairs of partners.
      rho_sq_bar[i] <- mean(rho_vals^2)
      r2_partner[i] <- rho_sq_bar[i] / nu_bar[i]^2
      benefit_kernel[i] <- nu_bar[i]^2 - rho_sq_bar[i]

      generalist_benefit[i] <- m * alpha * (1 - 1 / ki) * benefit_kernel[i]
    }
  }

  # Coevolving species track their moving adaptive peak. The distance to the
  # peak is therefore zero, and log GMF equals the average height of the peak.
  peak_height <- -var_selective_target + generalist_benefit
  ln_gmf_approx <- peak_height

  indirect <- compute_indirect_effects(A=A, T_mat=T_mat)
  R_terms <- compute_R_terms(A=A, T_mat=T_mat)

  out <- data.frame(
    ln_gmf_approx=ln_gmf_approx,
    gmf_approx=exp(ln_gmf_approx),
    peak_height=peak_height,
    peak_distance=rep(0, N),
    var_selective_target=var_selective_target,
    generalist_benefit=generalist_benefit,
    nu_bar=nu_bar,
    rho_sq_bar=rho_sq_bar,
    r2_partner=r2_partner,
    sigma_ii=sigma_ii,
    sigma_jj=sigma_jj,
    benefit_kernel=benefit_kernel,
    indirect_network=indirect$indirect_network,
    indirect_species=indirect$indirect_species,
    indirect_partners=indirect$indirect_partners,
    R_i=R_terms$R_i,
    R_bar=R_terms$R_bar,
    stringsAsFactors=FALSE
  )

  return(out)
}


#---------------------------------------------#
#-----------Figure layout helpers-------------#
#---------------------------------------------#

# Place the species of one guild along the vertical axis of the bipartite
# network insets of Figures 3 and 4. Highlighted species are moved toward
# the centre, and the remaining species fill the outer positions in
# decreasing order of their number of partners.

place_guild_nodes <- function(names, degree, highlighted,
                              y_min=0.06, y_max=0.94) {

  n <- length(names)
  slots <- seq(y_min, y_max, length.out=n)

  highlighted <- unique(highlighted[highlighted %in% names])

  center_order <- order(abs(slots - 0.50))
  central_slots <- sort(center_order[seq_along(highlighted)])

  remaining <- setdiff(names, highlighted)
  remaining <- remaining[order(degree[match(remaining, names)],
                               decreasing=TRUE)]

  ordered_names <- rep(NA_character_, n)

  if (length(highlighted) > 0) {
    ordered_names[central_slots] <- highlighted
  }

  ordered_names[is.na(ordered_names)] <- remaining

  out <- data.frame(species=ordered_names,
                    y=slots,
                    stringsAsFactors=FALSE)

  return(out)
}
