# ====================================================================== #
# main_fig3.R
#
# Figure 3 of the main text.
#
# This script builds Figure 3, which shows how coevolution reshapes the
# geometric mean fitness (GMF) of species through the height of their
# adaptive peaks. Panels a and b compare two focal plants with five
# mutualistic partners that achieve contrasting GMF. Panel c relates GMF to
# the correlation among partners, r_i^2, for species with two to six partners.
#
# Because coevolving species track their moving adaptive peaks, the
# landscapes of panels a-b are centered on the optimum trait value of the
# focal species at each time step. Centering improves visualization. It
# removes the horizontal movement of the peak and leaves the temporal
# variation in peak height, which is the quantity that separates the two
# focal species.
#
# Inputs:
#   results/coevolution_gmf.csv
#   results/network_subset.csv
#   src/empirical_networks/M_PL_061_46.csv
#   src/empirical_networks/M_PL_068.csv
#
# Outputs:
#   figures/main_fig3.svg
#   figures/main_fig3.png
# ====================================================================== #


#---------------------------------------------#
#---------Loading packages and functions------#
#---------------------------------------------#

library(here)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)
library(grid)
library(svglite)

source(here("src", "aux_functions.R"))
source(here("src", "main_functions.R"))

dir.create(here("figures"), showWarnings=FALSE, recursive=TRUE)


#---------------------------------------------#
#--------Selected species and parameters------#
#---------------------------------------------#

# The two focal plants have the same number of partners and coevolve within
# pollination networks of similar size, yet their numerical GMF differs
# strongly.
results_file <- here("results", "coevolution_gmf.csv")
network_subset_file <- here("results", "network_subset.csv")

selected_species <- data.frame(
  panel=c("a", "b"),
  network=c("M_PL_061_46", "M_PL_068"),
  focal_species=c("Soulamea terminaloides", "Macleania rupestris"),
  panel_title=c("Higher GMF\n(five partners)", "Lower GMF\n(five partners)"),
  color=c("#0072B2", "#D55E00"),
  stringsAsFactors=FALSE
)

# These values are the ones used by run_coevolution.R to produce
# results_file. That script assigned each network the seed seed + net_id,
# where net_id is the row index of the network in the network subset table.
alpha <- 0.2
m <- 0.5
phi <- 0.5
coevo_nsteps <- 5000
coevo_n_burnin <- 300
coevo_eps <- 1e-3
coevo_mean_theta <- 0
coevo_sd_theta <- 2.88
seed <- 42

# The panels display the fitness landscapes of an evenly spaced subset of
# time steps. Peak-height marks and numerical GMF use every time step.
n_curves_shown <- 180
trait_offset_seq <- seq(-3.55, 3.55, length.out=321)

# All averages are computed in log space and exponentiated afterwards, so the
# plotted points and curves are geometric means.

#---------------------------------------------#
#-----------Graphical parameters--------------#
#---------------------------------------------#

# Axis titles, panel tags and annotations are enlarged relative to the base
# size so that, once the figure is scaled to the page, its text matches the
# single-panel Figure 1. Tick labels keep the base size.
theme_fig3 <- theme_pubr(base_size=14) +
  theme(legend.position="none",
        axis.title=element_text(color="gray10", size=18),
        axis.text=element_text(color="gray20"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid=element_blank(),
        plot.tag=element_text(face="bold", size=20),
        plot.tag.position=c(0.01, 0.985),
        # The top margin keeps the panel tags clear of the vertical axis
        # titles, which span most of the height of each panel.
        plot.margin=margin(16, 5.5, 5.5, 5.5))


#---------------------------------------------#
#----Fitness landscapes under coevolution-----#
#---------------------------------------------#

# Compute the exact log-fitness landscape of one focal species across the
# trait values in zi_seq. The traits of its partners, z_j, and its
# environmental optimum, theta_i, remain fixed at their values at the
# selected time step.

compute_exact_log_landscape <- function(z_j, theta_i, zi_seq, alpha, m) {

  sq_dist <- outer(zi_seq, z_j, "-")^2

  # Evaluate the mutualistic term with a stable log-sum-exp form, vectorized
  # over the possible trait values.
  log_kernel <- -alpha * sq_dist
  max_kernel <- apply(log_kernel, 1, max)
  log_S <- max_kernel +
    log(rowMeans(exp(log_kernel - max_kernel)))

  log_fitness <- (m / alpha) * log_S -
    (1 - m) * (theta_i - zi_seq)^2

  return(log_fitness)
}


# Simulate again the coevolution of one focal species, with the seed that
# run_coevolution.R gave to its network, and record its fitness landscape at
# each time step. Each landscape is centered on the coevolved trait of the
# focal species, so its value at zero is the height of the peak at that time
# step. The temporal average of these peak heights is the numerical log GMF.
# Returns the centered landscapes, their geometric mean fitness landscape,
# the peak height at every time step and a summary of the focal species.

simulate_focal_coevo_history <- function(net, focal_species, panel,
                                         net_seed, alpha, m, phi,
                                         nsteps, n_burnin, eps,
                                         mean_theta, sd_theta,
                                         trait_offset_seq,
                                         n_curves_shown) {

  focal_id <- which(net$species$species == focal_species &
                      net$species$type == "plant")

  if (length(focal_id) != 1) {
    stop("Focal plant not uniquely identified in network ", net$network,
         ": ", focal_species)
  }

  partner_id <- which(net$A[focal_id, ] == 1)

  if (length(partner_id) != 5) {
    stop("The focal species must have exactly five partners: ",
         focal_species)
  }

  shown_time <- unique(round(seq(1, nsteps, length.out=n_curves_shown)))
  curves <- vector("list", length(shown_time))
  curve_id <- 1

  realized_log_fitness <- rep(NA_real_, nsteps)
  realized_trait <- rep(NA_real_, nsteps)
  centered_log_landscape_sum <- rep(0, length(trait_offset_seq))
  set.seed(net_seed)

  # Species traits start at the mean of the environmental optima. At each
  # later time step, traits start from the steady state of the previous one.
  z <- rep(mean_theta, nrow(net$A))

  for (t in seq_len(nsteps)) {

    theta <- rnorm(nrow(net$A), mean_theta, sd_theta)

    z <- compute_z_eq(A=net$A,
                      z=z,
                      theta=theta,
                      phi=phi,
                      m=m,
                      alpha=alpha,
                      n_burnin=n_burnin,
                      eps=eps)

    realized_log_fitness[t] <- compute_ln_w(
      z_i=z[focal_id],
      z_j=z[partner_id],
      theta_i=theta[focal_id],
      alpha=alpha,
      m=m
    )
    realized_trait[t] <- z[focal_id]

    # Evaluate the fitness landscape of the focal species centered on its
    # coevolved trait, delta = x - z_i^*(t), holding partner traits at their
    # steady state and theta_i at its value at this time step. At delta = 0,
    # the focal species has its coevolved trait.
    log_landscape <- compute_exact_log_landscape(
      z_j=z[partner_id],
      theta_i=theta[focal_id],
      zi_seq=z[focal_id] + trait_offset_seq,
      alpha=alpha,
      m=m
    )

    # Average the centered landscapes in log scale. Exponentiating this
    # temporal average later gives the geometric mean fitness landscape
    # displayed in the panel.
    centered_log_landscape_sum <- centered_log_landscape_sum + log_landscape

    # Keep the landscapes of an evenly spaced subset of time steps for the
    # pale curves. Every time step still contributes one peak height to the
    # numerical GMF.
    if (t %in% shown_time) {
      curves[[curve_id]] <- data.frame(
        panel=panel,
        time=t,
        trait_value=trait_offset_seq,
        log_fitness=log_landscape,
        stringsAsFactors=FALSE
      )

      curve_id <- curve_id + 1
    }
  }

  curves <- do.call(rbind, curves)
  mean_trait <- mean(realized_trait)

  # Every time step contributes one peak height. The temporal average of
  # these log peak heights is the numerical log GMF.
  peaks <- data.frame(
    panel=panel,
    time=seq_len(nsteps),
    trait_value=rep(0, nsteps),
    coevolved_trait=realized_trait,
    log_fitness=realized_log_fitness,
    stringsAsFactors=FALSE
  )

  avg_centered <- data.frame(
    panel=panel,
    trait_value=trait_offset_seq,
    log_fitness=centered_log_landscape_sum / nsteps,
    model="Centered GM landscape",
    stringsAsFactors=FALSE
  )

  summary <- data.frame(
    panel=panel,
    network=net$network,
    species=focal_species,
    degree=length(partner_id),
    mean_trait=mean_trait,
    log_gmf_exact=mean(realized_log_fitness),
    gmf_exact=exp(mean(realized_log_fitness)),
    stringsAsFactors=FALSE
  )

  return(list(curves=curves,
              avg_centered=avg_centered,
              peaks=peaks,
              summary=summary))
}


#---------------------------------------------#
#--------Bipartite network inset helper-------#
#---------------------------------------------#

# First and last vertical slot of the network inset, passed to
# place_guild_nodes() in aux_functions.R.
inset_y_min <- 0.06
inset_y_max <- 0.94


# Build the bipartite network inset of one focal plant. The focal plant, its
# five pollinators and their interactions carry the panel color, and all
# other species and interactions remain gray.

plot_network_context_inset <- function(net, focal_species, panel_color) {

  incidence <- net$incidence
  focal_row <- which(rownames(incidence) == focal_species)
  partner_col <- which(incidence[focal_row, ] == 1)
  partner_names <- colnames(incidence)[partner_col]

  plant_pos <- place_guild_nodes(
    names=rownames(incidence),
    degree=rowSums(incidence),
    highlighted=focal_species,
    y_min=inset_y_min,
    y_max=inset_y_max
  )

  animal_pos <- place_guild_nodes(
    names=colnames(incidence),
    degree=colSums(incidence),
    highlighted=partner_names,
    y_min=inset_y_min,
    y_max=inset_y_max
  )

  plant_pos$x <- 0
  animal_pos$x <- 1
  plant_pos$guild <- "plant"
  animal_pos$guild <- "animal"

  edge_id <- which(incidence == 1, arr.ind=TRUE)

  edges <- data.frame(
    plant=rownames(incidence)[edge_id[, 1]],
    animal=colnames(incidence)[edge_id[, 2]],
    stringsAsFactors=FALSE
  )

  edges$x <- 0
  edges$xend <- 1
  edges$y <- plant_pos$y[match(edges$plant, plant_pos$species)]
  edges$yend <- animal_pos$y[match(edges$animal, animal_pos$species)]
  edges$edge_class <- ifelse(edges$plant == focal_species,
                             "direct", "background")

  plant_pos$node_class <- ifelse(plant_pos$species == focal_species,
                                 "focal", "background")
  animal_pos$node_class <- ifelse(animal_pos$species %in% partner_names,
                                  "partner", "background")

  background_edges <- edges[edges$edge_class == "background", ]
  direct_edges <- edges[edges$edge_class == "direct", ]

  p <- ggplot() +
    geom_segment(data=background_edges,
                 aes(x=x, y=y, xend=xend, yend=yend),
                 color="gray82", linewidth=0.38, alpha=0.45) +
    geom_segment(data=direct_edges,
                 aes(x=x, y=y, xend=xend, yend=yend),
                 color=panel_color, linewidth=0.95, alpha=0.90) +
    geom_point(data=plant_pos[plant_pos$node_class == "background", ],
               aes(x=x, y=y), color="gray80", size=1.8) +
    geom_point(data=animal_pos[animal_pos$node_class == "background", ],
               aes(x=x, y=y), color="gray80", size=1.8) +
    geom_point(data=animal_pos[animal_pos$node_class == "partner", ],
               aes(x=x, y=y), color=panel_color, size=2.9) +
    geom_point(data=plant_pos[plant_pos$node_class == "focal", ],
               aes(x=x, y=y), color=panel_color, size=4.2) +
    annotate("text", x=0, y=1.09, label="Plants",
             size=4.0, color="gray35") +
    annotate("text", x=1, y=1.09, label="Pollinators",
             size=4.0, color="gray35") +
    coord_cartesian(xlim=c(-0.24, 1.24), ylim=c(0, 1.15),
                    clip="off") +
    theme_void(base_size=8) +
    theme(plot.background=element_rect(fill="white", color=NA),
          panel.background=element_rect(fill="white", color=NA),
          plot.margin=margin(2, 2, 2, 2))

  return(p)
}


#---------------------------------------------#
#-------Simulating the two focal species------#
#---------------------------------------------#

if (!file.exists(results_file) || !file.exists(network_subset_file)) {
  stop("Figure 3 requires the empirical result and network-subset tables.")
}

empirical_results <- utils::read.csv(results_file, stringsAsFactors=FALSE)
network_subset <- utils::read.csv(network_subset_file,
                                  stringsAsFactors=FALSE)
network_subset$net_id <- seq_len(nrow(network_subset))

selected_species <- selected_species |>
  left_join(network_subset[, c("network", "network_file", "net_id")],
            by="network") |>
  left_join(empirical_results[, c("network", "species", "degree",
                                   "ln_gmf_exact", "gmf_exact",
                                   "ln_gmf_approx")],
            by=c("network", "focal_species"="species"))

if (any(is.na(selected_species$net_id)) ||
    any(is.na(selected_species$gmf_exact))) {
  stop("Could not match the selected Figure 3 species to the empirical results.")
}

if (any(selected_species$degree != 5)) {
  stop("Both selected Figure 3 focal species must have degree five.")
}

# Simulate again the coevolution of each focal species with the
# network-specific seed used by run_coevolution.R.
simulated <- lapply(seq_len(nrow(selected_species)), function(i) {

  selected_i <- selected_species[i, ]

  net_i <- read_network(file=here("src", "empirical_networks",
                                  paste0(selected_i$network, ".csv")),
                        max_species=300)

  simulation_i <- simulate_focal_coevo_history(
    net=net_i,
    focal_species=selected_i$focal_species,
    panel=selected_i$panel,
    net_seed=seed + selected_i$net_id,
    alpha=alpha,
    m=m,
    phi=phi,
    nsteps=coevo_nsteps,
    n_burnin=coevo_n_burnin,
    eps=coevo_eps,
    mean_theta=coevo_mean_theta,
    sd_theta=coevo_sd_theta,
    trait_offset_seq=trait_offset_seq,
    n_curves_shown=n_curves_shown
  )

  # Verify that each simulated GMF matches its stored value before using the
  # landscapes in the figure.
  if (abs(simulation_i$summary$gmf_exact -
          selected_i$gmf_exact) > 1e-10) {
    stop("Reconstructed GMF does not match the empirical result for ",
         selected_i$focal_species)
  }

  inset_i <- plot_network_context_inset(
    net=net_i,
    focal_species=selected_i$focal_species,
    panel_color=selected_i$color
  )

  return(list(simulation=simulation_i, inset=inset_i))
})


#---------------------------------------------#
#-----Panels a-b: GMF of two focal plants-----#
#---------------------------------------------#

# Common limits keep the two panels directly comparable. The horizontal axis
# is centered on the peak of the focal species at every time step.
plot_x_limits <- c(-3.42, 3.42)
plot_y_limits <- c(0, 1.05)

# Factory for the panels on the height of adaptive peaks under coevolution.
# Pale curves depict the exact fitness landscapes centered on the coevolved
# peak of the focal species. Short marks depict the height of the peak at
# every time step. The solid curve and filled point depict the numerical
# geometric mean fitness landscape and GMF, and the dashed curve and open
# point their analytical counterparts.
make_coevo_peak_height_panel <- function(selected_i, simulated_i,
                                         show_y_title=TRUE,
                                         show_legend=TRUE) {

  curves_i <- simulated_i$simulation$curves
  avg_centered_i <- simulated_i$simulation$avg_centered
  peaks_i <- simulated_i$simulation$peaks
  summary_i <- simulated_i$simulation$summary

  # Short horizontal marks depict the height of the peak at each time step.
  peak_half_width <- 0.005 * diff(plot_x_limits)
  peaks_i$x <- peaks_i$trait_value - peak_half_width
  peaks_i$xend <- peaks_i$trait_value + peak_half_width
  peaks_i$fitness <- exp(peaks_i$log_fitness)

  # The analytical centered landscape is ln G(delta) = ln G(0) - delta^2,
  # with ln G(0) the analytical GMF. The numerical landscape averages the
  # centered exact landscapes pointwise in log scale.
  avg_centered_i$model <- "Numerical landscape"

  analytical_centered_i <- data.frame(
    panel=selected_i$panel,
    trait_value=avg_centered_i$trait_value,
    log_fitness=selected_i$ln_gmf_approx -
      avg_centered_i$trait_value^2,
    model="Analytical landscape",
    stringsAsFactors=FALSE
  )

  centered_landscapes_i <- rbind(avg_centered_i,
                                 analytical_centered_i)

  gmf_i <- data.frame(
    trait_value=c(0, 0),
    fitness=c(summary_i$gmf_exact,
              exp(selected_i$ln_gmf_approx)),
    model=c("Numerical GMF", "Analytical GMF"),
    stringsAsFactors=FALSE
  )

  # The two-line panel title is centered on this point, which keeps both
  # lines aligned with each other in the upper left of the panel.
  title_x <- plot_x_limits[1] + 0.17 * diff(plot_x_limits)
  title_y <- plot_y_limits[2] - 0.06 * diff(plot_y_limits)

  p <- ggplot() +
    geom_line(data=curves_i,
              aes(x=trait_value,
                  y=exp(log_fitness),
                  group=time),
              color=selected_i$color, linewidth=0.24, alpha=0.040) +
    geom_line(data=centered_landscapes_i,
              aes(x=trait_value,
                  y=exp(log_fitness),
                  linetype=model),
              color=selected_i$color, linewidth=1.18) +
    geom_segment(data=peaks_i,
                 aes(x=x, xend=xend,
                     y=fitness, yend=fitness),
                 color=selected_i$color,
                 linewidth=0.36, alpha=0.10) +
    geom_point(data=gmf_i,
               aes(x=trait_value, y=fitness, shape=model),
               color=selected_i$color, fill=selected_i$color,
               size=3.2, stroke=0.85) +
    annotate("text", x=title_x, y=title_y,
             label=selected_i$panel_title,
             hjust=0.5, vjust=1, lineheight=0.9, size=5.2,
             color=selected_i$color) +
    # Legend keys are enlarged so that the point shapes and line types stay
    # legible once the figure is scaled to the page.
    scale_shape_manual(values=c("Numerical GMF"=21,
                                "Analytical GMF"=1),
                       name=NULL,
                       guide=if (show_legend) {
                         guide_legend(override.aes=list(size=4.6))
                       } else "none") +
    scale_linetype_manual(values=c("Numerical landscape"="solid",
                                   "Analytical landscape"="longdash"),
                          name=NULL,
                          guide=if (show_legend) {
                            guide_legend(override.aes=list(linewidth=1.5))
                          } else "none") +
    scale_x_continuous(breaks=c(-3, -1.5, 0, 1.5, 3),
                       expand=expansion(mult=c(0, 0))) +
    scale_y_continuous(breaks=seq(0, 1, by=0.2)) +
    coord_cartesian(xlim=plot_x_limits, ylim=plot_y_limits,
                    clip="on") +
    labs(x="Trait relative to coevolved optimum",
         y=ifelse(show_y_title, "Population mean fitness", ""),
         tag=selected_i$panel) +
    theme_fig3

  p <- p +
    inset_element(simulated_i$inset,
                  left=0.54, bottom=0.60,
                  right=0.99, top=0.985,
                  align_to="panel")

  return(p)
}

p3.1 <- make_coevo_peak_height_panel(selected_i=selected_species[1, ],
                                      simulated_i=simulated[[1]],
                                      show_y_title=TRUE,
                                      show_legend=TRUE)

p3.2 <- make_coevo_peak_height_panel(selected_i=selected_species[2, ],
                                      simulated_i=simulated[[2]],
                                      show_y_title=TRUE,
                                      show_legend=FALSE)


#---------------------------------------------#
#-----Panel c: correlation among partners-----#
#---------------------------------------------#

# r2_partner is the correlation among partners, r_i^2 in the main text. It
# is the temporal correlation between the squared trait differences of a
# species with its partners, averaged across all pairs of partners, and is
# stored in the results table. Panel c shows species with two to six
# partners.
panel_c_data <- empirical_results |>
  filter(degree > 1,
         degree <= 6,
         is.finite(r2_partner),
         is.finite(nu_bar),
         nu_bar > 0,
         is.finite(gmf_exact)) |>
  mutate(degree_group=factor(degree, levels=2:6))

pal_degree <- c("2"="#1B9E77",
                "3"="#7570B3",
                "4"="#E7298A",
                "5"="#66A61E",
                "6"="#A6761D")

p3.3 <- ggplot(panel_c_data,
               aes(x=r2_partner,
                   y=gmf_exact,
                   color=degree_group)) +
  geom_point(size=1.05, alpha=0.38) +
  scale_color_manual(values=pal_degree,
                     name="Number of partners",
                     guide=guide_legend(override.aes=list(size=4.2,
                                                          alpha=1))) +
  scale_x_continuous(breaks=seq(0, 0.20, by=0.05),
                     expand=expansion(mult=c(0, 0))) +
  scale_y_continuous(breaks=seq(0, 0.30, by=0.05),
                     expand=expansion(mult=c(0, 0))) +
  coord_cartesian(xlim=c(0, 0.22), ylim=c(0, 0.30),
                  clip="on") +
  labs(x="Correlation among partners",
       y="Geometric mean fitness",
       tag="c") +
  theme_fig3


#---------------------------------------------#
#-------------Assembling figure---------------#
#---------------------------------------------#

# Panel c keeps its own legend inside the panel, in its empty upper right
# corner. The shared legend at the bottom of the figure then holds only the
# landscape and GMF keys of panels a-b, which fit in a single row.
p3.3 <- p3.3 +
  theme(legend.position="inside",
        legend.position.inside=c(0.76, 0.90),
        legend.direction="horizontal",
        legend.text=element_text(size=13),
        legend.title=element_text(size=13.5),
        legend.key.width=unit(0.9, "cm"),
        legend.key.height=unit(0.7, "cm"),
        legend.spacing.x=unit(0.10, "cm"))

# The base theme hides legends, so the top row re-enables them before the
# patchwork collects the keys of panels a-b into the shared bottom legend.
# The collected keys take their text and key sizes from these panels, so
# the sizes are set here rather than on the figure.
fig3 <- ((p3.1 | p3.2) &
           theme(legend.position="bottom",
                 legend.text=element_text(size=13),
                 legend.key.width=unit(1.2, "cm"),
                 legend.key.height=unit(0.7, "cm"),
                 legend.spacing.x=unit(0.35, "cm"))) /
  (wrap_plots(p3.3) + plot_layout(guides="keep"))

fig3 <- fig3 +
  plot_layout(heights=c(1, 0.98), guides="collect") +
  plot_annotation(theme=theme(legend.position="bottom",
                              legend.direction="horizontal",
                              legend.box="horizontal",
                              legend.box.just="center",
                              legend.spacing.x=unit(0.35, "cm")))

print(fig3)

fig3_svg <- here("figures", "main_fig3.svg")
fig3_png <- here("figures", "main_fig3.png")

ggsave(fig3_svg, fig3,
       height=9.4, width=11.2, device=svglite::svglite)

ggsave(fig3_png, fig3,
       height=9.4, width=11.2, dpi=320)

cat("Saved Figure 3 to:\n")
cat("  ", fig3_svg, "\n", sep="")
cat("  ", fig3_png, "\n", sep="")
