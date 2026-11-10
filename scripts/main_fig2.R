# ====================================================================== #
# main_fig2.R
#
# Figure 2 of the main text.
#
# This script builds the scenario with independent fluctuations, in which the
# trait of the focal species remains fixed across generations, whereas
# partner traits and the environmental optimum fluctuate independently.
# Panels a-b show the exact fitness landscapes at each generation and the
# analytical geometric mean fitness landscape of a specialist and a
# generalist. Panel c compares the numerical GMF with its analytical
# approximation across the number of mutualistic partners.
#
# The model is computed on the log-fitness scale. For visualization only, the
# figure plots exp(c * log fitness), where c reduces compression close to zero.
#
# Inputs:
#   src/empirical_networks/*.csv
#   figures/butterfly_1.png, figures/butterfly_2.png
#   figures/plant_1.png, ..., figures/plant_5.png
#
# Outputs:
#   figures/main_fig2.svg
#   figures/main_fig2.png
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
library(png)
library(svglite)

source(here("src", "aux_functions.R"))
source(here("src", "main_functions.R"))

dir.create(here("figures"), showWarnings=FALSE, recursive=TRUE)


#---------------------------------------------#
#----------Simulation configuration-----------#
#---------------------------------------------#

# The fixed focal trait isolates how changes in partner traits and in the
# environmental optimum reshape fitness across generations. Partner traits
# vary enough to expose the portfolio effect of generalists across the number
# of partners.
alpha <- 0.2
m <- 0.5
nsteps <- 500
nsteps_degree <- 1000

# The trait of the focal species remains one unit away from the mean of
# partner traits and environmental optima. Consequently, the population pays
# a persistent cost for its distance to the average position of the peak.
fixed_zi <- 1
mean_zj <- 0
sd_zj <- 1.5
mean_theta <- 0
sd_theta <- 1

focal_degrees <- c(1, 5)
focal_labels <- c("Specialist",
                  "Generalist")

# This transformation changes only the display scale to improve
# visualization. All averages and analytical predictions are computed from
# log fitness before transformation.
fitness_log_multiplier <- 0.5

pal_species <- c("Specialist"="#0072B2",
                 "Generalist"="#D55E00")

pal_degree_points <- c("Specialist"=pal_species[1],
                       "Generalist"=pal_species[2],
                       "Other species"="gray72")

# Use common fixed bounds in both landscape panels so the distance between
# the trait of the focal species and the average peak can be compared
# directly.
landscape_x_lim <- c(-3.5, 3.5)
landscape_y_lim <- c(0, 1.07)
degree_y_lim <- c(0.22, 0.52)

zi_seq <- seq(landscape_x_lim[1], landscape_x_lim[2], length.out=260)


#---------------------------------------------#
#-----------Graphical parameters--------------#
#---------------------------------------------#

# Axis titles, panel tags and annotations are enlarged relative to the base
# size so that, once the figure is scaled to the page, its text matches the
# single-panel Figure 1. Tick labels keep the base size.
theme_fig2 <- theme_pubr(base_size=14) +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        axis.title=element_text(color="gray10", size=18),
        axis.text=element_text(color="gray20"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid=element_blank(),
        legend.text=element_text(size=13),
        legend.key.width=unit(1.1, "cm"),
        legend.spacing.x=unit(0.16, "cm"),
        plot.tag=element_text(face="bold", size=20),
        plot.tag.position=c(0.01, 0.98))


# Convert log fitness to the bounded scale displayed in the figure.
plot_fitness <- function(log_fitness, multiplier) {

  return(exp(multiplier * log_fitness))
}

butterfly_1_img <- readPNG(here("figures", "butterfly_1.png"))
butterfly_2_img <- readPNG(here("figures", "butterfly_2.png"))

# Load five distinct plant silhouettes. The specialist uses the first icon,
# whereas the generalist inset assigns one icon to each mutualistic partner.
plant_imgs <- lapply(seq_len(5), function(i) {
  readPNG(here("figures", paste0("plant_", i, ".png")))
})


# Place one ecological icon within data coordinates.
add_png_icon <- function(p, img, x, y, width=0.34, height=0.075) {

  p <- p +
    annotation_custom(
      grob=rasterGrob(img, interpolate=TRUE),
      xmin=x - width / 2,
      xmax=x + width / 2,
      ymin=y - height / 2,
      ymax=y + height / 2
    )

  return(p)
}


# Add a network inset that identifies the number of partners represented in
# each landscape panel. Plant icons denote distinct mutualistic partners.
add_network_inset <- function(p, species_name) {

  if (species_name == "Specialist") {

    butterfly_pos <- data.frame(x=-0.28, y=0.93)
    plant_pos <- data.frame(x=0.36, y=0.93)
    butterfly_img <- butterfly_1_img
    plant_icons <- plant_imgs[1]

  } else {

    butterfly_pos <- data.frame(x=0.00, y=0.87)
    plant_pos <- data.frame(x=c(-1.20, -0.60, 0.00, 0.60, 1.20),
                             y=c(0.96, 0.985, 0.995, 0.985, 0.96))
    butterfly_img <- butterfly_2_img
    plant_icons <- plant_imgs
  }

  network_edges <- data.frame(x=butterfly_pos$x,
                              y=butterfly_pos$y,
                              xend=plant_pos$x,
                              yend=plant_pos$y)

  p <- p +
    geom_segment(data=network_edges,
                 aes(x=x, y=y, xend=xend, yend=yend),
                 inherit.aes=FALSE,
                 color="gray70", linewidth=0.35, alpha=0.65)

  p <- add_png_icon(p, butterfly_img,
                    x=butterfly_pos$x, y=butterfly_pos$y,
                    width=0.58, height=0.124)

  for (i in seq_len(nrow(plant_pos))) {
    p <- add_png_icon(p, plant_icons[[i]],
                      x=plant_pos$x[i], y=plant_pos$y[i],
                      width=0.44, height=0.110)
  }

  return(p)
}


#---------------------------------------------#
#----Selective targets across generations-----#
#---------------------------------------------#

# Draw the partner traits and the environmental optimum of one focal species
# at each generation. These selective targets set the exact fitness landscape
# of the focal species at each generation.
simulate_random_target_timeseries <- function(degree, species_label, nsteps,
                                              mean_zj, sd_zj,
                                              mean_theta, sd_theta) {

  if (degree <= 0) {
    stop("degree must be positive for the random time series.")
  }

  out <- vector("list", nsteps)

  for (t in seq_len(nsteps)) {

    z_j <- rnorm(degree, mean_zj, sd_zj)
    theta_i <- rnorm(1, mean_theta, sd_theta)

    out[[t]] <- data.frame(
      time=t,
      species_id=species_label,
      degree=degree,
      z_j=I(list(z_j)),
      theta_i=theta_i,
      stringsAsFactors=FALSE
    )
  }

  out <- do.call(rbind, out)

  return(out)
}


# Compute the exact fitness landscape of the focal species at each
# generation, given the partner traits and the environmental optimum drawn at
# that generation.
build_exact_landscape_timeseries <- function(time_series, zi_seq, alpha, m) {

  curves <- vector("list", nrow(time_series))

  for (i in seq_len(nrow(time_series))) {

    row_i <- time_series[i, ]
    z_j_i <- row_i$z_j[[1]]

    log_fitness_i <- vapply(zi_seq, function(z) {

      compute_ln_w(z_i=z,
                   z_j=z_j_i,
                   theta_i=row_i$theta_i,
                   alpha=alpha,
                   m=m)

    }, numeric(1))

    curves[[i]] <- data.frame(
      time=row_i$time,
      species_id=row_i$species_id,
      degree=row_i$degree,
      zi=zi_seq,
      log_fitness=log_fitness_i,
      stringsAsFactors=FALSE
    )
  }

  curves <- do.call(rbind, curves)

  return(curves)
}


#---------------------------------------------#
#------Simulating the two focal species-------#
#---------------------------------------------#

set.seed(22011)

focal_series <- lapply(seq_along(focal_degrees), function(i) {

  simulate_random_target_timeseries(
    degree=focal_degrees[i],
    species_label=focal_labels[i],
    nsteps=nsteps,
    mean_zj=mean_zj,
    sd_zj=sd_zj,
    mean_theta=mean_theta,
    sd_theta=sd_theta
  )
})

focal_series <- do.call(rbind, focal_series)
focal_series$species_id <- factor(focal_series$species_id,
                                  levels=focal_labels)


#---------------------------------------------#
#----Loading empirical numbers of partners----#
#---------------------------------------------#

# Panel c uses the numbers of partners observed in the empirical networks. In
# the scenario with independent fluctuations, the network affects fitness
# only through the number of partners of the focal species, so only this
# number is read from each network.
network_dir <- here("src", "empirical_networks")
max_species <- 300
max_degree <- 10

network_files <- sort(list.files(network_dir, pattern="\\.csv$"))

if (length(network_files) == 0) {
  stop("Could not find the empirical networks in ", network_dir)
}

species_degrees <- lapply(network_files, function(file) {

  net <- read_network(file=file.path("src", "empirical_networks", file),
                      max_species=max_species)

  if (!net$included) return(NULL)

  data.frame(network=net$network,
             degree=net$species$degree,
             stringsAsFactors=FALSE)
})

species_degrees <- do.call(rbind, species_degrees)

empirical_b <- species_degrees |>
  filter(degree <= max_degree) |>
  mutate(degree_group=case_when(
    degree == 1 ~ "Specialist",
    degree == 5 ~ "Generalist",
    TRUE ~ "Other species"
  ))

empirical_b$degree_group <- factor(empirical_b$degree_group,
                                   levels=c("Other species",
                                            "Specialist",
                                            "Generalist"))

# In this scenario, species with the same number of partners experience the
# same distributions of partner traits and environmental optima. Each
# combination of network and number of partners receives one independent
# simulation, and the spread of points shows the variation expected over a
# finite number of generations.
degree_replicates <- empirical_b |>
  distinct(network, degree)

set.seed(22012)

degree_replicates_sim <- lapply(seq_len(nrow(degree_replicates)),
                                function(i) {

  row_i <- degree_replicates[i, ]

  sim_i <- simulate_random_gmf(
    degree=row_i$degree,
    alpha=alpha,
    m=m,
    nsteps=nsteps_degree,
    mean_zi=fixed_zi,
    sd_zi=0,
    mean_zj=mean_zj,
    sd_zj=sd_zj,
    mean_theta=mean_theta,
    sd_theta=sd_theta
  )

  data.frame(network=row_i$network,
             degree=row_i$degree,
             ln_gmf_exact=sim_i$ln_gmf_exact,
             stringsAsFactors=FALSE)
})

degree_replicates_sim <- do.call(rbind, degree_replicates_sim)

empirical_b <- empirical_b |>
  left_join(degree_replicates_sim,
            by=c("network", "degree"))


mean_by_degree <- empirical_b |>
  group_by(degree) |>
  summarise(ln_gmf_mean=mean(ln_gmf_exact, na.rm=TRUE),
            .groups="drop")

# Setting sd_zi to zero in the analytical approximation gives the expected
# GMF when the trait of the focal species remains fixed across generations.
analytical_b <- lapply(sort(unique(empirical_b$degree)), function(k) {

  approx <- compute_random_gmf_terms(
    degree=k,
    alpha=alpha,
    m=m,
    mean_zi=fixed_zi,
    sd_zi=0,
    mean_zj=mean_zj,
    sd_zj=sd_zj,
    mean_theta=mean_theta,
    sd_theta=sd_theta
  )

  data.frame(degree=k,
             ln_gmf_analytical=approx$ln_gmf_approx,
             stringsAsFactors=FALSE)
})

analytical_b <- do.call(rbind, analytical_b)


#---------------------------------------------#
#---Panel c: GMF and the number of partners---#
#---------------------------------------------#

# The black line depicts the average numerical GMF for each number of
# partners. The dashed line depicts the analytical approximation for the same
# fixed trait of the focal species.
p2.3 <- ggplot() +
  geom_point(data=empirical_b,
             aes(x=degree,
                 y=plot_fitness(ln_gmf_exact,
                                fitness_log_multiplier),
                 color=degree_group),
             alpha=0.32, size=1.25,
             position=position_jitter(width=0.12, height=0)) +
  geom_line(data=mean_by_degree,
            aes(x=degree,
                y=plot_fitness(ln_gmf_mean, fitness_log_multiplier)),
            color="black", linewidth=1.05) +
  geom_point(data=mean_by_degree,
             aes(x=degree,
                 y=plot_fitness(ln_gmf_mean, fitness_log_multiplier)),
             color="black", size=2.1) +
  geom_line(data=analytical_b,
            aes(x=degree,
                y=plot_fitness(ln_gmf_analytical,
                               fitness_log_multiplier)),
            color="gray25", linetype="longdash", linewidth=1.05) +
  geom_point(data=analytical_b,
             aes(x=degree,
                 y=plot_fitness(ln_gmf_analytical,
                                fitness_log_multiplier)),
             color="gray25", size=2.1) +
  annotate("text", x=7.5,
           y=plot_fitness(mean_by_degree$ln_gmf_mean[
             mean_by_degree$degree == 8
           ], fitness_log_multiplier) - 0.026,
           label="Simulation", hjust=0, size=4.4, color="black") +
  annotate("text", x=5.3,
           y=plot_fitness(analytical_b$ln_gmf_analytical[
             analytical_b$degree == 5
           ], fitness_log_multiplier) + 0.030,
           label="Analytical",
           hjust=0, size=4.4, color="gray25") +
  scale_color_manual(values=pal_degree_points) +
  scale_x_continuous(breaks=1:10, limits=c(0.65, 10.35)) +
  scale_y_continuous(limits=degree_y_lim,
                     breaks=seq(0.25, 0.50, by=0.05)) +
  labs(x="Number of mutualistic partners",
       y="Geometric mean fitness",
       tag="c") +
  theme_fig2 +
  theme(legend.position="none")


#---------------------------------------------#
#-----Panels a-b: analytical GMF landscapes---#
#---------------------------------------------#

landscape_exact <- build_exact_landscape_timeseries(
  time_series=focal_series,
  zi_seq=zi_seq,
  alpha=alpha,
  m=m
)

# The solid curve is the analytical geometric mean fitness landscape. To draw
# it, compute_random_gmf_terms() evaluates the closed-form approximation
# (equation 3 of the main text) at each possible trait value of the focal
# species.

avg_approx <- lapply(seq_along(focal_degrees), function(i) {

  landscape_i <- lapply(zi_seq, function(zi_candidate) {

    terms_i <- compute_random_gmf_terms(
      degree=focal_degrees[i],
      alpha=alpha,
      m=m,
      mean_zi=zi_candidate,
      sd_zi=0,
      mean_zj=mean_zj,
      sd_zj=sd_zj,
      mean_theta=mean_theta,
      sd_theta=sd_theta
    )

    data.frame(
      species_id=focal_labels[i],
      degree=focal_degrees[i],
      zi=zi_candidate,
      log_fitness=terms_i$ln_gmf_approx,
      stringsAsFactors=FALSE
    )
  })

  return(do.call(rbind, landscape_i))
})

avg_approx <- do.call(rbind, avg_approx)
avg_approx$species_id <- factor(avg_approx$species_id,
                                levels=focal_labels)

# Analytical GMF at the fixed focal trait. The same value appears at the
# corresponding number of partners in panel c.
summary_landscape <- lapply(seq_along(focal_degrees), function(i) {

  terms_i <- compute_random_gmf_terms(
    degree=focal_degrees[i],
    alpha=alpha,
    m=m,
    mean_zi=fixed_zi,
    sd_zi=0,
    mean_zj=mean_zj,
    sd_zj=sd_zj,
    mean_theta=mean_theta,
    sd_theta=sd_theta
  )

  data.frame(
    species_id=focal_labels[i],
    degree=focal_degrees[i],
    zi_focal=fixed_zi,
    log_fitness_focal=terms_i$ln_gmf_approx,
    stringsAsFactors=FALSE
  )
})

summary_landscape <- do.call(rbind, summary_landscape)
summary_landscape$species_id <- factor(summary_landscape$species_id,
                                       levels=focal_labels)

# Plot one landscape panel. Faint curves depict the exact fitness landscapes
# at each generation, and the solid curve the analytical geometric mean
# fitness landscape.
make_landscape_panel <- function(species_name, tag_label,
                                 show_y_title=FALSE,
                                 show_gmf_label=FALSE,
                                 show_reference=FALSE) {

  landscape_i <- landscape_exact[
    landscape_exact$species_id == species_name, ]
  avg_approx_i <- avg_approx[avg_approx$species_id == species_name, ]
  summary_i <- summary_landscape[
    summary_landscape$species_id == species_name, ]

  reference_i <- avg_approx[
    avg_approx$species_id == "Specialist" &
      show_reference, ]
  reference_summary <- summary_landscape[
    summary_landscape$species_id == "Specialist" &
      show_reference, ]

  p <- ggplot() +
    geom_line(data=landscape_i,
              aes(x=zi,
                  y=plot_fitness(log_fitness, fitness_log_multiplier),
                  group=time, color=species_id),
              linewidth=0.24, alpha=0.045) +
    geom_line(data=reference_i,
              aes(x=zi,
                  y=plot_fitness(log_fitness, fitness_log_multiplier),
                  color=species_id),
              linewidth=0.70, alpha=0.55) +
    geom_line(data=avg_approx_i,
              aes(x=zi,
                  y=plot_fitness(log_fitness, fitness_log_multiplier),
                  color=species_id),
              linewidth=1.25) +
    geom_point(data=reference_summary,
               aes(x=zi_focal,
                   y=plot_fitness(log_fitness_focal,
                                  fitness_log_multiplier),
                   fill=species_id),
               shape=21, color=pal_species["Specialist"],
               alpha=0.55, size=2.3, stroke=0.45,
               show.legend=FALSE) +
    geom_point(data=summary_i,
               aes(x=zi_focal,
                   y=plot_fitness(log_fitness_focal,
                                  fitness_log_multiplier),
                   fill=species_id),
               shape=21, color="black", size=2.7, stroke=0.45,
               show.legend=FALSE) +
    annotate("text", x=-3.32, y=1.035,
             label=species_name,
             hjust=0, size=5.2, color=pal_species[species_name]) +
    scale_color_manual(values=pal_species) +
    scale_fill_manual(values=pal_species) +
    scale_x_continuous(breaks=c(-3, -1.5, 0, 1.5, 3),
                       expand=expansion(mult=c(0, 0))) +
    scale_y_continuous(breaks=seq(0, 1, by=0.2)) +
    coord_cartesian(xlim=landscape_x_lim,
                    ylim=landscape_y_lim,
                    clip="on") +
    labs(x="Trait value of focal species",
         y=ifelse(show_y_title, "Population mean fitness", ""),
         tag=tag_label) +
    theme_fig2 +
    theme(legend.position="none")

  if (show_gmf_label) {
    p <- p +
      annotate("text", x=summary_i$zi_focal + 0.12,
               y=plot_fitness(summary_i$log_fitness_focal,
                              fitness_log_multiplier) + 0.030,
               label="Analytical GMF",
               hjust=0, size=4.2, color="gray20")
  }

  p <- add_network_inset(p, species_name)

  return(p)
}

p2.1 <- make_landscape_panel(species_name="Specialist",
                             tag_label="a",
                             show_y_title=TRUE,
                             show_gmf_label=TRUE)

p2.2 <- make_landscape_panel(species_name="Generalist",
                             tag_label="b",
                             show_y_title=TRUE,
                             show_gmf_label=TRUE,
                             show_reference=TRUE)


#---------------------------------------------#
#-------------Assembling figure---------------#
#---------------------------------------------#

fig2 <- ((p2.1 | p2.2) / p2.3) +
  plot_layout(heights=c(1.12, 0.90))

print(fig2)

fig2_svg <- here("figures", "main_fig2.svg")
fig2_png <- here("figures", "main_fig2.png")

ggsave(fig2_svg, fig2,
       height=8.6, width=11.2, device=svglite::svglite)

ggsave(fig2_png, fig2,
       height=8.6, width=11.2, dpi=320)

cat("Saved Figure 2 to:\n")
cat("  ", fig2_svg, "\n", sep="")
cat("  ", fig2_png, "\n", sep="")
