# ====================================================================== #
# main_fig4.R
#
# Figure 4 of the main text.
#
# This script builds Figure 4, which shows how indirect evolutionary effects
# link network structure to the height of adaptive peaks. Networks with more
# pathways connecting species have a larger contribution of indirect
# effects. In turn, species coevolving within such networks experience a
# larger conflict of evolutionary effects and a higher correlation among
# their mutualistic partners.
#
# The three panels show that:
#
#   a) the conflict of evolutionary effects, S_i - m S_bar (stored as R_i
#      and R_bar), increases with the contribution of indirect effects;
#
#   b) the correlation among partners, r_i^2 (stored as r2_partner), also
#      increases with the contribution of indirect effects;
#
#   c) the contribution of indirect effects increases with the number of
#      pathways connecting species, 1 - 1/lambda, where lambda is the
#      leading eigenvalue of the adjacency matrix.
#
# The contribution of indirect effects appears on the horizontal axis of
# panels a and b and on the vertical axis of panel c.
#
# The conflict of evolutionary effects, S_i - m S_bar, is the variance of
# selective targets divided by the variance of the environmental optima,
# because E[V_i] = sd_theta^2 (S_i - m S_bar). The conflict increases when
# the coevolutionary pathways reaching the focal species are more diverse
# than the ones reaching its partners. In this case, the trait of the focal
# species is affected by many environmental optima, while each of its
# partners tracks only a few of them.
#
# Panels a and b describe species with three mutualistic partners. Fixing the
# number of partners isolates the effects of network structure from those of
# the number of partners. Both panels average the species with three partners
# within each network, so each point is one empirical network.
#
# Inputs:
#   results/coevolution_gmf.csv
#   src/empirical_networks/*.csv
#
# Outputs:
#   figures/main_fig4.svg
#   figures/main_fig4.png
#   results/fig4_indirecteffects_fits.csv
# ====================================================================== #


#---------------------------------------------#
#---------Loading packages and settings-------#
#---------------------------------------------#

library(here)
library(ggplot2)
library(ggpubr)
library(data.table)
library(patchwork)
library(svglite)

source(here("src", "aux_functions.R"))

dir.create(here("figures"), showWarnings=FALSE, recursive=TRUE)


#---------------------------------------------#
#-----------Data and figure settings----------#
#---------------------------------------------#

results_file <- here("results", "coevolution_gmf.csv")

# Panels a and b fix the number of partners at three, since species with
# three partners are present in most networks of the dataset.
focal_degree <- 3

# Colors of the three example networks whose structure is drawn in the small
# bipartite diagrams.
pal_inset <- c(very_low="#8F5A28",
               low="#C78342",
               high="#2A6F77")

# First and last vertical slot of an inset, passed to place_guild_nodes() in
# aux_functions.R.
inset_y_min <- 0.07
inset_y_max <- 0.93

if (!file.exists(results_file)) {
  stop("Figure 4 requires the coevolutionary result table: ", results_file)
}


#---------------------------------------------#
#---------Preparing network-level data--------#
#---------------------------------------------#

coevo <- fread(results_file)

m <- unique(coevo$m)

if (length(m) != 1) {
  stop("The result table must use a single value of m.")
}

# The conflict of evolutionary effects, S_i - m S_bar in the main text. R_i
# is the Simpson diversity of the coevolutionary pathways reaching species i,
# and R_bar is the same quantity averaged over its partners. Subtracting
# m R_bar measures how much more diverse the pathways reaching the focal
# species are than the ones reaching its partners, which increases the
# variance of its selective targets.
coevo[, conflict := R_i - m * R_bar]

# Panels a and b use the species with three partners whose T-matrix terms are
# finite. Their correlation among partners, r2_partner, is the temporal
# correlation between the squared trait differences of a species with its
# partners, averaged across all pairs of partners (r_i^2 in the main text).

focal_species <- coevo[
  degree == focal_degree &
    is.finite(r2_partner) &
    is.finite(nu_bar) &
    nu_bar > 0 &
    is.finite(conflict) &
    is.finite(indirect_network)
]

if (nrow(focal_species) == 0) {
  stop("No species with ", focal_degree, " partners carry finite T-matrix terms.")
}

# Species-level terms are averaged within each network, so each row describes
# one network.
degree_data <- focal_species[, .(
  r2_partner=mean(r2_partner),
  conflict=mean(conflict),
  indirect_network=unique(indirect_network),
  leading_eigenvalue=unique(leading_eigenvalue),
  n_species_degree=.N
), by=network]

# Panel c uses every analysed network, not only those containing species
# with three partners.
network_data <- unique(coevo[
  is.finite(indirect_network) & is.finite(leading_eigenvalue),
  .(network, indirect_network, leading_eigenvalue, n_species)
])

# Number of pathways connecting species, 1 - 1/lambda (Methods).
network_data[, lambda_linear := 1 - 1 / leading_eigenvalue]

cat("\nNetworks with species of degree", focal_degree, ":",
    nrow(degree_data), "\n")
cat("Networks in panel c:", nrow(network_data), "\n")
cat("Species represented in panels a and b:",
    sum(degree_data$n_species_degree), "\n\n")


#---------------------------------------------#
#----------Model summaries by panel-----------#
#---------------------------------------------#

# Simple linear regressions for the three panel relationships (Methods).
summarize_fit <- function(data, response, predictor, panel) {

  fit <- lm(reformulate(predictor, response), data=data)
  pred <- predict(fit)
  slope_ci <- confint(fit)[2, ]

  out <- data.table(
    panel=panel,
    response=response,
    predictor=predictor,
    n=nrow(data),
    pearson=cor(data[[response]], data[[predictor]]),
    spearman=cor(data[[response]], data[[predictor]], method="spearman"),
    adj_r2=summary(fit)$adj.r.squared,
    slope=coef(fit)[2],
    slope_low=slope_ci[1],
    slope_high=slope_ci[2],
    p_value=summary(fit)$coefficients[2, 4],
    rmse=sqrt(mean((data[[response]] - pred)^2))
  )

  return(out)
}

fit_summary <- rbind(
  summarize_fit(degree_data, "conflict", "indirect_network", "a"),
  summarize_fit(degree_data, "r2_partner", "indirect_network", "b"),
  summarize_fit(network_data, "indirect_network", "lambda_linear", "c")
)

fwrite(fit_summary, here("results", "fig4_indirecteffects_fits.csv"))


#---------------------------------------------#
#----------Network inset preparation----------#
#---------------------------------------------#

# Choose three example networks along the gradient of indirect effects. Two
# networks come from the low-contribution region and one from the
# high-contribution region, all close to the fitted line.
degree_data[, fitted_r2 := as.numeric(predict(
  lm(r2_partner ~ indirect_network)
))]

select_inset_example <- function(data, indirect_prob) {

  d <- copy(data)

  indirect_range <- diff(range(d$indirect_network))
  r2_range <- diff(range(d$r2_partner))
  n_range <- diff(range(d$n_species_degree))

  if (indirect_range == 0) indirect_range <- 1
  if (r2_range == 0) r2_range <- 1
  if (n_range == 0) n_range <- 1

  target_indirect <- as.numeric(quantile(d$indirect_network,
                                         probs=indirect_prob,
                                         names=FALSE))
  target_n <- as.numeric(quantile(d$n_species_degree,
                                  probs=0.50,
                                  names=FALSE))

  d[, example_score :=
      abs(indirect_network - target_indirect) / indirect_range +
      abs(r2_partner - fitted_r2) / r2_range +
      0.25 * abs(n_species_degree - target_n) / n_range]

  out <- d[which.min(example_score)]

  return(out)
}

very_low_example <- select_inset_example(degree_data, indirect_prob=0.05)
low_example <- select_inset_example(degree_data, indirect_prob=0.20)
high_example <- select_inset_example(degree_data, indirect_prob=0.80)

very_low_example[, inset_position := "very_low"]
low_example[, inset_position := "low"]
high_example[, inset_position := "high"]

inset_examples <- rbind(very_low_example, low_example, high_example)

inset_examples[, inset_position := factor(inset_position,
                                          levels=c("very_low",
                                                   "low",
                                                   "high"))]

# The same three networks are marked in every panel. Panel c needs their
# number of pathways connecting species.
inset_examples <- merge(
  inset_examples,
  network_data[, .(network, lambda_linear)],
  by="network",
  sort=FALSE
)


# Build the species and interactions of the inset for one example network,
# rescaled into the x and y units of the panel that carries them.
build_network_inset <- function(example, x0, x1, y0, y1) {

  net <- read_network(file=here("src", "empirical_networks",
                                paste0(example$network, ".csv")),
                      max_species=300)

  incidence <- net$incidence

  plant_degree <- rowSums(incidence)
  animal_degree <- colSums(incidence)

  plant_highlight <- rownames(incidence)[plant_degree == focal_degree]
  animal_highlight <- colnames(incidence)[animal_degree == focal_degree]

  plant_pos <- place_guild_nodes(
    names=rownames(incidence),
    degree=plant_degree,
    highlighted=plant_highlight,
    y_min=inset_y_min,
    y_max=inset_y_max
  )

  animal_pos <- place_guild_nodes(
    names=colnames(incidence),
    degree=animal_degree,
    highlighted=animal_highlight,
    y_min=inset_y_min,
    y_max=inset_y_max
  )

  plant_pos$x <- 0
  animal_pos$x <- 1

  plant_pos$node_class <- ifelse(plant_pos$species %in% plant_highlight,
                                 "highlight", "background")
  animal_pos$node_class <- ifelse(animal_pos$species %in% animal_highlight,
                                  "highlight", "background")

  node_data <- rbind(
    data.frame(species=plant_pos$species,
               x=plant_pos$x,
               y=plant_pos$y,
               node_class=plant_pos$node_class,
               stringsAsFactors=FALSE),
    data.frame(species=animal_pos$species,
               x=animal_pos$x,
               y=animal_pos$y,
               node_class=animal_pos$node_class,
               stringsAsFactors=FALSE)
  )

  edge_id <- which(incidence == 1, arr.ind=TRUE)

  edge_data <- data.frame(
    plant=rownames(incidence)[edge_id[, 1]],
    animal=colnames(incidence)[edge_id[, 2]],
    stringsAsFactors=FALSE
  )

  edge_data$x <- 0
  edge_data$xend <- 1
  edge_data$y <- plant_pos$y[match(edge_data$plant, plant_pos$species)]
  edge_data$yend <- animal_pos$y[match(edge_data$animal, animal_pos$species)]

  edge_data$edge_class <- ifelse(
    edge_data$plant %in% plant_highlight |
      edge_data$animal %in% animal_highlight,
    "highlight", "background"
  )

  node_data$x <- x0 + node_data$x * (x1 - x0)
  node_data$y <- y0 + node_data$y * (y1 - y0)
  edge_data$x <- x0 + edge_data$x * (x1 - x0)
  edge_data$xend <- x0 + edge_data$xend * (x1 - x0)
  edge_data$y <- y0 + edge_data$y * (y1 - y0)
  edge_data$yend <- y0 + edge_data$yend * (y1 - y0)

  node_data$inset_position <- example$inset_position
  edge_data$inset_position <- example$inset_position

  return(list(nodes=node_data, edges=edge_data))
}


#---------------------------------------------#
#-------------Panel axis limits---------------#
#---------------------------------------------#

pad_range <- function(values, lower=0.06, upper=0.06, floor_at_zero=FALSE) {

  limits <- range(values)
  span <- diff(limits)

  out <- c(limits[1] - lower * span, limits[2] + upper * span)

  if (floor_at_zero) out[1] <- max(0, out[1])

  return(out)
}

indirect_lim <- pad_range(degree_data$indirect_network)
r2_lim <- pad_range(degree_data$r2_partner, lower=0.10, upper=0.10,
                    floor_at_zero=TRUE)

# The conflict panel carries the three network diagrams, so it is given
# extra headroom above the point cloud.
conflict_lim <- pad_range(degree_data$conflict, lower=0.10, upper=0.34)

lambda_lim <- pad_range(network_data$lambda_linear)
indirect_c_lim <- pad_range(network_data$indirect_network,
                            lower=0.10, upper=0.10, floor_at_zero=TRUE)


#---------------------------------------------#
#---------------Panel factory-----------------#
#---------------------------------------------#

# Axis titles and panel tags are enlarged relative to the base size so that
# the saved figure matches the text of Figures 2 and 3. Tick labels keep the
# base size.
theme_fig <- theme_pubr(base_size=15) +
  theme(legend.position="none",
        axis.title=element_text(color="gray10", size=17),
        axis.text=element_text(color="gray20"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid=element_blank(),
        plot.margin=margin(t=20, r=6, b=6, l=6),
        plot.tag=element_text(face="bold", size=20),
        plot.tag.position="topleft",
        plot.tag.location="margin")


# Draw one scatter panel with its fitted line, its gray point cloud, and the
# three highlighted example networks. The network diagrams are drawn only on
# the panel that carries them.
plot_relationship <- function(data, x, y, x_lim, y_lim, x_label, y_label,
                              example_data, insets=NULL) {

  p <- ggplot(data, aes(x=.data[[x]], y=.data[[y]])) +
    geom_smooth(method="lm", formula=y ~ x, se=FALSE,
                color="black", linewidth=0.78, alpha=0.95) +
    geom_point(shape=21, size=2.45,
               fill="gray72", color="white",
               stroke=0.42, alpha=0.78)

  if (!is.null(insets)) {

    p <- p +
      geom_segment(data=insets$edges[edge_class == "background"],
                   aes(x=x, y=y, xend=xend, yend=yend),
                   inherit.aes=FALSE,
                   color="gray72", linewidth=0.14, alpha=0.36) +
      geom_segment(data=insets$edges[edge_class == "highlight"],
                   aes(x=x, y=y, xend=xend, yend=yend,
                       color=inset_position),
                   inherit.aes=FALSE,
                   linewidth=0.32, alpha=0.78) +
      geom_point(data=insets$nodes[node_class == "background"],
                 aes(x=x, y=y),
                 inherit.aes=FALSE,
                 color="gray62", size=0.55, alpha=0.68) +
      geom_point(data=insets$nodes[node_class == "highlight"],
                 aes(x=x, y=y, color=inset_position),
                 inherit.aes=FALSE,
                 size=1.35, alpha=0.95)
  }

  p <- p +
    geom_point(data=example_data,
               aes(x=.data[[x]], y=.data[[y]], fill=inset_position),
               inherit.aes=FALSE,
               shape=21, size=2.70,
               color="black",
               stroke=0.50, alpha=0.95) +
    scale_color_manual(values=pal_inset, guide="none") +
    scale_fill_manual(values=pal_inset, guide="none") +
    scale_x_continuous(expand=expansion(mult=c(0.02, 0.02))) +
    scale_y_continuous(expand=expansion(mult=c(0.02, 0.04))) +
    coord_cartesian(xlim=x_lim, ylim=y_lim, clip="on") +
    labs(x=x_label, y=y_label) +
    theme_fig

  return(p)
}


#---------------------------------------------#
#----------------Making panels----------------#
#---------------------------------------------#

# Placement of the network diagrams within panel a, in the units of its axes.
host_lim <- conflict_lim

inset_examples[, inset_x0 := indirect_lim[1] + 0.040 * diff(indirect_lim)]
inset_examples[, inset_x1 := indirect_lim[1] + 0.185 * diff(indirect_lim)]

inset_examples[inset_position == "high",
               inset_y0 := host_lim[2] - 0.150 * diff(host_lim)]
inset_examples[inset_position == "high",
               inset_y1 := host_lim[2] - 0.012 * diff(host_lim)]

inset_examples[inset_position == "low",
               inset_y0 := host_lim[2] - 0.315 * diff(host_lim)]
inset_examples[inset_position == "low",
               inset_y1 := host_lim[2] - 0.177 * diff(host_lim)]

inset_examples[inset_position == "very_low",
               inset_y0 := host_lim[2] - 0.480 * diff(host_lim)]
inset_examples[inset_position == "very_low",
               inset_y1 := host_lim[2] - 0.342 * diff(host_lim)]

inset_list <- lapply(seq_len(nrow(inset_examples)), function(i) {
  build_network_inset(inset_examples[i],
                      x0=inset_examples$inset_x0[i],
                      x1=inset_examples$inset_x1[i],
                      y0=inset_examples$inset_y0[i],
                      y1=inset_examples$inset_y1[i])
})

panel_a_insets <- list(
  edges=rbindlist(lapply(inset_list, `[[`, "edges")),
  nodes=rbindlist(lapply(inset_list, `[[`, "nodes"))
)

indirect_label <- "Contribution of indirect effects"

# Panel a carries the conflict of evolutionary effects, and with it the
# three network diagrams. Panel b repeats the same networks against the
# correlation among partners, so the two share an axis and can be read side
# by side.
p1.1 <- plot_relationship(
  data=degree_data,
  x="indirect_network",
  y="conflict",
  x_lim=indirect_lim,
  y_lim=conflict_lim,
  x_label=indirect_label,
  y_label="Conflicting evolutionary effects",
  example_data=inset_examples,
  insets=panel_a_insets
)

p1.2 <- plot_relationship(
  data=degree_data,
  x="indirect_network",
  y="r2_partner",
  x_lim=indirect_lim,
  y_lim=r2_lim,
  x_label=indirect_label,
  y_label="Correlation among partners",
  example_data=inset_examples
)

p1.3 <- plot_relationship(
  data=network_data,
  x="lambda_linear",
  y="indirect_network",
  x_lim=lambda_lim,
  y_lim=indirect_c_lim,
  x_label="Number of pathways connecting species",
  y_label=indirect_label,
  example_data=inset_examples
)


#---------------------------------------------#
#-------------Assembling figure---------------#
#---------------------------------------------#

fig <- (p1.1 | p1.2) / p1.3 +
  plot_layout(heights=c(1, 1)) +
  plot_annotation(tag_levels="a")

print(fig)


#---------------------------------------------#
#----------------Saving figure----------------#
#---------------------------------------------#

fig_svg <- here("figures", "main_fig4.svg")
fig_png <- here("figures", "main_fig4.png")

ggsave(fig_svg, fig,
       height=8.40, width=8.60, device=svglite::svglite)

ggsave(fig_png, fig,
       height=8.40, width=8.60, dpi=320)

cat("Saved Figure 4 to:\n")
cat("  ", fig_svg, "\n", sep="")
cat("  ", fig_png, "\n", sep="")

cat("\nExample networks:\n")
print(inset_examples[, .(inset_position, network, indirect_network,
                         r2_partner, conflict, n_species_degree)])

cat("\nPanel fits:\n")
print(fit_summary[, .(panel, response, predictor, n, pearson,
                      adj_r2, slope, p_value)])
