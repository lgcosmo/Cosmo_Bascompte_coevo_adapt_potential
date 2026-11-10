# ====================================================================== #
# main_fig1.R
#
# Figure 1 of the main text.
#
# The figure is a single ridgeline panel. Each ridge is the fitness
# landscape of the focal species at one generation, stacked from generation
# 1 at the top downwards, and the bottom ridge is the temporal average of
# these landscapes, the geometric mean fitness landscape. The figure shows
# that the peak of the landscape moves and changes in height from one
# generation to the next because the traits of mutualistic partners and the
# environmental optimum change.
#
# Each landscape is computed from the selective targets of its generation,
# which are drawn on the baseline of the ridge: plant icons mark the trait
# values of the two mutualistic partners, a green diamond marks the
# environmental optimum, and a butterfly marks the fixed trait of the focal
# species. The adaptive peak is the weighted average of these targets, and
# the height of the peak follows from their spread (equation 2 of the main
# text):
#
#   log w_i(z_i, t) = -V_i(t) + B_i(t) - [z_i - z_i,opt(t)]^2,
#
# where V_i is the variance of selective targets and B_i = (m alpha / 2) V_D
# is the contribution of the variance in trait matching, V_D.
#
# The vertical axis shows log population mean fitness. Each ridge runs from
# log fitness -2 at its baseline to 0 (w = 1) at its top, so that the same
# values can be read on every ridge. These values are illustrative.
#
# Inputs:
#   figures/butterfly_1.png
#   figures/plant_1.png, figures/plant_2.png
#
# Outputs:
#   figures/main_fig1.svg
#   figures/main_fig1.png
# ====================================================================== #


#---------------------------------------------#
#---------Loading packages and setup----------#
#---------------------------------------------#

library(here)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(grid)
library(png)
library(svglite)

dir.create(here("figures"), showWarnings=FALSE, recursive=TRUE)


#---------------------------------------------#
#-----------Graphical parameters--------------#
#---------------------------------------------#

# Proportion of the selective pressures that comes from mutualisms, and how
# strongly the benefits of mutualisms decay with trait mismatches.
m <- 0.5
alpha <- 0.16

# Horizontal range of the figure, and the trait values at which each
# landscape is evaluated.
x_limits <- c(0.55, 7.05)
zi_seq <- seq(x_limits[1], x_limits[2], length.out=700)

color_average <- "black"
color_target_env <- "#009E73"

# Guides for the focal species. The key label uses a slightly darker shade
# of the butterfly blue (#0072B2).
color_focal <- "gray25"
color_focal_label <- "#035B8F"

theme_fig1 <- theme_pubr(base_size=15) +
  theme(legend.position="none",
        legend.title=element_blank(),
        axis.title=element_text(color="gray10"),
        axis.text=element_text(color="gray20"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid=element_blank())

butterfly_img <- readPNG(here("figures", "butterfly_1.png"))

# One plant silhouette per mutualistic partner, so that the same partner can
# be followed across generations.
plant_imgs <- lapply(seq_len(2), function(i) {
  readPNG(here("figures", paste0("plant_", i, ".png")))
})


#---------------------------------------------#
#------Landscape helper functions-------------#
#---------------------------------------------#

# Compute a weighted variance for the selective targets experienced by the
# focal species. The environmental optimum has weight 1-m, while all partners
# together have weight m and each partner has weight m/k_i.
weighted_variance <- function(x, w) {

  w <- w / sum(w)
  x_bar <- sum(w * x)
  v <- sum(w * x^2) - x_bar^2

  return(v)
}


# Compute the portfolio effect of generalists at the peak of the approximate
# landscape, B_i = (m alpha / 2) V_D, where V_D is the variance in trait
# matching among partners.
compute_partner_benefit <- function(z_j, z_peak, m, alpha) {

  D <- (z_j - z_peak)^2
  V_D <- mean(D^2) - mean(D)^2
  B_i <- (m * alpha / 2) * V_D

  return(B_i)
}


# Compute the position and the height of the adaptive peak of one
# generation from its selective targets, together with the squared distance
# of the focal species to the peak.
compute_landscape_terms <- function(z_j, theta_i, z_i, m, alpha, time) {

  k_i <- length(z_j)
  targets <- c(theta_i, z_j)
  weights <- c(1 - m, rep(m / k_i, k_i))

  z_peak <- sum(weights * targets)
  V_i <- weighted_variance(x=targets, w=weights)
  B_i <- compute_partner_benefit(z_j=z_j, z_peak=z_peak, m=m, alpha=alpha)
  peak_height <- -V_i + B_i

  out <- data.frame(time=time,
                    z_i=z_i,
                    theta_i=theta_i,
                    z_peak=z_peak,
                    peak_height=peak_height,
                    distance=(z_i - z_peak)^2,
                    stringsAsFactors=FALSE)

  return(out)
}


# Place one ecological icon within data coordinates.
add_png_icon <- function(p, img, x, y, width, height) {

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


#---------------------------------------------#
#-------Selective targets per generation------#
#---------------------------------------------#

# The focal species keeps the same trait across generations, while the
# traits of its two partners and its environmental optimum change. The
# three generations were chosen so that both elements of the landscape
# change markedly. In generation 1, the environmental optimum and the
# partners pull in opposite directions, so the peak is low and sits at an
# intermediate distance from the focal species. In generation 2, the
# targets cluster together further to the right, so the peak is highest but
# also farthest from the focal species, which has its lowest fitness. In
# generation 3, the targets straddle the focal species more closely, so the
# peak is of intermediate height but the focal species sits right at it and
# has its highest fitness.
z_focal <- 2.40

targets <- list(
  list(time=1L, theta_i=2.00, z_j=c(3.10, 3.70)),
  list(time=2L, theta_i=3.05, z_j=c(2.78, 3.42)),
  list(time=3L, theta_i=1.80, z_j=c(2.80, 3.20))
)

# Curvature used to draw the landscapes. In the decomposition of log
# fitness, the squared distance to the peak has a coefficient of one, which
# gives very wide ridges in this conceptual figure. The squared distance to
# the peak is therefore scaled by this factor for display only. The position
# and the height of the peak of every generation do not depend on it.
curvature <- 2.2

terms <- lapply(targets, function(g) {
  compute_landscape_terms(z_j=g$z_j, theta_i=g$theta_i, z_i=z_focal,
                          m=m, alpha=alpha, time=g$time)
})
terms <- do.call(rbind, terms)

# Log fitness of the focal species at each generation, with the squared
# distance to the peak scaled by the display curvature.
terms$log_fitness_i <- terms$peak_height - curvature * terms$distance

# Long table with one row per partner and generation, used to draw the
# plant icons.
partners <- lapply(targets, function(g) {
  data.frame(time=g$time,
             partner=seq_along(g$z_j),
             z_j=g$z_j,
             stringsAsFactors=FALSE)
})
partners <- do.call(rbind, partners)

# Landscape of each generation across the sequence of focal trait values.
landscapes <- lapply(seq_len(nrow(terms)), function(i) {
  data.frame(time=terms$time[i],
             zi=zi_seq,
             log_fitness=terms$peak_height[i] -
               curvature * (zi_seq - terms$z_peak[i])^2,
             stringsAsFactors=FALSE)
})
landscapes <- do.call(rbind, landscapes)

# The geometric mean fitness landscape is the pointwise temporal average of
# the three landscapes on the log scale.
landscape_average <- landscapes %>%
  group_by(zi) %>%
  summarise(log_fitness=mean(log_fitness), .groups="drop")

log_gmf <- approx(x=landscape_average$zi, y=landscape_average$log_fitness,
                  xout=z_focal)$y


#---------------------------------------------#
#-----------Ridgeline coordinates-------------#
#---------------------------------------------#

# Each generation gets one horizontal baseline. The height of a ridge above
# its baseline is the log fitness of that generation, rescaled so that the
# baseline corresponds to log fitness z_floor (lower values are clipped)
# and the top of the ridge to log fitness 0, i.e., w = 1. Baselines are
# spaced so that the icons and the tick labels of one generation do not
# touch the ridge of the previous generation.
z_floor <- -2
ridge_height <- 1.40
baseline_step <- 1.75

# Ridges are numbered from the bottom: the temporal average first, then the
# generations in reverse order, so that generation 1 sits at the top.
ridge_rank <- rev(c(terms$time - 1, nrow(terms)))

ridge_baselines <- data.frame(
  time=c(terms$time, 0L),
  y_base=0.60 + baseline_step * ridge_rank,
  label=c(paste("Time =", terms$time), "Temporal average"),
  stringsAsFactors=FALSE
)

# Convert log fitness into the vertical position of a ridge.
ridge_y <- function(log_fitness, y_base) {

  y <- y_base + ridge_height * (pmax(log_fitness, z_floor) - z_floor) /
    abs(z_floor)

  return(y)
}

landscapes <- merge(landscapes,
                    ridge_baselines[ridge_baselines$time > 0,
                                    c("time", "y_base")],
                    by="time")
landscapes$y_ridge <- ridge_y(landscapes$log_fitness, landscapes$y_base)

landscape_average$time <- 0L
landscape_average$y_base <- ridge_baselines$y_base[ridge_baselines$time == 0L]
landscape_average$y_ridge <- ridge_y(landscape_average$log_fitness,
                                     landscape_average$y_base)

terms <- merge(terms,
               ridge_baselines[ridge_baselines$time > 0,
                               c("time", "y_base")],
               by="time")
terms$y_focal <- ridge_y(terms$log_fitness_i, terms$y_base)
terms$y_peak <- ridge_y(terms$peak_height, terms$y_base)

partners <- merge(partners,
                  ridge_baselines[ridge_baselines$time > 0,
                                  c("time", "y_base")],
                  by="time")

y_gmf <- ridge_y(log_gmf, landscape_average$y_base[1])

# Log fitness values shown on the vertical axis of every ridge, from the
# baseline (z_floor) to the top of the ridge (0).
log_ticks <- c(-2, -1, 0)
y_breaks <- as.vector(outer(ridge_y(log_ticks, 0), ridge_baselines$y_base,
                            "+"))
y_labels <- rep(format(log_ticks), times=nrow(ridge_baselines))

# All generations share the same light grey, so that only the black
# temporal average stands out.
pal_time <- rep("#B5B5B5", nrow(terms))
names(pal_time) <- as.character(terms$time)

# Vertical offsets of the icons above each baseline, and icon sizes in data
# units. The environmental optimum shares the icon row, so the target values
# above keep every diamond at least one icon width away from the plants and
# the butterfly.
icon_lift <- 0.23
icon_width <- 0.36
icon_height <- 0.32
point_size <- 3.2
diamond_size <- 4.0

y_limits <- c(0.25, 7.95)
x_label_right <- 7.00


#---------------------------------------------#
#-------------Building the figure-------------#
#---------------------------------------------#

p1 <- ggplot() +
  geom_hline(data=ridge_baselines,
             aes(yintercept=y_base),
             color="gray88", linewidth=0.35) +
  # Thin guides from the baseline to the apex of each peak, so that the
  # peak can be read against the positions of the selective targets.
  geom_segment(data=terms,
               aes(x=z_peak, xend=z_peak, y=y_base, yend=y_peak),
               color="gray55", linewidth=0.4, linetype="dotted") +
  geom_ribbon(data=landscapes,
              aes(x=zi, ymin=y_base, ymax=y_ridge,
                  group=time, fill=factor(time)),
              alpha=0.30) +
  geom_line(data=landscapes,
            aes(x=zi, y=y_ridge, group=time, color=factor(time)),
            linewidth=0.95) +
  geom_ribbon(data=landscape_average,
              aes(x=zi, ymin=y_base, ymax=y_ridge),
              fill="gray25", alpha=0.20) +
  geom_line(data=landscape_average,
            aes(x=zi, y=y_ridge),
            color=color_average, linewidth=1.25) +
  # Thin guide from the trait of the focal species to its mean fitness.
  geom_segment(data=terms,
               aes(x=z_i, xend=z_i, y=y_base + icon_lift + icon_height / 2,
                   yend=y_focal),
               color=color_focal, linewidth=0.4, linetype="dotted") +
  geom_point(data=terms,
             aes(x=z_i, y=y_focal, fill=factor(time)),
             shape=21, color="black", size=point_size, stroke=0.4) +
  # Environmental optimum of each generation.
  geom_point(data=terms,
             aes(x=theta_i, y=y_base + icon_lift),
             shape=23, color="black", fill=color_target_env,
             size=diamond_size, stroke=0.4) +
  # Geometric mean fitness on the temporal average landscape.
  annotate("point", x=z_focal, y=y_gmf,
           color="black", fill="white", shape=21,
           size=point_size + 0.4, stroke=0.7) +
  annotate("text", x=0.70, y=y_gmf + 0.58,
           label="Geometric mean fitness (log)",
           hjust=0, size=4.4, color="black") +
  annotate("segment", x=1.72, xend=z_focal - 0.08,
           y=y_gmf + 0.42, yend=y_gmf + 0.09,
           arrow=arrow(length=unit(0.10, "cm")),
           color="black", linewidth=0.55) +
  geom_text(data=ridge_baselines,
            aes(x=x_label_right, y=y_base + 0.14, label=label),
            hjust=1, size=4.0, color="gray20") +
  scale_color_manual(values=pal_time) +
  scale_fill_manual(values=pal_time) +
  scale_x_continuous(breaks=1:7,
                     limits=x_limits,
                     expand=expansion(mult=c(0, 0))) +
  scale_y_continuous(breaks=y_breaks,
                     labels=y_labels,
                     limits=y_limits,
                     expand=expansion(mult=c(0, 0))) +
  coord_cartesian(clip="off") +
  labs(x="Trait value of focal species",
       y="Population mean fitness (log)") +
  theme_fig1 +
  theme(axis.text.y=element_text(size=9.5),
        plot.margin=margin(5.5, 14, 5.5, 5.5))

# Distance to the peak and height of the peak, both labelled on the first
# generation. The distance arrow runs from the focal species to the peak at
# the height of its mean fitness. The height arrow stands to the right of
# the icons and is connected to the apex by a dotted guide. A third label
# points to the mean fitness of the focal species at this generation, in
# the same way as the geometric mean fitness label on the bottom ridge.
first_row <- terms[terms$time == 1L, ]
height_arrow_x <- max(partners$z_j[partners$time == 1L]) + icon_width

p1 <- p1 +
  annotate("text", x=0.70, y=first_row$y_focal + 0.62,
           label="Fitness~at~time~italic(t)~(log)", parse=TRUE,
           hjust=0, size=4.4, color="black") +
  annotate("segment", x=1.72, xend=first_row$z_i - 0.10,
           y=first_row$y_focal + 0.46, yend=first_row$y_focal + 0.10,
           arrow=arrow(length=unit(0.10, "cm")),
           color="black", linewidth=0.55) +
  annotate("segment",
           x=first_row$z_i + 0.08, xend=first_row$z_peak,
           y=first_row$y_focal, yend=first_row$y_focal,
           arrow=arrow(length=unit(0.12, "cm"), ends="both"),
           color="gray20", linewidth=0.5) +
  annotate("text", x=(first_row$z_i + first_row$z_peak) / 2 + 0.05,
           y=first_row$y_focal - 0.26,
           label="Distance\nto the peak",
           size=3.1, lineheight=0.85, color="gray20") +
  annotate("segment",
           x=first_row$z_peak, xend=height_arrow_x,
           y=first_row$y_peak, yend=first_row$y_peak,
           color="gray55", linewidth=0.4, linetype="dotted") +
  annotate("segment",
           x=height_arrow_x, xend=height_arrow_x,
           y=first_row$y_base + 0.02, yend=first_row$y_peak,
           arrow=arrow(length=unit(0.12, "cm"), ends="both"),
           color="gray20", linewidth=0.5) +
  annotate("text", x=height_arrow_x + 0.13,
           y=(first_row$y_base + first_row$y_peak) / 2,
           label="Peak height",
           hjust=0, size=3.1, color="gray20")

# Icons of the selective targets on the baseline of each generation.
for (i in seq_len(nrow(partners))) {
  p1 <- add_png_icon(p1, plant_imgs[[partners$partner[i]]],
                     x=partners$z_j[i],
                     y=partners$y_base[i] + icon_lift,
                     width=icon_width, height=icon_height)
}

for (i in seq_len(nrow(terms))) {
  p1 <- add_png_icon(p1, butterfly_img,
                     x=terms$z_i[i],
                     y=terms$y_base[i] + icon_lift,
                     width=icon_width * 1.15, height=icon_height)
}

# Key for the icons, placed in the upper right of the figure, where the
# landscape of the first generation is flat.
key_x <- 4.60
key_y <- c(7.70, 7.32, 6.94)

p1 <- add_png_icon(p1, butterfly_img, x=key_x, y=key_y[1],
                   width=icon_width * 1.15, height=icon_height)
p1 <- add_png_icon(p1, plant_imgs[[1]], x=key_x, y=key_y[2],
                   width=icon_width, height=icon_height)

p1 <- p1 +
  annotate("point", x=key_x, y=key_y[3],
           shape=23, color="black", fill=color_target_env,
           size=diamond_size, stroke=0.4) +
  annotate("text", x=key_x + 0.30, y=key_y,
           label=c("Focal species trait value",
                   "Mutualistic partners trait values",
                   "Environmental optimum"),
           hjust=0, size=4.0,
           color=c(color_focal_label, "black", color_target_env))


#---------------------------------------------#
#-------------Saving the figure---------------#
#---------------------------------------------#

print(p1)

ggsave(here("figures", "main_fig1.svg"), p1,
       height=7.0, width=7.4, device=svglite::svglite)

ggsave(here("figures", "main_fig1.png"), p1,
       height=7.0, width=7.4, dpi=320)

cat("Saved Figure 1 to:\n")
cat("  ", here("figures", "main_fig1.svg"), "\n")
cat("  ", here("figures", "main_fig1.png"), "\n")

# Fitness values behind each ridge, on the natural scale.
summary_terms <- terms[, c("time", "theta_i", "z_peak", "distance",
                           "peak_height", "log_fitness_i")]
summary_terms$w_peak <- exp(summary_terms$peak_height)
summary_terms$w_focal <- exp(summary_terms$log_fitness_i)
print(summary_terms)
cat("Geometric mean fitness of the focal species:", exp(log_gmf), "\n")
