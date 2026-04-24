#==========================================
# Grid spatial coverage
# sdm-spatial-coverage
# Fran Izquierdo — April 2026
#==========================================

# OBJECTIVE:
# Exploratory analysis of how observation points are distributed across the
# regular prediction grid used as the spatial domain in raster-based SDMs
# (GAM, GLM, BART, RF, MaxEnt...). The goal is to assess whether the grid
# resolution is appropriate for the available data.

# Key concern: too many observations per cell means the model cannot resolve
# fine-scale spatial variation.. Too few occupied cells suggests the grid 
# is finer than the data can support.

# Two grid resolutions (medium vs coarse) are compared to illustrate how
# cell size affects the observation-to-cell distribution.

# INPUTS (input/data/):
#   observations/df_subset_anchovy.rds  — X, Y in km (UTM Zone 33N), Year
#   grid/grid_30km.tif                  — medium grid (30 km cells)
#   grid/grid_70km.tif                  — coarse grid (70 km cells)

# Grid rasters must be:
#   - Projected (metres, e.g. EPSG:32633 UTM Zone 33N)
#   - Masked to the marine study domain (NA = outside domain)

# OUTPUTS (output/grid_coverage_output/):
#   grid_coverage_<label>.png     — histogram + CDF + map
#   metrics_<label>.csv           — coverage metrics
#   diagnostic_<label>.txt        — heuristic interpretation
#   grid_coverage_comparison.csv  — side-by-side table
#   grid_coverage_comparison.png  — comparison bar chart

# NOTE: Thresholds values and metrics are purely heuristics.

# NOTE: The anchovy occurrences and grids used here are toy examples
# intended to illustrate the workflow. They are not meant to represent a
# validated species distribution modelling setup.
#==========================================

library(sf)
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(rnaturalearth)

# 0. Set dirs ------------------------------------------------------------------

path_obs  <- "input/data/observations"
path_grid <- "input/data/grid"
out_dir   <- "output/grid_coverage_output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 1. Load observations ---------------------------------------------------------

df_obs    <- readRDS(file.path(path_obs, "df_subset_anchovy.rds"))
coords_df <- as.data.frame(df_obs[, c("X", "Y")])  # km

cat("Observations loaded:", nrow(coords_df), "\n")
cat("Years:", paste(sort(unique(df_obs$Year)), collapse = ", "), "\n")

# 2. Coastline -----------------------------------------------------------------

coast_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid() %>%
  st_transform(32633)  # metres, matches grid CRS

# 3. Grid coverage fn ----------------------------------------------------------

run_grid_coverage <- function(grid_rast, res_label, coords_df_km,
                               out_dir, coast_sf) {

  cat("\n=== Grid coverage:", res_label, "===\n")

  # Non-NA cells = prediction domain
  grid_sf <- as.polygons(grid_rast, dissolve = FALSE) %>%
    st_as_sf() %>%
    mutate(cell_id = row_number())

  n_cells <- nrow(grid_sf)

  # Domain area (km²)
  res_km     <- as.numeric(res(grid_rast)[1]) / 1000
  domain_km2 <- round(n_cells * res_km^2, 0)

  # Convert observation coords km → m (matching grid CRS)
  coords_m <- data.frame(X = coords_df_km$X * 1000,
                         Y = coords_df_km$Y * 1000)
  pres_sf  <- st_as_sf(coords_m, coords = c("X", "Y"), crs = st_crs(grid_sf))

  # Count observations per cell
  hits          <- st_intersects(grid_sf, pres_sf)
  grid_sf$n_pts <- lengths(hits)

  st_write(grid_sf,
           file.path(out_dir, paste0("cells_", res_label, ".gpkg")),
           delete_dsn = TRUE, quiet = TRUE)

  # Metrics
  pts     <- grid_sf$n_pts
  pts_occ <- pts[pts > 0]
  n_occ   <- sum(pts > 0)

  metrics <- list(
    n_cells     = n_cells,
    cells_occ   = n_occ,
    cells_empty = n_cells - n_occ,
    pct_occ     = n_occ / n_cells,
    domain_km2  = domain_km2,
    mean_occ    = ifelse(length(pts_occ) > 0, mean(pts_occ),                         NA),
    median_occ  = ifelse(length(pts_occ) > 0, median(pts_occ),                       NA),
    p95         = ifelse(length(pts_occ) > 0, as.numeric(quantile(pts_occ, 0.95)),    NA),
    max         = ifelse(length(pts)     > 0, max(pts),                               NA),
    ratio       = ifelse(n_occ           > 0, sum(pts) / n_occ,                       NA)
  )

  write.csv(as.data.frame(metrics),
            file.path(out_dir, paste0("metrics_", res_label, ".csv")),
            row.names = FALSE)

  diag_txt <- paste0(
    "=== GRID SPATIAL COVERAGE — ", res_label, " ===\n",
    "NOTE: Thresholds are heuristic and defined empirically by the authors.\n",
    "There are no universally accepted values — use these as a starting point.\n\n",
    "Resolution    : ", res_label, "\n",
    "Total cells   : ", n_cells, "\n",
    "Occupied      : ", n_occ,
    " (", round(100 * metrics$pct_occ, 2), "%)\n",
    "Domain area   : ", format(domain_km2, big.mark = ","), " km²\n\n",
    "Ratio (obs per occupied cell):\n",
    "  < 1.5  → grid too coarse\n",
    "  1.5–3  → good balance\n",
    "  > 3    → grid too fine\n",
    "  → Your value: ", round(metrics$ratio, 2), "\n\n",
    "P95 (95th percentile per occupied cell):\n",
    "  < 5    → well distributed\n",
    "  5–8    → moderate clustering\n",
    "  > 8    → high clustering\n",
    "  → Your value: ", round(metrics$p95, 2), "\n\n",
    "% Occupied cells:\n",
    "  < 10%  → strong extrapolation (most domain unobserved)\n",
    "  10–30% → moderate geographic coverage\n",
    "  > 30%  → good geographic coverage\n",
    "  → Your value: ", round(100 * metrics$pct_occ, 2), "%"
  )

  cat(diag_txt)
  writeLines(diag_txt,
             file.path(out_dir, paste0("diagnostic_", res_label, ".txt")))

  # Plots -----------------------------------------------------------------------

  pct_label <- round(100 * (n_occ / n_cells), 2)

  df_hist <- data.frame(n_pts = pts_occ) %>%
    mutate(category = cut(n_pts,
                          breaks = c(-Inf, 0, 3, 8, Inf),
                          labels = c("Empty", "Low (1-3)", "Moderate (4-8)", "High (>8)")))

  p_hist <- ggplot(df_hist, aes(x = n_pts, fill = category)) +
    geom_histogram(binwidth = 1, color = "white", boundary = 0.5) +
    scale_fill_manual(
      values = c("Low (1-3)"      = "#a1dab4",
                 "Moderate (4-8)" = "#41b6c4",
                 "High (>8)"      = "#225ea8"),
      name = "Density category") +
    scale_x_continuous(breaks = seq(0, max(df_hist$n_pts, na.rm = TRUE), by = 1)) +
    labs(title    = "Sampling points per cell",
         subtitle = paste0("Occupied cells: ", n_occ, " of ", n_cells,
                           " (", pct_label, "%)  |  Resolution: ", res_label),
         x = "Points per cell",
         y = "Count (cells)") +
    theme_light() +
    theme(legend.position  = "bottom",
          panel.grid.minor = element_blank(),
          plot.title       = element_text(face = "bold", size = 12))

  p_cdf <- ggplot(data.frame(n_pts = pts_occ), aes(x = n_pts)) +
    stat_ecdf(geom = "step", color = "aquamarine3", linewidth = 1) +
    geom_vline(xintercept = metrics$p95, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = seq(0, max(pts_occ, na.rm = TRUE), by = 1)) +
    scale_y_continuous(labels = scales::percent) +
    labs(title    = "Cumulative distribution",
         subtitle = paste0("P95 = ", round(metrics$p95, 2),
                           " pts  |  Resolution: ", res_label),
         x = "Points per cell",
         y = "Cumulative proportion (%)") +
    theme_light() +
    theme(panel.grid.minor = element_blank(),
          plot.title       = element_text(face = "bold", size = 12))

  bbox   <- st_bbox(grid_sf)
  xlim_p <- c(bbox["xmin"], bbox["xmax"])
  ylim_p <- c(bbox["ymin"], bbox["ymax"])

  p_map <- ggplot() +
    geom_sf(data  = grid_sf, aes(fill = n_pts),
            color = "grey80", linewidth = 0.05) +
    geom_sf(data  = coast_sf, fill = "grey85",
            color = "grey60", linewidth = 0.2) +
    geom_point(data  = coords_m, aes(x = X, y = Y),
               color = "firebrick", size = 0.5, alpha = 0.7) +
    scale_fill_viridis_c(option = "mako", direction = -1, name = "N. points") +
    coord_sf(xlim = xlim_p, ylim = ylim_p, datum = NULL) +
    labs(title    = "Spatial distribution",
         subtitle = paste0("Ratio: ", round(metrics$ratio, 2),
                           " | Domain: ", format(domain_km2, big.mark = ","), " km²",
                           " | Occupied: ", round(100 * metrics$pct_occ, 2), "%",
                           " | Cells: ", n_cells),
         x = "X (m)", y = "Y (m)") +
    theme_light() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 13))

  fig_final <- p_hist + p_cdf + p_map +
    plot_layout(design = "AACC\nBBCC") +
    plot_annotation(
      title    = paste0("Grid spatial coverage — ", res_label),
      subtitle = "Anchovy observations over prediction grid cells",
      theme = theme(
        plot.title    = element_text(size = 16, face = "bold",  hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15))
      )
    )

  print(fig_final)
  ggsave(file.path(out_dir, paste0("grid_coverage_", res_label, ".png")),
         fig_final, width = 16, height = 11, dpi = 300)
  cat("Plot saved:", res_label, "\n")

  invisible(list(metrics = metrics, grid_sf = grid_sf))
}

# 4. Load grids and run --------------------------------------------------------

# Load grid (spatial domain)
grid_30km <- rast(file.path(path_grid, "grid_30km.tif"))
grid_70km <- rast(file.path(path_grid, "grid_70km.tif"))

# Run grid coverage function
diag_30 <- run_grid_coverage(grid_30km, "30km", coords_df, out_dir, coast_sf)
diag_70 <- run_grid_coverage(grid_70km, "70km", coords_df, out_dir, coast_sf)

# 5. Comparison table ----------------------------------------------------------

comp_table <- data.frame(
  Resolution = c("30km", "70km"),
  Cells      = c(diag_30$metrics$n_cells,    diag_70$metrics$n_cells),
  Cells_occ  = c(diag_30$metrics$cells_occ,  diag_70$metrics$cells_occ),
  Pct_occ    = round(c(diag_30$metrics$pct_occ, diag_70$metrics$pct_occ) * 100, 2),
  Domain_km2 = c(diag_30$metrics$domain_km2, diag_70$metrics$domain_km2),
  Ratio      = round(c(diag_30$metrics$ratio, diag_70$metrics$ratio), 2),
  P95        = round(c(diag_30$metrics$p95,   diag_70$metrics$p95),   2),
  Max        = c(diag_30$metrics$max,         diag_70$metrics$max)
)

print(comp_table)
write.csv(comp_table,
          file.path(out_dir, "grid_coverage_comparison.csv"),
          row.names = FALSE)

# 6. Comparison plot -----------------------------------------------------------

comp_long <- pivot_longer(comp_table,
                          cols      = c(Pct_occ, Ratio, P95),
                          names_to  = "Metric",
                          values_to = "Value")
comp_long$Metric <- recode(comp_long$Metric,
                            Pct_occ = "% Occupied cells",
                            Ratio   = "Ratio (obs/occ. cell)",
                            P95     = "P95 (pts per cell)")
comp_long$Resolution <- factor(comp_long$Resolution, levels = c("30km", "70km"))

p_comp <- ggplot(comp_long, aes(x = Resolution, y = Value, fill = Resolution)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(values = c("30km" = "#41b6c4", "70km" = "#225ea8")) +
  labs(title    = "Grid resolution comparison",
       subtitle = "Medium grid (30km) vs coarse grid (70km) — how cell size affects observation coverage",
       x = NULL, y = "Value") +
  theme_light() +
  theme(plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        strip.text    = element_text(face = "bold"))

print(p_comp)
ggsave(file.path(out_dir, "grid_coverage_comparison.png"),
       p_comp, width = 10, height = 5, dpi = 300)

cat("\nScript complete. Outputs saved to:", out_dir, "\n")
