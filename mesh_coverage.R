#==========================================
# sdm-spatial-coverage
# mesh-coverage-analysis
# Fran Izquierdo — April 2026
#==========================================

# OBJECTIVE:
# Exploratory analysis of how observation points are distributed across the
# triangular mesh used as the spatial domain in INLA-based SDMs (sdmTMB,
# R-INLA). The goal is to assess whether the mesh resolution is appropriate
# for the available data.

# Key concern: too many observations per triangle means the spatial random
# field cannot resolve fine-scale variability. Too few occupied triangles
# suggests the mesh is finer than the data can support.

# Two mesh resolutions (fine vs coarse) are compared to illustrate how
# triangle size affects the observation-to-triangle distribution.

# INPUTS (input/data/):
#   observations/df_subset_anchovy.rds  — X, Y in km (UTM Zone 33N), Year
#   mesh/mesh_inla_35km.rds             — fine mesh   (MaxEdge ~ 37.5 km)
#   mesh/mesh_inla_75km.rds             — coarse mesh (MaxEdge ~ 75 km)

# OUTPUTS (output/mesh_coverage_output/):
#   mesh_coverage_<label>.png     — histogram + CDF + map
#   metrics_<label>.csv           — coverage metrics
#   diagnostic_<label>.txt        — heuristic interpretation
#   mesh_coverage_comparison.csv  — side-by-side table
#   mesh_coverage_comparison.png  — comparison bar chart

# NOTE: Thresholds values and metrics are purely heuristics.

# NOTE: The anchovy occurrences and grids used here are toy examples
# intended to illustrate the workflow. They are not meant to represent a
# validated species distribution modelling setup.
#==========================================

library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(rnaturalearth)
library(INLA)

# 0. Set dirs ------------------------------------------------------------------

path_obs  <- "input/data/observations"
path_mesh <- "input/data/mesh"
out_dir   <- "output/mesh_coverage_output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

crs_km <- "+proj=utm +zone=33 +units=km +ellps=WGS84 +datum=WGS84 +no_defs"

# 1. Load observations ---------------------------------------------------------

df_obs    <- readRDS(file.path(path_obs, "df_subset_anchovy.rds"))
coords_df <- as.data.frame(df_obs[, c("X", "Y")])

cat("Observations loaded:", nrow(coords_df), "\n")
cat("Years:", paste(sort(unique(df_obs$Year)), collapse = ", "), "\n")

# 2. Coastline -----------------------------------------------------------------

bbox_wgs84 <- st_bbox(c(xmin = 5, xmax = 18, ymin = 34, ymax = 47),
                      crs = st_crs(4326))
coast_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid() %>%
  st_crop(bbox_wgs84) %>%
  st_transform(crs_km)

# 3. Mesh coverage fn ----------------------------------------------------------

run_mesh_coverage <- function(mesh_inla, mesh_label, coords_df,
                               out_dir, coast_sf) {

  cat("\n=== Mesh coverage:", mesh_label, "===\n")

  n_nodes <- mesh_inla$n

  # Build triangle polygons (coordinates in km)
  tri_idx  <- mesh_inla$graph$tv
  node_xy  <- mesh_inla$loc[, 1:2]
  n_tri    <- nrow(tri_idx)

  tri_polys <- lapply(seq_len(n_tri), function(i) {
    pts <- node_xy[tri_idx[i, ], ]
    pts <- rbind(pts, pts[1, ])
    st_polygon(list(pts))
  })

  tri_sf <- st_sf(
    triangle_id = seq_len(n_tri),
    geometry    = st_sfc(tri_polys, crs = crs_km)
  )

  # Domain area from inner boundary (consistent across mesh resolutions).
  # All triangles include the outer extension buffer, which scales with max.edge
  # and would give different areas for fine vs coarse meshes.
  domain_km2 <- tryCatch({
    seg        <- mesh_inla$segm$int
    if (is.null(seg) || nrow(seg$idx) == 0) stop("no inner boundary")
    poly_nodes  <- seg$idx[, 1]
    poly_coords <- mesh_inla$loc[poly_nodes, 1:2]
    inner_sf    <- st_polygon(list(rbind(poly_coords, poly_coords[1, ]))) %>%
      st_sfc(crs = crs_km) %>% st_make_valid()
    round(as.numeric(st_area(inner_sf)), 0)
  }, error = function(e) {
    round(as.numeric(st_area(st_union(tri_sf))), 0)
  })

  # Count observations per triangle
  pres_sf      <- st_as_sf(coords_df, coords = c("X", "Y"), crs = crs_km)
  hits         <- st_intersects(tri_sf, pres_sf)
  tri_sf$n_pts <- lengths(hits)

  st_write(tri_sf,
           file.path(out_dir, paste0("triangles_", mesh_label, ".gpkg")),
           delete_dsn = TRUE, quiet = TRUE)

  # Metrics
  pts     <- tri_sf$n_pts
  pts_occ <- pts[pts > 0]
  n_occ   <- sum(pts > 0)

  metrics <- list(
    n_nodes    = n_nodes,
    n_tri      = n_tri,
    tri_occ    = n_occ,
    tri_empty  = n_tri - n_occ,
    pct_occ    = n_occ / n_tri,
    domain_km2 = domain_km2,
    mean_occ   = ifelse(length(pts_occ) > 0, mean(pts_occ),                         NA),
    median_occ = ifelse(length(pts_occ) > 0, median(pts_occ),                       NA),
    p95        = ifelse(length(pts_occ) > 0, as.numeric(quantile(pts_occ, 0.95)),    NA),
    max        = ifelse(length(pts)     > 0, max(pts),                               NA),
    ratio      = ifelse(n_occ           > 0, sum(pts) / n_occ,                       NA)
  )

  write.csv(as.data.frame(metrics),
            file.path(out_dir, paste0("metrics_", mesh_label, ".csv")),
            row.names = FALSE)

  diag_txt <- paste0(
    "=== MESH SPATIAL COVERAGE — ", mesh_label, " ===\n",
    "NOTE: Thresholds are heuristic and defined empirically by the authors.\n",
    "There are no universally accepted values — use these as a starting point.\n\n",
    "Mesh nodes      : ", n_nodes, "\n",
    "Total triangles : ", n_tri, "\n",
    "Occupied        : ", n_occ,
    " (", round(100 * metrics$pct_occ, 2), "%)\n",
    "Domain area     : ", format(domain_km2, big.mark = ","), " km²\n\n",
    "Ratio (obs per occupied triangle):\n",
    "  < 1.5  → mesh too coarse\n",
    "  1.5–3  → good balance\n",
    "  > 3    → mesh too fine\n",
    "  → Your value: ", round(metrics$ratio, 2), "\n\n",
    "P95 (95th percentile per occupied triangle):\n",
    "  < 5    → well distributed\n",
    "  5–8    → moderate clustering\n",
    "  > 8    → high clustering\n",
    "  → Your value: ", round(metrics$p95, 2), "\n\n",
    "% Occupied triangles:\n",
    "  < 10%  → many empty triangles\n",
    "  10–30% → reasonable\n",
    "  > 30%  → mesh relatively coarse\n",
    "  → Your value: ", round(100 * metrics$pct_occ, 2), "%"
  )

  cat(diag_txt)
  writeLines(diag_txt,
             file.path(out_dir, paste0("diagnostic_", mesh_label, ".txt")))

  # Plots -----------------------------------------------------------------------

  pct_label <- round(100 * (n_occ / n_tri), 2)

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
    labs(title    = "Sampling points per triangle",
         subtitle = paste0("Occupied triangles: ", n_occ, " of ", n_tri,
                           " (", pct_label, "%)  |  Nodes: ", n_nodes),
         x = "Points per triangle",
         y = "Count (triangles)") +
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
                           " pts  |  Mesh: ", mesh_label),
         x = "Points per triangle",
         y = "Cumulative proportion (%)") +
    theme_light() +
    theme(panel.grid.minor = element_blank(),
          plot.title       = element_text(face = "bold", size = 12))

  bbox   <- st_bbox(tri_sf)
  xlim_p <- c(bbox["xmin"], bbox["xmax"])
  ylim_p <- c(bbox["ymin"], bbox["ymax"])

  p_map <- ggplot() +
    geom_sf(data  = tri_sf, aes(fill = n_pts),
            color = "grey75", linewidth = 0.15) +
    geom_sf(data  = coast_sf, fill = "grey85",
            color = "grey60", linewidth = 0.2) +
    geom_point(data  = coords_df, aes(x = X, y = Y),
               color = "firebrick", size = 0.6, alpha = 0.8) +
    scale_fill_viridis_c(option = "mako", direction = -1, name = "N. points") +
    coord_sf(xlim = xlim_p, ylim = ylim_p, datum = NULL) +
    labs(title    = "Spatial distribution",
         subtitle = paste0("Ratio: ", round(metrics$ratio, 2),
                           " | Domain: ", format(domain_km2, big.mark = ","), " km²",
                           " | Occupied: ", round(100 * metrics$pct_occ, 2), "%",
                           " | Nodes: ", n_nodes),
         x = "X (km)", y = "Y (km)") +
    theme_light() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 13))

  fig_final <- p_hist + p_cdf + p_map +
    plot_layout(design = "AACC\nBBCC") +
    plot_annotation(
      title    = paste0("Mesh spatial coverage — ", mesh_label),
      subtitle = "Anchovy observations over INLA mesh triangles",
      theme = theme(
        plot.title    = element_text(size = 16, face = "bold",  hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 15))
      )
    )

  print(fig_final)
  ggsave(file.path(out_dir, paste0("mesh_coverage_", mesh_label, ".png")),
         fig_final, width = 16, height = 11, dpi = 300)
  cat("Plot saved:", mesh_label, "\n")

  invisible(list(metrics = metrics, tri_sf = tri_sf))
}

# 4. Load meshes and run -------------------------------------------------------

# Load mesh (spatial domain)
mesh_fine   <- readRDS(file.path(path_mesh, "mesh_inla_35km.rds"))
mesh_coarse <- readRDS(file.path(path_mesh, "mesh_inla_75km.rds"))

# Run mesh coverage function
diag_fine   <- run_mesh_coverage(mesh_fine,   "35km", coords_df, out_dir, coast_sf)
diag_coarse <- run_mesh_coverage(mesh_coarse, "75km", coords_df, out_dir, coast_sf)

# 5. Comparison table ----------------------------------------------------------

comp_table <- data.frame(
  Mesh       = c("35km", "75km"),
  Nodes      = c(diag_fine$metrics$n_nodes,    diag_coarse$metrics$n_nodes),
  Triangles  = c(diag_fine$metrics$n_tri,       diag_coarse$metrics$n_tri),
  Tri_occ    = c(diag_fine$metrics$tri_occ,     diag_coarse$metrics$tri_occ),
  Pct_occ    = round(c(diag_fine$metrics$pct_occ,  diag_coarse$metrics$pct_occ) * 100, 2),
  Domain_km2 = c(diag_fine$metrics$domain_km2, diag_coarse$metrics$domain_km2),
  Ratio      = round(c(diag_fine$metrics$ratio, diag_coarse$metrics$ratio), 2),
  P95        = round(c(diag_fine$metrics$p95,   diag_coarse$metrics$p95),   2),
  Max        = c(diag_fine$metrics$max,         diag_coarse$metrics$max)
)

print(comp_table)
write.csv(comp_table,
          file.path(out_dir, "mesh_coverage_comparison.csv"),
          row.names = FALSE)

# 6. Comparison plot -----------------------------------------------------------

comp_long <- pivot_longer(comp_table,
                          cols      = c(Pct_occ, Ratio, P95),
                          names_to  = "Metric",
                          values_to = "Value")
comp_long$Metric <- recode(comp_long$Metric,
                            Pct_occ = "% Occupied triangles",
                            Ratio   = "Ratio (obs/occ. triangle)",
                            P95     = "P95 (pts per triangle)")

p_comp <- ggplot(comp_long, aes(x = Mesh, y = Value, fill = Mesh)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(values = c("35km" = "#41b6c4", "75km" = "#225ea8")) +
  labs(title    = "Mesh resolution comparison",
       subtitle = "Fine mesh (35km) vs coarse mesh (75km) — how triangle size affects observation coverage",
       x = NULL, y = "Value") +
  theme_light() +
  theme(plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        strip.text    = element_text(face = "bold"))

print(p_comp)
ggsave(file.path(out_dir, "mesh_coverage_comparison.png"),
       p_comp, width = 10, height = 5, dpi = 300)

cat("\nScript complete. Outputs saved to:", out_dir, "\n")
