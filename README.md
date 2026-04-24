# sdm-spatial-coverage

Exploratory analysis of how well observation points cover a spatial prediction domain — for any species distribution model (SDM).

Two complementary scripts, each built around a single reusable function:

| Script | Core function | Domain | Use case |
|--------|--------------|--------|----------|
| `mesh_coverage.R` | `run_mesh_coverage()` | INLA triangular mesh | sdmTMB, R-INLA |
| `grid_coverage.R` | `run_grid_coverage()` | Regular raster grid | GAM, GLM, BART, RF, MaxEnt… |

Both answer the same question: **what fraction of the prediction domain contains at least one observation?**  
High coverage → mostly interpolation. Low coverage → mostly extrapolation.

> **Note**: The anchovy occurrences, meshes, and grids used here are toy examples to illustrate the workflow. They are not meant to represent a validated SDM setup.
> **UTM (km)**: observations (species location) and spatial domain are projected in km
> **Toy example**: European anchovy (*Engraulis encrasicolus*) presence records from the MEDITS bottom trawl survey (2014–2020) in the Central Mediterranean and Northern Adriatic Sea (UTM Zone 33N). See [Worked example](#worked-example) for details.

---

## How it works

The prediction domain is divided into spatial units — triangles (mesh) or cells (grid).  
For each unit, the script counts how many observations fall inside.  
Key metrics:

| Metric | Description |
|--------|-------------|
| **% Occupied** | Fraction of units with ≥ 1 observation |
| **Ratio** | Mean observations per occupied unit |
| **P95** | 95th percentile of observations per unit (robust clustering indicator) |
| **Domain area (km²)** | Total area of the prediction domain |

### Interpretation thresholds (heuristic)

| Metric | Range | Interpretation |
|--------|-------|----------------|
| Ratio | < 1.5 | Too coarse — loss of spatial detail |
| Ratio | 1.5–3 | Good balance |
| Ratio | > 3 | Too fine — many empty units |
| P95 | < 5 | Well distributed |
| P95 | 5–8 | Moderate clustering |
| P95 | > 8 | High clustering — check sampling bias |
| % Occupied | < 10% | Most domain unobserved — strong extrapolation |
| % Occupied | 10–30% | Moderate geographic coverage |
| % Occupied | > 30% | Good geographic coverage |

> These thresholds are exploratory guides, not hard rules.  
> Interpret metrics together, not in isolation.

---

## Repository structure

```
sdm-spatial-coverage/
├── input/
│   └── data/
│       ├── observations/          # occurrence or abundance records
│       │   └── df_subset_anchovy.rds
│       ├── mesh/                  # INLA mesh objects (.rds)
│       │   ├── mesh_inla_35km.rds
│       │   └── mesh_inla_75km.rds
│       └── grid/                  # prediction grids (.tif, masked to domain)
│           ├── grid_30km.tif
│           └── grid_70km.tif
├── output/
│   ├── mesh_coverage_output/
│   └── grid_coverage_output/
├── mesh_coverage.R
├── grid_coverage.R
└── README.md
```

---

## Input requirements

### Observations

An `.rds` file (R data frame) with at least:

| Column | Description |
|--------|-------------|
| `X` | Easting in **km** (UTM projected) |
| `Y` | Northing in **km** (UTM projected) |
| `Year` | Sampling year (used to filter the period of interest) |

### Mesh (mesh_coverage.R)

Standard INLA mesh object created with `INLA::inla.mesh.2d()`.  
Coordinates must be in km (sdmTMB convention).

Example (fine mesh):
```r
mesh_inla_35km <- inla.mesh.2d(
  loc      = coords_km,
  max.edge = c(37.5, 75),   # inner / outer edge (km)
  cutoff   = 12.5,
  offset   = c(50, 150)
)
saveRDS(mesh_inla_35km, "input/data/mesh/mesh_inla_35km.rds")
```

### Grid (grid_coverage.R)

A `.tif` raster file:
- Projected CRS in **metres** (e.g. EPSG:32633 UTM Zone 33N)
- Masked to the marine study domain — cells outside the domain must be `NA`
- Cell values are not used (the script only needs cell geometry)

Example (create grid from a domain polygon):
```r
library(terra)
domain_v  <- vect(domain_sf)           # your study area polygon (terra SpatVector)
template  <- rast(ext(domain_v), res = 15000, crs = "EPSG:32633")
values(template) <- 1
grid_15km <- mask(template, domain_v)
writeRaster(grid_15km, "input/data/grid/grid_15km.tif", overwrite = TRUE)
```

---

## Worked example

Species: **European anchovy** (*Engraulis encrasicolus*)  
Area: Central Mediterranean Sea and Northern Adriatic Sea  
Period: **2014–2020** (7 years)  
Source: MEDITS bottom trawl survey (presence records)  
Projection: UTM Zone 33N (EPSG:32633)

The example evaluates two mesh resolutions (fine vs coarse) and two grid resolutions (medium vs coarse) to illustrate how spatial coverage metrics change with the discretisation of the prediction domain.

### Meshes compared

| File | MaxEdge (inner) | Expected triangles |
|------|-----------------|--------------------|
| `mesh_inla_35km.rds` | ~37.5 km (fine) | ~500 |
| `mesh_inla_75km.rds` | ~75 km (coarse) | ~150 |

### Grids compared

| File | Cell size | Expected cells |
|------|-----------|----------------|
| `grid_30km.tif` | 30 km | medium |
| `grid_70km.tif` | 70 km | low |

---

## Outputs

Each script produces, per resolution:

- `mesh_coverage_<label>.png` / `grid_coverage_<label>.png`  
  Three-panel figure: histogram of points per unit · CDF with P95 marker · spatial map

- `metrics_<label>.csv` — all metrics in tabular form  
- `diagnostic_<label>.txt` — plain-text heuristic interpretation

Plus a cross-resolution comparison:

- `mesh_coverage_comparison.csv` / `grid_coverage_comparison.csv`  
- `mesh_coverage_comparison.png` / `grid_coverage_comparison.png`

---

## Dependencies

```r
# mesh_coverage.R
library(sf); library(ggplot2); library(dplyr); library(tidyr)
library(patchwork); library(rnaturalearth); library(INLA)

# grid_coverage.R
library(sf); library(terra); library(ggplot2); library(dplyr); library(tidyr)
library(patchwork); library(rnaturalearth)
```

Install INLA:
```r
install.packages("INLA",
  repos = c(getOption("repos"),
            INLA = "https://inla.r-inla-download.org/R/stable"))
```

---

## Notes

- **Mesh outer extension**: INLA meshes include a buffer extension around the study area. Triangles in the extension are always empty. This is expected and does not affect model performance — it inflates the total triangle count and lowers `% Occupied` compared to grid metrics.
- **This is not model validation**: coverage analysis evaluates geographic spread of observations relative to the prediction domain. It does not replace cross-validation or residual diagnostics.
- **Thresholds are exploratory**: the heuristics provided are based on practical experience and SPDE theory (Lindgren et al., 2011). Adapt them to your study system.

---

## References

Lindgren, F., Rue, H., & Lindström, J. (2011). An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. *Journal of the Royal Statistical Society B*, 73(4), 423–498.
