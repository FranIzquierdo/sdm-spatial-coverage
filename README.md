# sdm-spatial-coverage

Exploratory analysis of species observations (e.g., occurrences) over a selected spatial domain in an SDM context. 

See `/output/mesh_coverage_output/mesh_coverage_35km.png` and `sdm-spatial-coverage.pdf` for more info and to see the applied example.

## What it does

These scripts help to address whether species observations are well distributed across the prediction domain by computing some coverage metrics (heuristic) and comparative plots.

| Script | Spatial unit | Use case |
|--------|-------------|----------|
| `mesh_coverage.R` | Irregular triangles | `sdmTMB`, `R-INLA` |
| `grid_coverage.R` | Regular cells | `GAM`, `GLM`, `BART` |


## Directories

```text
├── input/
│   └── data/
│       ├── observations/   # occurrence or abundance records
│       ├── mesh/           # INLA mesh objects (.rds)
│       └── grid/           # prediction grids (.tif)
├── output/                 # plots and metrics
├── mesh_coverage.R         # triangle-based analysis
└── grid_coverage.R         # cell-based analysis
