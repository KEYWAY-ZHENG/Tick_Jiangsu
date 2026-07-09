# =============================================================================
# Code_tick_maxent.R - MaxEnt Modeling for Tick Distribution (Simplified Version)
# Corresponds to Paper Section 2.3, 2.4
# Function: Correlation screening, ENMeval parameter tuning
# Data: Directly reads final occurrence points from Tick.csv
# Run order: Step 1
# =============================================================================

# install.packages(c("terra", "ENMeval", "sf", "maxnet"))
library(terra)
library(ENMeval)
library(sf)
library(maxnet)

species <- "Hf"#Hl
paper_reproduce <- TRUE
res_aggregate <- if (paper_reproduce) 1 else 4
occ_file <- "Hf.csv"#Hl.csv
n_replicates <- 1
env_dir <- file.path("ASC_file", species)

# -------------------------------
# 1. Read species occurrence data
# -------------------------------
if (!file.exists(occ_file)) stop("Occurrence file not found: ", occ_file)
if (!dir.exists(env_dir)) stop("Env directory not found: ", env_dir)

ticks <- read.csv(occ_file)
ticks <- ticks[ticks$species == species, ]
if (nrow(ticks) == 0) stop("No records for species ", species)

occs <- data.frame(
  longitude = ticks$longitude,
  latitude = ticks$latitude
)

message("Loaded ", nrow(occs), " occurrence records for species: ", species)

# -------------------------------
# 2. Read environmental variables
# -------------------------------
asc_files <- list.files(env_dir, full.names = TRUE, pattern = "\\.asc$")
if (length(asc_files) == 0) stop("No .asc files in ", env_dir)

env <- rast(asc_files)

if (res_aggregate > 1) {
  env <- aggregate(env, fact = res_aggregate)
}

# -------------------------------
# 3. Environmental variable correlation screening
# -------------------------------
occs_vect <- vect(occs, geom = c("longitude", "latitude"), crs = "EPSG:4326")
env_vals <- extract(env, occs_vect, ID = FALSE)
env_vals <- env_vals[complete.cases(env_vals), ]

if (ncol(env_vals) >= 2) {
  
  cors <- cor(env_vals)
  
  high_cor <- which(abs(cors) >= 0.8 & upper.tri(cors), arr.ind = TRUE)
  
  exclude <- character(0)
  
  if (nrow(high_cor) > 0) {
    exclude <- unique(colnames(cors)[high_cor[, 2]])
  }
  
  if (length(exclude) > 0) {
    keep_vars <- setdiff(names(env), exclude)
    env <- env[[keep_vars]]
    message("Correlation screening: excluded ", paste(exclude, collapse = ", "))
  } else {
    message("No highly correlated variables found.")
  }
}

# -------------------------------
# 4. ENMeval parameter settings
# -------------------------------
tune.args <- list(
  rm = seq(0.5, 4, 0.5),
  fc = c("L", "LQ", "LQH", "LQHP")
)

partition.settings <- list(orientation = "lat_lon")

# -------------------------------
# 5. Run ENMeval (single replicate)
# -------------------------------
set.seed(1)

message("Running ENMeval (single replicate)")

eval <- ENMeval::ENMevaluate(
  occs = occs,
  envs = env,
  algorithm = "maxnet",
  tune.args = tune.args,
  partitions = "block",
  partition.settings = partition.settings,
  parallel = FALSE
)

# -------------------------------
# 6. Select best model
# -------------------------------
paper_rm <- if (species == "Hf") 2 else if (species == "Hl") 1.5 else NA
paper_fc <- if (species == "Hf") "LQHP" else if (species == "Hl") "LQ" else NA

if (paper_reproduce && !is.na(paper_rm) && !is.na(paper_fc)) {
  
  idx <- which(
    abs(eval@results$rm - paper_rm) < 0.01 &
      eval@results$fc == paper_fc
  )
  
  best_idx <- if (length(idx) > 0) as.integer(idx[1]) else which.min(eval@results$AICc)
  
} else {
  best_idx <- which.min(eval@results$AICc)
}

# -------------------------------
# 7. Get prediction results
# -------------------------------
pred_obj <- eval@predictions

if (inherits(pred_obj, "SpatRaster") && nlyr(pred_obj) > 1) {
  best_pred <- pred_obj[[best_idx]]
} else {
  best_pred <- pred_obj
}

# -------------------------------
# 8. Save result files
# -------------------------------
eval_results <- eval@results
best_mod <- eval@models[[best_idx]]

best_fc <- eval_results$fc[best_idx]
best_rm <- eval_results$rm[best_idx]

write.csv(
  eval_results,
  paste0("MaxEnt_", species, "_models.csv"),
  row.names = FALSE
)

#writeRaster(
#  best_pred,
#  paste0("MaxEnt_", species, "_prediction.tif"),
#  overwrite = TRUE
#)

#saveRDS(
#  best_mod,
#  paste0("MaxEnt_", species, "_best_model.rds")
#)

message("Best model: fc=", best_fc, ", rm=", best_rm)

# -------------------------------
# 9. Generate parameter documentation file
# -------------------------------
param_text <- c(
  "# MaxEnt software manual run parameters (obtained from ENMeval tuning)",
  paste0("# Species: ", species),
  paste0("# Generated time: ", Sys.time()),
  "",
  paste0("Regularization multiplier (RM) = ", best_rm),
  paste0("Feature class (FC) = ", best_fc),
  "",
  paste0("Occurrence file: Tick.csv"),
  paste0("Environmental variables directory: ASC_file/", species, "/"),
  paste0("Environmental variables: ", paste(names(env), collapse = ", ")),
  paste0("ENMeval replicates: ", n_replicates)
)

writeLines(
  param_text,
  paste0("MaxEnt_", species, "_MaxEnt_para.txt"),
  useBytes = TRUE
)

message("Saved: MaxEnt_", species, "_MaxEnt_para.txt")
