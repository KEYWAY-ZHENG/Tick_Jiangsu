# ============================================================
# Calculate high suitability area using Jenks natural breaks classification
# Classify MaxEnt logistic output (0-1) into 4 classes: unsuitable, low, moderate, high
# Output: high suitability pixel count, area (rounded to integer), and percentage change
# ============================================================

library(terra)
library(classInt)

out_dir <- "results"
if(!dir.exists(out_dir)) dir.create(out_dir)

total_area_km2 <- 107200   # Total area of Jiangsu Province (km²)

# -------------------------------
# 1. Species and ASC files (no thresholds - will use Jenks natural breaks)
# -------------------------------
species_files <- list(
  "Hf" = list(
    "current" = list(file = "MaxEnt_Hf_maxent_output/future2041_2060/Hf_avg.asc"),
    "future_2021_2040" = list(file = "MaxEnt_Hf_maxent_output/future2021_2040/Hf_Future_final_avg.asc"),
    "future_2041_2060" = list(file = "MaxEnt_Hf_maxent_output/future2041_2060/Hf_Future_final_avg.asc")
  ),
  "Hl" = list(
    "current" = list(file = "MaxEnt_Hl_maxent_output/future2021_2040/Hl_avg.asc"),
    "future_2021_2040" = list(file = "MaxEnt_Hl_maxent_output/future2021_2040/Hl_Future_final_avg.asc"),
    "future_2041_2060" = list(file = "MaxEnt_Hl_maxent_output/future2041_2060/Hl_Future_final_avg.asc")
  )
)

# -------------------------------
# 2. Function to calculate high suitability area using Jenks natural breaks
# -------------------------------
calc_high_area_jenks <- function(r, total_area_km2, jenks_n_max = 6000L, seed = 42){
  
  # Exclude invalid values, keep range 0~1
  v <- values(r)
  v <- v[!is.na(v) & is.finite(v) & v >= 0 & v <= 1]
  
  n_total_pixels <- length(v)
  
  if(n_total_pixels < 10L){
    return(list(n_high = 0, n_total_pixels = n_total_pixels, area_high = 0, breaks = NA))
  }
  
  # For large datasets, sample pixels for Jenks calculation (faster)
  # Then apply breaks to all pixels for area calculation
  if(length(v) > jenks_n_max){
    set.seed(seed)
    v_jenks <- v[sample.int(length(v), jenks_n_max)]
  } else {
    v_jenks <- v
  }
  
  # Calculate Jenks natural breaks into 4 classes
  ci <- classIntervals(v_jenks, n = 4, style = "jenks", warnLargeN = FALSE)
  breaks <- ci$brks
  
  # Class boundaries:
  # Class 1 (unsuitable): [0, breaks[2])
  # Class 2 (low): [breaks[2], breaks[3])
  # Class 3 (moderate): [breaks[3], breaks[4])
  # Class 4 (high): [breaks[4], 1]
  
  # Count high suitability pixels (>= breaks[4])
  n_high <- sum(v >= breaks[4])
  
  # Calculate area
  area_high <- round((n_high / n_total_pixels) * total_area_km2)  # Round to integer
  
  return(list(
    n_high = n_high, 
    n_total_pixels = n_total_pixels, 
    area_high = area_high,
    breaks = breaks
  ))
}

# -------------------------------
# 3. Loop through species and periods
# -------------------------------
for(species in names(species_files)){
  
  cat("====================================================\n")
  cat("Processing species:", species, "\n")
  
  future_files <- species_files[[species]]
  
  results <- data.frame(
    period = character(),
    high_pixel_count = numeric(),
    total_pixels = numeric(),
    high_area_km2 = numeric(),
    pct_change = numeric(),
    jenks_break_high = numeric(),
    stringsAsFactors = FALSE
  )
  
  curr_area <- NA
  
  for(period_name in names(future_files)){
    
    asc_file <- future_files[[period_name]]$file
    
    if(!file.exists(asc_file)){
      warning("File not found: ", asc_file)
      next
    }
    
    r <- rast(asc_file)
    stats <- calc_high_area_jenks(r, total_area_km2)
    
    if(period_name == "current") curr_area <- stats$area_high
    
    pct_chg <- ifelse(period_name == "current", 0,
                      100 * (stats$area_high - curr_area)/curr_area)
    
    jenks_high_thresh <- if(!any(is.na(stats$breaks))) round(stats$breaks[4], 6) else NA
    
    results <- rbind(results, data.frame(
      period = period_name,
      high_pixel_count = stats$n_high,
      total_pixels = stats$n_total_pixels,
      high_area_km2 = stats$area_high,
      pct_change = round(pct_chg, 1),  # Percentage rounded to 1 decimal place
      jenks_break_high = jenks_high_thresh
    ))
    
    # Console output
    cat("Completed:", period_name, "\n")
    if(!any(is.na(stats$breaks))){
      cat("Jenks breaks (4 classes):", paste(round(stats$breaks, 4), collapse = ", "), "\n")
      cat("High suitability threshold (breaks[4]):", round(stats$breaks[4], 6), "\n")
    }
    cat("High suitability area:", stats$area_high, "km²\n",
        "High suitability pixel count:", stats$n_high, "\n",
        "Total valid pixels:", stats$n_total_pixels, "\n",
        "Change relative to current:", round(pct_chg,1), "%\n\n")
  }
  
  # -------------------------------
  # Save summary table for this species
  # -------------------------------
  write.csv(
    results,
    file.path(out_dir, paste0("MaxEnt_", species, "_high_area_summary.csv")),
    row.names = FALSE
  )
  
  cat("Species", species, "high suitability statistics completed ✔\n")
}
