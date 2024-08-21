library(SingleCellExperiment)
library(slingshot)
library(dplyr)

# Load the data from the CSV file
data <- read.csv('C:/Users/gases/Desktop/112723_slingshot.csv')

# Create a SingleCellExperiment object with a dummy matrix for counts
sce <- SingleCellExperiment(
  assays = list(counts = matrix(1, ncol = nrow(data), nrow = 1)),
  colData = data
)

# Ensure PHATE coordinates are stored in reducedDims under the name 'PHATE'
reducedDims(sce)$PHATE <- as.matrix(data[, c("PHATE_1", "PHATE_2", "PHATE_3")])

# Running Slingshot with PHATE coordinates
sce <- slingshot(sce, clusterLabels = 'kmeans_clusters_k10', reducedDim = 'PHATE', start.clus = 3)

# Extracting pseudotime
pseudotimes <- slingPseudotime(SlingshotDataSet(sce), na = TRUE)

# Handle pseudotimes data
if (is.list(pseudotimes)) {
  for (i in seq_along(pseudotimes)) {
    lineage_name <- names(pseudotimes)[i]
    data[[paste0("Pseudotime_", lineage_name)]] <- pseudotimes[[i]]
  }
} else if (is.matrix(pseudotimes)) {
  for (lineage_name in colnames(pseudotimes)) {
    data[[paste0("Pseudotime_", lineage_name)]] <- pseudotimes[, lineage_name]
  }
}

# Save the updated dataframe with pseudotimes as a new CSV
write.csv(data, "C:/Users/gases/Desktop/updated_with_pseudotime.csv", row.names = FALSE)

# Extracting the Slingshot lineages (trajectories)
lineages <- slingCurves(SlingshotDataSet(sce))

# Prepare a dataframe to hold all the lineages for exporting
lineages_df <- do.call(rbind, lapply(seq_along(lineages), function(i) {
  curve_coords <- lineages[[i]]$s
  curve_df <- as.data.frame(curve_coords)
  curve_df$Lineage <- factor(i)
  curve_df$Index = seq_len(nrow(curve_df))
  return(curve_df)
}))

# Assign column names based on the dimensionality of the data (3D in this case)
names(lineages_df) <- c("PHATE_1", "PHATE_2", "PHATE_3", "Lineage", "Index")

# Save the trajectories as a CSV for use in Python
write.csv(lineages_df, "C:/Users/gases/Desktop/slingshot_trajectories.csv", row.names = FALSE)

