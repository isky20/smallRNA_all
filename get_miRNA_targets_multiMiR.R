# Load necessary libraries
library(readxl)     # For reading Excel files (not used directly here but included if needed)
library(multiMiR)   # For querying miRNA target data

# Define the function to retrieve miRNA targets
get_mirna_targets <- function(input_file, output_folder) {
  # Read CSV file containing miRNA names
  MM <- read.csv(input_file)
  mirna <- MM$miRNA.names

  # Set output directory
  setwd(output_folder)

  # Initialize a list to store results
  target <- list()

  # Loop through each miRNA and query multiMiR database
  for(i in mirna) {
    tryCatch({
      target[[i]] <- get_multimir(mirna = noquote(i), summary = TRUE)
    }, error = function(e) {
      message(paste("Error for miRNA:", i, "-", e$message))
    })
  }

  # Save each result as a CSV file
  for(k in seq_along(target)) {
    print(names(target)[k])
    write.csv(target[[k]]@data, paste0(names(target)[k], ".csv"))
  }
}

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure both input file and output folder are provided
if (length(args) != 2) {
  stop("Usage: Rscript get_mirna_targets.R <input_file> <output_folder>")
}

# Assign arguments
input_file <- args[1]
output_folder <- args[2]

# Run the function
get_mirna_targets(input_file, output_folder)
