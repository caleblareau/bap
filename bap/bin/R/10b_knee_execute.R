options(warn=-1)

suppressMessages(suppressWarnings(library(tools)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))

# I/O
args <- commandArgs(trailingOnly = FALSE)

# Parse parameters working backwards since I am keeping all of them
# to be able to source the knee call Rscript
nn <- length(args)
vals_file <- args[nn-2]
logTransform <- args[nn-1]
column_name <- args[nn]

# Grab knee call functions from figuring out where this script lives and then changing the file name
needle <- "--file="
match <- grep(needle, args)
source(normalizePath(gsub("10b_knee_execute.R", "00_knee_CL.R", gsub(needle, "", args[match]))))

# Import data supplied by the user
vec_values <- as.numeric(fread(vals_file)[[as.character(column_name)]])

# Call knee
estimated_knee_threshold <- get_density_threshold_CL(vec_values, "bead", as.logical(as.numeric(as.character(logTransform))))

# Write value to table
write.table(data.frame(estimated_knee_threshold), row.names = FALSE, col.names = FALSE, 
            quote = FALSE, sep = "", file = paste0(vals_file, "_kneeValue.txt"))


