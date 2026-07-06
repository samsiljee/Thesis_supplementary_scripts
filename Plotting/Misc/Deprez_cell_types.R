# Small script to get cell type distribution from Deprez paper

# Packages
library(dplyr) # Data manipulation

# Files and directories
input_dir <- paste0(getwd(), "/data/")
output_dir <- paste0(getwd(), "/../produced_plots/Misc/")
cells_file <- "meta.tsv"
epithelial_cells_file <- "epithelial_cells.txt"

# Load data
cells <- read.table(paste0(input_dir, cells_file), header = TRUE, sep = "\t")
epithelial_cells <- readLines(paste0(input_dir, epithelial_cells_file))

# Filter data and group by cell type for distribution comparison
cells %>%
    filter(Position == "Intermediate" & CellType %in% epithelial_cells) %>%
    group_by(CellType) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) %>%
    write.csv(paste0(output_dir, "Deprez cell type counts.csv"))
    
