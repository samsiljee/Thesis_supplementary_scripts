# Small script to plot densitometry results

# Packages
library(dplyr) # Data manipulation
library(ggplot2) # Plotting

# Files and directories
input_dir <- paste0(getwd(), "/data/")
output_dir <- paste0(getwd(), "/../produced_plots/Misc/")
densitometry_file <- "p53_KD_validation_WB.csv"
colours_file <- "colours.csv"

# Load data
densitometry <- read.csv(paste0(input_dir, densitometry_file))

# Custom colour scheme for plotting
custom_colours <- read.csv(paste0(input_dir, "../../", colours_file))$Colour
names(custom_colours) <- read.csv(paste0(input_dir, "../../", colours_file))$Name

# Make the plot
densitometry %>%
    mutate(Sample = paste(Donor, Condition)) %>%
    ggplot(aes(x = Sample, y = Densitometry.area, fill = Donor)) +
    geom_col() +
    labs(
        x = "Sample",
        y = "Densitometry area",
        title = "p53 knockdown validation"
    ) +
    theme_bw() +
    scale_fill_manual(name = "Donor", values = custom_colours) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Save the plot
ggsave(paste0(output_dir, "p53 KD validation densitometry.png"), width = 8, height = 5.5)
