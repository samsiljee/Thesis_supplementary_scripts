# Script to select a suitable peptide for peak width calculations
# Sam Siljee
# 8th October 2025

# PeptideGroups and PSMs data loaded from Final_proteomics_LFQ_QC.Rmd script

# Get a list of Peptides present in all runs
common_peptides <- PeptideGroups %>%
  group_by(`Annotated Sequence`) %>%
  summarize(n_runs = n_distinct(Run)) %>%
  filter(n_runs == n_distinct(PeptideGroups$Run)) %>%
  pull(`Annotated Sequence`)

# Filter PSMs
arranged_peptides <- PSMs %>%
  filter(`Annotated Sequence` %in% common_peptides) %>%
  group_by(Run, `Annotated Sequence`) %>%
  summarise(Abundance = sum(`Precursor Abundance`)) %>% # bias towards multiple PSMs
  group_by(`Annotated Sequence`) %>%
  summarise(Abundance = min(Abundance)) %>% # check a good intensity in the least abundant run
  arrange(desc(Abundance))

# Get top peptide
top_peptide <- as.character(arranged_peptides[1,1])

# Get mass
top_masses <- PSMs$`m/z [Da]`[PSMs$`Annotated Sequence` == top_peptide]
top_RTs <- PSMs$`RT [min]`[PSMs$`Annotated Sequence` == top_peptide]

# Plot distributions
hist(top_masses)
hist(top_RTs)

# Get top mass - using mode (rounded values)
top_mass <- names(sort(-table(round(top_masses, digits = 4))))[1] %>% as.numeric()

# Histogram of RTs with outliers removed
hist(top_RTs[top_RTs > 26 & top_RTs < 41])

# get ranges
RT_ranges <- range(top_RTs[top_RTs > 26 & top_RTs < 41])
RT_mean <- mean(top_RTs[top_RTs > 26 & top_RTs < 41])

# Use the peptide: [R].VATVSLPR.[S]
# RT around 39.1 minutes
# m/z: 421.7584
# This peptide belongs to pig trypsin