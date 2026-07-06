# Script to load functions used to plot graphs
# Sam Siljee
# 3rd October 2025

# Load required packages
library(dplyr) # Data manipulation
library(stringr) # String manipulation
library(rawrr) # Load Thermo .raw files
library(biomaRt) # Uniprot/gene ID mappings

# Load custom colour scheme for plotting - make a copy for names without underscores
custom_colours <- c(
  read.csv(paste0(sub("(Graphing/).*", "\\1", input_dir), "colours.csv"))$Colour,
  read.csv(paste0(sub("(Graphing/).*", "\\1", input_dir), "colours.csv"))$Colour
)
names(custom_colours) <- c(
  read.csv(paste0(sub("(Graphing/).*", "\\1", input_dir), "colours.csv"))$Name,
  str_replace(read.csv(paste0(sub("(Graphing/).*", "\\1", input_dir), "colours.csv"))$Name, "_", " ")
)

# Function to import the data given a data type key
import_txt_data <- function(type_key, input_dir, file_list, annotations) {
  # Filter files list to only include relevant files
  filtered_file_list <- grep(file_list, pattern = paste0("_", type_key, ".txt$"), value = TRUE)
  
  # initialise data.frame
  dat <- data.frame()
  
  # Loop through and import files - add index to connect with annotations
  for (i in filtered_file_list) {
    # Save variable to use as the Spectrum file column
    Spectrum_file <- str_remove_all(gsub(paste0("_", type_key, ".txt"), ".raw", gsub(input_dir, "", i)), "/")
    
    # read in the data and bind in
    dat <- bind_rows(
      dat,
      mutate(
        vroom(i),
        "Spectrum File" = Spectrum_file
      )
    )
  }
  
  # Merge with annotations
  dat <- merge(
    dat,
    mutate(annotations, "Spectrum File" = paste0(Run, ".raw")),
    by = "Spectrum File",
    all.x = TRUE
  )
  
  # return combined data.frame
  return(dat)
}

# Function to map gene symbols
map_gene_symbols <- function(string) {
  gene_symbol_list <- character()
  for (i in 1:length(string)) {
    gene_symbol_list <- gene_symbols_map$Gene_symbol[match(string, gene_symbols_map$Protein)]
    gene_symbol_list <- gene_symbol_list[!is.na(gene_symbol_list)]
    return(paste(gene_symbol_list, collapse = "; "))
  }
}

# Get ID mapping
if (file.exists(paste0(input_dir, "gene_map.csv"))) {
  gene_map <- read.csv(paste0(input_dir, "gene_map.csv"))
} else {
  # Set up the connection to the Ensembl database
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Get UniProt Swissprot IDs for all human gene
  all_uniprot_ids <- getBM(attributes = "uniprotswissprot", mart = ensembl)
  
  # Remove empty entries and duplicates
  all_uniprot_ids <- unique(all_uniprot_ids$uniprotswissprot[all_uniprot_ids$uniprotswissprot != ""])
  
  # Get gene symbols for all used uniprot IDs
  gene_map <- getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id", "uniprotswissprot"),
    filters = "uniprotswissprot",
    values = all_uniprot_ids,
    mart = ensembl
  ) %>%
    group_by(uniprotswissprot) %>% # remove duplicates of uniprot ID
    arrange(!is.na(ensembl_gene_id)) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  # Save table for future reference
  write.csv(
    gene_map,
    paste0(input_dir, "gene_map.csv"),
    row.names = FALSE
  )
}

# Function to format t-test results
interpret_t_test <- function(test, data, metric, conditions = c("p53 KD", "Control")) {
  interpretation <- paste0(
    "(p-value ",
    round(test$p.value, digits = 3),
    ", paired t-test, p53 KD mean ",
    round(mean(data[[metric]][data$Condition == conditions[1]]), digits = 3),
    ", control mean ",
    round(mean(data[[metric]][data$Condition == conditions[2]]), digits = 3),
    ")"
  )
  return(interpretation)
}

# Function to determine peak width from a chromatogram
calculate_peak_width <- function(raw_file, peak_mass, tol, min_RT, max_RT, peak_RT_tol = NA) {
  # Create an extracted ion chromatogram
  XIC <- rawrr::readChromatogram(rawfile = raw_file, mass = peak_mass, tol = tol, type = "xic")
  
  # Extract to a data.frame
  xic_data <- data.frame(
    RT = as.numeric(XIC[[1]]$times),
    Intensity = XIC[[1]]$intensities
  )
  
  # Filter by retention times
  xic_data <- filter(xic_data, RT > min_RT & RT < max_RT)
  
  # Determine the height of the highest peak
  max_peak <- max(xic_data$Intensity)
  
  # Get RT from highest peak to filter data more strictly, TMT runs because of low intensity peak have falsely wide peak widths
  peak_RT <- xic_data$RT[which.max(xic_data$Intensity)]
  
  # Filter again to within 20s either side - if peak RT tollerance is provided
  if(!is.na(peak_RT_tol)) {
    xic_data <- filter(xic_data, RT > peak_RT + peak_RT_tol & RT < peak_RT - peak_RT_tol) 
  }
  
  # Calculate full width at half maximum
  peak_range <- xic_data %>%
    filter(Intensity > (max_peak / 2)) %>%
    .$RT %>%
    range()
  
  # Return peak width
  return(peak_range[2] - peak_range[1])
}

# Function to make a heatmap
make_custom_heatmap <- function(input_matrix,
                                condition_annotations,
                                donor_annotations,
                                filter_list = NULL,
                                column_title,
                                include_NAs = FALSE,
                                na_col = "gray",
                                annot_col = custom_colours,
                                output_dir,
                                file_name,
                                width = 8,
                                height = 5.5,
                                units = "in",
                                res = 300) {
  # Create annotations for heatmap - conditions
  heatmap_annotations <- HeatmapAnnotation(
    Condition = condition_annotations,
    Donor = donor_annotations,
    col = list(
      Condition = annot_col[unique(condition_annotations)],
      Donor = annot_col[unique(donor_annotations)]
    )
  )
  
  # Create the input to heatmap - with filtering if filter list given
  ifelse(
    is.null(filter_list), # Check if filtering or not
    # Input - all proteins
    heatmap_input <- ({
      df <- input_matrix %>%
        t() %>%
        scale() %>%
        t()
      df
    }),
    # Input - filtered for differentially abundant proteins
    heatmap_input <- ({
      df <- input_matrix %>%
        t() %>%
        scale() %>%
        t()
      df[row.names(df) %in% filter_list, ]
    })
  )
  
  # Remove NA values if include_NAs false
  if (!include_NAs) {
    heatmap_input <- na.omit(heatmap_input)
  } else { # Otherwise calculate distances
    d <- dist(heatmap_input)
    d[is.na(d)] <- 10^50
  }
  
  # Create heatmap - and save as png
  png(
    paste0(output_dir, file_name),
    width = width,
    height = height,
    units = units,
    res = res
  )
  
  draw(
    if (include_NAs) { # Draw heatmap with NAs
      Heatmap(
        matrix = heatmap_input,
        row_title = "Proteins",
        clustering_distance_rows = d,
        clustering_method_rows = "single",
        na_col = na_col,
        column_title = column_title,
        show_row_dend = FALSE,
        show_column_dend = TRUE,
        column_names_gp = gpar(fontsize = 8),
        bottom_annotation = heatmap_annotations,
        show_row_names = FALSE,
        show_heatmap_legend = FALSE
      )
    } else { # Draw heatmap without NAs
      Heatmap(
        matrix = heatmap_input,
        row_title = "Proteins",
        column_title = column_title,
        show_row_dend = FALSE,
        show_column_dend = TRUE,
        column_names_gp = gpar(fontsize = 8),
        bottom_annotation = heatmap_annotations,
        show_row_names = FALSE,
        show_heatmap_legend = FALSE
      )
    }
  )
  
  dev.off()
}

# Function to format results file for IPA
format_for_IPA <- function(dat, uniprot_ids = NULL) {
  # Go from long to wide format
  # Get list of subgroups to iterate through
  subgroups <- unique(dat$Label)
  
  # Find proteins that are in uniprot_ids but not in the dataset
  missing_proteins <- setdiff(uniprot_ids, unique(dat$Protein))
  
  # Add missing proteins to the dataset with default values
  # Create a dataframe for missing proteins across all subgroups
  if (length(missing_proteins) > 0) {
    missing_data <- expand.grid(
      Protein = missing_proteins,
      Label = subgroups,
      stringsAsFactors = FALSE
    )
    # Add default values for missing proteins
    missing_data$log2FC <- NA
    missing_data$adj.pvalue <- NA
    
    # Append missing proteins to the original dataset
    dat <- bind_rows(dat, missing_data)
  }
  
  # Initialise a blank data.frame
  merged_dat <- data.frame(Protein = unique(dat$Protein))
  
  # Loop through subgroups, renaming columns and merging with dat
  for (i in subgroups) {
    subgroup_dat <- filter(dat, Label == i) %>%
      dplyr::select(c("Protein", "log2FC", "adj.pvalue"))
    # Rename subgroup specifically
    colnames(subgroup_dat) <- c("Protein", paste(i, "log2FC"), paste(i, "adj.pvalue"))
    # Merge with other data
    merged_dat <- merge(merged_dat, subgroup_dat, by = "Protein", all.x = TRUE)
  }
  
  # Make a new column with single Uniprot IDs, taking the first entry
  merged_dat <- mutate(
    merged_dat,
    `First Uniprot ID` = str_split_i(Protein, pattern = ";", i = 1),
    `First gene symbol` = gene_symbols_map$Gene_symbol[match(`First Uniprot ID`, gene_symbols_map$Protein)],
    `Trimmed first Uniprot ID` = str_remove(`First Uniprot ID`, "-\\d+$"),
    `Trimmed first gene symbol` = gene_symbols_map$Gene_symbol[match(`Trimmed first Uniprot ID`, gene_symbols_map$Protein)],
    .before = Protein
  )
  
  return(merged_dat)
}

# Function to read in names from FASTA file
read_FASTA_names <- function(FASTA_path) {
  # Read in the data
  FASTA_lines <- readLines(FASTA_path)
  
  # Filter out lines containing sequence only
  FASTA_lines <- FASTA_lines[grep(">", FASTA_lines)]
  
  # Return extracted Uniprot IDs
  return(sub("^[^|]*\\|([^|]*)\\|.*", "\\1", FASTA_lines))
}

# Function to parse modifications string from UniProt query to extract list of modifications
parse_phosphosites <- function(x) {
  # Return NA if no modifications data
  if(!grepl(pattern = "/note=Phospho", x)) {
    return(NA)
  } else {
    # Separate the modifications
    x <- unlist(str_split(x, "MOD_RES "))
    
    # Get residue number
    res <- str_split_i(x, pattern = "; /note=", i = 1)
    
    # Get amino acid code
    aa <- ifelse(
      grepl("/note=Phosphoserine", x), "S",
      ifelse(
        grepl("/note=Phosphothreonine", x), "T",
        ifelse(
          grepl("/note=Phosphotyrosine", x), "Y", NA)))
    
    # Identify phosphosites matching expected pattern
    phos <- (aa %in% c("S", "T", "Y")) & (!is.na(as.numeric(res)))
    
    return(paste(paste0(aa[phos], res[phos]), collapse = ","))
  }
}

# Function to parse phosphosite probabilities
parse_phos_confidence <- function(x) {
  # Return NA if missing or no phospho data
  if (is.na(x) || !grepl(pattern = "(Phospho): ", x, fixed = TRUE)) {
    return(list(phosphosite = NA, probability = NA))
  } else {
    # Separate the modifications
    x <- unlist(str_split(x, "; "))
    
    # Split each modification on "(Phospho): "
    parts <- strsplit(x, "(Phospho): ", fixed = TRUE)
    
    # Get residue numbers and probability scores
    res   <- sapply(parts, `[`, 1)
    probs <- as.numeric(sapply(parts, `[`, 2))
    
    # Find highest probability
    best_prob <- max(probs)
    
    # Concatenate all sites matching probability
    best_res  <- paste(res[probs == best_prob], collapse = ";")
    
    return(list(phosphosite = best_res, probability = best_prob))
  }
}

# Function to identify protein targets from enrichment analysis
identify_enrichment_targets <- function(
    pathway_name,
    domain = "biological_process",
    Dif_IDs = DEGs,
    gene_lookup = gene_map,
    GO_annotations = GO_annotation
){
  # Get the GO term ID for your pathway of interest first
  term_id <- GO_annotations[[domain]]$annotation$id[
    GO_annotations[[domain]]$annotation$name == pathway_name
  ]
  
  # Then look up the genes for that term ID
  pathway_genes <- names(GO_annotations[[domain]]$g[[term_id]]$weights)
  
  # return table
  return(filter(gene_lookup, ensembl_gene_id %in% intersect(Dif_IDs, pathway_genes)))
}
