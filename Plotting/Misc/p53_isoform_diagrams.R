# Script to make a diagram of p53 isoforms
# Sam Siljee
# 12/3/26

# ============================================================
# WHY MANUAL DATAFRAME CONSTRUCTION?
#
# Two issues affect the automatic approach (get_features +
# extract_transcripts):
#
#  1. SINGLE CHAIN: The UniProt REST API changed its JSON format
#     after ~2022. drawProteins::extract_transcripts() can no longer
#     reliably parse isoform VSP entries, so only the canonical
#     chain is returned.
#
#  2. NO DOMAINS: p53's functional regions are stored in UniProt as
#     type "REGION" (not "DOMAIN"), so draw_domains() finds nothing.
#
# SOLUTION: Build the drawProteins dataframe directly from published
# domain coordinates, and add the Δ160 isoforms (absent from UniProt)
# by truncating the Δ133 transcripts at residue 160.
#
# ISOFORM STRUCTURE (canonical numbering, isoform 1 = 393 aa):
#
#   Name            N-start   C-end   Length
#   p53α               1       393     393
#   p53β               1       331     331
#   p53γ               1       346     346
#   Δ40p53α           40       393     354
#   Δ40p53β           40       331     292
#   Δ40p53γ           40       346     307
#   Δ133p53α         133       393     261
#   Δ133p53β         133       331     199
#   Δ133p53γ         133       346     214
#   Δ160p53α         160       393     234
#   Δ160p53β         160       331     172
#   Δ160p53γ         160       346     187
#
# DOMAIN COORDINATES (canonical numbering, updated from literature):
#   TAD1:  1 –  42   Transactivation domain 1 (AD1)
#                      Wikipedia/Raj et al. 2016 (PMID 28007893): residues 1-42
#   TAD2: 43 –  63   Transactivation domain 2 (AD2)
#                      Wikipedia: residues 43-63
#   PRR:  64 –  92   Proline-rich region
#                      Consistent across sources (Walker & Levine 1996)
#   DBD:  94 – 292   DNA-binding domain
#                      Most consistent value; Joerger & Fersht 2008,
#                      Frontiers Oncol. 2020 (PMID 33194618)
#   NLS: 316 – 325   Nuclear localisation signal
#                      Wikipedia (updated from 305-322)
#   OD:  325 – 356   Oligomerisation / tetramerisation domain
#                      Jeffrey et al. 1995; Clore et al. 1995;
#                      ScienceDirect 2009 (PMID 19800327): 325-355;
#                      Wikipedia: dimerization interface 325-356
#   REG: 356 – 393   C-terminal regulatory domain
#                      PMC11270737 (2024): REG(356-393)
# ============================================================

# Setup
# Packages
library(drawProteins)
library(ggplot2)

# Set some directories
output_dir <- paste0(getwd(), "/Plotting/produced_plots/Misc/")

# ============================================================
# Define isoforms
# ============================================================

# Each isoform defined by N-terminal start and C-terminal end
# in canonical (full-length p53α) coordinates.
isoforms <- data.frame(
  name = c(
    "p53\u03b1 (full-length)",
    "p53\u03b2",
    "p53\u03b3",
    "\u039440p53\u03b1",
    "\u039440p53\u03b2",
    "\u039440p53\u03b3",
    "\u0394133p53\u03b1",
    "\u0394133p53\u03b2",
    "\u0394133p53\u03b3",
    "\u0394160p53\u03b1",
    "\u0394160p53\u03b2",
    "\u0394160p53\u03b3",
    "\u0394246p53\u03b1",
    "\u0394246p53\u03b2",
    "\u0394246p53\u03b3"
  ),
  nstart = c(1, 1, 1, 40, 40, 40, 133, 133, 133, 160, 160, 160, 246, 246, 246),
  cend = c(393, 331, 346, 393, 331, 346, 393, 331, 346, 393, 331, 346, 393, 331, 346),
  stringsAsFactors = FALSE
)
# Assign order in reverse so that p53α (row 1) gets the highest
# order value and therefore appears at the TOP of the plot.
# drawProteins places order=1 at the bottom, order=N at the top.
isoforms$order <- rev(seq_len(nrow(isoforms)))
n_chains <- nrow(isoforms)

# ============================================================
# Define domain coordinates (canonical numbering)
# ============================================================

domains <- data.frame(
  description = c("TAD1", "TAD2", "PXXP", "DBD", "NLS", "OD", "REG"),
  begin = c(1, 43, 64, 94, 316, 325, 356),
  end = c(42, 63, 92, 292, 325, 356, 393),
  stringsAsFactors = FALSE
)

# ============================================================
# Build drawProteins-compatible dataframe
# ============================================================
# All coordinates are kept in canonical (full-length p53α) numbering
# throughout, so truncated isoforms naturally align with the
# full-length protein on the x-axis.

build_df <- function(isoforms_df, domains_df) {
  rows <- list()

  for (i in seq_len(nrow(isoforms_df))) {
    iso <- isoforms_df[i, ]
    ns <- iso$nstart
    ce <- iso$cend
    ord <- iso$order
    nm <- iso$name

    # CHAIN row – use canonical coordinates so all isoforms align
    # on the same x-axis (truncated isoforms start further right).
    rows[[length(rows) + 1]] <- data.frame(
      type = "CHAIN",
      description = nm,
      begin = ns, # canonical start, not 1
      end = ce, # canonical end
      order = ord,
      entryName = nm,
      stringsAsFactors = FALSE
    )

    # DOMAIN rows – clipped to isoform boundaries, canonical coords
    for (j in seq_len(nrow(domains_df))) {
      d <- domains_df[j, ]
      clipped_beg <- max(d$begin, ns)
      clipped_end <- min(d$end, ce)

      if (clipped_beg > clipped_end) next # domain outside isoform

      rows[[length(rows) + 1]] <- data.frame(
        type = "DOMAIN",
        description = d$description,
        begin = clipped_beg, # canonical coords
        end = clipped_end,
        order = ord,
        entryName = nm,
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}

p53_plot_df <- build_df(isoforms, domains)

# ============================================================
# Draw the diagram
# ============================================================

# Canvas
p <- draw_canvas(p53_plot_df)

# Protein chains (backbone bars)
# label_chains = FALSE because canonical coordinates push the
# draw_chains labels far to the right; we add them manually below.
p <- draw_chains(p, p53_plot_df,
  label_chains = FALSE,
  label_size   = 3,
  fill         = "lightgrey",
  outline      = "grey40"
)

# Add isoform name labels to the left of the chains.
# The x-axis is expanded to -120 to make room, but breaks are
# manually set to only show positive protein coordinates.
p <- p +
  geom_text(
    data        = isoforms,
    mapping     = aes(x = 470, y = order, label = name),
    hjust       = 1,
    size        = 3,
    inherit.aes = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, 470),
    breaks = c(1, 50, 100, 150, 200, 250, 300, 350, 393)
  ) +
  coord_cartesian(clip = "off")

# Domains – drawn with geom_rect for full colour control
# (draw_domains() is not used because p53's UniProt annotations
#  use type "REGION" rather than "DOMAIN", returning nothing)
domain_rows <- p53_plot_df[p53_plot_df$type == "DOMAIN", ]

domain_colours <- c(
  "TAD1" = "#E63946",
  "TAD2" = "#F4A261",
  "PXXP" = "#F1C40F",
  "DBD"  = "#2A9D8F",
  "NLS"  = "#A8DADC",
  "OD"   = "#457B9D",
  "REG"  = "#9B59B6"
)

p <- p +
  geom_rect(
    data = domain_rows,
    mapping = aes(
      xmin  = begin,
      xmax  = end,
      ymin  = order - 0.25,
      ymax  = order + 0.25,
      fill  = description
    ),
    colour = "grey40",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    name   = "Domain",
    values = domain_colours,
    limits = c("TAD1", "TAD2", "PXXP", "DBD", "NLS", "OD", "REG")
  )

# Theme and labels
p <- p +
  labs(
    #title = "Human p53 / TP53 isoforms",
    x = "Amino acid position"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 8, colour = "grey30"),
    axis.title.x       = element_text(size = 10),
    axis.text.y        = element_blank(), # remove isoform order numbers
    axis.ticks.y       = element_blank(), # remove corresponding ticks
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "bottom",
    legend.title       = element_text(size = 9),
    legend.text        = element_text(size = 8)
  ) +
  guides(fill = guide_legend(nrow = 1))

# Save the plot
ggsave(paste0(output_dir, "p53 isoforms diagram.png"), plot = p, width = 8, height = 5.5)
