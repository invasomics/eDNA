# ============================================
# LITTLE BLUE PENGUIN DNA ANALYSIS 
# Follows: https://benjjneb.github.io/dada2/tutorial_1_8.html
# ============================================

# ---- Packages ----
suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
  library(readxl)
  library(ggplot2)
  library(vegan)
  library(DESeq2)
  library(ggeasy)
  library(patchwork)
})

# ---- Config ----
WD <- "~/Desktop/University/Waikato/FINAL DATA ANALYSIS FILES/Little Blue Data Analysis"
setwd(WD)
theme_set(theme_bw())

# Central palettes
pal <- list(
  season = c("spring" = "#1F78B4", "summer" = "#33A02C"),
  compass = c("east" = "#FF7F00", "west" = "#6A3D9A", "both" = "#6A3D9A"),
  terrain = c("Aquatic" = "#1F78B4", "Terrestrial" = "#33A02C", "Both" = "#B2B200"),
  location = c(
    "accommodation"       = "orange",
    "cable_bay"           = "forestgreen",
    "caretakers_cottage"  = "yellow",
    "gun_emplacement"     = "#D5006D",
    "north_point"         = "cornflowerblue",
    "nursery"             = "red",
    "wharf"               = "purple",
    "workshop_paddock"    = "darkblue",
    "PCR_Control"         = "black",
    "Control_Ex"          = "grey"
  )
)

# ---- Helpers ----
factor_clean <- function(x, levels = NULL) {
  x <- as.factor(x)
  if (!is.null(levels)) x <- factor(x, levels = levels)
  droplevels(x)
}

pie_plot <- function(df, fill_col, y = "Proportion", title = "", palette = NULL) {
  ggplot(df, aes(x = "", y = !!sym(y), fill = !!sym(fill_col))) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    (if (!is.null(palette)) scale_fill_manual(values = palette) else NULL) +
    labs(title = title, x = NULL, y = NULL, fill = tools::toTitleCase(fill_col)) +
    theme_void() + theme(plot.title = element_text(hjust = 0.5))
}

bar_prop_plot <- function(df, x, y = "Proportion", title = "", palette = NULL, y_label = "Proportion") {
  ggplot(df, aes(x = !!sym(x), y = !!sym(y), fill = !!sym(x))) +
    geom_col() +
    (if (!is.null(palette)) scale_fill_manual(values = palette) else NULL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    theme_minimal() +
    labs(title = title, x = tools::toTitleCase(x), y = y_label)
}

# ============================================
# 1) Load data
# ============================================
# --- Robust ASV count reader ---
read_asv_counts <- function(path) {
  df <- read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Ensure first column is ASV IDs
  if (!"ASV" %in% names(df)) names(df)[1] <- "ASV"
  
  # Fix sample column names (no NA/blank, unique, underscores not dots)
  cn <- names(df)
  missing <- which(is.na(cn) | cn == "")
  if (length(missing)) cn[missing] <- paste0("sample_", missing)
  cn <- make.names(cn, unique = TRUE)       # valid, unique
  cn <- gsub("\\.", "_", cn)                 # underscores
  names(df) <- cn
  
  # Move ASV IDs to rownames
  if (any(is.na(df$ASV) | df$ASV == "")) {
    stop("ASV column has empty/NA IDs; please fix the input file's first column.")
  }
  df$ASV <- make.unique(df$ASV)              # guard against duplicate IDs
  df <- tibble::column_to_rownames(df, "ASV")
  
  as.matrix(df)
}

asv_mat <- read_asv_counts("ASVs_counts_blues.tsv")

tax_mat <- read.csv("MergedTaxonomy_blues.csv", header = TRUE, stringsAsFactors = FALSE)
if ("X" %in% colnames(tax_mat)) {
  tax_mat <- tax_mat %>%
    filter(X != "" & !is.na(X)) %>%
    distinct(X, .keep_all = TRUE) %>%
    tibble::column_to_rownames("X")
}
tax_mat <- as.matrix(tax_mat)

samples_df <- read_excel("LB_penguinDNA_metadata_new.xlsx") %>%
  distinct(samples.name, .keep_all = TRUE) %>%
  tibble::column_to_rownames("samples.name")

# Align
shared_taxa    <- intersect(rownames(asv_mat), rownames(tax_mat))
shared_samples <- intersect(colnames(asv_mat), rownames(samples_df))

asv_mat     <- asv_mat[shared_taxa, shared_samples, drop = FALSE]
tax_mat     <- tax_mat[shared_taxa, , drop = FALSE]
samples_df  <- samples_df[shared_samples, , drop = FALSE]

ps <- phyloseq(
  otu_table(asv_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(samples_df)
)
ps

# ============================================
# 2) Host & Environmental Proportions 
# ============================================
# Host (Eudyptula minor)
host_asvs <- taxa_names(ps)[which(tax_table(ps)[, "Species"] == "Eudyptula_minor")]
if (!length(host_asvs)) {
  message("No exact host match; trying partial...")
  host_asvs <- taxa_names(ps)[which(grepl("Eudyptula", tax_table(ps)[, "Species"], TRUE))]
}
if (!length(host_asvs)) stop("No host ASVs found. Check tax_table Species labels.")

total_reads <- sample_sums(ps)
host_reads  <- colSums(otu_table(ps)[host_asvs, , drop = FALSE])
host_prop   <- host_reads / total_reads

host_df <- tibble(
  Sample          = names(total_reads),
  Total_Reads     = as.numeric(total_reads),
  Host_Reads      = as.numeric(host_reads[match(names(total_reads), names(host_reads))]),
  Host_Proportion = as.numeric(host_prop[match(names(total_reads), names(host_prop))])
)

# Host ASV share
asv_summary <- tibble(
  Total_ASVs      = ntaxa(ps),
  Host_ASVs       = length(host_asvs),
  NonHost_ASVs    = ntaxa(ps) - length(host_asvs),
  Prop_Host_ASVs  = length(host_asvs) / ntaxa(ps),
  Prop_NonHost_ASVs = 1 - (length(host_asvs) / ntaxa(ps))
)
asv_summary

# Diet species list (edit as needed)
diet_species <- c(
  "Hyporhamphus_ihi","Engraulis_japonicus","Stolephorus_sp.",
  "Forsterygion_flavonigrum","Thyrsites_atun","Pseudophycis_bachus"
)
diet_asvs <- taxa_names(ps)[which(tax_table(ps)[, "Species"] %in% diet_species)]
diet_reads <- if (length(diet_asvs)) colSums(otu_table(ps)[diet_asvs, , drop = FALSE]) else rep(0, nsamples(ps))
diet_prop  <- diet_reads / total_reads

host_df <- host_df %>%
  mutate(Diet_Reads = as.numeric(diet_reads[match(Sample, names(diet_reads))]),
         Diet_Proportion = as.numeric(diet_prop[match(Sample, names(diet_prop))]))
head(host_df)

# ============================================
# 3) Remove host/human etc., explore terrain. 
# ============================================
ps_nohomo <- subset_taxa(ps, !(Genus %in% c("Eudyptula","Gambusia","Chalinolobus","Homo","Poecilia","Pygoscelis")))
ps_nohomo <- subset_taxa(ps_nohomo, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# --- Terrain summaries from tax_table(ps) ---

# 0) Work on the host-filtered object you already made
stopifnot(exists("ps_nohomo"))

# 1) Totals per ASV across all samples
otu_df <- as(otu_table(ps_nohomo), "matrix") |> as.data.frame() |> tibble::rownames_to_column("ASV")
otu_totals <- otu_df |>
  dplyr::mutate(TotalReads = rowSums(dplyr::across(-ASV))) |>
  dplyr::select(ASV, TotalReads)

# 2) Join taxonomy (has Terrain) to totals
tax_df <- as.data.frame(tax_table(ps_nohomo)) |> tibble::rownames_to_column("ASV")
by_asv <- dplyr::left_join(tax_df, otu_totals, by = "ASV")

# --- A) ASV counts per Terrain (unweighted) ---
terrain_asv <- by_asv |>
  dplyr::filter(!is.na(Terrain), Terrain != "") |>
  dplyr::count(Terrain, name = "ASV_Count") |>
  dplyr::mutate(Prop_ASVs = ASV_Count / sum(ASV_Count))

terrain_asv

# --- B) Read-weighted proportions per Terrain ---
terrain_reads <- by_asv |>
  dplyr::filter(!is.na(Terrain), Terrain != "") |>
  dplyr::group_by(Terrain) |>
  dplyr::summarise(TotalReads = sum(TotalReads), .groups = "drop") |>
  dplyr::mutate(Prop_Reads = TotalReads / sum(TotalReads))

terrain_reads

# ============================================
# 4) Coverage / library sizes ###exploratory
# ============================================
sample_sum_df <- tibble(ReadCount = sample_sums(ps))
ggplot(sample_sum_df, aes(x = ReadCount)) +
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of Sample Sequencing Depth") + xlab("Read counts") +
  theme(axis.title.y = element_blank())
write.csv(sample_sums(ps), "ITS_Reads_Per_Sample.csv")

df <- as(sample_data(ps), "data.frame") %>%
  mutate(LibrarySize = sample_sums(ps),
         Index = rank(LibrarySize, ties.method = "first"))

# Library size by season
ggplot(df, aes(Index, LibrarySize, colour = season)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = pal$season) +
  scale_y_continuous(trans = "sqrt") +
  theme_minimal() + labs(x = NULL, y = "Library Size")

# Library size by location
ggplot(df, aes(Index, LibrarySize, colour = location)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = pal$location) +
  scale_y_continuous(trans = "sqrt") +
  theme_minimal() + labs(x = NULL, y = "Library Size")

# Library size by compass
ggplot(df, aes(Index, LibrarySize, colour = compass)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = pal$compass) +
  scale_y_continuous(trans = "sqrt") +
  theme_minimal() + labs(x = NULL, y = "Library Size")

# ============================================
# 5) Taxonomic summaries & bar plots (after host removal) ###generates data for figure 2A & B, S3, and table 2
# ============================================
plot_bar(ps_nohomo, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack")

# Terrain-specific subsets
ps_nohomo_terrestrial <- subset_taxa(ps_nohomo, Terrain == "Terrestrial")
ps_nohomo_aquatic     <- subset_taxa(ps_nohomo, Terrain == "Aquatic")
ps_nohomo_both        <- subset_taxa(ps_nohomo, Terrain == "Both")
ntaxa(ps_nohomo); ntaxa(ps_nohomo_terrestrial); ntaxa(ps_nohomo_aquatic); ntaxa(ps_nohomo_both)

# Abundance by Terrain (reads)
otu_tax <- as.data.frame(otu_table(ps_nohomo)) %>%
  rownames_to_column("ASV") %>%
  left_join(as.data.frame(tax_table(ps_nohomo)) %>% rownames_to_column("ASV"), by = "ASV")

terrain_abundance <- otu_tax %>%
  filter(!is.na(Terrain), Terrain != "") %>%
  group_by(Terrain) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") %>%
  mutate(TotalAbundance = rowSums(across(where(is.numeric))),
         Proportion = TotalAbundance / sum(TotalAbundance)) %>%
  select(Terrain, TotalAbundance, Proportion)

bar_prop_plot(terrain_abundance, "Terrain",
              title = "Proportion of Total Reads by Terrain", palette = pal$terrain, y_label = "Proportion of Reads")
pie_plot(terrain_abundance, "Terrain", title = "Reads by Terrain", palette = pal$terrain)

# Species richness by Terrain
otu_long <- as.data.frame(otu_table(ps_nohomo)) %>%
  rownames_to_column("ASV") %>%
  pivot_longer(-ASV, names_to = "SampleID", values_to = "Abundance")

tax_df    <- as.data.frame(tax_table(ps_nohomo)) %>% rownames_to_column("ASV")
samp_df   <- as.data.frame(sample_data(ps_nohomo)) %>% rownames_to_column("SampleID")

otu_tax_meta <- otu_long %>%
  left_join(tax_df, by = "ASV") %>%
  left_join(samp_df, by = "SampleID") %>%
  filter(!is.na(Terrain), Terrain != "", !is.na(Species), Species != "", !is.na(Phylum), Phylum != "", Abundance > 0)

species_richness <- otu_tax_meta %>% group_by(Terrain) %>% summarise(SpeciesRichness = n_distinct(Species), .groups = "drop")
ggplot(species_richness, aes(Terrain, SpeciesRichness, fill = Terrain)) +
  geom_col() + scale_fill_manual(values = pal$terrain) +
  theme_minimal() + labs(title = "Species Richness by Terrain", x = "Terrain", y = "Unique Species")

# Phylum proportions by Terrain
phylum_prop_sorted <- otu_tax_meta %>%
  group_by(Terrain, Phylum) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  group_by(Terrain) %>%
  mutate(Proportion = TotalAbundance / sum(TotalAbundance)) %>%
  arrange(Terrain, desc(Proportion)) %>% ungroup()

ggplot(phylum_prop_sorted, aes(Terrain, Proportion, fill = Phylum)) +
  geom_col(position = "stack") + theme_minimal() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Proportion of Phyla by Terrain", x = "Terrain", y = "Proportion of Total Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================
# 6) Normalisation (DESeq2 poscounts) on ps_nohomo
# ============================================
ps_nohomo <- prune_taxa(taxa_sums(ps_nohomo) > 0, ps_nohomo)
ps_nohomo <- prune_samples(sample_sums(ps_nohomo) > 0, ps_nohomo)

sdf <- as(sample_data(ps_nohomo), "data.frame")
sdf$location <- factor_clean(sdf$location)
sdf$season   <- factor_clean(sdf$season)
sample_data(ps_nohomo)$location <- sdf$location
sample_data(ps_nohomo)$season   <- sdf$season

diagdds <- phyloseq_to_deseq2(ps_nohomo, ~ location + season)
diagdds <- estimateSizeFactors(diagdds, type = "poscounts")
norm_counts <- counts(diagdds, normalized = TRUE)

ps_nohomo_norm <- phyloseq(
  otu_table(norm_counts, taxa_are_rows = TRUE),
  tax_table(ps_nohomo),
  sample_data(ps_nohomo)
)

# ============================================
# 7) Top 20 Species (relative, without host)
# ============================================
stopifnot("Species" %in% colnames(tax_table(ps_nohomo)))

# Species label per ASV
species_vec <- as.character(tax_table(ps_nohomo)[, "Species"])
names(species_vec) <- taxa_names(ps_nohomo)

# Total reads per ASV
asv_reads <- taxa_sums(ps_nohomo)

# Aggregate reads to species (drop NA/blank species)
top20_by_reads <- tibble(
  ASV     = names(asv_reads),
  Species = species_vec[names(asv_reads)],
  Reads   = as.numeric(asv_reads)
) %>%
  filter(!is.na(Species), Species != "") %>%
  group_by(Species) %>%
  summarise(TotalReads = sum(Reads), .groups = "drop") %>%
  mutate(PropOfAllReads = TotalReads / sum(TotalReads)) %>%
  arrange(desc(TotalReads)) %>%
  slice_head(n = 20)

top20_by_reads  # <- nice table

# Subset phyloseq to ASVs belonging to these species
keep_asvs_top20_reads <- taxa_names(ps_nohomo)[species_vec %in% top20_by_reads$Species]
ps_top20_species_reads <- prune_taxa(keep_asvs_top20_reads, ps_nohomo)
ps_top20_species_reads

# ============================================
# 8) Alpha diversity (Shannon) ### generates figure 3
# ============================================
# Compute Shannon diversity
alpha_div <- estimate_richness(ps_nohomo, measures = "Shannon") %>%
  tibble::rownames_to_column("samplename")

# Extract metadata
geo_df <- as(sample_data(ps_nohomo), "data.frame") %>%
  tibble::rownames_to_column("samplename")

# Merge diversity with metadata
geo_data <- geo_df %>%
  left_join(alpha_div, by = "samplename") %>%
  mutate(Shannon = as.numeric(Shannon))

# Build long frame for season + location
geo_data_long <- bind_rows(
  geo_data %>% select(samplename, Shannon, group = season)   %>% filter(!is.na(group)),
  geo_data %>% select(samplename, Shannon, group = location) %>% filter(!is.na(group))
)

# Define the desired order
status_groups <- c("spring", "summer")
location_groups <- c(
  "caretakers_cottage",
  "gun_emplacement",
  "north_point",
  "nursery",
  "wharf",
  "workshop_paddock",
  "accommodation",
  "cable_bay"
)

# Final x-axis order: status → custom locations
group_levels <- c(status_groups, location_groups)
geo_data_long$group <- factor(geo_data_long$group, levels = group_levels)

# Define color palettes
status_colors <- c("spring"="#1F78B4", "summer"="#33A02C")
location_colors <- c(
  "accommodation"="orange",
  "cable_bay"="forestgreen",
  "caretakers_cottage"="yellow",
  "gun_emplacement"="#D5006D",
  "north_point"="cornflowerblue",
  "nursery"="red",
  "wharf"="purple",
  "workshop_paddock"="darkblue"
)
combined_colors <- c(status_colors, location_colors)

# Final combined plot
ggplot(geo_data_long, aes(x = group, y = Shannon, color = group)) +
  geom_point(size = 4, alpha = 0.7, position = position_jitter(width = 0.2)) +
  scale_color_manual(values = combined_colors, na.translate = FALSE) +
  theme_minimal() +
  labs(x = "Season & Location", y = "Alpha Diversity (Shannon)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


# ============================================
# 9) Ordination (PCoA on Bray) — normalised. ###plots figure 4
# ============================================
sample_data(ps_nohomo_norm)$compass <- factor(sample_data(ps_nohomo_norm)$compass, levels = c("east","west","both"))
ord <- ordinate(ps_nohomo_norm, method = "MDS", distance = "bray", na.rm = TRUE)
eigs <- ord$values$Eigenvalues

plot_pcoa <- function(color_var, palette) {
  p12 <- plot_ordination(ps_nohomo_norm, ord, color = color_var) +
    geom_point(size = 4) +
    scale_color_manual(values = palette) +
    coord_fixed(sqrt(eigs[2] / eigs[1])) +
    xlim(-1,1) + ylim(-1,1) +
    ggtitle(paste0(color_var, "; Axis 1 vs 2")) +
    ggeasy::easy_center_title() + theme_bw() + theme(legend.position = "right") +
    stat_ellipse(aes_string(group = color_var))
  p13 <- plot_ordination(ps_nohomo_norm, ord, axes = c(1,3), color = color_var) +
    geom_point(size = 4) +
    scale_color_manual(values = palette) +
    coord_fixed(sqrt(eigs[3] / eigs[1])) +
    xlim(-1,1) + ylim(-1,1) +
    ggtitle(paste0(color_var, "; Axis 1 vs 3")) +
    ggeasy::easy_center_title() + theme_bw() + theme(legend.position = "right") +
    stat_ellipse(aes_string(group = color_var))
  list(p12 = p12, p13 = p13)
}

p_loc   <- plot_pcoa("location", pal$location)
p_season <- plot_pcoa("season", pal$season)
p_comp  <- plot_pcoa("compass", pal$compass)

# Compact combo (no legends)
(p_season$p12 + p_comp$p12 + p_loc$p12) + plot_layout(ncol = 2) & theme(legend.position = "none")

# ============================================
# 10) PERMANOVA + beta dispersion (normalised) ###generates figure S5 A, B, C
# ============================================
# Remove empty samples (still possible after normalisation)
table(sample_sums(ps_nohomo_norm) == 0)
ps_nohomo_norm_nonzero <- prune_samples(sample_sums(ps_nohomo_norm) > 0, ps_nohomo_norm)

# Distance matrix (Bray-Curtis)
dist.matrix <- phyloseq::distance(ps_nohomo_norm_nonzero, method = "bray")
samples_df <- data.frame(sample_data(ps_nohomo_norm_nonzero))

# PERMANOVA & Beta Dispersion: season
perm_season <- vegan::adonis2(dist.matrix ~ season, data = samples_df)
print(perm_season)

bd_season <- vegan::betadisper(dist.matrix, samples_df$season)
print(bd_season)

boxplot(bd_season, 
        col = status_colors,
        names = c("Spring", "Summer"),
        xlab = "season",
        ylab = "Distance to Centroid")

permutest(bd_season, pairwise = TRUE)

# PERMANOVA & Beta Dispersion: Location
# with gun emplacement
perm_location <- vegan::adonis2(dist.matrix ~ location, data = samples_df)
print(perm_location)

bd_location <- vegan::betadisper(dist.matrix, samples_df$location)
print(bd_location)

mycol_location <- c(
  "accommodation"       = "orange", 
  "cable_bay"           = "forestgreen", 
  "caretakers_cottage"  = "yellow", 
  "gun_emplacement"     = "#D5006D", 
  "north_point"         = "cornflowerblue", 
  "nursery"             = "red", 
  "wharf"               = "purple", 
  "workshop_paddock"    = "darkblue"
)

boxplot(bd_location,
        col = mycol_location,
        las = 2,
        xlab = "Location",
        ylab = "Distance to Centroid")

permutest(bd_location, pairwise = TRUE)

# without gun emplacement (n=2)
samples_df_filtered <- samples_df[samples_df$location != "gun_emplacement", ]
dist.matrix_full <- as.matrix(dist.matrix)
dist.matrix_subset <- dist.matrix_full[rownames(samples_df_filtered), rownames(samples_df_filtered)]
dist.matrix_filtered <- as.dist(dist.matrix_subset)

perm_location <- vegan::adonis2(dist.matrix_filtered ~ location, data = samples_df_filtered)
print(perm_location)

bd_location <- vegan::betadisper(dist.matrix_filtered, samples_df_filtered$location)
print(bd_location)

mycol_location <- c(
  "accommodation"       = "orange", 
  "cable_bay"           = "forestgreen", 
  "caretakers_cottage"  = "yellow", 
  "north_point"         = "cornflowerblue", 
  "nursery"             = "red", 
  "wharf"               = "purple", 
  "workshop_paddock"    = "darkblue"
)

boxplot(bd_location,
        col = mycol_location,
        las = 2,
        xlab = "Location",
        ylab = "Distance to Centroid")

permutest(bd_location, pairwise = TRUE)

# PERMANOVA & Beta Dispersion: compass
perm_compass <- vegan::adonis2(dist.matrix ~ compass, data = samples_df)
print(perm_compass)

bd_compass <- vegan::betadisper(dist.matrix, samples_df$compass)
print(bd_compass)

mycol_compass <- c(
  "east" = "#E31A1C",
  "west" = "#FF7F00",
  "both" = "#6A3D9A"
)

boxplot(bd_compass,
        col = mycol_compass,
        names = c("East", "West"),
        xlab = "compass",
        ylab = "Distance to Centroid")

permutest(bd_compass, pairwise = TRUE)

# ============================================
# 11) CLAM tests (normalised)  ###generates figure 5 A-D, tables S1 & S2
# ============================================
clam_or_stop <- function(ps_obj, group_var, groups_keep = NULL, alpha = 0.01, specialization = 2/3,
                         outfile_prefix = "specialist") {
  stopifnot(group_var %in% sample_variables(ps_obj))
  if (!is.null(groups_keep)) {
    ps_obj <- subset_samples(ps_obj, !!as.name(group_var) %in% groups_keep)
    ps_obj <- prune_samples(sample_sums(ps_obj) > 0, ps_obj)
  }
  mat <- t(as(otu_table(ps_obj), "matrix"))  # species in columns for clamtest
  grp <- as.data.frame(sample_data(ps_obj))[[group_var]]
  if (!requireNamespace("vegan", quietly = TRUE)) stop("vegan needed for clamtest")
  res <- vegan::clamtest(mat, grp, alpha = alpha, specialization = specialization)
  
  # Extract groups from names like "Specialist_spring"
  labs <- unique(res$Classes)
  tax_df <- as.data.frame(tax_table(ps_obj))
  for (lv in labs[grepl("^Specialist_", labs)]) {
    key <- sub("^Specialist_", "", lv)
    asvs <- res$Species[res$Classes == lv]
    out  <- tax_df[rownames(tax_df) %in% asvs, , drop = FALSE]
    write.csv(out, paste0(outfile_prefix, "_", key, "_normalised.csv"))
  }
  invisible(res)
}

# season specialists
season_res <- clam_or_stop(ps_nohomo_norm, "season", groups_keep = c("spring","summer"),
                           outfile_prefix = "specialist_season")
# compass specialists (east/west only)
compass_res <- clam_or_stop(ps_nohomo_norm, "compass", groups_keep = c("east","west"),
                            outfile_prefix = "specialist_compass")

# ---- Load & summarise a CLAM export (example: season spring/summer) ----
summarise_specialist_csv <- function(csv_path, title_prefix) {
  if (!file.exists(csv_path)) return(invisible(NULL))
  df <- read.csv(csv_path, stringsAsFactors = FALSE) %>%
    rownames_to_column("ASV") %>%
    mutate(across(c(Terrain, Phylum, Species), ~ as.character(.)))
  # Terrain (ASV count)
  terr_counts <- df %>% filter(!is.na(Terrain), Terrain != "") %>%
    count(Terrain) %>% mutate(Proportion = n / sum(n))
  print(terr_counts); print(
    pie_plot(terr_counts, "Terrain",
             title = paste(title_prefix, "ASVs by Terrain"), palette = pal$terrain)
  )
  
  # Phylum (ASV count)
  phy_counts <- df %>% filter(!is.na(Phylum), Phylum != "") %>% count(Phylum)
  print(phy_counts); print(
    ggplot(phy_counts, aes(reorder(Phylum, -n), n, fill = Phylum)) +
      geom_col() + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = paste(title_prefix, "ASVs by Phylum"), x = "Phylum", y = "ASV Count")
  )
  
  # ASV vs species richness per phylum
  asv_counts <- df %>% filter(!is.na(Phylum), Phylum != "") %>% count(Phylum, name = "ASV_Count")
  sp_counts  <- df %>% filter(!is.na(Phylum), Phylum != "", !is.na(Species), Species != "") %>%
    distinct(Phylum, Species) %>% count(Phylum, name = "Species_Count")
  summary_tbl <- full_join(asv_counts, sp_counts, by = "Phylum") %>% arrange(desc(ASV_Count))
  print(summary_tbl, n = Inf)
}

# Examples (run those that exist in your folder)
summarise_specialist_csv("specialist_season_spring_normalised.csv",  "Egg Laying Specialists (Spring)")
summarise_specialist_csv("specialist_season_summer_normalised.csv",  "Chick Rearing Specialists (Summer)")
summarise_specialist_csv("specialist_compass_east_normalised.csv",   "East Specialists")
summarise_specialist_csv("specialist_compass_west_normalised.csv",   "West Specialists")


