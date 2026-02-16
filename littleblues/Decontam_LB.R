# ============================================
# LITTLE BLUE PENGUIN DNA ANALYSIS
# Including decontam filtering and prokaryote/viral ASV removal
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
  library(decontam)
})

if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = TRUE)
}


# ---- Config ----
theme_set(theme_bw())

# Central palettes
pal <- list(
  season = c("egg_laying" = "#1F78B4", "chick_rearing" = "#33A02C"),
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

# ---- Helper Functions ----
clean_sample_names <- function(x) {
  x <- gsub("\\.", "_", x)
  x <- gsub("-", "_", x)
  x
}

read_asv_counts <- function(path) {
  df <- read.delim(path, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(nrow(df) > 0, ncol(df) > 1)
  names(df)[1] <- "ASV"
  stopifnot(!any(is.na(df$ASV) | df$ASV == ""))
  df$ASV <- make.unique(df$ASV)
  df <- tibble::column_to_rownames(df, "ASV")
  mat <- as.matrix(df)
  storage.mode(mat) <- "numeric"
  stopifnot(nrow(mat) > 0, ncol(mat) > 0)
  mat
}

# ---- Load Data ----
asv_mat <- read_asv_counts("Dec25_newASVs_counts.tsv")
colnames(asv_mat) <- clean_sample_names(colnames(asv_mat))

tax_df <- read.csv("Dec25_newTaxonomy.csv", stringsAsFactors = FALSE)
names(tax_df) <- make.names(names(tax_df), unique = TRUE)
tax_df <- tax_df %>% filter(ASV != "" & !is.na(ASV)) %>% distinct(ASV, .keep_all = TRUE) %>%
  tibble::column_to_rownames("ASV")
tax_mat <- as.matrix(tax_df)

samples_df <- readxl::read_excel("metadata.xlsx", sheet = 1) %>%
  dplyr::select(
    plate_no,
    samples.name,
    tube_number,
    sample_name,
    location,
    compass,
    collection_month
  ) %>%
  dplyr::mutate(
    samples.name = clean_sample_names(samples.name),
    compass = tolower(compass),
    location = tolower(location)
  ) %>%
  tibble::column_to_rownames("samples.name")

# ---- Create phyloseq object ----
ps <- phyloseq(
  otu_table(asv_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(samples_df)
)
ps

# ============================================
Host, Prokaryote/Virus & Decontam Filtering
# ============================================

#Define taxa to remove ----
host_genera <- c("Eudyptula", "Gambusia", "Chalinolobus", "Homo", "Poecilia", "Pygoscelis")
prok_viral_phyla <- c("Pseudomonadota", "Actinomycetota", "Bacteroidota",
                      "Bacillota", "Peploviricota", "Uroviricota")

#Remove host and unwanted microbial/viral taxa ----
ps_nohomo <- subset_taxa(
  ps,
  !(Genus %in% host_genera) &
    !(Phylum %in% prok_viral_phyla) &
    !is.na(Phylum) &
    !Phylum %in% c("", "uncharacterized")
)

ps_nohomo

#QC: identify samples with zero reads after taxonomic filtering ----
empty_after_tax_filter <- sample_names(ps_nohomo)[sample_sums(ps_nohomo) == 0]

cat("Samples with zero reads after taxonomic filtering (will be dropped):\n")
print(empty_after_tax_filter)
cat("Total:", length(empty_after_tax_filter), "\n")

#Explicitly remove empty samples (prevents decontam warnings) ----
ps_nohomo <- prune_samples(sample_sums(ps_nohomo) > 0, ps_nohomo)

#Mark negative controls for decontam ----
neg_samples <- c("Control_Ext", "Control_PCR")
sample_data(ps_nohomo)$is.neg <-
  sample_data(ps_nohomo)$sample_name %in% neg_samples

# Quick check
table(sample_data(ps_nohomo)$is.neg)

#Identify contaminants using decontam ----
contam_df <- isContaminant(
  ps_nohomo,
  method = "prevalence",
  neg = "is.neg"
)

# Inspect contaminant counts
table(contam_df$contaminant)

#Remove contaminant ASVs ----
contaminants <- rownames(contam_df)[contam_df$contaminant]
ps_clean <- prune_taxa(!taxa_names(ps_nohomo) %in% contaminants, ps_nohomo)

#Inspect cleaned object ----
ps_clean
ps_nohomo <- ps_clean

library(phyloseq)

#Count ASVs before and after filtering 
total_ASVs_before <- ntaxa(ps)       # total ASVs in original object
total_ASVs_after  <- ntaxa(ps_nohomo) # total ASVs after filtering
ASVs_lost <- total_ASVs_before - total_ASVs_after

#Count unique genera before and after 
genera_before <- length(unique(tax_table(ps)[, "Genus"]))
genera_after  <- length(unique(tax_table(ps_nohomo)[, "Genus"]))
genera_lost   <- genera_before - genera_after

#Count unique phyla before and after
phyla_before <- length(unique(tax_table(ps)[, "Phylum"]))
phyla_after  <- length(unique(tax_table(ps_nohomo)[, "Phylum"]))
phyla_lost   <- phyla_before - phyla_after

#Summary table 
summary_table <- data.frame(
  Level = c("ASVs", "Genera", "Phyla"),
  Before = c(total_ASVs_before, genera_before, phyla_before),
  After  = c(total_ASVs_after, genera_after, phyla_after),
  Lost   = c(ASVs_lost, genera_lost, phyla_lost)
)

print(summary_table)

# ============================================
# Host Proportions
# ============================================
host_asvs <- taxa_names(ps)[which(tax_table(ps)[, "Species"] == "Eudyptula_minor")]
if (!length(host_asvs)) {
  host_asvs <- taxa_names(ps)[which(grepl("Eudyptula", tax_table(ps)[, "Species"], TRUE))]
}
if (!length(host_asvs)) stop("No host ASVs found.")

total_reads <- sample_sums(ps)
host_reads  <- colSums(otu_table(ps)[host_asvs, , drop = FALSE])
host_prop   <- host_reads / total_reads
host_prop
host_df <- tibble(
  Sample          = names(total_reads),
  Total_Reads     = as.numeric(total_reads),
  Host_Reads      = as.numeric(host_reads[match(names(total_reads), names(host_reads))]),
  Host_Proportion = as.numeric(host_prop[match(names(total_reads), names(host_prop))])
)
ggplot(host_df, aes(x = Sample, y = Host_Proportion)) +
  geom_col(fill = "#1F78B4") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Proportion of Host (Eudyptula minor) Reads per Sample",
    x = "Sample",
    y = "Host Read Proportion"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
overall_total_reads <- sum(host_df$Total_Reads, na.rm = TRUE)
overall_host_reads  <- sum(host_df$Host_Reads,  na.rm = TRUE)

overall_host_prop <- overall_host_reads / overall_total_reads
print(host_df, n = Inf)
overall_host_prop
#Extract ASV table with taxonomy from ps_nohomo 
# Convert phyloseq to a tidy data frame
by_asv <- as.data.frame(otu_table(ps_nohomo)) |> 
  rownames_to_column(var = "ASV") |>
  pivot_longer(-ASV, names_to = "Sample", values_to = "Abundance") |>
  left_join(as.data.frame(tax_table(ps_nohomo)) |> 
              rownames_to_column(var = "ASV"), by = "ASV") |>
  group_by(ASV, Phylum) |>
  summarise(TotalReads = sum(Abundance), .groups = "drop")

#Proportion of ASVs per Phylum 
phyla_asv <- by_asv |>
  group_by(Phylum) |>
  summarise(ASV_Count = n(), .groups = "drop") |>
  mutate(Prop_ASVs = ASV_Count / sum(ASV_Count)) |>
  arrange(desc(Prop_ASVs))

#Proportion of reads per Phylum (raw counts)
phyla_reads <- by_asv |>
  group_by(Phylum) |>
  summarise(TotalReads = sum(TotalReads), .groups = "drop") |>
  mutate(Prop_Reads = TotalReads / sum(TotalReads)) |>
  arrange(desc(Prop_Reads))

#Print results
print("Proportion of ASVs per Phylum:")
print(phyla_asv)

print("Proportion of raw reads per Phylum:")
print(phyla_reads)

# ============================================
#Terrain summaries & Taxonomic summaries (after host removal)
# ============================================

#Convert OTU table to data frame
otu_df <- as(otu_table(ps_nohomo), "matrix") |>
  as.data.frame() |>
  rownames_to_column("ASV")

#Get total reads per ASV across all samples
otu_totals <- otu_df |>
  mutate(TotalReads = rowSums(across(-ASV))) |>
  select(ASV, TotalReads)

#Taxonomy table with ASV as column
tax_df <- as.data.frame(tax_table(ps_nohomo)) |>
  rownames_to_column("ASV")

# Combine taxonomy and read totals
by_asv <- left_join(tax_df, otu_totals, by = "ASV") |>
  filter(!is.na(Terrain), Terrain != "")

#Proportions of ASVs by terrain
terrain_asv <- by_asv |>
  group_by(Terrain) |>
  summarise(ASV_Count = n(), .groups = "drop") |>
  mutate(Prop_ASVs = ASV_Count / sum(ASV_Count))

#Proportions of reads by terrain
terrain_reads <- by_asv |>
  group_by(Terrain) |>
  summarise(TotalReads = sum(TotalReads), .groups = "drop") |>
  mutate(Prop_Reads = TotalReads / sum(TotalReads))

#Proportion of ASVs per phylum within each terrain
terrain_phylum_asv <- by_asv |>
  filter(!is.na(Phylum), Phylum != "") |>
  group_by(Terrain, Phylum) |>
  summarise(ASV_Count = n(), .groups = "drop_last") |>
  mutate(
    Terrain_Total_ASVs = sum(ASV_Count),
    Prop_ASVs = ASV_Count / Terrain_Total_ASVs
  ) |>
  ungroup()

#Print results
print("Proportion of ASVs per terrain:")
print(terrain_asv)

print("Proportion of reads per terrain:")
print(terrain_reads)

print("Proportion of ASVs per phylum within each terrain:")
print(terrain_phylum_asv)

#######################################################
# Top 20 Species Relative Abundance Analysis 
#######################################################
#Remove zero-count samples
ps_nohomo <- prune_samples(sample_sums(ps_nohomo) > 0, ps_nohomo)

#Keep only ASVs with valid Species assignment
ps_species_only <- subset_taxa(ps_nohomo, !is.na(Species) & Species != "")

#Aggregate ASVs to Species level
ps_species_glom <- tax_glom(ps_species_only, taxrank = "Species", NArm = TRUE)

#Transform to relative abundance per sample
ps_species_rel <- transform_sample_counts(ps_species_glom, function(x) x / sum(x))

#Get total abundance per species
species_totals <- taxa_sums(ps_species_rel)

#Select Top 20 species
top20_species <- names(sort(species_totals, decreasing = TRUE))[1:20]

#Prune phyloseq object to Top 20
ps_top20_species <- prune_taxa(top20_species, ps_species_rel)

#Create Top 20 table with Phylum
top20_species_df <- data.frame(
  Species = tax_table(ps_top20_species)[, "Species"],
  Phylum = tax_table(ps_top20_species)[, "Phylum"],
  Total_Reads = taxa_sums(ps_top20_species),
  Proportion = round(taxa_sums(ps_top20_species) / sum(taxa_sums(ps_top20_species)), 4)
) %>%
  arrange(desc(Proportion))

print(top20_species_df)

#Plot Top 20 species per sample colored by Phylum
plot_bar(ps_top20_species, x = "Sample", fill = "Phylum") +
  labs(
    title = "Top 20 Species by Relative Abundance (Exclude NA Species)",
    x = "Sample",
    y = "Relative Abundance"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "bottom")


# ============================================
Normalisation (DESeq2 poscounts)
# ============================================

#Report samples & taxa to be removed
zero_samples <- sample_names(ps_nohomo)[sample_sums(ps_nohomo) == 0]
zero_taxa    <- taxa_names(ps_nohomo)[taxa_sums(ps_nohomo) == 0]

cat("Samples removed (zero total reads):", length(zero_samples), "\n")
if (length(zero_samples) > 0) print(zero_samples)

cat("\nTaxa removed (zero total reads):", length(zero_taxa), "\n")
if (length(zero_taxa) > 0) print(zero_taxa)

#Prune empty samples & taxa
ps_nohomo <- prune_taxa(taxa_sums(ps_nohomo) > 0, ps_nohomo)
ps_nohomo <- prune_samples(sample_sums(ps_nohomo) > 0, ps_nohomo)

#Ensure factors for DESeq2 compatibility
sdf <- as(sample_data(ps_nohomo), "data.frame")
sdf$location <- factor(sdf$location)
sdf$season   <- factor(sdf$season)

sample_data(ps_nohomo)$location <- sdf$location
sample_data(ps_nohomo)$season   <- sdf$season

#DESeq2 size-factor normalisation 
diagdds <- phyloseq_to_deseq2(ps_nohomo, ~ 1)
diagdds <- estimateSizeFactors(diagdds, type = "poscounts")

norm_counts <- counts(diagdds, normalized = TRUE)

#Rebuild normalised phyloseq object 
ps_nohomo_norm <- phyloseq(
  otu_table(norm_counts, taxa_are_rows = TRUE),
  tax_table(ps_nohomo),
  sample_data(ps_nohomo)
)

ps_nohomo_norm
zero_samples
zero_taxa
# ============================================
# Alpha Diversity (Shannon)
# ============================================

#Calculate Shannon diversity
alpha_div <- estimate_richness(ps_nohomo, measures = "Shannon")
alpha_div$samplename <- rownames(alpha_div)

#Extract metadata and merge
geo_data <- data.frame(sample_data(ps_nohomo))
geo_data$samplename <- rownames(geo_data)
geo_data <- merge(geo_data, alpha_div, by = "samplename")

#Ensure Shannon is numeric
geo_data$Shannon <- as.numeric(geo_data$Shannon)

#Remove controls and NAs
geo_data <- geo_data %>%
  filter(
    !is.na(location),
    !location %in% "control",
    !is.na(compass),
    !compass %in% "control"
  )

#Clean & set factor levels
geo_data <- geo_data %>%
  mutate(
    location = factor(trimws(location),
                      levels = sort(c(
                        "accommodation", "cable_bay", "caretakers_cottage",
                        "gun_emplacement", "north_point", "nursery",
                        "wharf", "workshop_paddock"
                      ))),
    compass = factor(tolower(trimws(compass)),
                     levels = c("east", "west"))
  )

#Define color palettes
location_colors <- c(
  "accommodation"      = "orange",
  "cable_bay"          = "forestgreen",
  "caretakers_cottage" = "yellow",
  "gun_emplacement"    = "#D5006D",
  "north_point"        = "cornflowerblue",
  "nursery"            = "red",
  "wharf"              = "purple",
  "workshop_paddock"   = "darkblue"
)

compass_colors <- c("east" = "red", "west" = "orange")

#Function to plot a factor variable
plot_alpha <- function(df, xvar, xlab, colors) {
  ggplot(df, aes_string(x = xvar, y = "Shannon", color = xvar)) +
    geom_point(size = 4, alpha = 0.6, position = position_jitter(width = 0.15)) +
    stat_summary(fun = mean, geom = "point", shape = 95, size = 6, color = "lightgray") +
    scale_color_manual(values = colors) +
    theme_minimal() +
    labs(x = xlab, y = "Alpha Diversity (Shannon)") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}
#Plot by each factor
p_location <- plot_alpha(geo_data, "location", "Location", location_colors)
p_compass  <- plot_alpha(geo_data, "compass", "Compass Direction", compass_colors)

#Display plots combined
combined_plot <- p_location / p_compass +
  plot_annotation(
    title = "Alpha Diversity (Shannon) by Factor",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

combined_plot
p_location

# ============================================
#Statistical tests (Wilcoxon / Kruskal-Wallis)
# ============================================
alpha_test <- function(df, factor_var) {
  factor_levels <- unique(df[[factor_var]])
  
  if(length(factor_levels) == 2) {
    test_result <- wilcox.test(Shannon ~ get(factor_var), data = df)
    test_name <- "Wilcoxon rank-sum test"
  } else if(length(factor_levels) > 2) {
    test_result <- kruskal.test(Shannon ~ get(factor_var), data = df)
    test_name <- "Kruskal-Wallis test"
  } else {
    stop("Factor has fewer than 2 levels.")
  }
  
  cat("\n", test_name, "for", factor_var, ":\n")
  print(test_result)
  return(test_result)
}

alpha_test(geo_data, "location")
alpha_test(geo_data, "compass")

# ============================================
# Ordination (PCoA on Bray)
# ============================================
ord_nohomo <- ordinate(ps_nohomo_norm, method = "MDS", distance = "bray", na.rm = TRUE)
eigen_vals <- ord_nohomo$values$Eigenvalues

#PCoA plots: LOCATION
p_loc_1_2 <- plot_ordination(ps_nohomo_norm, ord_nohomo, color = "location") +
  geom_point(size = 4) +
  scale_color_manual(values = location_colors) +
  coord_fixed(sqrt(eigen_vals[2] / eigen_vals[1])) +
  xlim(-1, 1) + ylim(-1, 1) +
  ggtitle("Location; Axis 1 vs 2") +
  ggeasy::easy_center_title() +
  theme_bw() + theme(legend.position = "right") +
  stat_ellipse(aes(group = location))


#PCoA plots: COMPASS
p_compass_1_2 <- plot_ordination(ps_nohomo_norm, ord_nohomo, color = "compass") +
  geom_point(size = 4) +
  scale_color_manual(values = compass_colors) +
  coord_fixed(sqrt(eigen_vals[2] / eigen_vals[1])) +
  xlim(-1, 1) + ylim(-1, 1) +
  ggtitle("Compass; Axis 1 vs 2") +
  ggeasy::easy_center_title() +
  theme_bw() + theme(legend.position = "right") +
  stat_ellipse(aes(group = compass))

#PERMANOVA & Beta Dispersion
#Remove empty samples (after normalization)
ps_nohomo_norm_nonzero <- prune_samples(
  sample_sums(ps_nohomo_norm) > 0, 
  ps_nohomo_norm
)

dist.matrix <- phyloseq::distance(ps_nohomo_norm_nonzero, method = "bray")
samples_df <- data.frame(sample_data(ps_nohomo_norm_nonzero))

#Ensure factors are clean
samples_df <- samples_df %>%
  mutate(
    location = factor(trimws(location),
                      levels = sort(c(
                        "accommodation", "cable_bay", "caretakers_cottage",
                        "gun_emplacement", "north_point", "nursery",
                        "wharf", "workshop_paddock"
                      ))),
    compass = factor(tolower(trimws(compass)),
                     levels = c("east", "west")),
    collection_month = factor(trimws(collection_month),
                              levels = c("Oct", "Nov", "Dec"))
  )

month_colors <- c("Oct" = "#D89000", "Nov" = "#A3A500", "Dec" = "#00BF7D")

#--------------------------------------------------------
# PERMANOVA & Beta Dispersion: location
#--------------------------------------------------------
perm_location <- vegan::adonis2(dist.matrix ~ location, data = samples_df)
bd_location <- vegan::betadisper(dist.matrix, samples_df$location)
boxplot(bd_location, col = location_colors, las = 2,
        xlab = "Location", ylab = "Distance to Centroid")
permutest(bd_location, pairwise = TRUE)

# Extract R²
R2_location <- perm_location$R2[1]
cat("Location PERMANOVA R²:", R2_location, "\n")

#--------------------------------------------------------
# PERMANOVA & Beta Dispersion: compass
#--------------------------------------------------------
perm_compass <- vegan::adonis2(dist.matrix ~ compass, data = samples_df)
bd_compass <- vegan::betadisper(dist.matrix, samples_df$compass)
boxplot(bd_compass, col = compass_colors, names = levels(samples_df$compass),
        xlab = "Compass", ylab = "Distance to Centroid")
permutest(bd_compass, pairwise = TRUE)

# Extract R²
R2_compass <- perm_compass$R2[1]
cat("Compass PERMANOVA R²:", R2_compass, "\n")

#--------------------------------------------------------
# PERMANOVA & Beta Dispersion: collection_month
#--------------------------------------------------------
perm_month <- vegan::adonis2(dist.matrix ~ collection_month, data = samples_df)
bd_month <- vegan::betadisper(dist.matrix, samples_df$collection_month)
boxplot(bd_month, names = levels(samples_df$collection_month),
        col = month_colors,
        xlab = "Collection Month", ylab = "Distance to Centroid")
permutest(bd_month, pairwise = TRUE)

# Extract R²
R2_month <- perm_month$R2[1]
cat("Collection Month PERMANOVA R²:", R2_month, "\n")


# =============================================
# FULL CLAM WORKFLOW: Tests, Plots, Printed Stats, Charts
# =============================================

month_groups <- list(
  "Oct" = c("Nov", "Dec"),
  "Nov" = c("Oct", "Dec"),
  "Dec" = c("Oct", "Nov")
)

alpha_val <- 0.01       # significance threshold
spec_val  <- 2/3        # specialization threshold

for (focal_month in names(month_groups)) {
  
  other_months <- month_groups[[focal_month]]
  other_label <- paste(other_months, collapse = "_")   
  comparison_label <- paste0(focal_month, "_vs_", other_label)
  
  #Subset phyloseq object
  month_subset <- subset_samples(
    ps_nohomo_norm,
    collection_month %in% c(focal_month, other_months)
  )
  
  #Drop zero-sum ASVs
  month_subset <- prune_taxa(taxa_sums(month_subset) > 0, month_subset)
  
  #Extract OTU table
  otu_t <- t(otu_table(month_subset))
  
  #Grouping vector
  group_vector <- ifelse(
    sample_data(month_subset)$collection_month == focal_month,
    focal_month,
    other_label
  )
  
  #Run CLAM test
  month_clam <- clamtest(
    otu_t,
    group_vector,
    alpha = alpha_val,
    specialization = spec_val
  )
  cat("\n===== CLAM STATISTICS: ", comparison_label, " =====\n")
  cat("Alpha: ", alpha_val, "\n")
  cat("Specialization threshold: ", spec_val, "\n")
  cat("Total ASVs tested: ", nrow(month_clam), "\n")
  cat("Samples in focal month: ", sum(group_vector == focal_month), "\n")
  cat("Samples in other months: ", sum(group_vector == other_label), "\n")
  cat("Specialists (", focal_month, "): ",
      sum(month_clam$Classes == paste0("Specialist_", focal_month)), "\n", sep = "")
  cat("Specialists (", other_label, "): ",
      sum(month_clam$Classes == paste0("Specialist_", other_label)), "\n", sep = "")
  cat("Generalists: ", sum(month_clam$Classes == "Generalist"), "\n")
  cat("Too rare: ", sum(month_clam$Classes == "Too_rare"), "\n")
  cat("============================================\n\n")
}

  #Extract taxonomy
  tax_df <- as.data.frame(tax_table(month_subset))
  
  #Export specialists - focal month
  focal_class <- paste0("Specialist_", focal_month)
  focal_specialists <- month_clam[month_clam$Classes == focal_class, ]
  focal_asvs <- focal_specialists$Species
  
  tax_subset_focal <- tax_df[row.names(tax_df) %in% focal_asvs, ]
  tax_subset_focal$Specialist_Month <- focal_month
  tax_subset_focal$Comparison <- comparison_label
  
  csv_focal <- paste0("Specialist_", focal_month, "_vs_", other_label, ".csv")
  write.csv(tax_subset_focal, csv_focal, row.names = TRUE)
  cat("Saved specialist CSV:", csv_focal, "\n")
  
  #Export specialists, combined group
  other_class <- paste0("Specialist_", other_label)
  other_specialists <- month_clam[month_clam$Classes == other_class, ]
  other_asvs <- other_specialists$Species
  
  tax_subset_other <- tax_df[row.names(tax_df) %in% other_asvs, ]
  tax_subset_other$Specialist_Month <- other_label
  tax_subset_other$Comparison <- comparison_label
  
  csv_other <- paste0("Specialist_", other_label, "_vs_", focal_month, ".csv")
  write.csv(tax_subset_other, csv_other, row.names = TRUE)
  cat("Saved specialist CSV:", csv_other, "\n")
  
  #Export generalists and too rare
  write.csv(
    month_clam[month_clam$Classes == "Generalist", ],
    paste0("generalists_", comparison_label, ".csv"),
    row.names = TRUE
  )
  
  write.csv(
    month_clam[month_clam$Classes == "Too_rare", ],
    paste0("too_rare_", comparison_label, ".csv"),
    row.names = TRUE
  )
  
  cat("Completed:", comparison_label, "\n")
}

# ============================================================
# CLAM SPECIALISTS: 
# ============================================================
# Utility: coerce Terrain column
.fix_terrain <- function(x) {
  vapply(x, function(xx) {
    if (length(xx) == 0 || all(is.na(xx))) NA_character_
    else if (length(xx) > 1) paste(xx, collapse = ";")
    else as.character(xx)
  }, FUN.VALUE = character(1))
}
# Utility: draw plots and print summary
.print_month_outputs <- function(tax_subset, title_prefix) {
  # ---- Terrain Pie ----
  terrain_counts <- tax_subset %>%
    filter(!is.na(Terrain), Terrain != "") %>%
    count(Terrain) %>%
    mutate(Proportion = n / sum(n))
  
  print(
    ggplot(terrain_counts, aes(x = "", y = Proportion, fill = Terrain)) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar(theta = "y") +
      scale_fill_manual(values = c(
        "Aquatic" = "#1F78B4",
        "Terrestrial" = "#33A02C",
        "Both" = "#B2B200"
      )) +
      labs(title = paste0(title_prefix, " — Terrain"), fill = "Terrain") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5))
  )
  
  #Phylum Bar
  phylum_counts <- tax_subset %>%
    filter(!is.na(Phylum), Phylum != "") %>%
    count(Phylum)
  
  print(
    ggplot(phylum_counts, aes(x = reorder(Phylum, -n), y = n, fill = Phylum)) +
      geom_bar(stat = "identity", show.legend = FALSE) +
      labs(title = paste0(title_prefix, " — Phylum"),
           y = "Number of ASVs", x = "Phylum") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  
  #Summary Table (print) 
  asv_counts <- tax_subset %>%
    filter(!is.na(Phylum), Phylum != "") %>%
    count(Phylum, name = "ASV_Count")
  
  species_counts <- tax_subset %>%
    filter(!is.na(Phylum), Phylum != "", !is.na(Species), Species != "") %>%
    distinct(Phylum, Species) %>%
    count(Phylum, name = "Species_Count")
  
  phylum_summary <- full_join(asv_counts, species_counts, by = "Phylum") %>%
    arrange(desc(ASV_Count))
  
  cat("\n---------- PHYLUM SUMMARY:", title_prefix, "----------\n")
  print(phylum_summary, n = nrow(phylum_summary))
  cat("------------------------------------------------------\n")
}

# ============================================================
# OCTOBER
# ============================================================
csv_path <- "Specialist_Oct_vs_Nov_Dec.csv"  # <-- update if needed
stopifnot(file.exists(csv_path))

cat("\nReading CSV:", csv_path, "\n")
tax_oct <- read_csv(csv_path, show_col_types = FALSE)

# Ensure columns exist and are clean text
needed <- c("Terrain", "Phylum", "Species")
for (nm in needed) if (!nm %in% names(tax_oct)) tax_oct[[nm]] <- NA_character_
tax_oct <- tax_oct |>
  mutate(
    Terrain = trimws(as.character(Terrain)),
    Phylum  = trimws(as.character(Phylum)),
    Species = trimws(as.character(Species))
  )

#Terrain Pie
terrain_counts <- tax_oct |>
  dplyr::filter(!is.na(Terrain), Terrain != "") |>
  dplyr::count(Terrain, name = "n") |>
  dplyr::mutate(Proportion = n / sum(n))

#Print table
print(terrain_counts)

print(
  ggplot(terrain_counts, aes(x = "", y = Proportion, fill = Terrain)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c("Aquatic"="#1F78B4", "Terrestrial"="#33A02C", "Both"="#B2B200")) +
    labs(title = "Oct Specialists — Terrain", fill = "Terrain") +
    theme_void()
)

#Phylum Bar
phylum_counts <- tax_oct |>
  dplyr::filter(!is.na(Phylum), Phylum != "") |>
  dplyr::group_by(Phylum) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop")
print(phylum_counts)
print(
  ggplot(phylum_counts, aes(x = reorder(Phylum, -n), y = n, fill = Phylum)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(title = "Oct Specialists — Phylum", y = "Number of ASVs", x = "Phylum") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

#Phylum summary 
asv_counts <- tax_oct |>
  dplyr::filter(!is.na(Phylum), Phylum != "") |>
  dplyr::group_by(Phylum) |>
  dplyr::summarise(ASV_Count = dplyr::n(), .groups = "drop")

species_counts <- tax_oct |>
  dplyr::filter(!is.na(Phylum), Phylum != "",
                !is.na(Species), Species != "") |>
  dplyr::distinct(Phylum, Species) |>
  dplyr::group_by(Phylum) |>
  dplyr::summarise(Species_Count = dplyr::n(), .groups = "drop")


phylum_summary <- full_join(asv_counts, species_counts, by = "Phylum") |>
  arrange(desc(ASV_Count))

cat("\n---- PHYLUM SUMMARY: OCT ----\n")
print(phylum_summary, n = nrow(phylum_summary))

# ============================================================
# NOVEMBER
# ============================================================
csv_path <- "Specialist_Nov_vs_Oct_Dec.csv"  # <-- update if needed
stopifnot(file.exists(csv_path))

cat("\nReading CSV:", csv_path, "\n")
tax_nov <- readr::read_csv(csv_path, show_col_types = FALSE)

# Ensure columns exist and are clean text
needed <- c("Terrain", "Phylum", "Species")
for (nm in needed) if (!nm %in% names(tax_nov)) tax_nov[[nm]] <- NA_character_

tax_nov <- tax_nov |>
  dplyr::mutate(
    Terrain = trimws(as.character(Terrain)),
    Phylum  = trimws(as.character(Phylum)),
    Species = trimws(as.character(Species))
  )

#Terrain Pie
terrain_counts <- tax_nov |>
  dplyr::filter(!is.na(Terrain), Terrain != "") |>
  dplyr::count(Terrain, name = "n") |>
  dplyr::mutate(Proportion = n / sum(n))

print(terrain_counts)

print(
  ggplot(terrain_counts, aes(x = "", y = Proportion, fill = Terrain)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c("Aquatic"="#1F78B4",
                                 "Terrestrial"="#33A02C",
                                 "Both"="#B2B200")) +
    labs(title = "Nov Specialists — Terrain", fill = "Terrain") +
    theme_void()
)

#Phylum Bar
phylum_counts <- tax_nov |>
  dplyr::filter(!is.na(Phylum), Phylum != "") |>
  dplyr::group_by(Phylum) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

print(phylum_counts)

print(
  ggplot(phylum_counts, aes(x = reorder(Phylum, -n), y = n, fill = Phylum)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(title = "Nov Specialists — Phylum", y = "Number of ASVs", x = "Phylum") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

#Phylum Summary
asv_counts <- tax_nov |>
  dplyr::filter(!is.na(Phylum), Phylum != "") |>
  dplyr::group_by(Phylum) |>
  dplyr::summarise(ASV_Count = dplyr::n(), .groups = "drop")

species_counts <- tax_nov |>
  dplyr::filter(!is.na(Phylum), Phylum != "",
                !is.na(Species), Species != "") |>
  dplyr::distinct(Phylum, Species) |>
  dplyr::group_by(Phylum) |>
  dplyr::summarise(Species_Count = dplyr::n(), .groups = "drop")

phylum_summary <- dplyr::full_join(asv_counts, species_counts, by = "Phylum") |>
  dplyr::arrange(desc(ASV_Count))

cat("\n---- PHYLUM SUMMARY: NOV ----\n")
print(phylum_summary, n = nrow(phylum_summary))

# ============================================================
# DECEMBER
# ============================================================
csv_path <- "Specialist_Dec_vs_Oct_Nov.csv"  # <-- update if needed
stopifnot(file.exists(csv_path))

cat("\nReading CSV:", csv_path, "\n")
tax_dec <- readr::read_csv(csv_path, show_col_types = FALSE)

# Ensure columns exist and are clean text
needed <- c("Terrain", "Phylum", "Species")
for (nm in needed) if (!nm %in% names(tax_dec)) tax_dec[[nm]] <- NA_character_

tax_dec <- tax_dec |>
  dplyr::mutate(
    Terrain = trimws(as.character(Terrain)),
    Phylum  = trimws(as.character(Phylum)),
    Species = trimws(as.character(Species))
  )

#Terrain Pie
terrain_counts <- tax_dec |>
  dplyr::filter(!is.na(Terrain), Terrain != "") |>
  dplyr::count(Terrain, name = "n") |>
  dplyr::mutate(Proportion = n / sum(n))

print(terrain_counts)

print(
  ggplot(terrain_counts, aes(x = "", y = Proportion, fill = Terrain)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c("Aquatic"="#1F78B4",
                                 "Terrestrial"="#33A02C",
                                 "Both"="#B2B200")) +
    labs(title = "Dec Specialists — Terrain", fill = "Terrain") +
    theme_void()
)

#Phylum Bar
phylum_counts <- tax_dec |>
  dplyr::filter(!is.na(Phylum), Phylum != "") |>
  dplyr::group_by(Phylum) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

print(phylum_counts)

print(
  ggplot(phylum_counts, aes(x = reorder(Phylum, -n), y = n, fill = Phylum)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(title = "Dec Specialists — Phylum", y = "Number of ASVs", x = "Phylum") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

#Phylum summary
asv_counts <- tax_dec |>
  dplyr::filter(!is.na(Phylum), Phylum != "") |>
  dplyr::group_by(Phylum) |>
  dplyr::summarise(ASV_Count = dplyr::n(), .groups = "drop")

species_counts <- tax_dec |>
  dplyr::filter(!is.na(Phylum), Phylum != "",
                !is.na(Species), Species != "") |>
  dplyr::distinct(Phylum, Species) |>
  dplyr::group_by(Phylum) |>
  dplyr::summarise(Species_Count = dplyr::n(), .groups = "drop")

phylum_summary <- dplyr::full_join(asv_counts, species_counts, by = "Phylum") |>
  dplyr::arrange(desc(ASV_Count))

cat("\n---- PHYLUM SUMMARY: DEC ----\n")
print(phylum_summary, n = nrow(phylum_summary))
