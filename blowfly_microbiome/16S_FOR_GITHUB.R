
############################## 16S MICROBIOME PIPELINE ##########################

################################################################################
## LIBRARIES
################################################################################

# >>> CHANGE: deduplicated + grouped libraries for clarity
library(BiocManager)
library(phyloseq)
library(vegan)
library(microbiome)
library(microViz)
library(decontam)
library(MicrobiomeStat)
library(MiscMetabar)
library(DECIPHER)
library(phangorn)
library(ape)
library(tidyverse)
library(dplyr)
library(tidyr)
library(plyr)
library(readr)
library(readxl)
library(tibble)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(gridExtra)
library(ggeasy)
library(aplot)
library(ggplotify)
library(showtext)
library(car)
library(GUniFrac)

# Set seed
set.seed(666)

################################################################################
## READ DATA
################################################################################

setwd("/Users/GGPC/OneDrive/Desktop/microbiome_november/16S")

asv_mat <- read_tsv("ASVs_counts_16S_bf.tsv")
tax_mat <- read.csv("MergedTaxonomy_16S_bf.csv", header = TRUE)
samples_df <- read_excel("metadata_subset_16S.xlsx")

# Set rownames
asv_mat <- asv_mat %>% column_to_rownames("ASV")
tax_mat <- tax_mat %>% column_to_rownames("ASV")
samples_df <- samples_df %>% column_to_rownames("sample")

# Convert to matrices
asv_mat <- as.matrix(asv_mat)
tax_mat <- as.matrix(tax_mat)

# Build phyloseq
ps <- phyloseq(
  otu_table(asv_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(samples_df)
)

ps   # 8148 taxa, 40 samples

################################################################################
## FILTERING
################################################################################

#Filter for uncharacterised reads
table(tax_table(ps)[, "Phylum"], exclude = NULL)
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) 
ps0 #7993 taxa and 40 samples
ps0 <- ps 

#Filter out human contamination
ps_nohomo <- (subset_taxa(ps, Family != "Hominidae"))
ps_nohomo #7128 taxa; 40 samples

#Remove plant contaminations
ps_clean <- subset_taxa(ps_nohomo, Phylum != "Streptophyta")
ps_clean #7,127 taxa; 40 samples 

#Remove archaea contaminations
ps_clean <- subset_taxa(ps_clean, Phylum != "Methanobacteriota")
ps_clean #7108
ps_nohomo <- ps_clean

################################################################################
## DECONTAM
################################################################################

library(decontam)
sample_data(ps_nohomo)$is.neg <- sample_names(ps_nohomo) %in% c("blank-16S-1", "blank-16S-2", "BLANK-16S-3", "BLANK-16S-4", 
                                                                "PCRNEG-STEP1-16S-1", "PCRNEG-STEP1-16S-3", 
                                                                "PCRNEG-STEP2-16S-2", "PCRNEG-STEP2-16S-4")
contamdf <- isContaminant(ps_nohomo, method="prevalence", neg="is.neg")
table(contamdf$contaminant) #false = 7102; true = 25
contaminants <- rownames(contamdf)[contamdf$contaminant == TRUE]
ps_clean <- prune_taxa(!taxa_names(ps_nohomo) %in% contaminants, ps_nohomo)
ps_clean #7,083 taxa; 40 samples

#remove taxa with 0 reads 
ps_nozero <- prune_taxa(taxa_sums(ps_nohomo) > 0, ps_nohomo)
ps_nohomo <- ps_nozero
ps_nohomo #2732

################################################################################
## BALANCE SAMPLES (8 PER SPECIES); also removes controls
################################################################################

sample_df <- as.data.frame(sample_data(ps_nohomo))

species_list <- split(rownames(sample_df), sample_df$BLAST_ID)

species_subs <- lapply(species_list, function(x) {
  if (length(x) >= 8) sample(x, 8) else NULL
})

selected_samples <- unlist(species_subs)

ps_balanced <- prune_samples(selected_samples, ps_nohomo)

# >>> CHANGE: explicit control removal
ps_final <- subset_samples(ps_balanced, BLAST_ID != "control")

ps_final #2732 taxa; 32 samples


################################################################################
## NORMALISATION (TSS)
################################################################################
source('mStat_convert_phyloseq_to_data_obj.R')
source('mStat_normalize_data.R')
source('mStat_validate_data.R')
source('update_data_obj_count.R')
source('mStat_calculate_alpha_diversity.R')
source('generate_alpha_test_single.R')

#convert phyloseq to data object
data.obj <- mStat_convert_phyloseq_to_data_obj(ps_final)

#normalise
norm.data.obj <- mStat_normalize_data(data.obj, "TSS")$data.obj.norm 

#validate
mStat_validate_data(norm.data.obj)

################################################################################
## RELATIVE ABUNDANCE PLOT -- FIGURE 1 
################################################################################

# Filter out unknown/ambiguous classes
ps_class <- tax_glom(ps_final, taxrank = "Class")
ps_class <- subset_taxa(ps_class, !grepl("Incertae_sedis|unassigned|Unknown", Class))
ps_class
# Step 2: Keep top 15 most abundant classes
top15 <- names(sort(taxa_sums(ps_class), decreasing = TRUE))[1:15]
ps_top15_class <- prune_taxa(top15, ps_class)

#make relative 
ps_top15_class <- transform_sample_counts(ps_top15_class, function(x) x / sum(x))

# Step 3: Melt to dataframe
df <- psmelt(ps_top15_class)

#Sort out the colours

myPal <- tax_palette(
  data = ps_class, rank = "Class", n = 15, pal = "greenArmytage",
  add = c(Other = "white")
)
tax_palette_plot(myPal)

# Common theme
base_theme <- theme_minimal(base_size = 11, base_family = "Montserrat") +
  theme(
    text = element_text(color = "black"),  # <-- this line makes all text black
    axis.text.x = element_text(angle = 90, size = 5, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none",
    plot.margin = margin(5, -5, 5, -5),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
  )
# Step 5: Subset by species
df_chilli  <- subset(df, BLAST_ID == "Calliphora_hilli")
df_stygia  <- subset(df, BLAST_ID == "Calliphora_stygia")
df_chryso  <- subset(df, BLAST_ID == "Chrysomya_rufifacies")
df_lucilia <- subset(df, BLAST_ID == "Lucilia_spp")

# Remove y-axis text for all but first plot
no_y_axis <- theme(
  axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.y = element_blank()
)

# Plot 1 — keep y-axis
p1 <- ggplot(df_chilli, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  labs(title = "C. hilli", y = "Relative abundance", x = NULL) +
  scale_fill_manual(values = myPal) +
  base_theme + theme(plot.title = element_text(face = "bold.italic"))

# Plot 2
p2 <- ggplot(df_stygia, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  labs(title = "C. stygia", x = NULL) +
  scale_fill_manual(values = myPal) +
  base_theme + no_y_axis + theme(plot.title = element_text(face = "bold.italic"))

# Plot 3
p3 <- ggplot(df_chryso, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  labs(title = "Ch. rufifacies", x = NULL) +
  scale_fill_manual(values = myPal) +
  base_theme + no_y_axis + theme(plot.title = element_text(face = "bold.italic"))

# Plot 4 — no "Sample" x-axis label
p4 <- ggplot(df_lucilia, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  labs(title = "Lucilia spp.", x = NULL) +
  scale_fill_manual(values = myPal) +
  base_theme + no_y_axis + theme(plot.title = element_text(face = "bold.italic"))

# Combine with adjusted widths & shared legend
final_plot <- (p1 + p2 + p3 + p4) +
  plot_layout(guides = "collect", widths = c(1, 1, 1, 1)) & 
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", family = "Montserrat", size = 13),
    legend.text = element_text(size = 12),
    legend.box = "vertical",
    legend.title.align = 0.5
  ) &
  guides(fill = guide_legend(title = "Class", title.position = "top"))

# Show plot
final_plot

################################################################################
## ALPHA DIVERSITY -- FIGURE 2
################################################################################

source("mStat_calculate_alpha_diversity.R")
source("generate_alpha_test_single.R")

alpha.res <- mStat_calculate_alpha_diversity(
  data.obj$feature.tab,
  alpha.name = "shannon"
)

alpha.df <- cbind(data.obj$meta.dat, alpha.res)

##PLOT 
# Clean species labels
alpha.df$Species_clean <- gsub("_", " ", alpha.df$BLAST_ID)
alpha.df$Species_clean <- dplyr::recode(alpha.df$Species_clean,
                                        "Calliphora hilli" = "C. hilli",
                                        "Calliphora stygia" = "C. stygia",
                                        "Chrysomya rufifacies" = "Ch. rufifacies",
                                        "Lucilia spp." = "Lucilia spp"
)
alpha.df$Species_clean <- factor(alpha.df$Species_clean,
                                 levels = c("C. hilli", "C. stygia", 
                                            "Ch. rufifacies", "Lucilia spp"))


# Calculate mean and SE per species
summary_df <- alpha.df %>%
  group_by(Species_clean) %>%
  summarise(
    mean = mean(shannon),
    se = sd(shannon) / sqrt(n())
  )

#plot
ggplot(alpha.df, aes(x = Species_clean, y = shannon, colour = Species_clean, fill = Species_clean)) +
  geom_boxplot(alpha = 0.5, width = 0.5, outlier.shape = NA, colour = "black") +
  geom_jitter(width = 0.15, size = 4, alpha = 0.7) +
  scale_colour_manual(values = c(
    "C. hilli" = "#009B77", #pantone coty 2013 - emerald
    "C. stygia" = "#DD4124",#pantone coty 2012 - tangerine tango
    "Ch. rufifacies" = "#AFAB23",#pantone coty 2009 - mimosa
    "Lucilia spp" = "#4169E1" #pantone coty 2001 - fuchsia rose
  )) +
  scale_fill_manual(values = c(
    "C. hilli" = "#009B77",#pantone coty 2013 - emerald
    "C. stygia" = "#DD4124",#pantone coty 2012 - tangerine tango
    "Ch. rufifacies" = "#AFAB23",#pantone coty 2009 - mimosa
    "Lucilia spp" = "#4169E1"#pantone coty 2001 - fuchsia rose
  )) +
  theme_minimal(base_size = 15, base_family = "Montserrat") +
  labs(x = NULL, y = "Alpha diversity (Shannon)") +
  theme(
    legend.position = "none",
    axis.title.y = element_text(margin = margin(r = 18), size = 16),
    axis.text.x = element_text(size = 15, colour = "black", face = "italic"),
    panel.grid.major.x = element_blank()
  ) + scale_x_discrete(labels = c(
    "C. hilli"     = expression(italic("C. hilli")),
    "C. stygia"    = expression(italic("C. stygia")),
    "Ch. rufifacies" = expression(italic("Ch. rufifacies")),
    "Lucilia spp"          = expression(italic("Lucilia")~spp.)
  ))

#Test for differences between groups
anova_model <- aov(shannon ~ BLAST_ID + location, data = alpha.df)

# View summary results:
summary(anova_model)
TukeyHSD(anova_model)

################################################################################
## BETA DIVERSITY (BRAY–CURTIS)
################################################################################

#Beta diversity (dissimilarities between groups) based on Bray-Curtis differences
library(aplot)
library(ggplotify)
library(dplyr)
library(ggplot2)
library(tibble)
source('mStat_calculate_beta_diversity.R')
source('generate_beta_ordination_single.R')
source('mStat_calculate_PC.R')
source('mStat_get_palette.R')
source('mStat_get_theme.R')

#generate beta ordination based on bray-curtis dissimilarities
beta.ord <- generate_beta_ordination_single(
  data.obj   = norm.data.obj,
  dist.obj   = NULL,
  pc.obj     = NULL,
  subject.var = NULL,
  time.var   = NULL,
  t.level    = NULL,
  group.var  = "BLAST_ID",
  strata.var = NULL,
  dist.name  = "BC",
  base.size  = 18,
  theme.choice = "bw",
  palette    = NULL,
  pdf        = FALSE
)
beta.ord

#Let's generate a nicer looking one
#calculate distances
beta.bc <- mStat_calculate_beta_diversity(
  data.obj = norm.data.obj,
  dist.name = "BC"
)

bc.dist <- beta.bc$BC
pcoa_16s <- pcoa(bc.dist) #run the pcoa
pcoa_coords <- as.data.frame(pcoa_16s$vectors) #extract the coordinates
pcoa_df <- pcoa_coords %>%
  rownames_to_column("Sample") %>%
  left_join(norm.data.obj$meta.dat %>% 
              rownames_to_column("Sample"),
            by = "Sample")

#get % variance explained for axis labels
var1 <- round(pcoa_16s$values$Relative_eig[1] * 100, 2)
var2 <- round(pcoa_16s$values$Relative_eig[2] * 100, 2)

#plot!

# ---- Colour palette ----
species_palette <- c(
  "Calliphora_hilli"       = "#009B77",
  "Calliphora_stygia"      = "#DD4124",
  "Chrysomya_rufifacies"   = "#AFAB23",
  "Lucilia_spp"            = "#4169E1"
)

# ---- Main plot ----
ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, colour = BLAST_ID)) +
  geom_point(size = 6, alpha = 1) +
  stat_ellipse(size = 1, alpha = 1) +
  
  # Colour scale
  scale_color_manual(
    values = species_palette,
    labels = c(
      expression(italic("C. hilli")),
      expression(italic("C. stygia")),
      expression(italic("Ch. rufifacies")),
      expression(italic("Lucilia")~spp.)
    )
  ) +
  
  # FIXED ASPECT RATIO (same as phyloseq)
  coord_fixed(sqrt(pcoa_16s$values$Relative_eig[1] /
                     pcoa_16s$values$Relative_eig[2])) +
  
  # Axis labels with % variance
  labs(
    x = paste0("Axis 2 (",
               round(100 * pcoa_16s$values$Relative_eig[1], 1), "%)"),
    y = paste0("Axis 3 (",
               round(100 * pcoa_16s$values$Relative_eig[2], 1), "%)"),
    colour = NULL
  ) +
  
  # Theme
  theme(
    legend.position = "right",
    legend.text = element_text(size = 13),
    axis.title = element_text(size = 14,
                              margin = margin(t = 10, r = 10)),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  )

#PERMANOVA
source('generate_beta_test_single.R')
library(GUniFrac)

generate_beta_test_single(
  data.obj = norm.data.obj,
  dist.obj = NULL,
  time.var = NULL,
  t.level = NULL,
  group.var = "BLAST_ID", 
  adj.vars = "location",
  dist.name = c('BC') 
)

#beta disper looks at the significance pairwisely
bd <- betadisper(bc.dist, pcoa_df$BLAST_ID)
bd


#PERMUTATION TEST
permutest(bd, pairwise=T)

# Extract distances to centroids
bd_df <- data.frame(
  Species = bd$group,
  DistanceToCentroid = bd$distances
)
# Clean species names and fix typo
bd_df$Species <- gsub("_", " ", bd$group)
bd_df$Species <- gsub("Chrysomya rufifacies", "Ch. rufifacies", bd_df$Species)
bd_df$Species <- gsub("Calliphora stygia", "C. stygia", bd_df$Species)
bd_df$Species <- gsub("Calliphora hilli", "C. hilli", bd_df$Species)

# Make it a factor with correct order
bd_df$Species <- factor(
  bd_df$Species,
  levels = c("C. hilli", "C. stygia", "Ch. rufifacies", "Lucilia spp")
)

# Define custom color palette
my_colors <- c(
  "C. hilli" = "#fb8500",
  "C. stygia" = "#ffb703",
  "Ch. rufifacies" = "#219ebc",
  "Lucilia spp" = "#023047"
)

ggplot(bd_df, aes(x = Species, y = DistanceToCentroid)) +
  # Boxplot with original fill colours
  geom_boxplot(aes(fill = Species), outlier.shape = NA, alpha = 0.6, width = 0.6) +
  
  # Jittered points with dark grey outline
  geom_jitter(aes(fill = Species),
              color = "grey10",       # dark grey outline
              width = 0.2, size = 2.5,
              shape = 21, stroke = 0.8) +
  
  # Colour palettes
  scale_fill_manual(values = my_colors) +         # for fill (box and point)
  
  # Theme and axis formatting
  theme_minimal(base_family = "Montserrat", base_size = 15) +
  labs(x = NULL, y = "Distance to centroid") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14, colour = "black", face = "italic"),
    axis.title.y = element_text(size = 16, margin = margin(r = 18)),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  )

################################################################################
############################## CLAM TESTS #######################################
######################## ASV + CLASS LEVEL ANALYSES #############################
################################################################################

library(vegan)
library(phyloseq)
library(dplyr)

# Split species
hilli   <- subset_samples(ps_final, BLAST_ID == "Calliphora_hilli")
stygia  <- subset_samples(ps_final, BLAST_ID == "Calliphora_stygia")
chryso  <- subset_samples(ps_final, BLAST_ID == "Chrysomya_rufifacies")
luci    <- subset_samples(ps_final, BLAST_ID == "Lucilia_spp")

no_hilli  <- subset_samples(ps_final, BLAST_ID != "Calliphora_hilli")
no_stygia <- subset_samples(ps_final, BLAST_ID != "Calliphora_stygia")
no_chryso <- subset_samples(ps_final, BLAST_ID != "Chrysomya_rufifacies")
no_luci   <- subset_samples(ps_final, BLAST_ID != "Lucilia_spp")

sample_data(no_hilli)$BLAST_ID  <- "Other"
sample_data(no_stygia)$BLAST_ID <- "Other"
sample_data(no_chryso)$BLAST_ID <- "Other"
sample_data(no_luci)$BLAST_ID   <- "Other"


################################################################################
## ASV-LEVEL CLAM TESTS
################################################################################

# >>> CHANGE: explicit list to store results consistently
clam_asv_results <- list()
specialist_counts <- data.frame(
  Species = character(),
  N_specialist_ASVs = integer()
)

run_asv_clam <- function(focal, others, species_name) {
  
  ps_pair <- merge_phyloseq(focal, others)
  
  otu <- as(otu_table(ps_pair), "matrix")
  if (taxa_are_rows(ps_pair)) otu <- t(otu)
  otu <- otu[, colSums(otu) > 0]
  
  groups <- as.character(sample_data(ps_pair)$BLAST_ID)
  
  clam_res <- clamtest(
    otu,
    groups,
    alpha = 0.01,
    specialization = 2/3
  )
  
  # Count specialist ASVs
  class_label <- paste0("Specialist_", species_name)
  n_spec <- sum(clam_res$Classes == class_label)
  
  # Store
  clam_asv_results[[species_name]] <<- clam_res
  
  specialist_counts <<- rbind(
    specialist_counts,
    data.frame(
      Species = species_name,
      N_specialist_ASVs = n_spec
    )
  )
  
  # Plot
  plot(
    clam_res,
    main = paste("CLAM (ASV level):", species_name, "vs others")
  )
  
  return(clam_res)
}

################################################################################
## RUN ASV-LEVEL CLAM FOR EACH SPECIES
################################################################################

# Subsets already created earlier:
# hilli, stygia, chryso, luci
# no_hilli, no_stygia, no_chryso, no_luci

clam_hilli_asv <- run_asv_clam(
  hilli, no_hilli, "Calliphora_hilli"
)

clam_stygia_asv <- run_asv_clam(
  stygia, no_stygia, "Calliphora_stygia"
)

clam_chryso_asv <- run_asv_clam(
  chryso, no_chryso, "Chrysomya_rufifacies"
)

clam_luci_asv <- run_asv_clam(
  luci, no_luci, "Lucilia_spp"
)

################################################################################
## ASV SPECIALIST COUNTS (REPORT THESE)
################################################################################

specialist_counts

write.csv(
  specialist_counts,
  "ASV_SPECIALIST_COUNTS_PER_SPECIES.csv",
  row.names = FALSE
)

################################################################################
## OPTIONAL: EXTRACT ASV LISTS PER SPECIES
################################################################################

extract_asv_specialists <- function(clam_res, species_name) {
  
  class_label <- paste0("Specialist_", species_name)
  
  clam_res %>%
    as.data.frame() %>%
    filter(Classes == class_label) %>%
    pull(Species)
}

ASVs_hilli   <- extract_asv_specialists(clam_hilli_asv, "Calliphora_hilli")
ASVs_stygia  <- extract_asv_specialists(clam_stygia_asv, "Calliphora_stygia")
ASVs_chryso  <- extract_asv_specialists(clam_chryso_asv, "Chrysomya_rufifacies")
ASVs_luci    <- extract_asv_specialists(clam_luci_asv, "Lucilia_spp")

tax_df <- as.data.frame(tax_table(ps_final))  # Extract taxonomy table
tax_df$ASV <- rownames(tax_df)               # Add ASV IDs as a column

get_specialist_species <- function(asv_vector, tax_df) {
  tax_df %>%
    filter(ASV %in% asv_vector) %>%   # keep only specialist ASVs
    group_by(Species) %>%              # group by fungal species
    summarise(n_ASVs = n()) %>%       # count ASVs per species
    arrange(desc(n_ASVs))
}

species_hilli  <- get_specialist_species(ASVs_hilli, tax_df)
species_stygia <- get_specialist_species(ASVs_stygia, tax_df)
species_chryso <- get_specialist_species(ASVs_chryso, tax_df)
species_luci   <- get_specialist_species(ASVs_luci, tax_df)

data.frame(
  Host = c("C. hilli", "C. stygia", "Ch. rufifacies", "Lucilia spp."),
  n_species = c(
    nrow(species_hilli),
    nrow(species_stygia),
    nrow(species_chryso),
    nrow(species_luci)
  )
)

write.csv(species_hilli, "specialists_hilli.csv") # to figure out how many ASVS are associated with species
write.csv(species_hilli, "specialists_stygia.csv")
write.csv(species_hilli, "specialists_chryso.csv")
write.csv(species_hilli, "specialists_luci.csv")

################################################################################
## CLASS-LEVEL CLAM TESTS (COLLAPSED)
################################################################################

run_clam_class_specialists <- function(physeq_obj,
                                       species_col = "BLAST_ID",
                                       taxrank = "Class",
                                       alpha = 0.01,
                                       specialization = 2/3,
                                       outfile_prefix = "OUTPUT") {
  
  ps_class <- tax_glom(physeq_obj, taxrank = taxrank, NArm = TRUE)
  
  otu <- as(otu_table(ps_class), "matrix")
  if (taxa_are_rows(ps_class)) otu <- t(otu)
  otu <- otu[, colSums(otu) > 0, drop = FALSE]
  
  groups <- as.factor(sample_data(ps_class)[[species_col]])
  
  clam_res <- clamtest(
    otu,
    groups,
    alpha = alpha,
    specialization = specialization
  )
  
  spec <- as.data.frame(clam_res) %>%
    filter(grepl("^Specialist_", Classes)) %>%
    mutate(Host_Species = sub("^Specialist_", "", Classes))
  
  tax_df <- as.data.frame(tax_table(ps_class))
  tax_df$Taxon <- rownames(tax_df)
  
  spec_out <- spec %>%
    left_join(tax_df, by = c("Species" = "Taxon"))
  
  write.csv(
    spec_out,
    paste0(outfile_prefix, "_CLASS_SPECIALISTS.csv"),
    row.names = FALSE
  )
  
  plot(
    clam_res,
    main = paste("CLAM (Class level):", outfile_prefix)
  )
  
  return(list(clam = clam_res, specialists = spec_out))
}

################################################################################
## RUN CLASS-LEVEL CLAM TESTS
################################################################################

hilli_class   <- run_clam_class_specialists(
  merge_phyloseq(hilli, no_hilli),
  outfile_prefix = "Calliphora_hilli"
)

stygia_class  <- run_clam_class_specialists(
  merge_phyloseq(stygia, no_stygia),
  outfile_prefix = "Calliphora_stygia"
)

chryso_class  <- run_clam_class_specialists(
  merge_phyloseq(chryso, no_chryso),
  outfile_prefix = "Chrysomya_rufifacies"
)

lucilia_class <- run_clam_class_specialists(
  merge_phyloseq(luci, no_luci),
  outfile_prefix = "Lucilia_spp"
)

################################################################################
## COMBINED CLASS-LEVEL TABLE
################################################################################

all_class_specialists <- bind_rows(
  hilli_class$specialists,
  stygia_class$specialists,
  chryso_class$specialists,
  lucilia_class$specialists
)

write.csv(
  all_class_specialists,
  "ALL_SPECIES_CLASS_SPECIALISTS.csv",
  row.names = FALSE
)

################################################################################
## SUPPLEMENTARY FIGURES
################################################################################

# SACs
library(MiscMetabar)
p <- accu_plot(ps_final, "BLAST_ID", add_nb_seq = TRUE, by.fact = TRUE, step = 10)
p + scale_color_manual(values = c("Lucilia_spp"          = "cornflowerblue",
                                  "Calliphora_stygia"    = "brown1",
                                  "Calliphora_hilli"     = "darksalmon",
                                  "Chrysomya_rufifacies" = "green4"))

p + guides(fill = "none")+
  theme(
    # AXIS TITLES
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    
    # LEGEND TEXT (italic + moderate spacing)
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "italic"),
    legend.key.height = unit(1, "lines"),
    
    # REMOVE GREY BACKGROUND
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank()
  )+
  scale_colour_discrete(labels = ~ gsub("_", " ", .x))

# Library size

#plot library size 
# Look for skew in coverage across sample types
DATA.2 <- ps_final  

df = as.data.frame(sample_data(DATA.2))
df$LibrarySize = sample_sums(DATA.2)
df = df[order(df$LibrarySize),]
df$Index = seq(nrow(df))
#Plot ordered library size & coloured by species  
level_order2 <- c('Calliphora_stygia', 'Calliphora_hilli', 'Chrysomya_rufifacies', 'Lucilia_spp') 

p <- ggplot(data=df, aes(x=Index, y=LibrarySize, colour = BLAST_ID)) +  
  geom_point(size = 4) +
  facet_wrap(~ factor(BLAST_ID, level = level_order2)) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  scale_colour_manual(
    values = c(
      "Lucilia_spp"          = "cornflowerblue",
      "Calliphora_stygia"    = "brown1",
      "Calliphora_hilli"     = "darksalmon",
      "Chrysomya_rufifacies" = "green4"
    ),
    labels = function(x) gsub("_", " ", x)
  ) +
  scale_y_continuous(trans = "sqrt") +
  theme(
    # AXIS TITLES
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    
    # LEGEND TEXT (italic + moderate spacing)
    legend.title = element_blank(),
    legend.text = element_text(size = 12, face = "italic"),
    legend.key.height = unit(0.9, "lines"),
    
    # REMOVE GREY BACKGROUND
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank()
  )


p

# Sample coverage
#Create a sample sum table to look at coverage metrics
# Extract read counts per sample
sample_sum_df <- data.frame(
  SampleID = names(sample_sums(ps_final)),
  Reads = sample_sums(ps_final)
)

# Extract metadata (e.g., species) from phyloseq object
metadata <- data.frame(sample_data(ps_final))

# Combine read counts with species info
sample_sum_df <- cbind(sample_sum_df, metadata[rownames(metadata) %in% sample_sum_df$SampleID, ])


# Plot histogram
ggplot(sample_sum_df, aes(x = Reads)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# Save to CSV
write.csv(sample_sum_df, "ITS_Reads_Per_Sample.csv", row.names = FALSE)