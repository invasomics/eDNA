#=========================================================
# 0. SETUP
#=========================================================

setwd("C:/add_here")

init_packages <- function() {
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install(
    c("dada2","DECIPHER","decontam"),
    ask = FALSE,
    update = FALSE
  )
  
  pkgs <- c(
    "phyloseq","vegan","ggplot2","dplyr","tidyr","tibble",
    "readr","readxl","ape","microbiome","lme4",
    "broom","ggpubr","ComplexUpset","betapart"
  )
  
  invisible(lapply(pkgs, require, character.only = TRUE))
  
  theme_set(theme_bw())
}

init_packages()

library(phyloseq)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(plyr)
library(dplyr)
library(readr)
library(readxl)
library(vegan)
library(car)
library(tidyverse)
library(ape)
library(microbiome)
library(gridExtra)
library(ggpubr)
library(ggeasy)
library(MiscMetabar)
library(decontam)
library(lme4)
library(broom)
library(janitor)
library(ComplexUpset)
library(betapart)

#=========================================================
# 1. LOAD DATA + BUILD PHYLOSEQ
#=========================================================

asv_mat <- read_tsv("ASVs_counts_Ant16S.tsv") %>%
  column_to_rownames("ASV") %>%
  as.matrix()

tax_mat <- read.csv("MergedTaxonomy.csv", header = TRUE) %>%
  column_to_rownames("ASV") %>%
  as.matrix()

samples_df <- read_excel("AntarcticBacteria_metadata.xlsx") %>%
  column_to_rownames("sample")

# Match ASVs
common_asvs <- intersect(rownames(asv_mat), rownames(tax_mat))

asv_mat <- asv_mat[common_asvs, , drop = FALSE]
tax_mat <- tax_mat[common_asvs, , drop = FALSE]

# Match samples
common_samples <- intersect(colnames(asv_mat), rownames(samples_df))

asv_mat <- asv_mat[, common_samples, drop = FALSE]
samples_df <- samples_df[common_samples, , drop = FALSE]

# Prefix sample names
colnames(asv_mat) <- paste0("S", colnames(asv_mat))
rownames(samples_df) <- paste0("S", rownames(samples_df))

ASV <- otu_table(asv_mat, taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax_mat))
SAM <- sample_data(samples_df)

ps <- phyloseq(ASV, TAX, SAM)

meta <- as(sample_data(ps), "data.frame") %>%
  mutate(
    sample_name = rownames(.),
    
    swab_type = case_when(
      sample_type == "Soil" ~ "Soil",
      sample_type %in% c("Boot", "Tent", "Box", "Helo") ~ sample_type,
      TRUE ~ NA_character_
    ),
    
    is_control = sample_type %in% c("Ex_Control", "PCR_Control"),
    
    day = factor(day, levels = c(1, 2, 3)),
    
    travel_order = factor(travel_order)
  )

sample_data(ps) <- sample_data(meta)

ps

#=========================================================
# 2. CLEANING MODULE (CONTROLS + DECONTAM)
#=========================================================

clean_ps <- function(ps) {
  
  # Remove empty taxa
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  # Identify contaminants
  contam_df <- isContaminant(
    ps,
    method = "prevalence",
    neg = "is_control"
  )
  
  # Safety check
  if (!"contaminant" %in% colnames(contam_df)) {
    stop("decontam failed: no 'contaminant' column found")
  }
  
  # Summarise contamination
  n_contam <- sum(contam_df$contaminant, na.rm = TRUE)
  
  message("Contaminants identified: ", n_contam)
  
  # Remove contaminants
  ps_clean <- prune_taxa(!contam_df$contaminant, ps)
  
  # Remove zero-count taxa again
  ps_clean <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)
  
  list(
    ps_clean = ps_clean,
    contam_table = contam_df,
    n_contaminants = n_contam
  )
}

cleaned <- clean_ps(ps)

ps_clean <- cleaned$ps_clean
contam_table <- cleaned$contam_table

ps_clean

#=========================================================
# 2.1 ANALYTICAL DATASETS (POST-DECONTAM)
#=========================================================

#---------------------------------------------------------
# Remove controls first
#---------------------------------------------------------

ps_no_controls <- subset_samples(
  ps_clean,
  !sample_type %in% c(
    "Ex_Control",
    "PCR_Control"
  )
)

ps_no_controls <- prune_samples(
  sample_sums(ps_no_controls) > 0,
  ps_no_controls
)

ps_no_controls <- prune_taxa(
  taxa_sums(ps_no_controls) > 0,
  ps_no_controls
)

#=========================================================
# MICROBIAL TRANSFER DATASET
# Excludes Box + Tent
# Retains Soil, Boot, Helo
#=========================================================

ps_microbial_dataset <- subset_samples(
  ps_no_controls,
  
  !sample_type %in% c(
    "Box",
    "Tent"
  ) &
    
    !is.na(sample_type)
)

ps_microbial_dataset <- prune_samples(
  sample_sums(ps_microbial_dataset) > 0,
  ps_microbial_dataset
)

ps_microbial_dataset <- prune_taxa(
  taxa_sums(ps_microbial_dataset) > 0,
  ps_microbial_dataset
)

#=========================================================
# DECONTAMINATION DATASET
# Excludes Boot + Helo
# Retains Soil, Box, Tent
#=========================================================

ps_decontamination_dataset <- subset_samples(
  ps_no_controls,
  
  !sample_type %in% c(
    "Boot",
    "Helo"
  ) &
    
    !is.na(sample_type)
)

ps_decontamination_dataset <- prune_samples(
  sample_sums(ps_decontamination_dataset) > 0,
  ps_decontamination_dataset
)

ps_decontamination_dataset <- prune_taxa(
  taxa_sums(ps_decontamination_dataset) > 0,
  ps_decontamination_dataset
)

#---------------------------------------------------------
# Additional metadata subsets
#---------------------------------------------------------

ps_location <- subset_samples(
  ps_no_controls,
  !is.na(location)
)

ps_location <- prune_samples(
  sample_sums(ps_location) > 0,
  ps_location
)

ps_location <- prune_taxa(
  taxa_sums(ps_location) > 0,
  ps_location
)

ps_sample_condition <- subset_samples(
  ps_no_controls,
  !is.na(sample_condition)
)

ps_sample_condition <- prune_samples(
  sample_sums(ps_sample_condition) > 0,
  ps_sample_condition
)

ps_sample_condition <- prune_taxa(
  taxa_sums(ps_sample_condition) > 0,
  ps_sample_condition
)

ps_day <- subset_samples(
  ps_no_controls,
  !is.na(day)
)

ps_day <- prune_samples(
  sample_sums(ps_day) > 0,
  ps_day
)

ps_day <- prune_taxa(
  taxa_sums(ps_day) > 0,
  ps_day
)

ps_travel_order <- subset_samples(
  ps_no_controls,
  !is.na(travel_order)
)

ps_travel_order <- prune_samples(
  sample_sums(ps_travel_order) > 0,
  ps_travel_order
)

ps_travel_order <- prune_taxa(
  taxa_sums(ps_travel_order) > 0,
  ps_travel_order
)

#---------------------------------------------------------
# Check objects
#---------------------------------------------------------

ps_no_controls

ps_microbial_dataset
ps_decontamination_dataset

ps_location
ps_sample_condition
ps_day
ps_travel_order

#=========================================================
# 3. FILTERING + NORMALISATION
#=========================================================

filter_ps <- function(ps) {
  
  # Remove controls
  ps <- subset_samples(
    ps,
    !sample_type %in% c("Ex_Control", "PCR_Control")
  )
  
  # Remove missing metadata
  ps <- prune_samples(
    !is.na(sample_data(ps)$sample_condition) &
      !is.na(sample_data(ps)$location),
    ps
  )
  
  # Remove empty samples
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  
  # Remove empty taxa
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  return(ps)
}

#---------------------------------------------------------
# Filter datasets
#---------------------------------------------------------

ps_filt <- filter_ps(ps_clean)

ps_microbial_dataset_filt <- filter_ps(
  ps_microbial_dataset
)

ps_decontamination_dataset_filt <- filter_ps(
  ps_decontamination_dataset
)

#---------------------------------------------------------
# Rarefaction
#---------------------------------------------------------

set.seed(7)

ps_rare <- rarefy_even_depth(
  ps_filt,
  sample.size = min(sample_sums(ps_filt)),
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

ps_microbial_dataset_rare <- rarefy_even_depth(
  ps_microbial_dataset_filt,
  sample.size = min(sample_sums(ps_microbial_dataset_filt)),
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

ps_decontamination_dataset_rare <- rarefy_even_depth(
  ps_decontamination_dataset_filt,
  sample.size = min(sample_sums(ps_decontamination_dataset_filt)),
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

#---------------------------------------------------------
# Check objects
#---------------------------------------------------------

ps_rare

ps_microbial_dataset_rare
ps_decontamination_dataset_rare

#=========================================================
# 4. ALPHA DIVERSITY
#=========================================================

alpha_div <- estimate_richness(
  ps_filt,
  measures = c("Observed", "Shannon")
)

#---------------------------------------------------------
# MASTER alpha diversity dataframe
# KEEP ALL SAMPLES
#---------------------------------------------------------

alpha_df <- cbind(
  as(sample_data(ps_filt), "data.frame"),
  alpha_div
)

# Restore metadata factors
alpha_df <- alpha_df %>%
  mutate(
    location = factor(location),
    day = factor(day),
    travel_order = factor(travel_order),
    sample_condition = factor(sample_condition),
    sample_type = factor(sample_type)
  )

#=========================================================
# STAGE-SPECIFIC DATAFRAME
# ONLY used for decontamination analyses
#=========================================================

alpha_stage_df <- alpha_df %>%
  mutate(
    Stage = case_when(
      sample_condition == "BeforeDeployment" ~ "Before",
      sample_condition == "AfterDeployment" ~ "Post",
      sample_condition == "AfterPostFieldCleaning" ~ "Clean",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Stage)) %>%
  droplevels()

alpha_stage_df <- alpha_stage_df %>%
  mutate(pair_id = sample_name)

alpha_stage_df$Stage <- factor(
  alpha_stage_df$Stage,
  levels = c("Before", "Post", "Clean")
)

#---------------------------------------------------------
# STAGE-SPECIFIC DATAFRAME
# ONLY used for microbial transfer analyses
#---------------------------------------------------------

alpha_micro <- estimate_richness(
  ps_microbial_dataset,
  measures = c("Observed", "Shannon")
)

alpha_microbial_dataset <- cbind(
  as(sample_data(ps_microbial_dataset), "data.frame"),
  alpha_micro
)

# Restore metadata factors
alpha_microbial_dataset <- alpha_microbial_dataset %>%
  mutate(
    location = factor(location),
    day = factor(day),
    travel_order = factor(travel_order),
    sample_condition = factor(sample_condition),
    sample_type = factor(
      sample_type,
      levels = c("Soil", "Boot", "Helo")
    )
  )

#=========================================================
# 4.0 ALPHA DIVERSITY — MICROBIAL DATASET
#=========================================================

# Set factor order
alpha_microbial_dataset$sample_type <- factor(
  alpha_microbial_dataset$sample_type,
  levels = c("Boot", "Helo", "Soil")
)

# Observed richness
ggplot(alpha_microbial_dataset,
       aes(x = sample_type,
           y = Observed,
           fill = sample_type,
           colour = sample_type)) +
  
  geom_boxplot(
    width = 0.5,
    alpha = 0.4,
    outlier.shape = NA
  ) +
  
  geom_point(
    size = 2,
    alpha = 0.9,
    position = position_jitter(width = 0.15)
  ) +
  
  theme_classic(base_size = 12) +
  
  labs(
    x = NULL,
    y = "ASV richness",
    title = "Observed ASV richness across sample types"
  )

# Shannon diversity
ggplot(alpha_microbial_dataset,
       aes(x = sample_type,
           y = Shannon,
           fill = sample_type,
           colour = sample_type)) +
  
  geom_boxplot(
    width = 0.5,
    alpha = 0.4,
    outlier.shape = NA
  ) +
  
  geom_point(
    size = 2,
    alpha = 0.9,
    position = position_jitter(width = 0.15)
  ) +
  
  theme_classic(base_size = 12) +
  
  labs(
    x = NULL,
    y = "Shannon diversity",
    title = "Shannon diversity across sample types"
  )

# ANOVA
aov_obs_sampletype <- aov(
  Observed ~ sample_type,
  data = alpha_microbial_dataset
)

summary(aov_obs_sampletype)

aov_shannon_sampletype <- aov(
  Shannon ~ sample_type,
  data = alpha_microbial_dataset
)

summary(aov_shannon_sampletype)

# Tukey HSD
TukeyHSD(aov_obs_sampletype)
TukeyHSD(aov_shannon_sampletype)

#=========================================================
# 4.1 ALPHA DIVERSITY — DECONTAMINATION STATUS × SAMPLE TYPE
#=========================================================

# Observed richness
ggplot(alpha_stage_df, aes(x = Stage, y = Observed)) +
  
  geom_boxplot(
    width = 0.5,
    alpha = 0.4,
    outlier.shape = NA
  ) +
  
  geom_point(
    size = 2,
    alpha = 0.9
  ) +
  
  facet_wrap(~sample_type, nrow = 1) +
  
  theme_classic(base_size = 12) +
  
  labs(
    x = NULL,
    y = "ASV richness",
    title = "Observed ASV richness across decontamination status"
  )

# Shannon diversity
ggplot(alpha_stage_df, aes(x = Stage, y = Shannon)) +
  
  geom_boxplot(
    width = 0.5,
    alpha = 0.4,
    outlier.shape = NA
  ) +
  
  geom_point(
    size = 2,
    alpha = 0.9
  ) +
  
  facet_wrap(~sample_type, nrow = 1) +
  
  theme_classic(base_size = 12) +
  
  labs(
    x = NULL,
    y = "Shannon diversity",
    title = "Shannon diversity across decontamination status"
  )

# ANOVA
aov_obs_stage <- aov(
  Observed ~ Stage + sample_type,
  data = alpha_stage_df
)

summary(aov_obs_stage)

aov_shannon_stage <- aov(
  Shannon ~ Stage + sample_type,
  data = alpha_stage_df
)

summary(aov_shannon_stage)

# Tukey HSD
TukeyHSD(aov_obs_stage)
TukeyHSD(aov_shannon_stage)

#=========================================================
# 4.2 ALPHA DIVERSITY — LOCATION
#=========================================================

alpha_location_df <- alpha_df %>%
  filter(!is.na(location))

# Observed richness
ggplot(alpha_location_df,
       aes(x = location, y = Observed)) +
  
  geom_boxplot(alpha = 0.4) +
  
  geom_point(
    size = 2,
    alpha = 0.8,
    position = position_jitter(width = 0.2)
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  labs(
    x = "Location",
    y = "ASV richness",
    title = "Observed ASV richness across locations"
  )

# Shannon diversity
ggplot(alpha_location_df,
       aes(x = location, y = Shannon)) +
  
  geom_boxplot(alpha = 0.4) +
  
  geom_point(
    size = 2,
    alpha = 0.8,
    position = position_jitter(width = 0.2)
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  labs(
    x = "Location",
    y = "Shannon diversity",
    title = "Shannon diversity across locations"
  )

# ANOVA
aov_obs_location <- aov(
  Observed ~ location,
  data = alpha_location_df
)

summary(aov_obs_location)

aov_shannon_location <- aov(
  Shannon ~ location,
  data = alpha_location_df
)

summary(aov_shannon_location)

# Tukey HSD
TukeyHSD(aov_obs_location)
TukeyHSD(aov_shannon_location)

#=========================================================
# 4.3 ALPHA DIVERSITY — DAY
#=========================================================

alpha_day_df <- alpha_df %>%
  filter(!is.na(day))

# Observed richness
ggplot(alpha_day_df,
       aes(x = factor(day), y = Observed)) +
  
  geom_boxplot(alpha = 0.4) +
  
  geom_point(
    size = 2,
    alpha = 0.8
  ) +
  
  theme_classic(base_size = 12) +
  
  labs(
    x = "Day",
    y = "ASV richness",
    title = "Observed ASV richness across sampling days"
  )

# Shannon diversity
ggplot(alpha_day_df,
       aes(x = factor(day), y = Shannon)) +
  
  geom_boxplot(alpha = 0.4) +
  
  geom_point(
    size = 2,
    alpha = 0.8
  ) +
  
  theme_classic(base_size = 12) +
  
  labs(
    x = "Day",
    y = "Shannon diversity",
    title = "Shannon diversity across sampling days"
  )

# ANOVA
aov_obs_day <- aov(
  Observed ~ day,
  data = alpha_day_df
)

summary(aov_obs_day)

aov_shannon_day <- aov(
  Shannon ~ day,
  data = alpha_day_df
)

summary(aov_shannon_day)

# Tukey HSD
TukeyHSD(aov_obs_day)
TukeyHSD(aov_shannon_day)

#=========================================================
# 4.4 ALPHA DIVERSITY — SAMPLE CONDITION
#=========================================================

alpha_condition_df <- alpha_df %>%
  filter(!is.na(sample_condition))

# Observed richness
ggplot(alpha_condition_df,
       aes(x = sample_condition, y = Observed)) +
  
  geom_boxplot(alpha = 0.4) +
  
  geom_point(
    size = 2,
    alpha = 0.8
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  labs(
    x = "Sample condition",
    y = "ASV richness",
    title = "Observed ASV richness across sample conditions"
  )

# Shannon diversity
ggplot(alpha_condition_df,
       aes(x = sample_condition, y = Shannon)) +
  
  geom_boxplot(alpha = 0.4) +
  
  geom_point(
    size = 2,
    alpha = 0.8
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  labs(
    x = "Sample condition",
    y = "Shannon diversity",
    title = "Shannon diversity across sample conditions"
  )

# ANOVA
aov_obs_condition <- aov(
  Observed ~ sample_condition,
  data = alpha_condition_df
)

summary(aov_obs_condition)

aov_shannon_condition <- aov(
  Shannon ~ sample_condition,
  data = alpha_condition_df
)

summary(aov_shannon_condition)

# Tukey HSD
TukeyHSD(aov_obs_condition)
TukeyHSD(aov_shannon_condition)

#=========================================================
# 4.5 ALPHA DIVERSITY — TRAVEL ORDER
#=========================================================

alpha_travel_df <- alpha_df %>%
  filter(!is.na(travel_order))

# Observed richness
ggplot(alpha_travel_df,
       aes(x = factor(travel_order), y = Observed)) +
  
  geom_boxplot(alpha = 0.4) +
  
  geom_point(
    size = 2,
    alpha = 0.8
  ) +
  
  theme_classic(base_size = 12) +
  
  labs(
    x = "Travel order",
    y = "ASV richness",
    title = "Observed ASV richness across travel order"
  )

# Shannon diversity
ggplot(alpha_travel_df,
       aes(x = factor(travel_order), y = Shannon)) +
  
  geom_boxplot(alpha = 0.4) +
  
  geom_point(
    size = 2,
    alpha = 0.8
  ) +
  
  theme_classic(base_size = 12) +
  
  labs(
    x = "Travel order",
    y = "Shannon diversity",
    title = "Shannon diversity across travel order"
  )

# ANOVA
aov_obs_travel <- aov(
  Observed ~ travel_order,
  data = alpha_travel_df
)

summary(aov_obs_travel)

aov_shannon_travel <- aov(
  Shannon ~ travel_order,
  data = alpha_travel_df
)

summary(aov_shannon_travel)

# Tukey HSD
TukeyHSD(aov_obs_travel)
TukeyHSD(aov_shannon_travel)

#=========================================================
# 4.6 SHANNON DIVERSITY FIGURES
# Location, Day, Travel Order, and Sample Type
#=========================================================

library(patchwork)

#---------------------------------------------------------
# LOCATION
#---------------------------------------------------------

p_location_shannon <- ggplot(
  alpha_location_df,
  aes(
    x = location,
    y = Shannon,
    colour = location
  )
) +
  
  geom_boxplot(
    alpha = 0.4,
    outlier.shape = NA,
    colour = "black"
  ) +
  
  geom_point(
    size = 2.5,
    alpha = 0.9,
    position = position_jitter(width = 0.15)
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  labs(
    x = "Location",
    y = "Shannon diversity",
    title = "Location"
  ) +
  
  guides(colour = "none")

#---------------------------------------------------------
# DAY
#---------------------------------------------------------

p_day_shannon <- ggplot(
  alpha_day_df,
  aes(
    x = factor(day),
    y = Shannon,
    colour = factor(day)
  )
) +
  
  geom_boxplot(
    alpha = 0.4,
    outlier.shape = NA,
    colour = "black"
  ) +
  
  geom_point(
    size = 2.5,
    alpha = 0.9,
    position = position_jitter(width = 0.15)
  ) +
  
  theme_classic(base_size = 12) +
  
  labs(
    x = "Sampling day",
    y = "Shannon diversity",
    title = "Sampling day"
  ) +
  
  guides(colour = "none")

#---------------------------------------------------------
# TRAVEL ORDER
#---------------------------------------------------------

# Ensure travel order is plotted sequentially
alpha_travel_df$travel_order <- factor(
  alpha_travel_df$travel_order,
  levels = 1:18
)

p_travel_shannon <- ggplot(
  alpha_travel_df,
  aes(
    x = travel_order,
    y = Shannon,
    colour = as.numeric(travel_order)
  )
) +
  
  geom_boxplot(
    alpha = 0.4,
    outlier.shape = NA,
    colour = "black"
  ) +
  
  geom_point(
    size = 2.5,
    alpha = 0.9,
    position = position_jitter(width = 0.15)
  ) +
  
  scale_colour_viridis_c(option = "viridis") +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 8
    )
  ) +
  
  labs(
    x = "Travel order",
    y = "Shannon diversity",
    title = "Travel order"
  ) +
  
  guides(colour = "none")

#---------------------------------------------------------
# SAMPLE TYPE
# Microbial dataset only
#---------------------------------------------------------

alpha_type_df <- alpha_df %>%
  filter(
    sample_type %in% c(
      "Boot",
      "Helo",
      "Soil"
    )
  )

# Set plotting order
alpha_type_df$sample_type <- factor(
  alpha_type_df$sample_type,
  levels = c("Boot", "Helo", "Soil")
)

p_type_shannon <- ggplot(
  alpha_type_df,
  aes(
    x = sample_type,
    y = Shannon,
    colour = sample_type
  )
) +
  
  geom_boxplot(
    alpha = 0.4,
    outlier.shape = NA,
    colour = "black"
  ) +
  
  geom_point(
    size = 2.5,
    alpha = 0.9,
    position = position_jitter(width = 0.15)
  ) +
  
  theme_classic(base_size = 12) +
  
  labs(
    x = "Sample type",
    y = "Shannon diversity",
    title = "Sample type"
  ) +
  
  guides(colour = "none")

#---------------------------------------------------------
# COMBINE FIGURES
#---------------------------------------------------------

shannon_combined <- (
  p_location_shannon |
    p_day_shannon
) / (
  p_travel_shannon |
    p_type_shannon
)

shannon_combined

#=========================================================
# 5. BETA DIVERSITY
#=========================================================

#---------------------------------------------------------
# 5.1 Rarefaction
#---------------------------------------------------------

set.seed(7)

ps_location_rare <- rarefy_even_depth(
  ps_location,
  sample.size = min(sample_sums(ps_location)),
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

ps_sample_condition_rare <- rarefy_even_depth(
  ps_sample_condition,
  sample.size = min(sample_sums(ps_sample_condition)),
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

ps_day_rare <- rarefy_even_depth(
  ps_day,
  sample.size = min(sample_sums(ps_day)),
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

ps_travel_order_rare <- rarefy_even_depth(
  ps_travel_order,
  sample.size = min(sample_sums(ps_travel_order)),
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

ps_soil_boot_helo_rare <- rarefy_even_depth(
  ps_soil_boot_helo,
  sample.size = min(sample_sums(ps_soil_boot_helo)),
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

ps_soil_box_tent_rare <- rarefy_even_depth(
  ps_soil_box_tent,
  sample.size = min(sample_sums(ps_soil_box_tent)),
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

#---------------------------------------------------------
# 5.2 MDS PLOTS
#---------------------------------------------------------

make_mds_plot <- function(ps, colour_var, title) {
  
  ord <- ordinate(
    ps,
    method = "MDS",
    distance = "bray"
  )
  
  eigs <- ord$values$Eigenvalues
  
  p <- plot_ordination(
    ps,
    ord,
    color = colour_var
  ) +
    
    geom_point(
      size = 3,
      alpha = 0.9
    ) +
    
    theme_classic(base_size = 12) +
    
    coord_fixed(
      sqrt(eigs[2] / eigs[1])
    ) +
    
    labs(
      color = colour_var,
      title = title
    )
  
  # Add ellipses for categorical variables
  if (colour_var != "travel_order") {
    
    p <- p +
      stat_ellipse(
        type = "t",
        linewidth = 0.6
      )
  }
  
  # Continuous scale for travel order
  if (colour_var == "travel_order") {
    
    p <- p +
      scale_color_viridis_c(
        option = "viridis",
        breaks = 1:18,
        limits = c(1, 18)
      )
  }
  
  return(p)
}

#---------------------------------------------------------
# Ensure travel order is numeric
#---------------------------------------------------------

sample_data(ps_travel_order_rare)$travel_order <- as.numeric(
  as.character(sample_data(ps_travel_order_rare)$travel_order)
)

#---------------------------------------------------------
# Create combined grouping for soil/box/tent condition
#---------------------------------------------------------

sample_data(ps_soil_box_tent_rare)$soil_box_tent_condition <- with(
  data.frame(sample_data(ps_soil_box_tent_rare)),
  
  ifelse(
    sample_type == "soil",
    "Soil",
    paste(sample_type, sample_condition, sep = "_")
  )
)

#---------------------------------------------------------
# Generate plots
#---------------------------------------------------------

mds_location <- make_mds_plot(
  ps_location_rare,
  "location",
  "MDS: Location"
)

mds_sample_condition <- make_mds_plot(
  ps_sample_condition_rare,
  "sample_condition",
  "MDS: Sample condition"
)

mds_day <- make_mds_plot(
  ps_day_rare,
  "day",
  "MDS: Day"
)

mds_travel_order <- make_mds_plot(
  ps_travel_order_rare,
  "travel_order",
  "MDS: Travel order"
)

mds_soil_boot_helo <- make_mds_plot(
  ps_soil_boot_helo_rare,
  "sample_type",
  "MDS: Soil, Boot, Helo"
)

# Custom soil / box / tent plot
ord_soil_box_tent <- ordinate(
  ps_soil_box_tent_rare,
  method = "MDS",
  distance = "bray"
)

eigs_soil_box_tent <- ord_soil_box_tent$values$Eigenvalues

mds_soil_box_tent <- plot_ordination(
  ps_soil_box_tent_rare,
  ord_soil_box_tent,
  color = "soil_box_tent_condition"
) +
  
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  
  # Ellipses by sample type
  stat_ellipse(
    aes(group = sample_type, color = sample_type),
    type = "t",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  theme_classic(base_size = 12) +
  
  coord_fixed(
    sqrt(eigs_soil_box_tent[2] / eigs_soil_box_tent[1])
  ) +
  
  labs(
    color = "Condition",
    title = "MDS: Soil, Box, Tent by Condition"
  )

#---------------------------------------------------------
# Print plots
#---------------------------------------------------------

mds_location
mds_sample_condition
mds_day
mds_travel_order
mds_soil_boot_helo
mds_soil_box_tent

#---------------------------------------------------------
# 5.3 PERMANOVA
#---------------------------------------------------------

run_permanova <- function(ps, variable) {
  
  dist <- phyloseq::distance(
    ps,
    method = "bray"
  )
  
  meta <- as(
    sample_data(ps),
    "data.frame"
  )
  
  form <- as.formula(
    paste("dist ~", variable)
  )
  
  fit <- adonis2(
    form,
    data = meta,
    permutations = 999
  )
  
  print(fit)
  
  return(fit)
}

permanova_location <- run_permanova(
  ps_location_rare,
  "location"
)

permanova_sample_condition <- run_permanova(
  ps_sample_condition_rare,
  "sample_condition"
)

permanova_day <- run_permanova(
  ps_day_rare,
  "day"
)

permanova_travel_order <- run_permanova(
  ps_travel_order_rare,
  "travel_order"
)

permanova_soil_boot_helo <- run_permanova(
  ps_soil_boot_helo_rare,
  "sample_type"
)

permanova_soil_box_tent <- run_permanova(
  ps_soil_box_tent_rare,
  "sample_type"
)

#---------------------------------------------------------
# 5.4 BETADISPER
#---------------------------------------------------------

run_betadisper <- function(ps, variable) {
  
  dist <- phyloseq::distance(
    ps,
    method = "bray"
  )
  
  meta <- as(
    sample_data(ps),
    "data.frame"
  )
  
  bd <- betadisper(
    dist,
    meta[[variable]]
  )
  
  perm <- permutest(
    bd,
    permutations = 999
  )
  
  print(perm)
  
  list(
    model = bd,
    test = perm
  )
}

betadisp_location <- run_betadisper(
  ps_location_rare,
  "location"
)

betadisp_sample_condition <- run_betadisper(
  ps_sample_condition_rare,
  "sample_condition"
)

betadisp_day <- run_betadisper(
  ps_day_rare,
  "day"
)

betadisp_travel_order <- run_betadisper(
  ps_travel_order_rare,
  "travel_order"
)

betadisp_soil_boot_helo <- run_betadisper(
  ps_soil_boot_helo_rare,
  "sample_type"
)

betadisp_soil_box_tent <- run_betadisper(
  ps_soil_box_tent_rare,
  "sample_type"
)

#---------------------------------------------------------
# 5.5 DISPERSION PLOTS
#---------------------------------------------------------

par(mfrow = c(2,3))

boxplot(
  betadisp_location$model,
  main = "Dispersion: location",
  las = 2
)

boxplot(
  betadisp_sample_condition$model,
  main = "Dispersion: sample_condition",
  las = 2
)

boxplot(
  betadisp_day$model,
  main = "Dispersion: day",
  las = 2
)

boxplot(
  betadisp_travel_order$model,
  main = "Dispersion: travel_order",
  las = 2
)

boxplot(
  betadisp_soil_boot_helo$model,
  main = "Dispersion: Soil/Boot/Helo",
  las = 2
)

boxplot(
  betadisp_soil_box_tent$model,
  main = "Dispersion: Soil/Box/Tent",
  las = 2
)
#=========================================================
# 6. ASV SHARING STRUCTURE
#=========================================================

ps_upset <- subset_samples(
  ps_clean,
  
  sample_type %in% c("Tent", "Box") &
    
    sample_condition %in% c(
      "BeforeDeployment",
      "AfterDeployment",
      "AfterPostFieldCleaning"
    )
)

ps_upset <- prune_taxa(
  taxa_sums(ps_upset) > 0,
  ps_upset
)

otu <- otu_table(ps_upset)

if (taxa_are_rows(otu)) {
  otu <- t(otu)
}

otu <- as.matrix(otu)

otu[otu > 0] <- 1

otu <- as.data.frame(otu)

meta <- as(sample_data(ps_upset), "data.frame")

meta$sample <- rownames(meta)

otu$sample <- rownames(otu)

otu_long <- otu %>%
  
  pivot_longer(
    -sample,
    names_to = "ASV",
    values_to = "presence"
  ) %>%
  
  filter(presence > 0) %>%
  
  left_join(
    meta[, c(
      "sample",
      "sample_condition",
      "sample_type"
    )],
    by = "sample"
  ) %>%
  
  distinct(
    sample,
    ASV,
    sample_condition,
    sample_type
  ) %>%
  
  mutate(
    group = paste(
      sample_condition,
      sample_type,
      sep = "_"
    )
  )

upset_matrix <- otu_long %>%
  
  mutate(value = 1) %>%
  
  group_by(ASV, group) %>%
  
  summarise(
    value = 1,
    .groups = "drop"
  ) %>%
  
  pivot_wider(
    names_from = group,
    values_from = value,
    values_fill = 0
  )

upset_matrix <- as.data.frame(upset_matrix)

rownames(upset_matrix) <- upset_matrix$ASV

upset_matrix$ASV <- NULL

upset_matrix[] <- lapply(
  upset_matrix,
  as.numeric
)

upset(
  upset_matrix,
  intersect = colnames(upset_matrix),
  name = "ASVs",
  min_size = 5,
  
  base_annotations = list(
    "Intersection size" = intersection_size()
  ),
  
  set_sizes = upset_set_size()
)

#=========================================================
# 7. DIRECTIONAL ASV CHANGE
#=========================================================

otu_pa <- otu_table(ps_clean)

if (taxa_are_rows(otu_pa)) {
  otu_pa <- t(otu_pa)
}

otu_pa <- as.matrix(otu_pa)

otu_pa[otu_pa > 0] <- 1

otu_pa <- as.data.frame(otu_pa)

meta <- as(sample_data(ps_clean), "data.frame")

meta$sample <- rownames(meta)

otu_pa$sample <- rownames(otu_pa)

otu_long <- otu_pa %>%
  
  pivot_longer(
    -sample,
    names_to = "ASV",
    values_to = "presence"
  ) %>%
  
  filter(presence > 0) %>%
  
  left_join(
    meta[, c(
      "sample",
      "sample_type",
      "sample_condition"
    )],
    by = "sample"
  ) %>%
  
  distinct(
    sample,
    ASV,
    sample_type,
    sample_condition
  )

asv_stage <- otu_long %>%
  
  mutate(
    stage = case_when(
      sample_condition == "BeforeDeployment" ~ "Before",
      sample_condition == "AfterDeployment" ~ "Post",
      sample_condition == "AfterPostFieldCleaning" ~ "Clean",
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(!is.na(stage)) %>%
  
  distinct(ASV, sample_type, stage)

compare_stages <- function(df, type, stage1, stage2) {
  
  a <- df %>%
    filter(
      sample_type == type,
      stage == stage1
    ) %>%
    pull(ASV) %>%
    unique()
  
  b <- df %>%
    filter(
      sample_type == type,
      stage == stage2
    ) %>%
    pull(ASV) %>%
    unique()
  
  data.frame(
    sample_type = type,
    comparison = paste(stage1, "vs", stage2),
    gained = length(setdiff(b, a)),
    lost = length(setdiff(a, b)),
    shared = length(intersect(a, b))
  )
}

box_dir <- compare_stages(
  asv_stage,
  "Box",
  "Before",
  "Post"
)

tent_dir <- bind_rows(
  
  compare_stages(
    asv_stage,
    "Tent",
    "Before",
    "Post"
  ),
  
  compare_stages(
    asv_stage,
    "Tent",
    "Before",
    "Clean"
  ),
  
  compare_stages(
    asv_stage,
    "Tent",
    "Post",
    "Clean"
  )
)

directional_df <- bind_rows(
  box_dir,
  tent_dir
)

directional_long <- directional_df %>%
  
  pivot_longer(
    cols = c(gained, lost, shared),
    names_to = "state",
    values_to = "count"
  )

directional_long$comparison <- factor(
  directional_long$comparison,
  levels = c(
    "Before vs Post",
    "Before vs Clean",
    "Post vs Clean"
  )
)

ggplot(
  directional_long,
  aes(
    x = comparison,
    y = count,
    fill = state
  )
) +
  
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.6
  ) +
  
  facet_wrap(
    ~sample_type,
    scales = "free_x"
  ) +
  
  scale_fill_manual(
    values = c(
      gained = "#F2C14E",
      lost = "#E67E22",
      shared = "#4E79A7"
    )
  ) +
  
  theme_classic(base_size = 12) +
  
  labs(
    x = NULL,
    y = "Number of ASVs",
    fill = "ASV fate",
    title = "Directional ASV change across sampling stages"
  )

#=========================================================
# 8. NESTEDNESS / TURNOVER PARTITIONING
#=========================================================

otu_pa <- otu_table(ps_clean)

if (taxa_are_rows(otu_pa)) {
  otu_pa <- t(otu_pa)
}

otu_pa <- as.matrix(otu_pa)

otu_pa <- (otu_pa > 0) * 1

meta <- as(sample_data(ps_clean), "data.frame")

meta$sample <- rownames(meta)

meta <- meta %>%
  
  mutate(
    Stage = case_when(
      sample_condition == "BeforeDeployment" ~ "Before",
      sample_condition == "AfterDeployment" ~ "Post",
      sample_condition == "AfterPostFieldCleaning" ~ "Clean",
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(
    sample_type %in% c("Box", "Tent"),
    !is.na(Stage)
  )

otu_pa <- otu_pa[
  meta$sample,
  ,
  drop = FALSE
]

meta$Stage <- factor(
  meta$Stage,
  levels = c("Before", "Post", "Clean")
)

get_bp <- function(mat, meta_sub, s1, s2) {
  
  m1 <- mat[
    meta_sub$Stage == s1,
    ,
    drop = FALSE
  ]
  
  m2 <- mat[
    meta_sub$Stage == s2,
    ,
    drop = FALSE
  ]
  
  mat_pair <- rbind(m1, m2)
  
  mat_pair <- mat_pair[
    rowSums(mat_pair) > 0,
    ,
    drop = FALSE
  ]
  
  if (nrow(mat_pair) < 2 || ncol(mat_pair) < 1)
    return(NULL)
  
  bp <- betapart.core(mat_pair)
  
  bp <- beta.pair(
    bp,
    index.family = "jaccard"
  )
  
  data.frame(
    turnover = mean(bp$beta.jtu, na.rm = TRUE),
    nestedness = mean(bp$beta.jne, na.rm = TRUE),
    total = mean(bp$beta.jac, na.rm = TRUE)
  )
}

meta_box <- meta %>%
  filter(sample_type == "Box")

meta_tent <- meta %>%
  filter(sample_type == "Tent")

mat_box <- otu_pa[
  meta_box$sample,
  ,
  drop = FALSE
]

mat_tent <- otu_pa[
  meta_tent$sample,
  ,
  drop = FALSE
]

bp_box_pre_post <- get_bp(
  mat_box,
  meta_box,
  "Before",
  "Post"
)

bp_tent_pre_post <- get_bp(
  mat_tent,
  meta_tent,
  "Before",
  "Post"
)

bp_tent_post_clean <- get_bp(
  mat_tent,
  meta_tent,
  "Post",
  "Clean"
)

bp_tent_pre_clean <- get_bp(
  mat_tent,
  meta_tent,
  "Before",
  "Clean"
)

panelC <- bind_rows(
  
  data.frame(
    group = "Box",
    comparison = "Pre vs Post",
    nestedness = bp_box_pre_post$nestedness,
    turnover = bp_box_pre_post$turnover
  ),
  
  data.frame(
    group = "Tent",
    comparison = "Pre vs Post",
    nestedness = bp_tent_pre_post$nestedness,
    turnover = bp_tent_pre_post$turnover
  ),
  
  data.frame(
    group = "Tent",
    comparison = "Pre vs Clean",
    nestedness = bp_tent_pre_clean$nestedness,
    turnover = bp_tent_pre_clean$turnover
  ),
  
  data.frame(
    group = "Tent",
    comparison = "Post vs Clean",
    nestedness = bp_tent_post_clean$nestedness,
    turnover = bp_tent_post_clean$turnover
  )
)

panelC_long <- panelC %>%
  
  pivot_longer(
    c(nestedness, turnover),
    names_to = "component",
    values_to = "value"
  )

panelC_long$comparison <- factor(
  panelC_long$comparison,
  levels = c(
    "Pre vs Post",
    "Pre vs Clean",
    "Post vs Clean"
  )
)

ggplot(
  panelC_long,
  aes(
    x = comparison,
    y = value,
    fill = component
  )
) +
  
  geom_col(position = "dodge") +
  
  facet_wrap(~group) +
  
  scale_fill_manual(
    values = c(
      nestedness = "#4E79A7",
      turnover = "#E15759"
    )
  ) +
  
  theme_classic(base_size = 12) +
  
  labs(
    x = NULL,
    y = "Jaccard dissimilarity",
    title = "Nestedness vs turnover across deployment stages"
  )