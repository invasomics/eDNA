#=========================================================
# 0. SETUP
#=========================================================

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
library(remotes)

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
# Additional metadata subsets on the microbial dataset
#---------------------------------------------------------

ps_location <- subset_samples(
  ps_microbial_dataset,
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

ps_day <- subset_samples(
  ps_microbial_dataset,
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
  ps_microbial_dataset,
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

ps_type <- subset_samples(
  ps_microbial_dataset,
  !is.na(sample_type)
)

ps_type <- prune_samples(
  sample_sums(ps_type) > 0,
  ps_type
)

ps_type <- prune_taxa(
  taxa_sums(ps_type) > 0,
  ps_type
)

#---------------------------------------------------------
# Check objects
#---------------------------------------------------------

ps_no_controls

ps_microbial_dataset
ps_decontamination_dataset

ps_location
ps_day
ps_travel_order
ps_type

#=========================================================
# 3. NORMALISATION/RAREFACTION
#=========================================================

set.seed(7)

ps_microbial_dataset_rare <- rarefy_even_depth(
  ps_microbial_dataset,
  sample.size = min(sample_sums(ps_microbial_dataset)),
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

ps_decontamination_dataset_rare <- rarefy_even_depth(
  ps_decontamination_dataset,
  sample.size = min(sample_sums(ps_decontamination_dataset)),
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

ps_location_rare <- rarefy_even_depth(
  ps_location,
  sample.size = min(sample_sums(ps_location)),
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

ps_type_rare <- rarefy_even_depth(
  ps_type,
  sample.size = min(sample_sums(ps_type)),
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = FALSE
)

#---------------------------------------------------------
# Check objects
#---------------------------------------------------------

ps_microbial_dataset_rare
ps_decontamination_dataset_rare

ps_location_rare
ps_day_rare
ps_travel_order_rare
ps_type_rare

#=========================================================
# 4. ALPHA DIVERSITY
#=========================================================

#---------------------------------------------------------
# 4.1 SET UP FOR MICROBIAL DATASET
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

#---------------------------------------------------------
# 4.2 ALPHA DIVERSITY BY SAMPLE TYPE
#---------------------------------------------------------

# Set factor order
alpha_microbial_dataset$sample_type <- factor(
  alpha_microbial_dataset$sample_type,
  levels = c("Boot", "Helo", "Soil")
)

# Observed richness by sample type
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
# 4.3 ALPHA DIVERSITY — LOCATION
#=========================================================

alpha_location_df <- alpha_microbial_dataset %>%
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
# 4.4 ALPHA DIVERSITY — DAY
#=========================================================

alpha_day_df <- alpha_microbial_dataset %>%
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
# 4.5 ALPHA DIVERSITY — TRAVEL ORDER
#=========================================================

alpha_travel_df <- alpha_location_df %>%
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
alpha_location_df$location <- factor(
  alpha_location_df$location,
  levels = c(
    "Scott_Base_Station",
    "Minna_Bluff",
    "Mier_Valley",
    "Mier_Ridge",
    "Commonwealth_Glacier",
    "Botany_Bay_Granite_Harbour",
    "Haystack_Ridge",
    "Victoria_Valley",
    "Marble_Point",
    "Don_Juan_Pond",
    "Bull_Pass"
  )
)

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
  
  scale_x_discrete(
    limits = c(
      "Scott_Base_Station",
      "Minna_Bluff",
      "Mier_Valley",
      "Mier_Ridge",
      "Commonwealth_Glacier",
      "Botany_Bay_Granite_Harbour",
      "Haystack_Ridge",
      "Victoria_Valley",
      "Marble_Point",
      "Don_Juan_Pond",
      "Bull_Pass"
    ),
    
    labels = c(
      "Scott_Base_Station" = "SBS",
      "Minna_Bluff" = "MB",
      "Mier_Valley" = "MV",
      "Mier_Ridge" = "MR",
      "Commonwealth_Glacier" = "CG",
      "Botany_Bay_Granite_Harbour" = "BBGH",
      "Haystack_Ridge" = "HR",
      "Victoria_Valley" = "VV",
      "Marble_Point" = "MP",
      "Don_Juan_Pond" = "DJP",
      "Bull_Pass" = "BP"
    )
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  labs(
    x = "",
    y = "Alpha diversity (Shannon)"
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
  
  scale_x_discrete(
    limits = c(
      "1",
      "2",
      "3"
    ),
    
    labels = c(
      "1" = "Day 1",
      "2" = "Day 2",
      "3" = "Day 3"
    ),
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  labs(
    x = "",
    y = "",
    title = ""
  ) +
  
  guides(colour = "none")

#---------------------------------------------------------
# TRAVEL ORDER
#---------------------------------------------------------

# Ensure travel order is plotted sequentially
alpha_travel_df$travel_order <- factor(
  alpha_travel_df$travel_order,
  levels = 1:14
)

alpha_travel_df <- alpha_travel_df |>
  dplyr::filter(!is.na(travel_order))

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
  
  scale_x_discrete(
    limits = c(
      "2",
      "3",
      "4",
      "5",
      "6",
      "7",
      "8",
      "9",
      "10",
      "11",
      "12",
      "13",
      "14"
    ),
    
    labels = c(
      "2" = "1a: SBS",
      "3" = "1b: MB",
      "4" = "1c: MV",
      "5" = "1d: MR",
      "6" = "1e: CG",
      "7" = "2a SBS",
      "8" = "2b: BBGH",
      "9" = "2c: HR",
      "10" = "2d: VV",
      "11" = "2e: MP",
      "12" = "3a: SBS",
      "13" = "3b: DJP",
      "14" = "3c: BP"
    )
  ) +
  
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
    )
  ) +
  
  labs(
    x = "",
    y = "Alpha diversity (Shannon)"
  ) +
  
  guides(colour = "none")

#---------------------------------------------------------
# SAMPLE TYPE
#---------------------------------------------------------

alpha_type_df <- alpha_microbial_dataset %>%
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
  
  scale_x_discrete(
    limits = c(
      "Boot",
      "Helo",
      "Soil"
    ),
    
    labels = c(
      "Boot" = "Boot",
      "Helo" = "Helicopter",
      "Soil" = "Soil"
    )
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
    )
  ) +
  
  labs(
    x = "",
    y = "",
    title = ""
  ) +
  
  guides(colour = "none")

#---------------------------------------------------------
# COMBINE FIGURES - FIGURE 2
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
# 4.7 DECONTAMINATION STAGE
#=========================================================

# Sample decontamination - sample type

alpha_decontam <- estimate_richness(
  ps_decontamination_dataset,
  measures = c("Observed", "Shannon")
)

alpha_decontam_dataset <- cbind(
  as(sample_data(ps_decontamination_dataset), "data.frame"),
  alpha_decontam
)

# Restore metadata factors
alpha_decontam_dataset <- alpha_decontam_dataset %>%
  mutate(
    Type = sample_type
    )

alpha_decontam_dataset <- alpha_decontam_dataset %>%
  mutate(pair_id = sample_name)

# Set factor order
alpha_decontam_dataset$Type <- factor(
  alpha_decontam_dataset$Type,
  levels = c("Tent", "Box", "Soil")
)

# Observed richness
ggplot(alpha_decontam_dataset, aes(x = Type, y = Observed, colour = Type)) +
  
  geom_boxplot(
    width = 0.5,
    alpha = 0.4,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8
  ) +
  
  theme_classic(base_size = 12) +
  
  labs(
    x = NULL,
    y = "ASV richness",
    title = "Observed ASV richness"
  )

# Shannon diversity
ggplot(alpha_decontam_dataset, aes(x = Type, y = Shannon, colour = Type)) +
  
  geom_boxplot(
    width = 0.5,
    alpha = 0.4,
    outlier.shape = NA
  ) +
  
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8
  ) +
  
  theme_classic(base_size = 12) +
  
  labs(
    x = NULL,
    y = "Shannon diversity",
    title = "Shannon diversity"
  )

# ANOVA
aov_obs_type <- aov(
  Observed ~ Type,
  data = alpha_decontam_dataset
)

summary(aov_obs_type)

aov_shannon_type <- aov(
  Shannon ~ Type,
  data = alpha_decontam_dataset
)

summary(aov_shannon_type)

# Tukey HSD
TukeyHSD(aov_obs_type)
TukeyHSD(aov_shannon_type)

# Sample decontamination stage

alpha_decontam <- estimate_richness(
  ps_decontamination_dataset,
  measures = c("Observed", "Shannon")
)

alpha_decontam_dataset <- cbind(
  as(sample_data(ps_decontamination_dataset), "data.frame"),
  alpha_decontam
)

# Restore metadata factors
alpha_decontam_dataset <- alpha_decontam_dataset %>%
  mutate(
    location = factor(location),
    day = factor(day),
    travel_order = factor(travel_order),
    sample_condition = factor(sample_condition),
    sample_type = factor(
      sample_type,
      levels = c("BeforeDeployment", "AfterDeployment", "AfterPostFieldCleaning")
    )
  )

alpha_decontam_dataset <- alpha_decontam_dataset %>%
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

alpha_decontam_dataset <- alpha_decontam_dataset %>%
  mutate(pair_id = sample_name)

# Set factor order
alpha_decontam_dataset$Stage <- factor(
  alpha_decontam_dataset$Stage,
  levels = c("Before", "Post", "Clean")
)

# Observed richness
ggplot(alpha_decontam_dataset, aes(x = Stage, y = Observed)) +
  
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
ggplot(alpha_decontam_dataset, aes(x = Stage, y = Shannon)) +
  
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
  Observed ~ Stage,
  data = alpha_decontam_dataset
)

summary(aov_obs_stage)

aov_shannon_stage <- aov(
  Shannon ~ Stage,
  data = alpha_decontam_dataset
)

summary(aov_shannon_stage)

# Tukey HSD
TukeyHSD(aov_obs_stage)
TukeyHSD(aov_shannon_stage)

#=========================================================
# 5. BETA DIVERSITY
#=========================================================

#---------------------------------------------------------
# 5.1 MDS PLOTS
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
    
    theme(
      legend.position = "none"
    ) +
    
    labs(
      x = paste0(
        "MDS1 (",
        round(eigs[1] / sum(eigs) * 100, 1),
        "%)"
      ),
      y = paste0(
        "MDS2 (",
        round(eigs[2] / sum(eigs) * 100, 1),
        "%)"
      ),
      title = NULL
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
        breaks = 1:14,
        limits = c(1, 14)
      )
  }
  
  p
}

#---------------------------------------------------------
# Ensure travel order is numeric for plotting
#---------------------------------------------------------

sample_data(ps_travel_order_rare)$travel_order <- as.numeric(
  as.character(sample_data(ps_travel_order_rare)$travel_order)
)

#---------------------------------------------------------
# Create combined grouping for soil/box/tent condition
#---------------------------------------------------------

sample_data(ps_decontamination_dataset_rare)$soil_box_tent_condition <- with(
  data.frame(sample_data(ps_decontamination_dataset_rare)),
  
  ifelse(
    sample_type == "soil",
    "Soil",
    paste(sample_type, sample_condition, sep = "_")
  )
)

#---------------------------------------------------------
# Generate plots
#---------------------------------------------------------

ord_global <- ordinate(
  ps_microbial_dataset_rare,
  method = "MDS",
  distance = "bray"
)

eigs_global <- ord_global$values$Eigenvalues

mds_location <- plot_ordination(
  ps_location_rare,
  ord_global,
  color = "location"
) +
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  
  # Ellipses by sample type
  stat_ellipse(
    aes(group = location, color = location),
    type = "t",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "none",
    plot.title = element_blank()
  ) +
  
  coord_fixed(
    sqrt(eigs_global[2] / eigs_global[1])
  )

mds_day <- plot_ordination(
  ps_day_rare,
  ord_global,
  color = "day"
) +
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  
  # Ellipses by sample type
  stat_ellipse(
    aes(group = day, color = day),
    type = "t",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "none",
    plot.title = element_blank()
  ) +
  
  coord_fixed(
    sqrt(eigs_global[2] / eigs_global[1])
  )

mds_travel_order <- plot_ordination(
  ps_travel_order_rare,
  ord_global,
  color = "travel_order"
) +
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "bottom"
  )
  
  coord_fixed(
    sqrt(eigs_global[2] / eigs_global[1])
  )

mds_sample_type <- plot_ordination(
  ps_type_rare,
  ord_global,
  color = "sample_type"
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
  
  theme(
    legend.position = "none",
    plot.title = element_blank()
  ) +
  
  coord_fixed(
    sqrt(eigs_global[2] / eigs_global[1])
  )


# Custom soil / box / tent plot
ord_decontam <- ordinate(
  ps_decontamination_dataset_rare,
  method = "MDS",
  distance = "bray"
)

eigs_decontam <- ord_decontam$values$Eigenvalues

mds_decontam <- plot_ordination(
  ps_decontamination_dataset_rare,
  ord_decontam,
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
  
  theme(
    legend.position = "bottom"
  ) +
  
  coord_fixed(
    sqrt(eigs_decontam[2] / eigs_decontam[1])
  ) +
  
  labs(
    x = paste0(
      "MDS1 (",
      round(eigs_decontam[1] / sum(eigs_decontam) * 100, 1),
      "%)"
    ),
    y = paste0(
      "MDS2 (",
      round(eigs_decontam[2] / sum(eigs_decontam) * 100, 1),
      "%)"
    ),
    title = NULL
  )

#---------------------------------------------------------
# Print plots - FIGURE 3
#---------------------------------------------------------

mds_combined <- (
  mds_location |
  mds_day
) / (
  mds_travel_order |
  mds_sample_type
)

mds_combined

mds_decontam

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

permanova_sample_type <- run_permanova(
  ps_microbial_dataset_rare,
  "sample_type"
)
permanova_day <- run_permanova(
  ps_day_rare,
  "day"
)

permanova_travel_order <- run_permanova(
  ps_travel_order_rare,
  "travel_order"
)

permanova_decontam <- run_permanova(
  ps_decontamination_dataset_rare,
  "sample_type"
)

permanova_decontam <- run_permanova(
  ps_decontamination_dataset_rare,
  "soil_box_tent_condition"
)


#---------------------------------------------------------
# 5.4 BETADISPER - only for categorical variables
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

betadisp_sample_type <- run_betadisper(
  ps_microbial_dataset_rare,
  "sample_type"
)

betadisp_day <- run_betadisper(
  ps_day_rare,
  "day"
)

betadisp_decontam <- run_betadisper(
  ps_decontamination_dataset_rare,
  "sample_type"
)

betadisp_decontam2 <- run_betadisper(
  ps_decontamination_dataset_rare,
  "soil_box_tent_condition"
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
  betadisp_sample_type$model,
  main = "Dispersion: sample_type",
  las = 2
)

boxplot(
  betadisp_day$model,
  main = "Dispersion: day",
  las = 2
)

boxplot(
  betadisp_decontam$model,
  main = "Dispersion: sample type",
  las = 2
)

boxplot(
  betadisp_decontam2$model,
  main = "Dispersion: stage",
  las = 2
)

#=========================================================
# 6.1 DISTANCE-BASED LINEAR MODELLING (dbRDA) 
#=========================================================

dist_bray <- phyloseq::distance(
  ps_microbial_dataset_rare,
  method = "bray"
)

meta_db <- as(
  sample_data(ps_microbial_dataset_rare),
  "data.frame"
)

# Ensure factors
meta_db$location <- factor(meta_db$location)
meta_db$day <- factor(meta_db$day)
meta_db$sample_type <- factor(meta_db$sample_type)
meta_db$travel_order <- as.numeric(
  as.character(meta_db$travel_order)
)

# Full model

db_full <- capscale(
  dist_bray ~
    location +
    day +
    sample_type +
    travel_order,
  data = meta_db
)

anova(db_full, by = "margin", permutations = 999)


#=========================================================
# 6.2. FORMAL HIERARCHICAL VARIANCE PARTITIONING
#=========================================================

varpart_res <- varpart(
  dist_bray,
  ~ location,
  ~ day,
  ~ sample_type,
  ~ travel_order,
  data = meta_db
)

varpart_res
plot(varpart_res)

#=========================================================
# 6.3. TWO SENSITIVITY ANALYSES FOR ORDINATION ANALYSES
#=========================================================

# Remove rare ASVs before ordination and rerun PERMANOVA/dbRDA.

ps_prev <- filter_taxa(
  ps_microbial_dataset_rare,
  function(x) sum(x > 0) >= 5,
  prune = TRUE
)

ntaxa(ps_prev)

meta_prev <- data.frame(
  sample_data(ps_prev)
)

meta_prev$travel_order <- as.numeric(
  as.character(meta_prev$travel_order)
)

bray_prev <- phyloseq::distance(
  ps_prev,
  method = "bray"
)

adonis_prev <- adonis2(
  bray_prev ~
    sample_type +
    location +
    travel_order,
  data = meta_prev,
  permutations = 999
)

adonis_prev

dbrda_prev <- capscale(
  bray_prev ~
    sample_type +
    location +
    travel_order,
  data = meta_prev
)

anova(
  dbrda_prev,
  permutations = 999
)

anova(
  dbrda_prev,
  by = "margin",
  permutations = 999
)

plot(
  dbrda_prev,
  display = c("sites")
)

scores_prev <- scores(
  dbrda_prev,
  display = "sites"
)

scores_prev <- as.data.frame(
  scores_prev
)

scores_prev$sample_type <- meta_prev$sample_type

ggplot(
  scores_prev,
  aes(
    CAP1,
    CAP2,
    colour = sample_type
  )
) +
  geom_point(size = 3) +
  theme_classic()

# CLR + Euclidean distance

otu_mat <- as(
  otu_table(ps_microbial_dataset),
  "matrix"
)

if (taxa_are_rows(ps_microbial_dataset)) {
  otu_mat <- t(otu_mat)
}

otu_mat_clr <- otu_mat + 1

library(compositions)

clr_mat <- clr(
  otu_mat_clr
)

aitchison_dist <- dist(
  clr_mat,
  method = "euclidean"
)

meta_clr <- data.frame(
  sample_data(ps_microbial_dataset)
)

meta_clr$travel_order <- as.numeric(
  as.character(meta_clr$travel_order)
)

meta_clr <- meta_clr[
  match(
    rownames(clr_mat),
    rownames(meta_clr)
  ),
]

adonis_clr <- adonis2(
  aitchison_dist ~
    sample_type +
    location +
    travel_order,
  data = meta_clr,
  permutations = 999
)

adonis_clr

dbrda_clr <- capscale(
  aitchison_dist ~
    sample_type +
    location +
    travel_order,
  data = meta_clr
)

anova(
  dbrda_clr,
  permutations = 999
)

anova(
  dbrda_clr,
  by = "margin",
  permutations = 999
)

scores_clr <- scores(
  dbrda_clr,
  display = "sites"
)

scores_clr <- as.data.frame(
  scores_clr
)

scores_clr$sample_type <-
  meta_clr$sample_type

ggplot(
  scores_clr,
  aes(
    CAP1,
    CAP2,
    colour = sample_type
  )
) +
  geom_point(size = 3) +
  theme_classic()

#=========================================================
# 7. NEW MICROBIAL TRANSFER ANALYSES - Indicator ASVs
#=========================================================

# 7.1 Route-specific Indicator Carry-over

# How many taxa at a site were already characteristic of upstream locations?

# Metadata
meta <- data.frame(
  sample_data(ps_microbial_dataset)
)

meta$SampleID <- rownames(meta)

# OTU matrix
otu_mat <- as(
  otu_table(ps_microbial_dataset),
  "matrix"
)

if (taxa_are_rows(ps_microbial_dataset)) {
  otu_mat <- t(otu_mat)
}

# Filter rare ASVs (keep only those present in at least 5 samples)
keep_asvs <- colSums(
  otu_mat > 0
) >= 5

otu_filt <- otu_mat[
  ,
  keep_asvs
]

# Define routes
routes <- list(
  
  Day1 = c(
    "Scott_Base_Station",
    "Minna_Bluff",
    "Mier_Valley",
    "Mier_Ridge",
    "Commonwealth_Glacier"
  ),
  
  Day2 = c(
    "Scott_Base_Station",
    "Botany_Bay_Granite_Harbour",
    "Haystack_Ridge",
    "Victoria_Valley",
    "Marble_Point"
  ),
  
  Day3 = c(
    "Scott_Base_Station",
    "Don_Juan_Pond",
    "Bull_Pass"
  )
  
)

# Calculate upstream indicator carry-over
route_results <- list()

counter <- 1

for(route_name in names(routes)) {
  
  route <- routes[[route_name]]
  
  for(i in 2:length(route)) {
    
    current_location <- route[i]
    
    upstream_locations <- route[1:(i - 1)]
    
    upstream_samples <- meta$SampleID[
      meta$location %in%
        upstream_locations
    ]
    
    current_samples <- meta$SampleID[
      meta$location ==
        current_location
    ]
    
    upstream_asvs <- colnames(
      otu_filt
    )[
      colSums(
        otu_filt[
          upstream_samples,
          ,
          drop = FALSE
        ] > 0
      ) >= 2
    ]
    
    for(s in current_samples) {
      
      sample_asvs <- names(
        otu_filt[s, ]
      )[
        otu_filt[s, ] > 0
      ]
      
      carry_prop <-
        length(
          intersect(
            sample_asvs,
            upstream_asvs
          )
        ) /
        max(
          length(sample_asvs),
          1
        )
      
      route_results[[counter]] <-
        data.frame(
          SampleID = s,
          Route = route_name,
          Carryover = carry_prop
        )
      
      counter <- counter + 1
      
    }
  }
}

# Combine
route_df <- do.call(
  rbind,
  route_results
)

route_df <- left_join(
  route_df,
  meta,
  by = "SampleID"
)

# Plot
ggplot(
  route_df,
  aes(
    sample_type,
    Carryover,
    fill = sample_type
  )
) +
  geom_boxplot() +
  facet_wrap(~Route) +
  scale_fill_viridis_d(option = "viridis") +
  theme_classic() +
  labs(
    x = NULL,
    y = "Carry-over"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )

# Model
m4b <- lmer(
  Carryover ~
    sample_type +
    Route +
    (1|location),
  data = route_df
)

summary(m4b)

# 7.4.3 Revised Scott Base-Specific Taxa Persistence

# Metadata
meta <- data.frame(
  sample_data(ps_microbial_dataset)
)

meta$SampleID <- rownames(meta)

meta$travel_order <- as.numeric(
  as.character(meta$travel_order)
)

# Create OTU matrix
otu_mat <- as(
  otu_table(ps_microbial_dataset),
  "matrix"
)

if (taxa_are_rows(ps_microbial_dataset)) {
  otu_mat <- t(otu_mat)
}

# Convert to relative abundances
rel_abund <- sweep(
  otu_mat,
  1,
  rowSums(otu_mat),
  "/"
)

rel_abund[is.na(rel_abund)] <- 0

#Identify the Scott Base soil sample
ps_scott_soil <- subset_samples(
  ps_microbial_dataset,
  location == "Scott_Base_Station" &
    sample_type == "Soil"
)

scott_sample <- sample_names(
  ps_scott_soil
)

scott_sample

# Define Scott Base-specific taxa
# Taxa abundant at Scott Base
# Keep taxa representing at least 0.5% relative abundance.
scott_abund <- rel_abund[
  scott_sample,
]

candidate_taxa <- names(
  scott_abund
)[
  scott_abund > 0.005
]

# Taxa rare elsewhere
# Keep taxa occurring in less than 20% of all other samples
other_samples <- setdiff(
  rownames(rel_abund),
  scott_sample
)

other_prevalence <- colMeans(
  rel_abund[
    other_samples,
    ,
    drop = FALSE
  ] > 0
)

scott_specific_taxa <- candidate_taxa[
  other_prevalence[
    candidate_taxa
  ] < 0.20
]

# Check how many taxa remain
length(scott_specific_taxa)

# Calculate persistence
# Count how many Scott Base-specific taxa occur in each sample.
indicator_df <- data.frame(
  SampleID = rownames(rel_abund)
)

indicator_df$ScottBaseSpecificTaxa <- apply(
  rel_abund,
  1,
  function(x){
    
    present_taxa <- names(x)[x > 0]
    
    length(
      intersect(
        present_taxa,
        scott_specific_taxa
      )
    )
    
  }
)

# Calculate proportional persistence
indicator_df$ScottBaseSpecificProportion <- apply(
  rel_abund,
  1,
  function(x){
    
    present_taxa <- names(x)[x > 0]
    
    length(
      intersect(
        present_taxa,
        scott_specific_taxa
      )
    ) /
      max(
        length(scott_specific_taxa),
        1
      )
    
  }
)

# Add metadata
indicator_df <- left_join(
  indicator_df,
  meta,
  by = "SampleID"
)

# Remove the Scott Base soil sample (otherwise the source sample trivially contains all source taxa)
indicator_plot_df <- indicator_df %>%
  filter(
    SampleID != scott_sample
  )

# Plot
ggplot(
  indicator_plot_df,
  aes(
    travel_order,
    ScottBaseSpecificProportion,
    colour = sample_type
  )
) +
  geom_point(size = 3) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  theme_classic()

# Mixed model
m4a_revised <- lmer(
  ScottBaseSpecificProportion ~
    travel_order +
    sample_type +
    (1 | location),
  data = indicator_plot_df
)

summary(m4a_revised)

# Vector-only analysis
# Do Scott Base-specific taxa persist on operational vectors as sites become progressively more distant from the hub?
indicator_vectors <- indicator_plot_df %>%
  filter(
    sample_type %in%
      c(
        "Boot",
        "Helo"
      )
  )

# Model
m4a_vectors <- lmer(
  ScottBaseSpecificProportion ~
    travel_order +
    sample_type +
    (1 | location),
  data = indicator_vectors
)

summary(m4a_vectors)

#=========================================================
# 8. ASV SHARING STRUCTURE
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
# 9. DIRECTIONAL ASV CHANGE
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
# 10. NESTEDNESS / TURNOVER PARTITIONING
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

