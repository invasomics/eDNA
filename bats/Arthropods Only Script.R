# Following: https://benjjneb.github.io/dada2/tutorial_1_8.html

# Load in the required libraries; some may need to be installed first;
library(BiocManager)
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
library(RColorBrewer)

########################################### 
#### PHYLOSEQ
# Ste the environment theme
theme_set(theme_bw())

# Read in the ASV counts file:
asv_mat<- read_tsv("ASVs_counts_Bats_bespokeCOI.tsv")   # Add ASV to column A header

# Read in the taxonomy file
# It first required some reformatting and other steps run outside of R; Ang provided the working file
tax_mat<- read.csv("FinalMergedTaxonomy_bats.csv", header=T)
samples_df <- read_excel("batDNA_metadata_F.xlsx")  # This is the metadata file; first remove controls plus any individuals that had columns summing 
# to zero or it will cause problems later

#Define the row names from the asv column
asv_mat <- asv_mat %>%
  tibble::column_to_rownames("ASV") 
#Sane for the two other matrices
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample") 

# Transform into matrices asv and tax tables (sample table can be left as data frame)
asv_mat <- as.matrix(asv_mat)
tax_mat <- as.matrix(tax_mat)

# Construct a phyloseq object
ASV <- otu_table(asv_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
samples <- sample_data(samples_df)

ps <- phyloseq(ASV, TAX, samples)
ps
# Get the raw read counts for each sample (columns of the OTU table)
raw_counts <- colSums(ASV)

# Calculate the range of raw read counts
range(raw_counts)
########################################### 

# DATA EXPLORATION
sample_names(ps)
rank_names(ps)
sample_variables(ps)

# Check your phyla & filter for uncharacterised reads
table(tax_table(ps)[, "Phylum"], exclude = NULL)
#subset to diet phyla - just arthropods
ps_arthro <- subset_taxa(ps, !is.na(Phylum) & Phylum %in% c(
  "Arthropoda"
))
ps_arthro 
ps0<-ps_arthro
ps0
# Compute prevalence of each feature, store as data.frame THESE NEXT 3 STEPS NOT REALLY USED
prevdf = apply(X = otu_table(ps0), 
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2), 
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame 
prevdf = data.frame(Prevalence = prevdf, 
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))

#Visualise prevalence of each phyla
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
plyr::ddply(prevdf, "Order", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
plyr::ddply(prevdf, "Class", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
plyr::ddply(prevdf, "Family", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
plyr::ddply(prevdf, "Genus", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
plyr::ddply(prevdf, "Species", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})


#Check ps object & if everything looks OK, rename to ps
ps0 
ps <- ps0
ps
# Create a sample sum table to look at coverage metrics
sample_sums(ps)
sample_sum_df <- data.frame(sum=sample_sums(ps))
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
write.csv(sample_sums(ps),"LibraryReads_Per_Sample.csv") #outputs a coverage csv file

# Look for skew in coverage across sample types
set.seed(711)
DATA.2 <- ps  

df = as.data.frame(sample_data(DATA.2))
df$LibrarySize = sample_sums(DATA.2)
df = df[order(df$LibrarySize),]
df$Index = seq(nrow(df))

#Plot ordered library size & coloured by 'collection method'  # REPLACE COLLECTION METHOD WITH YOUR VARIABLE/S OF INTEREST
level_order2 <- c('Female_Lactating', 'Female_Pregnant', 'Male')  # REPLACE WITH YOUR NAMED VARIABLES OF INTEREST IN A GIVEN COLUMN FROM METADATA FILE

ggplot(data=df, aes(x=Index, y=LibrarySize, colour= reproductive_status))+  # REPLACE COLLECTION METHOD WITH YOUR COLUMN HEADING FOR VARIABLE/S OF INTEREST
  geom_point(size =2.5)+
  facet_wrap(~ factor(reproductive_status, level = level_order2)) +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_color_manual(values = c(  "Female_Lactating" = "#E69F00",   # orange-ish, friendly to colorblind
                                  "Female_Pregnant" = "#0072B2",    # blue
                                  "Male" = "#D55E00" )) +
  scale_y_continuous(trans='sqrt') 

#Plot ordered library size & coloured by 'collection method'  # REPLACE COLLECTION METHOD WITH YOUR VARIABLE/S OF INTEREST
level_order3 <- c('Pre', 'Post')

ggplot(data = df, aes(x = Index, y = LibrarySize, colour = harvest_status)) +
  geom_point(size = 2.5) +
  facet_wrap(
    ~ factor(harvest_status, levels = level_order3), 
    scales = "free_x",
    strip.position = "top",
    ncol = 2,
    nrow = 2
  ) +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  ) +
  scale_color_manual(values = c(
    "Pre"  = "#4CAF50",   # Mid Green
    "Post" = "#654321"    # Dark Brown
  )) +
  scale_y_continuous(trans = 'sqrt') +
  labs(
    title = "Library Size by Harvest Status",
    y = "Library Size (sqrt)",
    x = "Sample Index"
  )



## Exploring the taxonomic diversity
#Gives counts of phyla etc, may or may not be useful
table(phyloseq::tax_table(ps)[, "Phylum"])
table(phyloseq::tax_table(ps)[, "Class"])  # could be useful if looking for a specific group!
table(phyloseq::tax_table(ps)[, "Order"])
table(phyloseq::tax_table(ps)[, "Family"])
table(phyloseq::tax_table(ps)[, "Genus"])
table(phyloseq::tax_table(ps)[, "Species"])

# Get the Species column from the taxonomy table
species_column <- phyloseq::tax_table(ps)[, "Species"]
# Count how many unique entries (species) are present
length(unique(species_column))
# Get the count of each species
species_count <- table(species_column)
# Display the result
species_count


#Basic bar graphs based on taxonomic group
plot_bar(ps, fill = "Class") + 
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")
plot_bar(ps, fill = "Order") + 
  geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")
# THESE GENERALLY GET LESS USEFUL AS YOU CHANGE PHYLUM TO CLASS, ORDER, ETC.

#pie chart of proportions 
# Pie chart of proportions by Phylum
# Extract the phylum-level taxonomy
phylum_table <- table(tax_table(ps)[, "Phylum"])
# Convert to a data frame for plotting
phylum_df <- as.data.frame(phylum_table)
colnames(phylum_df) <- c("Phylum", "Count")
# Compute proportions and create new labels for the legend
phylum_df <- phylum_df %>%
  mutate(Proportion = Count / sum(Count),
         LegendLabel = paste0(Phylum, " (", round(Proportion * 100, 1), "%)"))
# Plot pie chart with percentage-included legend for Phylum
ggplot(phylum_df, aes(x = "", y = Proportion, fill = LegendLabel)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Proportion of Taxa by Phylum") +
  theme_void() +
  theme(legend.position = "right")
# Pie chart of proportions by Family
# Extract the family-level taxonomy
family_table <- table(tax_table(ps)[, "Family"])

# Convert to a data frame for plotting
family_df <- as.data.frame(family_table)
colnames(family_df) <- c("Family", "Count")

# Compute proportions and create new labels for the legend
family_df <- family_df %>%
  mutate(Proportion = Count / sum(Count),
         LegendLabel = paste0(Family, " (", round(Proportion * 100, 1), "%)"))

# Plot pie chart with percentage-included legend for Family
ggplot(family_df, aes(x = "", y = Proportion, fill = LegendLabel)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Proportion of Taxa by Family") +
  theme_void() +
  theme(legend.position = "right")

# Extract the order-level taxonomy
order_table <- table(tax_table(ps)[, "Order"])
# Convert to a data frame for plotting
order_df <- as.data.frame(order_table)
colnames(order_df) <- c("Order", "Count")
# Compute proportions and create new labels for the legend
order_df <- order_df %>%
  mutate(Proportion = Count / sum(Count),
         LegendLabel = paste0(Order, " (", round(Proportion * 100, 1), "%)"))
# Plot pie chart with percentage-included legend
ggplot(order_df, aes(x = "", y = Proportion, fill = LegendLabel)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Proportion of Taxa by Order") +
  theme_void() +
  theme(legend.position = "right")
# Extract the class-level taxonomy
class_table <- table(tax_table(ps)[, "Class"])
# Convert to a data frame for plotting
class_df <- as.data.frame(class_table)
colnames(class_df) <- c("Class", "Count")
# Compute proportions and create new labels for the legend
class_df <- class_df %>%
  mutate(Proportion = Count / sum(Count),
         LegendLabel = paste0(Class, " (", round(Proportion * 100, 1), "%)"))

# Plot pie chart with percentage-included legend for Class
ggplot(class_df, aes(x = "", y = Proportion, fill = LegendLabel)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Proportion of Taxa by Class") +
  theme_void() +
  theme(legend.position = "right")


# IF SOME SAMPLES HAVE REALLY HIGH ABUNDANCE COMPARED TO OTHERS, MAY LIKE TO RAREFY:
# this is just part of data exploration for now, to understand the data
ps_rare1 <- rarefy_even_depth(ps, sample.size=300,rngseed=7,replace=FALSE,trimOTUs=TRUE,verbose=TRUE) # If removed taxon
# Try different values of sample.size
ps_rare1 # run this value to see how many taxa and samples are left when you reduced the dataset down to the abundance set by sample.size
# Replot and see what it looks like
plot_bar(ps_rare1, fill = "Class") + 
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")
# Plot will now show abundance (at x reads) across samples, with x = the value you set for sample.size

## LET'S LOOK AT THE TOP 20 HITS IN TERMS OF 20 relative abundance for unrarefied data w/o humans:
ps_rel_abund <- phyloseq::transform_sample_counts(ps, function(x){x / sum(x)})
ps_rel_abund
head(otu_table(ps_rel_abund))

top20 <- names(sort(taxa_sums(ps_rel_abund), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps_rel_abund, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
ps.top20
tax_table(ps.top20)


# If the above throws an error it's because samples in the metadata spreadsheet have 0 ASV counts for all ASVs and need to be removed from metadata file
plot_bar(ps.top20, x="Sample", fill="Phylum") + facet_wrap(~reproductive_status, scales="free_x") + labs(x="Sample names",y="Relative abundance")
plot_bar(ps.top20, x="Sample", fill="Class") + facet_wrap(~reproductive_status, scales="free_x") + labs(x="Sample names",y="Relative abundance")
plot_bar(ps.top20, x="Sample", fill="Order") + facet_wrap(~reproductive_status, scales="free_x") + labs(x="Sample names",y="Relative abundance")
plot_bar(ps.top20, x="Sample", fill="Family") + facet_wrap(~reproductive_status, scales="free_x") + labs(x="Sample names",y="Relative abundance")
plot_bar(ps.top20, x="Sample", fill="Genus") + facet_wrap(~reproductive_status, scales="free_x") + labs(x="Sample names",y="Relative abundance")
plot_bar(ps.top20, x="Sample", fill="Species") + facet_wrap(~reproductive_status, scales="free_x") + labs(x="Sample Names",y="Relative Abundance")

plot_bar(ps.top20, x="Sample", fill="Phylum") + facet_wrap(~harvest_status, scales="free_x") + labs(x="Sample names",y="Relative abundance")
plot_bar(ps.top20, x="Sample", fill="Class") + facet_wrap(~harvest_status, scales="free_x") + labs(x="Sample names",y="Relative abundance")
plot_bar(ps.top20, x="Sample", fill="Order") + facet_wrap(~harvest_status, scales="free_x") + labs(x="Sample names",y="Relative abundance")
plot_bar(ps.top20, x="Sample", fill="Family") + facet_wrap(~harvest_status, scales="free_x") + labs(x="Sample names",y="Relative abundance")
plot_bar(ps.top20, x="Sample", fill="Genus") + facet_wrap(~harvest_status, scales="free_x") + labs(x="Sample names",y="Relative abundance")
plot_bar(ps.top20, x="Sample", fill="Species") + facet_wrap(~harvest_status, scales="free_x") + labs(x="Sample names",y="Relative abundance")

#overall top species
plot_bar(ps.top20, x = "Sample", fill = "Species") +
  labs(title = "Top 20 Most Abundant Species (All Samples)",
       x = "Sample Names", y = "Relative Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal")

########################################### 
#RAREFACTION CURVES # not normalised
mycol <- c(
  "Pre"  = "#4CAF50",   # Mid Green
  "Post" = "#654321",
  "Female_Lactating" = "#E69F00",   # orange-ish, friendly to colorblind
  "Female_Pregnant" = "#0072B2",    # blue
  "Male" = "#D55E00"
)
rarecurve(t(asv_mat), step=100, col=mycol, lwd=2, ylab="ASVs", label=F, color=T, xlim=c(0,60000)) 
legend("topright", levels(mycol), lty=1, title="Reproductive Status", lwd=2) 
# legend doesn't work; need to add manually?

#############################################
### SACs    # not normalised

#install.packages("remotes")
#remotes::install_github("adrientaudiere/MiscMetabar")
install.packages("MiscMetabar")
library("MiscMetabar")

p_repro <- accu_plot(ps, "reproductive_status", add_nb_seq = TRUE, by.fact = TRUE, step = 10)
p_repro +
  scale_color_manual(values = c("Female_Lactating" = "cornflowerblue", "Female_Pregnant" = "brown1", "Male" = 'green'))

p_harvest <- accu_plot(ps, "harvest_status", add_nb_seq = TRUE, by.fact = TRUE, step = 10)
p_harvest +
  scale_color_manual(values = c("Pre"  = "#4CAF50",   # Mid Green
                                "Post" = "#000000"))  # Black
# see https://rdrr.io/github/adrientaudiere/MiscMetabar/man/accu_plot.html

########################################### 
#ALPHA DIVERSITY -- FULL DATASET, requires raw read counts/abundances for richness measures
########################################### 
#Visualise alpha-diversity:
# Looking for obvious systematic difference in alpha-diversity between e.g., collection methods
alpha_div <- estimate_richness(ps, measures = "Shannon")  #ignore the warning message about singletons (they're removed deliberately!)
alpha_div

alpha_div$samplename <- rownames(alpha_div)
geo_data <- data.frame(sample_data(ps))
geo_data$samplename = rownames(geo_data)
geo_data <- merge(geo_data, alpha_div, by = "samplename")

ggplot(geo_data, aes(x = reproductive_status, y = Shannon, color = reproductive_status)) +
  geom_point(size = 4, alpha = 0.5) +
  geom_smooth(method = "glm", se = FALSE, color = "lightgray", linetype = "dashed") +
  scale_color_manual(values = c(
    "Female_Lactating" = "#E69F00",   # orange-ish, friendly to colorblind
    "Female_Pregnant" = "#0072B2",    # blue
    "Male" = "#D55E00"                # reddish-brown, distinguishable from blue/orange
  )) +        # Blue
  scale_x_discrete(labels = c("Female_Lactating" = "Female (Lactating)",
                              "Female_Pregnant" = "Female (Pregnant)",
                              "Male" = "Male")) +
  theme_minimal() +
  labs(title = "",
       x = "Reproductive Status",
       y = "Alpha Diversity (Shannon)") +
  theme(legend.position = "none")


# Reorder the factor levels of harvest_status to ensure 'Pre' comes before 'Post'
# Reorder the factor levels
geo_data$harvest_status <- factor(geo_data$harvest_status, levels = c("Pre", "Post"))

# Plot
ggplot(geo_data, aes(x = harvest_status, y = Shannon, color = harvest_status)) +
  geom_point(size = 4, alpha = 0.5) +
  geom_smooth(method = "glm", se = FALSE, color = "lightgray", linetype = "dashed") +
  scale_color_manual(values = c("Pre"  = "#4CAF50",   # Mid Green
                                "Post" = "saddlebrown")) + # Black
  theme_minimal() +
  labs(title = "",
       x = "Harvest Status",
       y = "Alpha Diversity (Shannon)",
       color = "Harvest Status")  # Legend title


################################################################ 
#ORDINATION  -- NORMALISED DATA 
################################################################
ps_rare1 

#Generating and visualising the PCoA
vst_pcoa <- ordinate(ps_rare1, method="MDS", distance="bray")  # do we need the na.rm=T??

eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(ps_rare1, vst_pcoa, color = "reproductive_status") + 
  geom_point(size = 5) +
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) +
  theme_bw() + 
  theme(legend.position = "right") +
  scale_color_manual(
    values = c(
      "Female_Lactating" = "#E69F00",   # orange-ish
      "Female_Pregnant" = "#0072B2",    # blue
      "Male" = "#D55E00" , # reddish-brown
      "Extraction_Control" = 'red',
      "PCR_Control" = 'black'
    ),
    labels = c(
      "Female_Lactating" = "Female (Lactating)", 
      "Female_Pregnant" = "Female (Pregnant)", 
      "Male" = "Male"
    )
  ) +
  labs(color = "Reproductive Status") +  # legend title
  stat_ellipse()


# Ensure the harvest_status variable is a factor with the desired level order
geo_data$harvest_status <- factor(geo_data$harvest_status, levels = c("Post", "Pre"))
plot_ordination(ps_rare1, vst_pcoa, color="harvest_status") + 
  geom_point(size = 5) +
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) +
  theme_bw() + 
  theme(legend.position = "right") +
  scale_color_manual(
    values = c("Pre"  = "#4CAF50",   # Mid Green
               "Post" = "#7F7F7F", #grey
              "PCR_Control" = 'black', # Black
              "Extraction_Control" = 'red'),
    breaks = c("Pre", "Post", "PCR_Control", "Extraction_Control")        # Legend order: Pre first
  ) +
  labs(color = "Harvest Status") +
  stat_ellipse()

### non normalised looks very different to normalised
vst_pcoa1 <- ordinate(ps, method = "MDS", distance = "bray")  # MDS = PCoA with Bray-Curtis

# Extract eigenvalues for axis scaling
eigen_vals <- vst_pcoa1$values$Eigenvalues  # Used to fix aspect ratio of axes

# Access the sample data from the phyloseq object
sample_data_ps <- sample_data(ps_rare1)

# Combine reproductive status and harvest status into a new factor for color
sample_data_ps$combined_status <- sample_data_ps$reproductive_status

# Ensure combined_status is a factor
sample_data_ps$combined_status <- factor(sample_data_ps$combined_status, 
                                         levels = c("Female_Lactating", "Female_Pregnant", "Male", "PCR_Control", "Extraction_Control"))

# Ensure harvest_status is a factor
sample_data_ps$harvest_status <- factor(sample_data_ps$harvest_status, levels = c("Pre", "Post", "PCR_Control", "Extraction_Control"))

# Update the sample_data in the phyloseq object with the new combined_status
sample_data(ps_rare1) <- sample_data_ps

# Generate and visualize the PCoA
vst_pcoa <- ordinate(ps_rare1, method="MDS", distance="bray")

# Extract eigenvalues for scaling
eigen_vals <- vst_pcoa$values$Eigenvalues 

# Plot PCoA with reproductive status as color and harvest status as shape
plot_ordination(ps_rare1, vst_pcoa, color = "combined_status", shape = "harvest_status") + 
  geom_point(size = 5, stroke = 1.2) +
  coord_fixed(sqrt(eigen_vals[2] / eigen_vals[1])) +
  theme_bw() + 
  theme(legend.position = "right") +
  
  # Color mapping including controls
  scale_color_manual(
    values = c(
      "Female_Lactating" = "#E69F00",
      "Female_Pregnant" = "#0072B2",
      "Male" = "#D55E00",
      "PCR_Control" = "black", 
      "Extraction_Control" = "grey"
    ),
    na.translate = FALSE
  ) +

  # Shape mapping including controls
  scale_shape_manual(
    values = c(
      "Pre" = 16,              # solid circle
      "Post" = 17,             # solid triangle
      "PCR_Control" = 15,      # square
      "Extraction_Control" = 18 # diamond
    ),
    na.translate = FALSE
  ) +

  # Legend titles
  labs(
    color = "Reproductive Status", 
    shape = "Harvest Status"
  )


## Permanovas
# First, exclude control samples from your metadata
samples_df_nocontrol <- samples_df[!samples_df$reproductive_status %in% c("PCR_Control", "Extraction_Control"), ]

# Then subset the distance matrix to only those samples
dist_nocontrol <- dist.matrix[rownames(samples_df_nocontrol), rownames(samples_df_nocontrol)]

# Run betadisper only on non-control samples
bd <- betadisper(dist_nocontrol, samples_df_nocontrol$reproductive_status)

# Set proper factor levels to control box order (optional but good practice)
bd$group <- factor(bd$group, levels = c("Female_Lactating", "Female_Pregnant", "Male"))

# Color mapping
mycol_repro <- c(
  "Female_Lactating" = "#E69F00",   # orange
  "Female_Pregnant" = "#0072B2",    # blue
  "Male" = "#D55E00"                # reddish-brown
)

# Boxplot
boxplot(bd,
        col = mycol_repro,
        names = c("Female (Lactating)", "Female (Pregnant)", "Male"),
        xlab = "Reproductive Status",
        ylab = "Distance to Centroid")

# Permutation test for dispersion differences
permutest(bd, pairwise = TRUE)


# harvest status
dist.matrix2 <- phyloseq::distance(ps, method = 'bray')
# If no samples were lost during rarefaction:
ps.perm2<- adonis2(dist.matrix2 ~ samples$harvest_status, data = samples_df)
# OR, if some samples were lost:
samples_df2 <- samples_df[!rownames(samples_df) %in% c("UOW-24-01-028", "UOW-24-01-079"), ]
ps.perm2 <- adonis2(dist.matrix2 ~ samples_df2$harvest_status, data = samples_df2)
ps.perm2
 # Beta dispersion
bd2 <- betadisper(dist.matrix2, samples$harvest_status)
# OR
bd2 <- betadisper(dist.matrix2, samples_df2$harvest_status)
bd2
# Reorder factor levels for bd2 grouping
bd2$group <- factor(bd2$group, levels = c("Pre", "Post"))
# Custom colors for Pre and Post
mycol_harvest<- c("Pre" = "#4CAF50", "Post" = "saddlebrown")
# Plot dispersion
boxplot(bd2,
        col = mycol_harvest,
        xlab = "Harvest Status",
        ylab = "Distance to Centroid")
# Test significance of dispersion differences
permutest(bd2, pairwise = TRUE)

#### distances analysis 
library(phyloseq)
library(tidyverse)
library(vegan)
library(rstatix)  
library(patchwork)

# === ALPHA DIVERSITY vs DISTANCE (Spearman correlations + plots) ===

# Calculate Shannon diversity per sample
alpha_diversity <- estimate_richness(ps, measures = "Shannon") %>%
  rownames_to_column("SampleID")

# Extract metadata with distances
metadata_df <- data.frame(sample_data(ps)) %>%
  rownames_to_column("SampleID") %>%
  select(SampleID, Distance_to_road, Distance_to_water, Distance_to_waterpoint) %>%
  mutate(across(starts_with("Distance_to_"), ~ as.numeric(as.character(.))))

# Merge alpha diversity and metadata
alpha_meta <- left_join(alpha_diversity, metadata_df, by = "SampleID")

# Spearman correlations
cor_road <- cor.test(alpha_meta$Shannon, alpha_meta$Distance_to_road, method = "spearman", use = "complete.obs")
cor_water <- cor.test(alpha_meta$Shannon, alpha_meta$Distance_to_water, method = "spearman", use = "complete.obs")
cor_waterpoint <- cor.test(alpha_meta$Shannon, alpha_meta$Distance_to_waterpoint, method = "spearman", use = "complete.obs")

cat("Spearman correlations:\n")
print(cor_road)
print(cor_water)
print(cor_waterpoint)

# Step 1: Shannon diversity (already done)
alpha_diversity <- estimate_richness(ps, measures = "Shannon") %>%
  rownames_to_column("SampleID")

# Step 2: Extract metadata (including reproductive_status)
metadata_df <- data.frame(sample_data(ps)) %>%
  rownames_to_column("SampleID") %>%
  select(SampleID, reproductive_status, Distance_to_road, Distance_to_water, Distance_to_waterpoint)

# Step 3: Join metadata to alpha diversity
alpha_meta <- left_join(alpha_diversity, metadata_df, by = "SampleID")

library(ggplot2)
library(dplyr)
library(patchwork)  # for plot_layout()

# Extract metadata from phyloseq object as data frame
meta_df <- data.frame(sample_data(ps))

# Add SampleID as a column (rownames are sample IDs)
meta_df <- meta_df %>% mutate(SampleID = rownames(meta_df))

# Merge with alpha_meta by SampleID, keeping harvest_status and reproductive_status from meta_df
alpha_meta <- alpha_meta %>%
  left_join(meta_df %>% select(SampleID, harvest_status, reproductive_status), by = "SampleID")
# Use the meta_df columns (assuming those are correct)
alpha_meta <- alpha_meta %>%
  mutate(
    reproductive_status = reproductive_status.y,
    harvest_status = harvest_status.y
  )
alpha_meta <- alpha_meta %>%
  mutate(
    Distance_to_road = as.numeric(Distance_to_road),
    Distance_to_water = as.numeric(Distance_to_water),
    Distance_to_waterpoint = as.numeric(Distance_to_waterpoint)
  )
alpha_meta <- alpha_meta %>%
  mutate(
    Distance_to_road = as.numeric(Distance_to_road),
    Distance_to_water = as.numeric(Distance_to_water),
    Distance_to_waterpoint = as.numeric(Distance_to_waterpoint)
  )


# Define filled shapes for harvest status
harvest_shapes <- c(
  "Pre" = 21,   # filled circle
  "Post" = 24   # filled triangle up
)

# Define fill colors for reproductive status
repro_colors <- c(
  "Female_Lactating" = "#E69F00",   # orange-ish
  "Female_Pregnant" = "#0072B2",    # blue
  "Male" = "#D55E00"                 # reddish-brown
)

# p1: Distance to road vs Shannon Diversity
p1 <- ggplot(alpha_meta, aes(Distance_to_road, Shannon)) +
  geom_point(aes(shape = harvest_status, fill = reproductive_status), 
             alpha = 0.6, size = 3, color = "black") +  # black border
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  scale_shape_manual(values = harvest_shapes) +
  scale_fill_manual(values = repro_colors) +
  labs(x = "Distance to Road (m)", y = "Shannon Diversity", title = "Shannon vs Distance to Road") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# p2: Distance to water vs Shannon Diversity
p2 <- ggplot(alpha_meta, aes(Distance_to_water, Shannon)) +
  geom_point(aes(shape = harvest_status, fill = reproductive_status), 
             alpha = 0.6, size = 3, color = "black") +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  scale_shape_manual(values = harvest_shapes) +
  scale_fill_manual(values = repro_colors) +
  labs(x = "Distance to Water (m)", y = NULL, title = "Shannon vs Distance to Water") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

# p3: Distance to waterpoint vs Shannon Diversity
p3 <- ggplot(alpha_meta, aes(Distance_to_waterpoint, Shannon)) +
  geom_point(aes(shape = harvest_status, fill = reproductive_status), 
             alpha = 0.6, size = 3, color = "black") +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  scale_shape_manual(values = harvest_shapes) +
  scale_fill_manual(values = repro_colors) +
  labs(x = "Distance to Waterpoint (m)", y = NULL, title = "Shannon vs Distance to Waterpoint") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  )

# Combine the plots side-by-side
combined_plot <- (p1 + p2 + p3) + plot_layout(ncol = 3)

combined_plot

# === BETA DIVERSITY and PERMANOVA (with continuous distances) ===

# Transform to relative abundance and calculate Bray-Curtis distance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
otu_mat <- as.matrix(otu_table(ps_rel))
bray_dist <- vegdist(t(otu_mat), method = "bray")

# Prepare metadata and filter samples with complete distance data
meta_df <- data.frame(sample_data(ps_rel)) %>%
  mutate(across(c(Distance_to_road, Distance_to_water, Distance_to_waterpoint),
                ~ suppressWarnings(as.numeric(as.character(.)))))

meta_df_clean <- meta_df %>% filter(complete.cases(select(., Distance_to_road, Distance_to_water, Distance_to_waterpoint)))

# Subset distance matrix to samples with complete metadata
valid_samples <- rownames(meta_df_clean)
bray_dist_mat <- as.matrix(bray_dist)
bray_dist_clean <- as.dist(bray_dist_mat[valid_samples, valid_samples])

# PERMANOVA tests for each continuous distance variable
adonis_road <- adonis2(bray_dist_clean ~ Distance_to_road, data = meta_df_clean, permutations = 999)
adonis_water <- adonis2(bray_dist_clean ~ Distance_to_water, data = meta_df_clean, permutations = 999)
adonis_waterpoint <- adonis2(bray_dist_clean ~ Distance_to_waterpoint, data = meta_df_clean, permutations = 999)

cat("\nPERMANOVA Results:\n")
print(adonis_road)
print(adonis_water)
print(adonis_waterpoint)

library(phyloseq)
library(tidyverse)
library(vegan)
library(patchwork)

# ---- Beta dispersion test ----
# For beta dispersion, you typically need a grouping factor.
# Since your distances are continuous, let's categorize them into bins (e.g. tertiles) for testing.

meta_df_clean <- meta_df_clean %>%
  mutate(
    Road_cat = ntile(Distance_to_road, 3) %>% factor(labels = c("Low", "Medium", "High")),
    Water_cat = ntile(Distance_to_water, 3) %>% factor(labels = c("Low", "Medium", "High")),
    Waterpoint_cat = ntile(Distance_to_waterpoint, 3) %>% factor(labels = c("Low", "Medium", "High"))
  )


# Beta dispersion (homogeneity of multivariate dispersion) test for Distance_to_road categories
betadisp_road <- betadisper(bray_dist_clean, meta_df_clean$Road_cat)
permutest_road <- permutest(betadisp_road)

# Beta dispersion test for Distance_to_water categories
betadisp_water <- betadisper(bray_dist_clean, meta_df_clean$Water_cat)
permutest_water <- permutest(betadisp_water)

# Beta dispersion test for Distance_to_waterpoint categories
betadisp_waterpoint <- betadisper(bray_dist_clean, meta_df_clean$Waterpoint_cat)
permutest_waterpoint <- permutest(betadisp_waterpoint)

cat("Beta dispersion test results:\n")
print(permutest_road)
print(permutest_water)
print(permutest_waterpoint)
 
#############################################
#### CLAM TESTS -- normalised
#################################################
# First subset:
#removes female lactating from this 
no_lact <- subset_samples(ps, sample_data(ps)$reproductive_status != "Female_Lactating")
plot_bar(no_lact, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
no_lact_rare <- rarefy_even_depth(no_lact, sample.size=1000, rngseed=7, replace=FALSE, trimOTUs=TRUE, verbose=TRUE)
no_lact_clam <- clamtest(t(otu_table(no_lact_rare)), sample_data(no_lact_rare)$reproductive_status, alpha=0.01, specialization=2/3)
no_lact_general <- no_lact_clam[no_lact_clam$Classes == "Generalist",]
no_lact_rare <- no_lact_clam[no_lact_clam$Classes == "Too_rare",]
no_lact_specialist <- no_lact_clam[no_lact_clam$Classes == "Specialist_presence",]
no_lact_spec_abs <- no_lact_clam[no_lact_clam$Classes == "Specialist_absence",]
# Plot the no_lact_clam data and adjust labels
plot(no_lact_clam)
summary(no_lact_clam)
###
# Subset to remove Pregnant samples
no_preg <- subset_samples(ps, sample_data(ps)$reproductive_status != "Female_Pregnant")
plot_bar(no_preg, fill = "Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack")
no_preg_rare <- rarefy_even_depth(no_preg, sample.size = 1000, rngseed = 7, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
no_preg_clam <- clamtest(t(otu_table(no_preg_rare)), sample_data(no_preg_rare)$reproductive_status, alpha = 0.01, specialization = 2/3)
no_preg_general <- no_preg_clam[no_preg_clam$Classes == "Generalist", ]
no_preg_rare <- no_preg_clam[no_preg_clam$Classes == "Too_rare", ]
no_preg_specialist <- no_preg_clam[no_preg_clam$Classes == "Specialist_presence", ]
no_preg_spec_abs <- no_preg_clam[no_preg_clam$Classes == "Specialist_absence", ]
plot(no_preg_clam)
summary(no_preg_clam)
###
# Subset to remove Male samples
no_male <- subset_samples(ps, sample_data(ps)$reproductive_status != "Male")
plot_bar(no_male, fill = "Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack")
no_male_rare <- rarefy_even_depth(no_male, sample.size = 1000, rngseed = 7, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
no_male_clam <- clamtest(t(otu_table(no_male_rare)), sample_data(no_male_rare)$reproductive_status, alpha = 0.01, specialization = 2/3)
no_male_general <- no_male_clam[no_male_clam$Classes == "Generalist", ]
no_male_rare <- no_male_clam[no_male_clam$Classes == "Too_rare", ]
no_male_specialist <- no_male_clam[no_male_clam$Classes == "Specialist_presence", ]
no_male_spec_abs <- no_male_clam[no_male_clam$Classes == "Specialist_absence", ]
plot(no_male_clam)
summary(no_male_clam)
###
# Subset to just Pre and Post harvest_status
harvest_subset <- subset_samples(ps, harvest_status %in% c("Pre", "Post"))
table(sample_data(harvest_subset)$harvest_status)
harvest_rare <- rarefy_even_depth(harvest_subset, sample.size = 1000, rngseed = 7, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
harvest_clam <- clamtest(
  t(otu_table(harvest_rare)),
  sample_data(harvest_rare)$harvest_status,
  alpha = 0.01,
  specialization = 2/3)
harvest_general <- harvest_clam[harvest_clam$Classes == "Generalist",]
harvest_rare_taxa <- harvest_clam[harvest_clam$Classes == "Too_rare",]
harvest_specialist <- harvest_clam[harvest_clam$Classes == "Specialist_presence",]
harvest_spec_abs <- harvest_clam[harvest_clam$Classes == "Specialist_absence",]
plot(harvest_clam)
summary(harvest_clam)

# Specialist taxa = significantly different between the two habitats
# what are the specialist taxa??
# The plot method returns values invisibly and produces a bivariate scatterplot of species total abundances in the two habitats. 
# Symbols and boundary lines are shown for species groups.
# head(webs_clam) shows the ASV assignment to 'generalist, too rare, etc.
#########################
#no lact Male Specialist 
no_lact_specialist_male <- no_lact_clam[no_lact_clam$Classes == "Specialist_Male",]
ASV_no_lact_specialist_male <- no_lact_specialist_male$Species
no_lact_tax_table <- as.data.frame(tax_table(no_lact))
subset_tax_specialist_male <- no_lact_tax_table[rownames(no_lact_tax_table) %in% ASV_no_lact_specialist_male, ]
write.csv(subset_tax_specialist_male, file = "specialist_male_species_no_lact.csv")
head(subset_tax_specialist_male)
#no lact Preg Specialist 
no_lact_specialist_female_pregnant <- no_lact_clam[no_lact_clam$Classes == "Specialist_Female_Pregnant",]
ASV_no_lact_specialist_female_pregnant <- no_lact_specialist_female_pregnant$Species
no_lact_tax_table <- as.data.frame(tax_table(no_lact))
subset_tax_specialist_female_pregnant <- no_lact_tax_table[rownames(no_lact_tax_table) %in% ASV_no_lact_specialist_female_pregnant, ]
write.csv(subset_tax_specialist_female_pregnant, file = "specialist_female_pregnant_species_no_lact.csv")
head(subset_tax_specialist_female_pregnant)

#no preg Male Specialist
no_preg_specialist_male <- no_preg_clam[no_preg_clam$Classes == "Specialist_Male", ]
ASV_no_preg_specialist_male <- no_preg_specialist_male$Species
no_preg_tax_table <- as.data.frame(tax_table(no_preg))
subset_tax_specialist_male <- no_preg_tax_table[rownames(no_preg_tax_table) %in% ASV_no_preg_specialist_male, ]
write.csv(subset_tax_specialist_male, file = "specialist_male_species_no_preg.csv")
head(subset_tax_specialist_male)
#no preg Lact Specialist
no_preg_female_lactating <- no_preg_clam[no_preg_clam$Classes == "Specialist_Female_Lactating", ]
ASV_no_preg_female_lactating <- no_preg_female_lactating$Species
no_preg_tax_table <- as.data.frame(tax_table(no_preg))
subset_tax_female_lactating <- no_preg_tax_table[rownames(no_preg_tax_table) %in% ASV_no_preg_female_lactating, ]
write.csv(subset_tax_female_lactating, file = "specialist_female_lactating_species_no_preg.csv")
head(subset_tax_female_lactating)

#no male Preg Specialist 
no_male_female_pregnant <- no_male_clam[no_male_clam$Classes == "Specialist_Female_Pregnant", ]
ASV_no_male_female_pregnant <- no_male_female_pregnant$Species
no_male_tax_table <- as.data.frame(tax_table(no_male))
subset_tax_female_pregnant <- no_male_tax_table[rownames(no_male_tax_table) %in% ASV_no_male_female_pregnant, ]
write.csv(subset_tax_female_pregnant, file = "specialist_female_pregnant_species_no_male.csv")
head(subset_tax_female_pregnant)
#no male Lact Specialist
no_male_female_lactating <- no_male_clam[no_male_clam$Classes == "Specialist_Female_Lactating", ]
ASV_no_male_female_lactating <- no_male_female_lactating$Species
no_male_tax_table <- as.data.frame(tax_table(no_male))
subset_tax_female_lactating <- no_male_tax_table[rownames(no_male_tax_table) %in% ASV_no_male_female_lactating, ]
write.csv(subset_tax_female_lactating, file = "specialist_female_lactating_species_no_male.csv")
head(subset_tax_female_lactating)

#harvest specialist PRE
pre_harvest <- harvest_clam[harvest_clam$Classes == "Specialist_Pre", ]
ASV_pre_harvest <- pre_harvest$Species
pre_harvest_tax_table <- as.data.frame(tax_table(harvest_subset))
subset_tax_pre_harvest <- pre_harvest_tax_table[rownames(pre_harvest_tax_table) %in% ASV_pre_harvest, ]
write.csv(subset_tax_pre_harvest, file = "specialist_pre_harvest_species.csv")
head(subset_tax_pre_harvest)
#harvest specialist POST 
post_harvest <- harvest_clam[harvest_clam$Classes == "Specialist_Post", ]
ASV_post_harvest <- post_harvest$Species
post_harvest_tax_table <- as.data.frame(tax_table(harvest_subset))
subset_tax_post_harvest <- post_harvest_tax_table[rownames(post_harvest_tax_table) %in% ASV_post_harvest, ]
write.csv(subset_tax_post_harvest, file = "specialist_post_harvest_species.csv")
head(subset_tax_post_harvest)

library(ggplot2)
library(dplyr)

# Read the CSV file
data <- read.csv("Specialist Pregnant.csv", stringsAsFactors = FALSE)

# Count the number of species per Biostatus category
biostatus_df <- data %>%
  count(Biostatus) %>%
  mutate(Percentage = n / sum(n) * 100,
         Label = paste0(Biostatus, " (", round(Percentage, 1), "%)"))

# Create a pie chart
ggplot(biostatus_df, aes(x = "", y = n, fill = Label)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  labs(title = "Biostatus of Species", fill = "Biostatus")

# Function to generate pie chart for a given file
plot_biostatus_pie <- function(filename, title) {
  data <- read.csv(filename, stringsAsFactors = FALSE)
  
  biostatus_df <- data %>%
    count(Biostatus) %>%
    mutate(Percentage = n / sum(n) * 100,
           Label = paste0(Biostatus, " (", round(Percentage, 1), "%)"))
  
  ggplot(biostatus_df, aes(x = "", y = n, fill = Label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    theme_void() +
    labs(title = title, fill = "Biostatus")
}

# Generate pie charts
plot_biostatus_pie("Specialist Pregnant.csv", "Biostatus of Specialist Pregnant Species")
plot_biostatus_pie("Specialist Lactating.csv", "Biostatus of Specialist Lactating Species")
plot_biostatus_pie("Specialist Pre.csv", "Biostatus of Specialist Pre Species")
plot_biostatus_pie("Specialist Post.csv", "Biostatus of Specialist Post Species")
######
#pie charts for order for specialist species
# Load required packages
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# Function to plot pie chart for Order column with % in legend
plot_order_pie <- function(file, title) {
  # Read CSV
  data <- read.csv(file)
  
  # Check if 'Order' column exists
  if (!'Order' %in% names(data)) {
    stop(paste("'Order' column not found in", file))
  }
  
  # Calculate counts and percentages
  order_counts <- data %>%
    count(Order) %>%
    mutate(Percent = n / sum(n) * 100,
           Label = paste0(Order, " (", round(Percent, 1), "%)"))
  
  # Create pie chart
  ggplot(order_counts, aes(x = "", y = n, fill = Label)) +
    geom_col(width = 1) +
    coord_polar(theta = "y") +
    theme_void() +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_brewer(palette = "Set3", name = "Order")
}

# Generate pie charts
plot_order_pie("Specialist Pregnant.csv", "Order of Specialist Pregnant Species")
plot_order_pie("Specialist Lactating.csv", "Order of Specialist Lactating Species")
plot_order_pie("Specialist Pre.csv", "Order of Specialist Pre Species")
plot_order_pie("Specialist Post.csv", "Order of Specialist Post Species")

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tibble)
library(scales)  # for percent formatting

# Function to plot pie chart of proportions for a given phyloseq object and taxonomic rank
plot_tax_pie <- function(ps_obj, tax_rank, plot_title) {
  tax_df <- as.data.frame(tax_table(ps_obj))
  counts <- taxa_sums(ps_obj)
  
  df <- tax_df %>%
    rownames_to_column("ASV") %>%
    mutate(Count = counts[ASV]) %>%
    filter(!is.na(.data[[tax_rank]])) %>%
    group_by_at(tax_rank) %>%
    summarise(TotalCount = sum(Count)) %>%
    ungroup() %>%
    mutate(Proportion = TotalCount / sum(TotalCount),
           Label = paste0(.data[[tax_rank]], " (", percent(Proportion, accuracy = 0.1), ")"))
  
  df[[tax_rank]] <- factor(df[[tax_rank]], levels = df[[tax_rank]][order(df$Proportion, decreasing = TRUE)])
  df$Label <- factor(df$Label, levels = df$Label[order(df$Proportion, decreasing = TRUE)])
  
  ggplot(df, aes(x = "", y = Proportion, fill = Label)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    labs(title = plot_title, fill = tax_rank) +
    theme_void() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5))
}

# Helper function to subset phyloseq object by sample_data variable/value
subset_ps <- function(ps_obj, var, val) {
  samples <- sample_names(ps_obj)[sample_data(ps_obj)[[var]] == val]
  ps_sub <- prune_samples(samples, ps_obj)
  prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
}

# Create subsets
ps_pre <- subset_ps(ps, "harvest_status", "Pre")
ps_post <- subset_ps(ps, "harvest_status", "Post")
ps_female_lact <- subset_ps(ps, "reproductive_status", "Female_Lactating")
ps_female_preg <- subset_ps(ps, "reproductive_status", "Female_Pregnant")
ps_male <- subset_ps(ps, "reproductive_status", "Male")

# Plot all pies together (prints plots sequentially)
print(plot_tax_pie(ps_pre, "Class", "Pre - Class Proportions"))
print(plot_tax_pie(ps_pre, "Order", "Pre - Order Proportions"))

print(plot_tax_pie(ps_post, "Class", "Post - Class Proportions"))
print(plot_tax_pie(ps_post, "Order", "Post - Order Proportions"))

print(plot_tax_pie(ps_female_lact, "Class", "Female Lactating - Class Proportions"))
print(plot_tax_pie(ps_female_lact, "Order", "Female Lactating - Order Proportions"))

print(plot_tax_pie(ps_female_preg, "Class", "Female Pregnant - Class Proportions"))
print(plot_tax_pie(ps_female_preg, "Order", "Female Pregnant - Order Proportions"))

print(plot_tax_pie(ps_male, "Class", "Male - Class Proportions"))
print(plot_tax_pie(ps_male, "Order", "Male - Order Proportions"))



library(dplyr)
library(tibble)

# Function to compute order proportions for a given phyloseq object and label
compute_order_proportions <- function(ps_obj, group_label) {
  tax_df <- as.data.frame(tax_table(ps_obj))
  counts <- taxa_sums(ps_obj)
  
  df <- tax_df %>%
    rownames_to_column("ASV") %>%
    mutate(Count = counts[ASV]) %>%
    filter(!is.na(Order)) %>%
    group_by(Order) %>%
    summarise(TotalCount = sum(Count), .groups = "drop") %>%
    mutate(Proportion = TotalCount / sum(TotalCount),
           Group = group_label)
  
  return(df)
}

# Subset your data
subset_ps <- function(ps_obj, var, val) {
  samples <- sample_names(ps_obj)[sample_data(ps_obj)[[var]] == val]
  ps_sub <- prune_samples(samples, ps_obj)
  prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
}

# Create subset phyloseq objects
ps_pre <- subset_ps(ps, "harvest_status", "Pre")
ps_post <- subset_ps(ps, "harvest_status", "Post")
ps_female_lact <- subset_ps(ps, "reproductive_status", "Female_Lactating")
ps_female_preg <- subset_ps(ps, "reproductive_status", "Female_Pregnant")
ps_male <- subset_ps(ps, "reproductive_status", "Male")

# Generate summary tables for each
df_pre <- compute_order_proportions(ps_pre, "Pre")
df_post <- compute_order_proportions(ps_post, "Post")
df_female_lact <- compute_order_proportions(ps_female_lact, "Female_Lactating")
df_female_preg <- compute_order_proportions(ps_female_preg, "Female_Pregnant")
df_male <- compute_order_proportions(ps_male, "Male")

# Combine all into one table
order_summary_table <- bind_rows(df_pre, df_post, df_female_lact, df_female_preg, df_male)

# Format proportions as percentages for display
order_summary_table$Proportion <- round(order_summary_table$Proportion * 100, 2)

# Optional: view the table
print(order_summary_table)

# Optional: wide-format version (orders as rows, groups as columns)
order_summary_wide <- order_summary_table %>%
  select(Order, Group, Proportion) %>%
  tidyr::pivot_wider(names_from = Group, values_from = Proportion, values_fill = 0)

# Print all rows (assuming ~29 unique orders)
print(order_summary_wide, n = 29)
