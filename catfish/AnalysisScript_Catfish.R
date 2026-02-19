### Statistical analysis script in R v 4.4.2

#########################################
### Figure 2 Detection rates
#########################################
# Read the data
detections <- read.csv("Detections_F.csv")

# Ensure factor levels for Assay and Result are in the desired order
detections$Assay <- factor(detections$Assay, levels = c("All", "WA", "RV", "WV"))
detections$Result <- factor(detections$Result, levels = c("False negative", "Positive"))


# Create the plot
detection.plot <- ggplot(detections, aes(x = Assay, y = Total, fill = Result)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  xlab("") + ylab("") +
  theme(text = element_text(size = 20))

detection.plot
ggsave("Figure2.jpeg")

#########################################
## Table 3
#########################################
# Read in the data for Rotoiti
compositionsR <- read.csv("Compositions_R.csv")
library("Hmisc")
compositions_corR <- rcorr(as.matrix(compositionsR[,2:13], type = "spearman"))
print(compositions_corR)
# Read in the data for Oranga
compositionsO <- read.csv("Compositions_O.csv")
compositions_corO <- rcorr(as.matrix(compositionsO[,2:9], type = "spearman"))
print(compositions_corO)

#########################################
## Figure 3a
#########################################
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

df <- compositionsR

# ---- Reshape ----
long_df <- df %>%
  pivot_longer(cols = -Location, names_to = "Variable", values_to = "Value") %>%
  mutate(scale_group = case_when(
    Variable == "Historic" ~ "Historic",
    Variable == "BPUE" ~ "BPUE",
    Variable == "CPUE" ~ "CPUE",
    TRUE ~ "Other"
  ))

# ---- Explicit order for Other facets ----
other_order <- c("WAFine", "WACoarse", "WADacron",
                 "RVFine", "RVCoarse", "RVDacron",
                 "WVFine", "WVCoarse", "WVDacron")

long_df$Variable <- factor(long_df$Variable,
                           levels = c("Historic", "CPUE", "BPUE", other_order))

# ---- Custom labels ----
facet_labels <- c(
  "Historic" = "Historic",
  "BPUE" = "BPUE",
  "CPUE" = "CPUE",
  "WAFine" = "Catfish_specific (WA) fine",
  "WACoarse" = "Catfish_specific (WA) coarse",
  "WADacron" = "Catfish_specific (WA) dacron",
  "RVFine" = "Multispecies (RV) fine",
  "RVCoarse" = "Multispecies (RV) coarse",
  "RVDacron" = "Multispecies (RV) dacron",
  "WVFine" = "Multispecies (WV) fine",
  "WVCoarse" = "Multispecies (WV) coarse",
  "WVDacron" = "Multispecies (WV) dacron"
)

# ---- Individual panels ----
p_hist <- ggplot(filter(long_df, scale_group == "Historic"), aes(Location, Value)) +
  geom_col(fill = "forest green") +
  facet_wrap(~Variable, ncol = 1, labeller = labeller(Variable = facet_labels)) +
  ylim(0, 50) +
  labs(y = "Abundance") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p_cpue <- ggplot(filter(long_df, scale_group == "CPUE"), aes(Location, Value)) +
  geom_col(fill = "forest green") +
  facet_wrap(~Variable, ncol = 1, labeller = labeller(Variable = facet_labels)) +
  ylim(0, 2.5) +
  labs(y = "No. of individuals/hr") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p_bpue <- ggplot(filter(long_df, scale_group == "BPUE"), aes(Location, Value)) +
  geom_col(fill = "forest green") +
  facet_wrap(~Variable, ncol = 1, labeller = labeller(Variable = facet_labels)) +
  ylim(0, 12) +
  labs(y = "g/hr") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p_other <- ggplot(filter(long_df, scale_group == "Other"), aes(Location, Value)) +
  geom_col(fill = "steelblue") +
  facet_wrap(~Variable, ncol = 3, labeller = labeller(Variable = facet_labels)) +
  ylim(0, 0.0125) +
  labs(y = "eDNA catfish read composition", x = "Location") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- Combine with patchwork ----
final_plot <- (p_hist | p_cpue | p_bpue) / p_other +
  plot_layout(heights = c(0.6, 1.4)) # top row smaller, bottom row larger

# ---- Display ----
final_plot

# ---- Optional: Save ----
ggsave("Figure3a.jpeg", final_plot, width = 14, height = 10, dpi = 600)

#########################################
## Figure 3b
#########################################
df <- compositionsO

# ---- Reshape ----
long_df <- df %>%
  pivot_longer(cols = -Location, names_to = "Variable", values_to = "Value") %>%
  mutate(scale_group = case_when(
  Variable == "BPUE" ~ "BPUE",
  Variable == "CPUE" ~ "CPUE",
  TRUE ~ "Other"
))

# ---- Explicit order for Other facets ----
other_order <- c("WAFine", "WACoarse", "RVFine", "RVCoarse", "WVFine", "WVCoarse")

long_df$Variable <- factor(long_df$Variable,
                           levels = c("CPUE", "BPUE", other_order))

# ---- Custom labels ----
facet_labels <- c(
  "BPUE" = "BPUE",
  "CPUE" = "CPUE",
  "WAFine" = "Catfish_specific (WA) fine",
  "WACoarse" = "Catfish_specific (WA) coarse",
  "RVFine" = "Multispecies (RV) fine",
  "RVCoarse" = "Multispecies (RV) coarse",
  "WVFine" = "Multispecies (WV) fine",
  "WVCoarse" = "Multispecies (WV) coarse"
)

# ---- CPUE panel styled like facet strip ----
p_cpue <- ggplot(filter(long_df, scale_group == "CPUE"), aes(Location, Value)) +
  geom_col(fill = "forest green") +
  facet_wrap(~Variable, ncol = 1, labeller = labeller(Variable = facet_labels)) +
  ylim(0, 10) +
  labs(y = "No. of individuals/hr") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p_bpue <- ggplot(filter(long_df, scale_group == "BPUE"), aes(Location, Value)) +
  geom_col(fill = "forest green") +
  facet_wrap(~Variable, ncol = 1, labeller = labeller(Variable = facet_labels)) +
  ylim(0, 40) +
  labs(y = "g/hr") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p_other <- ggplot(filter(long_df, scale_group == "Other"), aes(Location, Value)) +
  geom_col(fill = "steelblue") +
  facet_wrap(~Variable, ncol = 2, labeller = labeller(Variable = facet_labels)) +
  ylim(0, 0.030) +
  labs(y = "eDNA catfish read composition", x = "Location") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- Combine with patchwork ----
final_plot <- (p_cpue | p_bpue) / p_other +
  plot_layout(heights = c(0.6, 1.4)) # top row smaller, bottom row larger

# ---- Display ----
final_plot

# ---- Optional: Save ----
ggsave("Figure3b.jpeg", final_plot, width = 14, height = 10, dpi = 600)
