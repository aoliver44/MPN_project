###########################
####### Phylogenize #######
###########################

#path_to_data <- ""
#setwd(path_to_data)
source("../scripts/gen_basic_env.R")

library(ggplot2)
library(tidyverse)
library(ggrepel)

#### Phylogenize plots####
spec_healthy <- read.csv("enr-table.csv")
spec_healthy <- spec_healthy[grepl("strong", spec_healthy$Gene_significance),]
ggplot(data = spec_healthy) +
  aes(x = abs(log10(q_value)), y = log10(Enrichment_odds_ratio)) +
  theme_classic(base_line_size = 1) +
  geom_point(pch = 21, aes(fill = Phylum), colour="black", size = 4) +
  geom_text_repel(data = subset(spec_healthy, q_value < .05), aes(label = Subsystem), size = 3, vjust = -1.5) +
  geom_vline(xintercept = abs(log10(.05)), linetype = "dashed", color = "red") +
  labs(x = "q-value (absolute log10)", y = "Odds ratio (log10)", title = "")
