##################################
######### Beta Diversity #########
##################################

#path_to_data <- ""
#setwd(path_to_data)
source("../scripts/gen_basic_env.R")

library(vegan)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(nlme)
library(janitor)

set.seed(999)

# Sequence Data
OTU_table <- otu_merged %>% filter(., SUBSTATUS != "")
OTU_table <- OTU_table %>% select(., 1,32:NCOL(OTU_table)) %>% column_to_rownames(., var = "Sample") %>%
  clean_names()

# Beta-diversity analysis, NMDS and Permanovas
beta.mds <- metaMDS(OTU_table, distance="bray", k=2)

stressplot(beta.mds)

sites <- as.data.frame(scores(beta.mds, display = "sites"))
species <- as.data.frame(scores(beta.mds, display = "species"))

nmds.sites <- merge(sites, metadata, by.x = "row.names", by.y = "Sample")
nmds.sites <- nmds.sites[order(nmds.sites$Timepoint),]
nmds.sites$Individual <- as.factor(nmds.sites$Individual)
cb_7 <- c("#d22154",
          "#007f36",
          "#670066",
          "#fdce5d",
          "#ccaaff",
          "#a65800",
          "#b25566")

ggplot(data = nmds.sites, aes(NMDS1, NMDS2, color = SUBSTATUS)) + 
  geom_point(aes(shape = SUBSTATUS), alpha = 0.7, size = 2) +
  stat_ellipse()


  
# statistics: permanova + pairwise
# Thank you Pedro Martinez Arbizu
source("parwise.adonis.r")
permanova_data <- merge(metadata, midas, by.x = "X.NAME", by.y = "row.names")

permanova_ind <- adonis(permanova_data[,11:NCOL(permanova_data)] ~ as.factor(Treatment), data = permanova_data, permutations = 999, parallel = 4, method = "bray")
permanova_ind

pairwise.adonis(x = permanova_data[,11:NCOL(permanova_data)], factors = permanova_data$Treatment, p.adjust.m = "BH")


temp_d_r <- subset(permanova_data, permanova_data$SAMPLETYPE == "Donors" | permanova_data$Treatment == "Post_FMT")
pairwise.adonis(x = temp_d_r[,11:NCOL(temp_d_r)], factors = temp_d_r$Treatment, p.adjust.m = "BH")

# comparisons of bray distances
tmp <- merge(metadata, midas, by.x = "X.NAME", by.y = "row.names") 
tmp <- subset(tmp, tmp$Individual_level != "3_D")
tmp <- tmp %>% select(., 1, 11:NCOL(permanova_data))
rownames(tmp) <- tmp$X.NAME
tmp$X.NAME <- NULL

tmp_bray_dist <- vegdist(tmp, method = "bray")

# match up the order of the metadata based on the distance matrix
metadata_delta <- metadata[match(rownames(as.data.frame(as.matrix(tmp_bray_dist))), metadata$X.NAME), ]

# take mean distances based on individual
ind_bcs <- meandist(tmp_bray_dist, grouping = metadata_delta$Specific_donor)
# get the diagonal of the matrix...basically the avg within an individual
bc_distances_diag <- as.data.frame(diag(as.matrix(ind_bcs)))

tmp2 <- as.data.frame(as.matrix(tmp_diag))
diag(tmp_diag) <- NA

bc_distances_off_diag 

