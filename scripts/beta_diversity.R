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

#### BRAY DIST COMPARISONS ####
OTU_table <- otu_merged %>% select(., 1,32:NCOL(OTU_table)) %>% column_to_rownames(., var = "Sample") %>%
  clean_names()
# type
bray_dist <- vegdist(OTU_table, method = "bray")
bray_dist_merge <- merge(metadata, as.data.frame(as.matrix(bray_dist)), by.x = "Sample", by.y = "row.names")

# match up the order of the metadata based on the distance matrix
metadata_ordered <- metadata[match(rownames(as.data.frame(as.matrix(bray_dist))), bray_dist_merge$Sample), ]

# take mean distances based on individual
ind_bcs <- meandist(bray_dist, grouping = metadata_ordered$Individual)
# get the diagonal of the matrix...basically the avg within an factor
within_bcs <- as.data.frame(diag(as.matrix(ind_bcs)))

# calculate the between BC distances
tmp_diag <- ind_bcs
diag(tmp_diag) <- NA
tmp_diag <- as.data.frame(as.matrix(tmp_diag))
tmp_melt <- reshape2::melt(tmp_diag)

between_bcs <- tmp_melt %>% group_by(., variable) %>% drop_na(.) %>% summarise(., mean = mean(value))

# Group within and between together and merge with metadata

bc_distances <- merge(within_bcs, between_bcs, by.x = "row.names", by.y = "variable")
bc_distances <- bc_distances %>% drop_na(.)
colnames(bc_distances) <- c("factor", "within", "between")
bc_distances_melted <- reshape2::melt(bc_distances, id.vars = "factor")
bc_distances_melted <- merge(bc_distances_melted, metadata, by.x = "factor", by.y = "Individual", all.x = T)
bc_distances_melted <- filter(bc_distances_melted, variable == "within")
bc_distances_melted <- bc_distances_melted[!duplicated(bc_distances_melted$factor), ]

ggplot(data = subset(bc_distances_melted, bc_distances_melted$variable == "within"), aes(x = STATUS, y = value, fill = STATUS)) +
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(width = 0.06, alpha = 0.6) +
  geom_line(aes(group = HOUSE), alpha = 0.2, colour = "black") +
  scale_fill_manual(values = c("#003F5C", "#D45087")) +
  labs(x = "", y = "Average Bray-Curtis\nDissimilarity") + ggpubr::stat_compare_means()


########## COMPARE BRAY DISTANCE ######################

library(vegan)
library(reshape2)

# create distance matrix with Bray-Curtis distances
set.seed(999)
OTU_table <- otu_merged %>% select(., 1,32:NCOL(OTU_table)) %>% column_to_rownames(., var = "Sample") %>%
  clean_names()
# type
bray_dist <- vegdist(OTU_table, method = "bray")

# melt into long form-- thanks OTUsummary package
matrixConvert <- function(triMatrix, colname = c("sp1","sp2","dist")){
  m <- as.matrix(triMatrix)
  m2 <- reshape2::melt(m)[reshape2::melt(upper.tri(m))$value,]
  names(m2) <- colname
  invisible(m2)
}
type_distance_melted <- matrixConvert(bray_dist)

# make a column for the comparisons
bray_dist_metadata <- metadata %>% select(., Sample, HOUSE, Individual, STATUS)
type_distance_melted_merged <- merge(type_distance_melted, bray_dist_metadata, by.x = "sp1", by.y = "Sample")

type_distance_melted_merged <- type_distance_melted_merged %>% mutate(., cohabitation = ifelse(as.factor(HOUSE) == "", "between", "within"))

ggplot(type_distance_melted_merged) + aes(x = type, y = mean_dist, fill = type) +
  geom_boxplot(outlier.alpha = 0) + geom_point(position = position_jitter(width = 0.1)) +
  scale_fill_manual(values = c("#A6611A", "#018571"))

wilcox.test(type_distance_melted_merged$mean_dist ~ type_distance_melted_merged$type)



