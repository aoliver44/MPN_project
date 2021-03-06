########################
###### CYTOKINES #######
########################

#path_to_data <- ""
#setwd(path_to_data)
source("../scripts/gen_basic_env.R")

library(ggplot2)
library(pheatmap)
library(tidyverse)
library(viridis)
library(ggpubr)

metadata <- read.csv("metadata_2.14.19.txt", check.names = F, sep = "\t")

cytokines <- read.csv("Cytokines.txt", sep = "\t")
# remove outliers (cytokine readings wayy too high)
cytokines <- cytokines %>% filter(., Sample != "124", Sample != "65") %>% replace(is.na(.), 0)
cytokine_rfp <- select(cytokines, Status, 3:14)
cytokines <- cytokines %>% filter(., Sample %in% metadata$Individual)
#cytokines[, 3:14] <- apply(cytokines[, 3:14], 1, scale)

color_row <- as.vector(cytokines$Status)
color_row <- str_replace(color_row, "HEALTHY", "blue")
color_row <- str_replace(color_row, "MPN", "orange")

# Generate heatmap of cytokine abundance
pheatmap(as.matrix(cytokines[, 3:14]),annotation_names_row = cytokines$Sample)
heatmap(as.matrix(cytokines[, 3:14]), scale = "column", RowSideColors = color_row, col=viridis(60), margins = c(3.5,18), labRow = F)
      

# Run Random Forest on cytokines
library(rfPermute)

RFP <- rfPermute(Status ~ ., 
                 data = cytokine_rfp, proximity = TRUE, 
                 importance = TRUE, parallel = 6, na.action = na.omit)

# plot plots
proximityPlot(RFP)
var_imp <- varImpPlot(RFP, type = 1)
plotConfMat(RFP)
tmp_plot <- impHeatmap(RFP, alpha = 0.05, n = 20, ranks = F)
tmp_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 0))
plot(RFP, alpha = 0.05)

# plot raw cytokine boxplots to support RF 
library(reshape2)
tmp <- as.data.frame(RFP$importance)
tmp <- tmp %>% arrange(., desc(MeanDecreaseAccuracy))
cytokine_rfp_melt <- melt(cytokine_rfp)
cytokine_rfp_melt$variable <- factor(cytokine_rfp_melt$variable, levels= rownames(tmp))

ggplot(data = cytokine_rfp_melt, aes(x = Status, y = log2(value))) +
  geom_boxplot(aes(fill = Status), outlier.colour = "white") + 
  geom_jitter(width = 0.2, aes(alpha = 0.3, shape = Status)) +
  facet_grid(. ~ variable) +
  scale_fill_manual(values=c("turquoise3", "tan3")) + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA), 
        #axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position="none",
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing = unit(0, "lines")) #+
  stat_compare_means(method = "kruskal.test")

