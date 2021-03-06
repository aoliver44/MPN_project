###########################
##### ALPHA DIVERSITY #####
###########################

#path_to_data <- ""
#setwd(path_to_data)
source("../scripts/gen_basic_env.R")

library(vegan)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(nlme)

# Generate alpha diversity data
OTU_table <- otu_merged %>% select(., 1,32:NCOL(otu_merged)) %>% column_to_rownames(., var = "Sample")
richness <- as.data.frame(specnumber(OTU_table))
shannon <- diversity(OTU_table)
Evenness <- as.data.frame(shannon/log(richness))

alpha_tmp <- merge(metadata, richness, by.y = "row.names", by.x = "Sample")
alpha_data <- merge(Evenness, alpha_tmp, by.x = "row.names", by.y = "Sample")
names(alpha_data)[2] <- "Evenness"
names(alpha_data)[33] <- "Richness"

richness <- ggplot(data = alpha_data, aes(x = STATUS, y = Richness, fill = STATUS)) +
  geom_boxplot(outlier.size = 0, alpha = 0.3) + 
  geom_point(pch = 21, position = position_jitterdodge()) +
  labs(x = '',
       y = 'Richness') +
  annotate("text", x=1, y=270, label= "Anova: p > 0.05", size = 3) + 
  theme_bw() + scale_fill_manual(values=c("turquoise3", "tan3")) + theme(legend.position = "none") 

evenness <- ggplot(data = alpha_data, aes(x = STATUS, y = Evenness, fill = STATUS)) +
  geom_boxplot(outlier.size = 0, alpha = 0.3) + 
  geom_point(pch = 21, position = position_jitterdodge()) +
  labs(x = '',
       y = 'Evenness') +
  annotate("text", x=2, y=0.52, label= "Kruskal-Wallis: p > 0.05", size = 3) +
  theme_bw() + scale_fill_manual(values=c("turquoise3", "tan3")) + theme(legend.position = "none") 

plot_grid(richness, evenness, ncol = 1, labels = NULL)
# stats

# Stats (averaged across individual to get rid of repeated measures)
# Also because evenness is not normally distributed so gotta use a kruskal
# test, which doesnt allow the repeated measures "Error" term.

stats_data <- alpha_data %>% select(., Individual, Richness, Evenness) %>% group_by(., Individual) %>%
  summarise_at(., c("Richness", "Evenness"), mean)
stats_data <- merge(stats_data, metadata, by.x = "row.names", by.y = "Individual")
stats_data <- stats_data[!duplicated(stats_data$Individual),]

summary(aov(Richness ~ STATUS, data = stats_data))
kruskal.test(Evenness ~ STATUS, data = stats_data)

####################### SUBSTATUS ####################### 

Substatus_alpha <- subset(alpha_data, alpha_data$SUBSTATUS != "")
Substatus_rich <- ggplot(data = Substatus_alpha, aes(x = SUBSTATUS, y = Richness, fill = SUBSTATUS)) +
  geom_boxplot(outlier.size = 0, alpha = 0.3) + 
  geom_point(pch = 21, position = position_jitterdodge()) +
  labs(x = '',
       y = 'Richness') +
  annotate("text", x=1, y=270, label= "Anova: p > 0.05", size = 3) + 
  theme_bw() + scale_fill_manual(values=c("tan4", "tan1")) + theme(legend.position = "none") 

Substatus_even <- ggplot(data = Substatus_alpha, aes(x = SUBSTATUS, y = Evenness, fill = SUBSTATUS)) +
  geom_boxplot(outlier.size = 0, alpha = 0.3) + 
  geom_point(pch = 21, position = position_jitterdodge()) +
  labs(x = '',
       y = 'Evenness') +
  annotate("text", x=2, y=0.52, label= "Kruskal-Wallis: p > 0.05", size = 3) +
  theme_bw() + scale_fill_manual(values=c("tan4", "tan1")) + theme(legend.position = "none") 

Substatus_stats <- subset(stats_data, stats_data$SUBSTATUS != "")

summary(aov(Richness ~ SUBSTATUS, data = Substatus_stats))
kruskal.test(Evenness ~ SUBSTATUS, data = Substatus_stats)

plot_grid(Substatus_rich, Substatus_even, ncol = 1, labels = NULL)
