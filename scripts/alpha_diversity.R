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
library(rstatix)

# Generate alpha diversity data
OTU_table <- otu_merged %>% select(., 1,32:NCOL(otu_merged)) %>% column_to_rownames(., var = "Sample")
richness <- as.data.frame(specnumber(OTU_table))
shannon <- diversity(OTU_table)
Evenness <- as.data.frame(shannon/log(richness))

alpha_tmp <- merge(metadata, richness, by.y = "row.names", by.x = "Sample")
alpha_data <- merge(Evenness, alpha_tmp, by.x = "row.names", by.y = "Sample")
names(alpha_data)[2] <- "Evenness"
names(alpha_data)[33] <- "Richness"

alpha_data_mean <- alpha_data %>% select(., STATUS, Individual, Richness, Evenness) %>%
  group_by(., STATUS, Individual) %>%
  summarise_all(., funs(mean))

richness <- ggplot(data = alpha_data_mean, aes(x = STATUS, y = Richness, fill = STATUS)) +
  geom_boxplot(outlier.size = 0, alpha = 0.3) + 
  geom_point(pch = 21, position = position_jitterdodge()) +
  labs(x = '',
       y = 'Richness') +
  annotate("text", x=1, y=270, label= "Anova: p > 0.05", size = 3) + 
  theme_bw() + scale_fill_manual(values=c("turquoise3", "tan3")) + theme(legend.position = "none") 

evenness <- ggplot(data = alpha_data_mean, aes(x = STATUS, y = Evenness, fill = STATUS)) +
  geom_boxplot(outlier.size = 0, alpha = 0.3) + 
  geom_point(pch = 21, position = position_jitterdodge()) +
  labs(x = '',
       y = 'Evenness') +
  annotate("text", x=2, y=0.52, label= "Kruskal-Wallis: p > 0.05", size = 3) +
  theme_bw() + scale_fill_manual(values=c("turquoise3", "tan3")) + theme(legend.position = "none") 

plot_grid(richness, evenness, ncol = 1, labels = NULL)

shapiro.test(alpha_data_mean$Richness)
levene_test(formula = alpha_data_mean$Richness ~ alpha_data_mean$STATUS, data = alpha_data_mean)
summary(aov(Richness ~ STATUS, data = alpha_data_mean))
kruskal.test(Evenness ~ STATUS, data = alpha_data_mean)


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

