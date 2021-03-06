
######################
### TAXA BAR PLOTS ###
######################

#path_to_data <- ""
#setwd(path_to_data)
source("../scripts/gen_basic_env.R")

library(tidyverse)
library(reshape2)
library(cowplot)

midas <- t(otu_with_taxa)
midas <- rownames_to_column(as.data.frame(midas))

# Split taxonomy into different sections L1-L7
midas_melt <- melt(midas, id.vars = "rowname") %>% separate(., col = rowname, into = c("L1","L2","L3","L4","L5","L6","L7"), sep = "; ", remove = T, extra = "drop")

# Make sure the rel abundance or counts are numeric
midas_melt$value <- as.numeric(midas_melt$value)

# Summarize by L5 genus (or any other taxa group)...
# if you change make sure you change L5 and prefix (^g_)
midas_summarize <- midas_melt %>% group_by(., L5) %>% filter(str_detect(L5, "f_")) %>% summarise(., top_bacteria = sum(value)) %>% arrange(., desc(top_bacteria)) %>% slice(., 1:11)

# group the main players together into a list
high_abundance <- split(midas_summarize$L5, 1:NROW(midas_summarize))

# change everything that is not a main player into a other catagory
midas_melt$L5[midas_melt$L5 %in% high_abundance != "TRUE"] <- "other"


stephen_12 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')

sarah_color <- c("#7F0A57", "#A64685", "#CD9ABB", "#0B447A", "#3F77AC", "#4176AA", "#74A9DD", "#007976", "#39A9AB", "#71CFC5", "#72D3C6", "#007947", "#3BAA78")

julio_color <- c("#003f5c", "#2f4b7c", "#665191", "#a05195", "#d45087", "#f95d6a", "#ff7c43", "#ffa600", "#7f0a57", "#cd9abb", "#39a9ab", "#71cfc5", "#007947", "#bebebe")
# Plot

plot_data <- merge(metadata, midas_melt, by.x = "Sample", by.y = "variable")
plot_data <- plot_data %>% group_by(., Individual, STATUS, L5) %>% summarise(., sum_abund = sum(value))

MPN <- ggplot(data = subset(plot_data, plot_data$STATUS == "MPN"), 
       aes(x = Individual, weight = sum_abund, fill = L5)) +
  geom_bar(position = position_fill()) +
  theme_bw(base_size = 16) + 
  facet_grid(. ~  Individual, scales = "free") + 
  scale_fill_manual(values = julio_color) +
  theme(panel.spacing = unit(0.1, "lines")) +   
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
  labs(x = '',
       y = '') +
  ggtitle("MPN") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(plot.title = element_text(hjust = 0.5))

HEALTHY <- ggplot(data = subset(plot_data, plot_data$STATUS == "HEALTHY"), 
                  aes(x = Individual, weight = sum_abund, fill = L5)) +
  geom_bar(position = position_fill()) +
  theme_bw(base_size = 16) + 
  facet_grid(. ~  Individual, scales = "free") + 
  scale_fill_manual(values = julio_color) +
  theme(panel.spacing = unit(0.1, "lines")) +
  theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
  labs(x = '',
       y = '') +
  ggtitle("Healthy") +
  scale_y_continuous(labels = scales::percent_format()) + theme(legend.position = "blank") +
  theme(plot.title = element_text(hjust = 0.5))

cowplot::plot_grid(MPN, HEALTHY, nrow = 2, labels = NULL)
