#######################
#### RANDOM FOREST ####
#######################

#path_to_data <- ""
#setwd(path_to_data)
source("../scripts/gen_basic_env.R")

set.seed(seed = 999)


library(ggplot2)
library(tidyverse)
library(janitor)
library(rfPermute)
library(nlme)

# run random forest on health statys
rfp_data_raw <- otu_merged %>% select(., 6, 32:NCOL(otu_merged)) %>% 
  select_if(negate(function(col) is.numeric(col) && sum(col) < 10)) %>% 
  clean_names() 

rfp_data_raw_mean <- rfp_data_raw %>%
  group_by(., individual, status) %>%
  summarise_all(., funs(mean)) %>% ungroup() %>%select(., -individual)

RFP <- rfPermute(as.factor(status) ~ ., 
                 data = rfp_data_raw, proximity = TRUE, 
                 importance = TRUE, parallel = 6, na.action = na.omit)

# plot plots
proximity <- proximityPlot(RFP)
var_imp <- varImpPlot(RFP, type = 1,n.var = 10)
plotConfMat(RFP)
tmp_plot <- impHeatmap(RFP, alpha = 0.05, n = 10, ranks = F)
tmp_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 0))
plot(RFP, alpha = 0.05)

# pull out top 30 taxa from RF
varimp.df <- as.data.frame(var_imp)
varimp.df$taxa <- NA
varimp.df$taxa <- rownames(varimp.df) 
top_20 <- varimp.df %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  slice(1:30)

# average out and check normality again
otu_average <- otu_merged %>% select(., Individual, STATUS, SUBSTATUS, 32:NCOL(otu_merged)) %>% 
  group_by(., Individual, STATUS, SUBSTATUS) %>% summarise_at(., vars(starts_with("k_")), mean)

# checking the raw data
ggplot(data = otu_merged, aes(x = STATUS, y = `k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Veillonellaceae;g__Phascolarctobacterium;s__Phascolarctobacterium`, fill = STATUS)) +
  geom_boxplot(outlier.size = 0, alpha = 0.3) + 
  geom_point(pch = 21, position = position_jitterdodge()) +
  labs(x = '',
       y = 'Phascolarctobacterium\n(OTU 64) reads') +
  annotate("text", x=2, y=60, label= "Kruskal-Wallis: p < 0.05", size = 3) +
  theme_bw() + scale_fill_manual(values=c("turquoise3", "tan3")) + theme(legend.position = "none") 
##################################################
# medication:
# k_bacteria_p_firmicutes_c_clostridia_o_clostridiales_f_lachnospiraceae_g_ruminococcus_s_gnavus_1
# k_bacteria_p_firmicutes_c_clostridia_o_clostridiales_f_clostridiales_g_clostridiales_s_clostridiales_26
# k_bacteria_p_firmicutes_c_clostridia_o_clostridiales_f_ruminococcaceae_g_ruminococcaceae_s_ruminococcaceae_4
# k_bacteria_p_firmicutes_c_clostridia_o_clostridiales_f_lachnospiraceae_g_coprococcus_s_coprococcus_8
# 
# ggplot(data = rfp_data_raw, aes(x = medication, y = `k_bacteria_p_firmicutes_c_clostridia_o_clostridiales_f_lachnospiraceae_g_coprococcus_s_coprococcus_8`, fill = medication)) +
#   geom_boxplot(outlier.size = 0, alpha = 0.3) + 
#   geom_point(pch = 21, position = position_jitterdodge()) +
#   labs(x = '',
#        y = 'ruminococcus_s_gnavus') +
#   #annotate("text", x=2, y=60, label= "Kruskal-Wallis: p < 0.05", size = 3) +
#   stat_compare_means() +
#   theme_bw() + scale_fill_manual(values=c("turquoise3", "tan3", "red")) + theme(legend.position = "none") 
###################################################

# linear mixed effects model with individual as repeated measure, data not normally distributed
tmp <- otu_merged
tmp$bacteria <- tmp$`k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Veillonellaceae;g__Phascolarctobacterium;s__Phascolarctobacterium`
anova(lme(bacteria ~ STATUS, random = ~ 1 | as.factor(Individual), data = tmp))

# Not normally distributed for stats reasons
shapiro.test(otu_merged$`k__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Veillonellaceae; g__Phascolarctobacterium; s__`)

# Still not normally distributed, but got rid of repeated measures problem
shapiro.test(otu_average$`k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Veillonellaceae;g__Phascolarctobacterium;s__Phascolarctobacterium`)

# KW test for non-parametric anova
kruskal.test(otu_average$`k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Veillonellaceae;g__Phascolarctobacterium;s__Phascolarctobacterium` ~ STATUS, data = otu_average)

# run random forest on health substatus
rfp_data_raw <- otu_merged %>% select(., 7, 32:NCOL(otu_merged)) %>% filter(., SUBSTATUS == "MF" | SUBSTATUS == "PV-ET") %>%
  clean_names()
names(rfp_data_raw)[2:NCOL(rfp_data_raw)] <- sapply(strsplit(names(rfp_data_raw)[2:NCOL(rfp_data_raw)], '_g_'), getElement, 2);
names(rfp_data_raw)[2:NCOL(rfp_data_raw)] <- paste0( "g_", names(rfp_data_raw)[2:NCOL(rfp_data_raw)])
names(rfp_data_raw) <- make.unique(colnames(rfp_data_raw))
rfp_data_raw$substatus <- droplevels(rfp_data_raw$substatus)

RFP <- rfPermute(as.factor(substatus) ~ ., 
                 data = rfp_data_raw, proximity = TRUE, 
                 importance = TRUE, parallel = 6, na.action = na.omit)

# plot plots
reprtree:::plot.getTree(RFP)
proximityPlot(RFP)
var_imp <- varImpPlot(RFP, type = 1,n.var = 30)
plotConfMat(RFP)
tmp_plot <- impHeatmap(RFP, alpha = 0.05, n = 20, ranks = F)
tmp_plot + theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 0))
plot(RFP, alpha = 0.05)
