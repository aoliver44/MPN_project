###########################
####### Correlations ######
###########################


#path_to_data <- ""
#setwd(path_to_data)
source("../scripts/gen_basic_env.R")

library(ggplot2)
library(corrr)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(rstatix)
library(vegan)
library(zoo)

cytokines <- read.csv("Cytokines.txt", sep = "\t")

midas <- t(otu_with_taxa)
midas <- rownames_to_column(as.data.frame(midas))

# Split taxonomy into different sections L1-L7
midas_melt <- melt(midas, id.vars = "rowname") %>% separate(., col = rowname, into = c("L1","L2","L3","L4","L5","L6","L7"), sep = ";", remove = T, extra = "drop")

plot_data <- merge(metadata, midas_melt, by.x = "Sample", by.y = "variable")

# Make sure the rel abundance or counts are numeric
plot_data$value <- as.numeric(plot_data$value)
midas_summarize <- plot_data %>% 
  group_by(., Individual) %>%
  mutate(., count = n_distinct(Sample)) %>%
  group_by(., Individual, L6, Sample, count) %>% 
  filter(str_detect(L6, "g_")) %>% 
  group_by(., L6, Individual) %>% summarise(., mean_count = mean(value))

midas_summarize <- dcast(Individual ~ L6, data = midas_summarize)

MPN <- metadata %>% filter(., STATUS == "MPN") %>% select(., Individual)
mpn_ind <- as.vector(unique(MPN$Individual))

midas_summarize <- midas_summarize %>% filter(., Individual %in% mpn_ind)

corr_data <- merge(midas_summarize, cytokines, by.x = "Individual", by.y = "Sample")
corr_data <- corr_data %>% select(., -Status) %>% column_to_rownames(var = "Individual")
cytokine_names <- names(cytokines)[3:14]

corr_data_tnf <- corr_data %>% select(1:145,157)
cor_test <- corr_data_tnf %>%
  cor_test(method = "spearman", vars = "TNFa")
cor_test <- cor_test %>% filter(., cor != "")
tnf_p_values1 <- as.data.frame(p.adjust(cor_test$p, method = "fdr"))
tnf_p_values2 <- cbind(tnf_p_values1, cor_test) %>% arrange(., cor)
colnames(tnf_p_values2)[1] <- "padjust"
  
ggplot(data = cor_test) + aes(var1, var2, fill = cor) + geom_tile() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(data = corr_data) + aes(x = g__Parabacteroides, y = TNFa) + geom_point()

# RDA
community <- corr_data %>% select(., 1:95)
env <- corr_data %>% select(., 96:107)

rda_attempt <- cca(x = community, y = env, )

# dblm
cytokines_mpn <- cytokines %>% filter(., Status == "MPN") 
mpn_w_cytokine <- as.vector(unique(cytokines_mpn$Sample))
cytokines_mpn <- cytokines_mpn %>% column_to_rownames(., var = "Sample") %>% select(., 2:NCOL(.)) %>% replace(is.na(.), 0)
mpn_otu <- subset(otu_merged, otu_merged$Individual %in% mpn_w_cytokine)
rownames(mpn_otu) <- NULL
mpn_otu <- mpn_otu %>% 
  select(., Individual, 32:NCOL(.)) %>% 
  group_by(., Individual) %>%
  summarise_all(funs(mean)) %>%
  column_to_rownames(., var = "Individual") 
cytokines_mpn <- cytokines_mpn[rownames(cytokines_mpn) %in% rownames(mpn_otu), ]
mpn_dist <- vegdist(mpn_otu[,32:ncol(mpn_otu)])

# since axis length 1 is greater than 4,we should use cca over rda
decorana(veg = cytokines_mpn)
# check what distance metric you should use...?
rankindex(grad = mpn_otu, veg = cytokines_mpn, indices = c("euc", "man", "gow","bra", "kul"), method = "spearman", stepacross = FALSE)

# determining which variables to use in the model
ordistep(cca(mpn_otu ~ TNFa + IFNg + IFNa2 + IL8 + IP10 + GRO + IL1b + IL1a + IL17a + IL1RA, data = cytokines_mpn), direction="forward", pstep=1000, R2scop=TRUE)
# make sure they dont have inflated variance. I think under 10 is good.
vif.cca(tmp)

# run the model
tmp <- cca(mpn_otu ~ TNFa + IFNg + IFNa2 + IL8 + IP10 + GRO + IL1b + IL1a + IL17a + IL1RA, data = cytokines_mpn, distance = "bray", add = T)
# check the model
anova.cca(tmp2, by="terms", permu=200)
# plot the model
ggord(tmp2)
plot(tmp2)
text(tmp2, 'species')


