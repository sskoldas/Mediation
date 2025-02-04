setwd("~/Desktop/Mediation/")
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
counts <- read.delim("./Dataset/otutab_raw.txt", row.names = 1)
taxa <- read.delim("./Dataset/OTU_tax.txt", row.names=1)
taxa[taxa == ""] <- NA
taxa$Species <- ifelse(grepl("uncultured", taxa$Species, ignore.case = TRUE), NA, taxa$Species)
names(taxa)[1] <- "Kingdom"
taxa$Genus<- gsub("-", "_", taxa$Genus)

#handling meta data
metadata <- read.delim("./Dataset/meta_filtered.tsv", row.names=1)
metadata$host_fev1[metadata$host_fev1 == "missing"] <- NA
metadata$host_fev1 <- as.numeric(metadata$host_fev1)
metadata <- metadata[complete.cases(metadata[, "host_clinical_state"]), ] 

table(metadata$host_clinical_state)
metadata$host_clinical_state <- factor(metadata$host_clinical_state, levels = c("Baseline", "Exacerbation", "Treatment","Recovery"))
meta <- metadata

# functions-----------------------
process_fev1_abund <- function(counts, meta, taxa, 
                               taxrank = NULL,
                               lib.size.cut.off = 1000){
  
  ps <- phyloseq(otu_table(counts,taxa_are_rows = TRUE),
                 sample_data(meta),
                 phyloseq::tax_table(as.matrix(taxa)))
  ps <- subset_taxa(ps,!is.na(Kingdom) & Kingdom != "" & Kingdom != "Archaea" & Kingdom != "Eukaryota")
  ps.agg <- tax_glom(ps, taxrank = taxrank)
  ps.agg <- prune_samples(sample_sums(ps.agg) >= lib.size.cut.off, ps.agg)
  agg.meta <- data.frame(sample_data(ps.agg))
  agg.taxa <- data.frame(tax_table(ps.agg))
  agg.otu <- data.frame(t(otu_table(ps.agg)))
  agg.otu <- agg.otu[, colSums(agg.otu != 0) > 0]
  merge <- merge(agg.meta, agg.otu, by = 0)
  return(merge)
}

process_for_top5_taxa <- function(counts, meta, taxa, 
                                  taxrank = NULL,
                                  lib.size.cut.off = 1000){
  
  ps <- phyloseq(otu_table(counts,taxa_are_rows = TRUE),
                 sample_data(meta),
                 phyloseq::tax_table(as.matrix(taxa)))
  ps <- subset_taxa(ps,!is.na(Kingdom) & Kingdom != "" & Kingdom != "Archaea" & Kingdom != "Eukaryota")
  ps.agg <- tax_glom(ps, taxrank = taxrank)
  ps.agg <- prune_samples(sample_sums(ps.agg) >= lib.size.cut.off, ps.agg)
  agg.meta <- data.frame(sample_data(ps.agg))
  agg.taxa <- data.frame(tax_table(ps.agg))
  agg.otu <- data.frame(t(otu_table(ps.agg)))
  agg.otu <- agg.otu[, colSums(agg.otu != 0) > 0]
  ps <- phyloseq(otu_table(t(agg.otu), taxa_are_rows = TRUE),
                 sample_data(agg.meta),
                 phyloseq::tax_table(as.matrix(agg.taxa)))
  return(ps)
}


process_for_beta <- function(counts, meta, taxa, 
                             taxrank = NULL,
                             lib.size.cut.off = 1000){
  
  ps <- phyloseq(otu_table(counts,taxa_are_rows = TRUE),
                 sample_data(meta),
                 phyloseq::tax_table(as.matrix(taxa)))
  ps <- subset_taxa(ps,!is.na(Kingdom) & Kingdom != "" & Kingdom != "Archaea" & Kingdom != "Eukaryota")
  ps.agg <- tax_glom(ps, taxrank = taxrank)
  ps.agg <- prune_samples(sample_sums(ps.agg) >= lib.size.cut.off, ps.agg)
  
  agg.meta <- data.frame(sample_data(ps.agg))
  agg.taxa <- data.frame(tax_table(ps.agg))
  agg.otu <- data.frame(t(otu_table(ps.agg)))
  agg.otu <- agg.otu[, colSums(agg.otu != 0) > 0]
  otu.total <- vegan::decostand(agg.otu, MARGIN = 1, method = "total")
  ps <- phyloseq(otu_table(t(otu.total), taxa_are_rows = TRUE),
                 sample_data(agg.meta),
                 phyloseq::tax_table(as.matrix(agg.taxa)))
  return(ps)
}

##################################################################
# Clinical States
##################################################################
clinical_states_comparisons <- list(
  c("Baseline", "Exacerbation"),
  c("Baseline", "Treatment"),
  c("Baseline", "Recovery"),
  c("Exacerbation", "Treatment"),
  c("Exacerbation", "Recovery"),
  c("Treatment", "Recovery")
)

# FEV1 Distribution Between-------------------
data_fev1_abund <- process_fev1_abund(counts, meta, taxa, 
                                            taxrank = "Genus", 
                                            lib.size.cut.off = 1000)
# plot
ggplot(data_fev1_abund, aes(y = host_fev1, x = host_clinical_state, fill = host_clinical_state)) +
  geom_boxplot(outlier.shape = NA,color = "#393939") +
  labs(x = "", y = "FEV1 %") +
  scale_fill_manual(values = c(
    "Baseline" = "gold",
    "Exacerbation" = "#ff684c",
    "Treatment" = "#51b364",
    "Recovery" = "#5fa2ce"
  )) +
  ylim(15,200) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(colour = "black", size = 13),
    axis.text.x = element_text(angle = 30, vjust = 0.6, size = 13),
    strip.text = element_text(size = 13, colour = "black"),
    axis.title.y = element_text(size = 13, colour = "black")
  ) +
  stat_compare_means(
    comparisons = clinical_states_comparisons,
    method = "wilcox.test",
    label = "p.signif",  
    hide.ns = FALSE, 
    tip.length = 0.01,
    vjust = 0.1
  )

# Top 5  Dominant Taxa---------------------------
ps_top5_alpha <- process_for_top5_taxa(counts, meta, taxa, 
                                                  taxrank = "Genus", 
                                                  lib.size.cut.off = 1000)

meta.tab.filt <- data.frame(sample_data(ps_top5_alpha))
meta.tab.filt <- meta.tab.filt %>% rownames_to_column(var = "ID")

genus.melt <- psmelt(ps_top5_alpha)
data <- genus.melt %>% select(host_clinical_state, Abundance, Genus)

aggregated_data <- data %>%
  group_by(host_clinical_state, Genus) %>%
  summarise(TotalCount = sum(Abundance), .groups = "drop")

top_bacteria <- aggregated_data %>%
  group_by(host_clinical_state) %>%
  slice_max(order_by = TotalCount, n = 5, with_ties = FALSE) %>%
  ungroup()

aggregated_data <- aggregated_data %>%
  group_by(host_clinical_state) %>%
  mutate(Genus = if_else(
    Genus %in% top_bacteria$Genus[top_bacteria$host_clinical_state == unique(host_clinical_state)],
    as.character(Genus),
    "Others"
  )) %>%
  ungroup()

aggregated_data <- aggregated_data %>%
  group_by(host_clinical_state, Genus) %>%
  summarise(TotalCount = sum(TotalCount), .groups = "drop")

aggregated_data <- aggregated_data %>%
  group_by(host_clinical_state) %>%
  mutate(Abundance = TotalCount / sum(TotalCount)) %>%
  ungroup()

aggregated_data$host_clinical_state <- factor(
  aggregated_data$host_clinical_state,
  levels = c("Recovery", "Treatment", "Exacerbation", "Baseline")
)

aggregated_data$Abundance <- aggregated_data$Abundance * 100

# plot
ggplot(aggregated_data, aes(fill = Genus, y = host_clinical_state, x = Abundance)) +
  geom_bar(stat = "identity") +
  ggthemes::scale_fill_tableau(palette = "Tableau 10", direction = -1) +
  labs(y = "", x = "Relative Abundance %", fill = "Genus") +
  theme_light() +
  theme(
    plot.margin = unit(c(0, 1.1, 0, 0), "cm"),
    axis.line = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.y = element_text(size = 13, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 13, angle = 0, vjust = 0.5, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    plot.title = element_blank(),
    legend.position = 'bottom',
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10.5),
    strip.text.x = element_text(size = 13, color = "black", face = "bold"),
    strip.background = element_rect(fill = NA, color = "black", linewidth = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )


# Alpha Diversities -----------------------
genus_alpha_div <- estimate_richness(ps_top5_alpha, measures = c("Observed",  "Shannon", "Simpson"))
genus_alpha_div <- round(genus_alpha_div, digits = 3)
colnames(genus_alpha_div) <- c("Observed Genera", "Shannon Index", "Simpson Index")
genus_alpha_div.meta <- genus_alpha_div %>% rownames_to_column(var = "ID") %>% 
  left_join(., meta.tab.filt %>% dplyr::select(ID, host_clinical_state), by = "ID")
genus_alpha_div.meta$host_clinical_state <- factor(genus_alpha_div.meta$host_clinical_state, 
                                                   levels = c("Baseline","Exacerbation", "Treatment", "Recovery"))

a <- reshape2::melt(genus_alpha_div.meta)

colors <- c(
  "Baseline" = "gold",
  "Exacerbation" = "#ff684c",
  "Treatment" = "#51b364",
  "Recovery" = "#5fa2ce"
)

# Observed index ----------------
observed <- a %>% filter(variable == "Observed Genera")

median_values_observed <- observed %>%
  group_by(host_clinical_state) %>%
  summarise(median_value = median(value, na.rm = TRUE)) %>%
  mutate(label = paste0("median=", round(median_value, 2)))

# plot
ggplot(observed, aes(x = host_clinical_state, y = value, fill = host_clinical_state)) +
  geom_jitter(aes(color = host_clinical_state), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),size = 3, alpha = 0.2) +
  geom_violin(fill = NA, color = "#393939", alpha = 0.4, width = 0.4, position = position_dodge(width = 0.5)) +
  geom_boxplot(fill = NA, color = "#393939", width = 0.15, position = position_dodge(width = 0.4), outlier.shape = NA) +
  geom_label(data = median_values_observed, aes(x = host_clinical_state, y = median_value, label = label),fill = "white",color = "black",label.size = 0.3,size = 3,nudge_x = 0.3,vjust = 0.5) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(y = "Observed Genera", x = " ", fill = "Genus") +
  stat_compare_means(comparisons = clinical_states_comparisons,
                     method = "wilcox.test",label = "p.signif",  
                     hide.ns = FALSE, tip.length = 0.01,vjust = 0.2) +
  theme_light() +
  theme(plot.margin = unit(c(0, 1.1, 0, 0), "cm"),
        axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 13, color = "black"),
        axis.title.x = element_text(size = 13, color = "black"),
        axis.text.x = element_text(size = 13, angle = 30, vjust = 0.6, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_blank(),
        strip.text.x = element_text(size = 13, color = "black", face = "bold"),
        strip.background = element_rect(fill = NA, color = "black", linewidth = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )



# shannon index-----------------
shannon <- a %>% filter(variable == "Shannon Index")

median_values_shannon <- shannon %>%
  group_by(host_clinical_state) %>%
  summarise(median_value = median(value, na.rm = TRUE)) %>%
  mutate(label = paste0("median=", round(median_value, 2)))

# plot
ggplot(shannon, aes(x = host_clinical_state, y = value, fill = host_clinical_state)) +
  geom_jitter(aes(color = host_clinical_state), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),size = 3, alpha = 0.2) +
  geom_violin(fill = NA, color = "#393939", alpha = 0.4, width = 0.4, position = position_dodge(width = 0.5)) +
  geom_boxplot(fill = NA, color = "#393939", width = 0.15, position = position_dodge(width = 0.4), outlier.shape = NA) +
  geom_label(data = median_values_shannon, aes(x = host_clinical_state, y = median_value, label = label),fill = "white",color = "black",label.size = 0.3,size = 3,nudge_x = 0.3,vjust = 0.5) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(y = "Shannon Index", x = " ", fill = "Genus") +
  stat_compare_means(comparisons = clinical_states_comparisons,
                     method = "wilcox.test",label = "p.signif",  
                     hide.ns = FALSE, tip.length = 0.01,vjust = 0.2) +
  theme_light() +
  theme(plot.margin = unit(c(0, 1.1, 0, 0), "cm"),
        axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 13, color = "black"),
        axis.title.x = element_text(size = 13, color = "black"),
        axis.text.x = element_text(size = 13, angle = 30, vjust = 0.6, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_blank(),
        strip.text.x = element_text(size = 13, color = "black", face = "bold"),
        strip.background = element_rect(fill = NA, color = "black", linewidth = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )


# Simpson index ------------------
simpson <- a %>% filter(variable == "Simpson Index")

median_values_simpson <- simpson %>%
  group_by(host_clinical_state) %>%
  summarise(median_value = median(value, na.rm = TRUE)) %>%
  mutate(label = paste0("median=", round(median_value, 2)))

# plot
ggplot(simpson, aes(x = host_clinical_state, y = value, fill = host_clinical_state)) +
  geom_jitter(aes(color = host_clinical_state), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),size = 3, alpha = 0.2) +
  geom_violin(fill = NA, color = "#393939", alpha = 0.4, width = 0.4, position = position_dodge(width = 0.5)) +
  geom_boxplot(fill = NA, color = "#393939", width = 0.15, position = position_dodge(width = 0.4), outlier.shape = NA) +
  geom_label(data = median_values_simpson, aes(x = host_clinical_state, y = median_value, label = label),fill = "white",color = "black",label.size = 0.3,size = 3,nudge_x = 0.3,vjust = 0.5) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(y = "Simpson Index", x = " ", fill = "Genus") +
  stat_compare_means(comparisons = clinical_states_comparisons,
                     method = "wilcox.test",label = "p.signif",  
                     hide.ns = FALSE, tip.length = 0.01,vjust = 0.2) +
  theme_light() +
  theme(plot.margin = unit(c(0, 1.1, 0, 0), "cm"),
        axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 13, color = "black"),
        axis.title.x = element_text(size = 13, color = "black"),
        axis.text.x = element_text(size = 13, angle = 30, vjust = 0.6, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_blank(),
        strip.text.x = element_text(size = 13, color = "black", face = "bold"),
        strip.background = element_rect(fill = NA, color = "black", linewidth = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )


# Beta Diversity ---------------------------------------------------------------
processed_for_beta <- process_for_beta(counts, meta, taxa, 
                                       taxrank = "Genus", 
                                       lib.size.cut.off = 1000)

genus_bray_dist = phyloseq::distance(processed_for_beta, method = "bray")
genus_ordination = ordinate(processed_for_beta, method = "PCoA", distance = genus_bray_dist)
genus_axis1.2 <- as.data.frame(genus_ordination$vectors[,1:2])
genus_axis1.2$rel_eigen <- genus_ordination$values$Relative_eig
colnames(genus_axis1.2) <- c("PC1", "PC2", "rel_eigen")
genus_axis1.2 <- genus_axis1.2 %>% rownames_to_column("ID")
genus_axis1.2$meta <- meta.tab.filt$host_clinical_state[match(genus_axis1.2$ID, meta.tab.filt$ID)]
genus_axis1.2$meta <- factor(genus_axis1.2$meta, levels = c("Baseline","Exacerbation", "Treatment", "Recovery"))
genus_axis1.2 <- drop_na(genus_axis1.2)


meta_colors <- c("gold","#ff684c","#51b364","#5fa2ce")
names(meta_colors) <- levels(genus_axis1.2$meta)

# plot
genus_axis1.2 %>%
  ggplot(aes(x = PC1, y = PC2, colour = meta)) +
  geom_point(size = 2.5, alpha = 0.6) + 
  stat_ellipse(linewidth = 0.8) + 
  scale_color_manual(name = "Clinical Status", labels = c("Baseline","Exacerbation", "Treatment", "Recovery"), values = meta_colors) +
  labs(x=paste0("PC1 (", round(genus_axis1.2$rel_eigen[1] * 100, digits = 2), "%)"), 
       y=paste0("PC2 (", round(genus_axis1.2$rel_eigen[2] * 100, digits = 2), "%)"),
       title = "Genus - PCoA (Bray-Curtis)") + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = 13, color = "black"),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        plot.title = element_text(size = 13),
        legend.text = element_text(size = 13, colour = "black"),
        legend.title = element_blank(),
        strip.text = element_text(size=13, colour = "black"),
        legend.position='right') 

# Pairwise PERMANOVA -----------------------
dist = phyloseq::distance(processed_for_beta, method = "bray")
ordination = ordinate(processed_for_beta, method = "PCoA", distance = dist)
meta.beta <- data.frame(sample_data(processed_for_beta))
cbn <- combn(x=unique(meta.beta$host_clinical_state), m = 2)
p <- c()
f_stat <- c()

set.seed(123)
for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(processed_for_beta, host_clinical_state %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  
  permanova_pairwise <- adonis2(
    phyloseq::distance(ps.subs, method = "bray") ~ host_clinical_state, 
    data = metadata_sub, permutations = 9999
  )
  print(permanova_pairwise)
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
  f_stat <- c(f_stat, permanova_pairwise$F[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), F_statistic = f_stat, p = p, p.adj = p.adj)[1:6,]
p.table <- as.data.frame(p.table)
p.table$p <- round(p.table$p, 3)
p.table$p.adj <- round(p.table$p.adj, 3)
p.table$F_statistic <- round(p.table$F_statistic, 3)
p.table <- as.data.frame(p.table)

# Comparison of Stretococcus and Gemella----------------------
ps <- process_for_beta(counts, meta, taxa, 
                      taxrank = "Genus", 
                      lib.size.cut.off = 1000)

meta.tab.filt <- data.frame(sample_data(ps))
meta.tab.filt <- meta.tab.filt %>% rownames_to_column(var = "ID")
otu.tab.filt <- data.frame(t(otu_table(ps)))
tax.tab.filt <- data.frame(tax_table(ps))
colnames(otu.tab.filt) <- tax.tab.filt$Genus[match(colnames(otu.tab.filt), rownames(tax.tab.filt))]
colnames(otu.tab.filt) <- make.names(colnames(otu.tab.filt), unique = TRUE)
otu.tab.filt <- otu.tab.filt %>% rownames_to_column(var = "ID")

data <- merge(otu.tab.filt, meta.tab.filt, by="ID")
sum(is.na(data))
data$host_clinical_state <- as.factor(data$host_clinical_state)
data[2:197] <- vegan::decostand(data[2:197], MARGIN = 1, method = "total") 
data$host_clinical_state <- factor(data$host_clinical_state, levels = c("Baseline","Exacerbation","Treatment","Recovery"))
star <- data %>% select(Streptococcus, Gemella, host_clinical_state)
star.long <- reshape2::melt(star)
star.long$value <- star.long$value * 100


clinical_states_comparisons <- list(
  c("Baseline", "Exacerbation"),
  c("Baseline", "Treatment"),
  c("Baseline", "Recovery"),
  c("Exacerbation", "Treatment"),
  c("Exacerbation", "Recovery"),
  c("Treatment", "Recovery")
)

# plot
ggplot(star.long, aes(fill = host_clinical_state, y = value, x = host_clinical_state)) +
  geom_boxplot(stat = "boxplot") +
  facet_wrap(~variable) +
  scale_fill_manual(values = c(
    "Baseline" = "gold",
    "Exacerbation" = "#ff684c",
    "Treatment" = "#51b364",
    "Recovery" = "#5fa2ce"
  )) +
  ylim(0,130) +
  labs(y = "Relative Abundance %", x = " ", fill = "Genus") +
  stat_compare_means(
    comparisons = clinical_states_comparisons,
    method = "wilcox.test",
    label = "p.signif",  
    hide.ns = FALSE, 
    tip.length = 0.01,
    vjust = 0.1
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(colour = "black", size = 13),
    axis.text.x = element_text(angle = 30, vjust = 0.6, size = 13),
    strip.text = element_text(size = 13, colour = "black"),
    strip.background = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.title.y = element_text(size = 13, colour = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) 



##################################################################
# Disease Agressiveness
##################################################################

# FEV1 Distribution------------------
disease_aggressiv_comparisons <- list(c("Mild", "Moderate_Severe"))
data_fev1_abund.2 <- process_fev1_abund(counts, meta, taxa, 
                                      taxrank = "Genus", 
                                      lib.size.cut.off = 1000)
data_fev1_abund.2 <- na.omit(data_fev1_abund.2, cols = "host_disease_aggressiveness")
data_fev1_abund.2 <- data_fev1_abund.2 %>% mutate(host_disease_aggressiveness = ifelse(host_disease_aggressiveness == "Moderate/Severe", "Moderate_Severe", host_disease_aggressiveness))

data_fev1_abund.2$host_disease_aggressiveness <- factor(data_fev1_abund.2$host_disease_aggressiveness, levels = c("Mild","Moderate_Severe"))

# plot
ggplot(data_fev1_abund.2, aes(y = host_fev1, x = host_disease_aggressiveness, fill = host_disease_aggressiveness)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "", y = "FEV1 %") +
  scale_fill_manual(values = c("Mild" = "#F6A245","Moderate_Severe" = "#7169A4")) +
  scale_x_discrete(labels = c("Mild" = "Mild","Moderate_Severe" = "Moderate/Severe")) +
  ylim(0, 115) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(colour = "black", size = 13),
    axis.text.x = element_text(angle = 30, vjust = 0.6, size = 13),
    strip.text = element_text(size = 13, colour = "black"),
    axis.title.y = element_text(size = 13, colour = "black")
  ) +
  stat_compare_means(
    comparisons = disease_aggressiv_comparisons,
    method = "wilcox.test",
    label = "p.signif",  
    hide.ns = FALSE, 
    tip.length = 0.01,
    vjust = 0.1
  )



# Top 5  Dominant Taxa in Disease agressiveness ----------------------
meta2 <- na.omit(meta, cols = "host_disease_aggressiveness")
ps_top5_alpha_disease_aggressiveness <- process_for_top5_taxa(counts, meta2, taxa, 
                                    taxrank = "Genus", 
                                    lib.size.cut.off = 1000)

meta.tab.filt.2 <- data.frame(sample_data(ps_top5_alpha_disease_aggressiveness))
meta.tab.filt.2 <- meta.tab.filt.2 %>% rownames_to_column(var = "ID")

genus.melt <- psmelt(ps_top5_alpha_disease_aggressiveness)
data <- genus.melt %>% select(host_disease_aggressiveness, Abundance, Genus)

aggregated_data <- data %>%
  group_by(host_disease_aggressiveness, Genus) %>%
  summarise(TotalCount = sum(Abundance), .groups = "drop")

top_bacteria <- aggregated_data %>%
  group_by(host_disease_aggressiveness) %>%
  slice_max(order_by = TotalCount, n = 5, with_ties = FALSE) %>%
  ungroup()

aggregated_data <- aggregated_data %>%
  group_by(host_disease_aggressiveness) %>%
  mutate(Genus = if_else(
    Genus %in% top_bacteria$Genus[top_bacteria$host_disease_aggressiveness == unique(host_disease_aggressiveness)],
    as.character(Genus),
    "Others"
  )) %>%
  ungroup()

aggregated_data <- aggregated_data %>%
  group_by(host_disease_aggressiveness, Genus) %>%
  summarise(TotalCount = sum(TotalCount), .groups = "drop")

aggregated_data <- aggregated_data %>%
  group_by(host_disease_aggressiveness) %>%
  mutate(Abundance = TotalCount / sum(TotalCount)) %>%
  ungroup()

aggregated_data$host_disease_aggressiveness <- factor(
  aggregated_data$host_disease_aggressiveness,
  levels = c("Mild", "Moderate/Severe")
)

aggregated_data$Abundance <- aggregated_data$Abundance * 100

# plot
ggplot(aggregated_data, aes(fill = Genus, y = host_disease_aggressiveness, x = Abundance)) +
  geom_bar(stat = "identity") +
  ggthemes::scale_fill_tableau(palette = "Tableau 10", direction = -1) +
  labs(y = "", x = "Relative Abundance %", fill = "Genus") +
  theme_light() +
  theme(
    plot.margin = unit(c(0, 1.1, 0, 0), "cm"),
    axis.line = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.y = element_text(size = 13, color = "black"),
    axis.title.x = element_text(size = 13, color = "black"),
    axis.text.x = element_text(size = 13, angle = 0, vjust = 0.5, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    plot.title = element_blank(),
    legend.position = 'bottom',
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 10.5),
    strip.text.x = element_text(size = 13, color = "black", face = "bold"),
    strip.background = element_rect(fill = NA, color = "black", linewidth = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )


# Alpha Diversities -----------------------
genus_alpha_div <- estimate_richness(ps_top5_alpha_disease_aggressiveness, measures = c("Observed",  "Shannon", "Simpson"))
genus_alpha_div <- round(genus_alpha_div, digits = 3)
colnames(genus_alpha_div) <- c("Observed Genera", "Shannon Index", "Simpson Index")
genus_alpha_div.meta <- genus_alpha_div %>% rownames_to_column(var = "ID") %>% 
  left_join(., meta.tab.filt.2 %>% dplyr::select(ID, host_disease_aggressiveness), by = "ID")
genus_alpha_div.meta <- na.omit(genus_alpha_div.meta, cols = "host_disease_aggressiveness")
genus_alpha_div.meta <- genus_alpha_div.meta %>% mutate(host_disease_aggressiveness = ifelse(host_disease_aggressiveness == "Moderate/Severe", "Moderate_Severe", host_disease_aggressiveness))
genus_alpha_div.meta$host_disease_aggressiveness <- factor(genus_alpha_div.meta$host_disease_aggressiveness, 
                                                   levels = c("Mild","Moderate_Severe"))

b <- reshape2::melt(genus_alpha_div.meta)

colors <- c("Mild" = "#F6A245","Moderate_Severe" = "#7169A4")

# Observed index ----------------
observed <- b %>% filter(variable == "Observed Genera")

median_values_observed <- observed %>%
  group_by(host_disease_aggressiveness) %>%
  summarise(median_value = median(value, na.rm = TRUE)) %>%
  mutate(label = paste0("median=", round(median_value, 2)))

# plot
ggplot(observed, aes(x = host_disease_aggressiveness, y = value, fill = host_disease_aggressiveness)) +
  geom_jitter(aes(color = host_disease_aggressiveness), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),size = 3, alpha = 0.2) +
  geom_violin(fill = NA, color = "#393939", alpha = 0.4, width = 0.4, position = position_dodge(width = 0.5)) +
  geom_boxplot(fill = NA, color = "#393939", width = 0.15, position = position_dodge(width = 0.4), outlier.shape = NA) +
  geom_label(data = median_values_observed, aes(x = host_disease_aggressiveness, y = median_value, label = label),fill = "white",color = "black",label.size = 0.3,size = 3,nudge_x = 0.3,vjust = 0.5) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  scale_x_discrete(labels = c("Mild" = "Mild","Moderate_Severe" = "Moderate/Severe")) +
  ylim(0,60) +
  labs(y = "Observed Genera", x = " ", fill = "Genus") +
  stat_compare_means(comparisons = disease_aggressiv_comparisons,
                     method = "wilcox.test",label = "p.signif",  
                     hide.ns = FALSE, tip.length = 0.01,vjust = 0.2) +
  theme_light() +
  theme(plot.margin = unit(c(0, 1.1, 0, 0), "cm"),
        axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 13, color = "black"),
        axis.title.x = element_text(size = 13, color = "black"),
        axis.text.x = element_text(size = 13, angle = 30, vjust = 0.6, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_blank(),
        strip.text.x = element_text(size = 13, color = "black", face = "bold"),
        strip.background = element_rect(fill = NA, color = "black", linewidth = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )



# shannon index-----------------
shannon <- b %>% filter(variable == "Shannon Index")

median_values_shannon <- shannon %>%
  group_by(host_disease_aggressiveness) %>%
  summarise(median_value = median(value, na.rm = TRUE)) %>%
  mutate(label = paste0("median=", round(median_value, 2)))

# plot
ggplot(shannon, aes(x = host_disease_aggressiveness, y = value, fill = host_disease_aggressiveness)) +
  geom_jitter(aes(color = host_disease_aggressiveness), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),size = 3, alpha = 0.2) +
  geom_violin(fill = NA, color = "#393939", alpha = 0.4, width = 0.4, position = position_dodge(width = 0.5)) +
  geom_boxplot(fill = NA, color = "#393939", width = 0.15, position = position_dodge(width = 0.4), outlier.shape = NA) +
  geom_label(data = median_values_shannon, aes(x = host_disease_aggressiveness, y = median_value, label = label),fill = "white",color = "black",label.size = 0.3,size = 3,nudge_x = 0.3,vjust = 0.5) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(y = "Shannon Index", x = " ", fill = "Genus") +
  ylim(0,2.85) +
  stat_compare_means(comparisons = disease_aggressiv_comparisons,
                     method = "wilcox.test",label = "p.signif",  
                     hide.ns = FALSE, tip.length = 0.01,vjust = 0.2) +
  theme_light() +
  theme(plot.margin = unit(c(0, 1.1, 0, 0), "cm"),
        axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 13, color = "black"),
        axis.title.x = element_text(size = 13, color = "black"),
        axis.text.x = element_text(size = 13, angle = 30, vjust = 0.6, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_blank(),
        strip.text.x = element_text(size = 13, color = "black", face = "bold"),
        strip.background = element_rect(fill = NA, color = "black", linewidth = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )



# Simpson index ------------------
simpson <- b %>% filter(variable == "Simpson Index")

median_values_simpson <- simpson %>%
  group_by(host_disease_aggressiveness) %>%
  summarise(median_value = median(value, na.rm = TRUE)) %>%
  mutate(label = paste0("median=", round(median_value, 2)))

# plot
ggplot(simpson, aes(x = host_disease_aggressiveness, y = value, fill = host_disease_aggressiveness)) +
  geom_jitter(aes(color = host_disease_aggressiveness), position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),size = 3, alpha = 0.2) +
  geom_violin(fill = NA, color = "#393939", alpha = 0.4, width = 0.4, position = position_dodge(width = 0.5)) +
  geom_boxplot(fill = NA, color = "#393939", width = 0.15, position = position_dodge(width = 0.4), outlier.shape = NA) +
  geom_label(data = median_values_simpson, aes(x = host_disease_aggressiveness, y = median_value, label = label),fill = "white",color = "black",label.size = 0.3,size = 3,nudge_x = 0.3,vjust = 0.5) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(y = "Simpson Index", x = " ", fill = "Genus") +
  ylim(0,1.0) +
  stat_compare_means(comparisons = disease_aggressiv_comparisons,
                     method = "wilcox.test",label = "p.signif",  
                     hide.ns = FALSE, tip.length = 0.01,vjust = 0.2) +
  theme_light() +
  theme(plot.margin = unit(c(0, 1.1, 0, 0), "cm"),
        axis.line = element_line(colour = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 13, color = "black"),
        axis.title.x = element_text(size = 13, color = "black"),
        axis.text.x = element_text(size = 13, angle = 30, vjust = 0.6, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        plot.title = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_blank(),
        strip.text.x = element_text(size = 13, color = "black", face = "bold"),
        strip.background = element_rect(fill = NA, color = "black", linewidth = 1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )


# Beta Diversity ---------------------------------------------------------------
processed_for_beta <- process_for_beta(counts, meta2, taxa, 
                                       taxrank = "Genus", 
                                       lib.size.cut.off = 1000)

genus_bray_dist = phyloseq::distance(processed_for_beta, method = "bray")
genus_ordination = ordinate(processed_for_beta, method = "PCoA", distance = genus_bray_dist)
genus_axis1.2 <- as.data.frame(genus_ordination$vectors[,1:2])
genus_axis1.2$rel_eigen <- genus_ordination$values$Relative_eig
colnames(genus_axis1.2) <- c("PC1", "PC2", "rel_eigen")
genus_axis1.2 <- genus_axis1.2 %>% rownames_to_column("ID")
genus_axis1.2$meta2 <- meta.tab.filt.2$host_disease_aggressiveness[match(genus_axis1.2$ID, meta.tab.filt.2$ID)]
genus_axis1.2$meta2 <- factor(genus_axis1.2$meta2, levels = c("Mild", "Moderate/Severe"))
genus_axis1.2 <- drop_na(genus_axis1.2)


# Pairwise PERMANOVA -----------------------
dist = phyloseq::distance(processed_for_beta, method = "bray")
ordination = ordinate(processed_for_beta, method = "PCoA", distance = dist)
meta.beta <- data.frame(sample_data(processed_for_beta))
cbn <- combn(x=unique(meta.beta$host_disease_aggressiveness), m = 2)
p <- c()

set.seed(123)
for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(processed_for_beta, host_disease_aggressiveness %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = "bray") ~ host_disease_aggressiveness, 
                                data = metadata_sub, permutations = 9999)
  print(permanova_pairwise)
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)[1:2,]
p.table <- as.data.frame(p.table)
p.table$p <- round(p.table$p, 3)
p.table$p.adj <- round(p.table$p.adj, 3)


meta_colors <- c("#F6A245","#7169A4")
names(meta_colors) <- levels(genus_axis1.2$meta2)

# plot
genus_axis1.2 %>%
  ggplot(aes(x = PC1, y = PC2, colour = meta2)) +
  geom_point(size = 2.5, alpha = 0.6) + 
  stat_ellipse(linewidth = 0.8) + 
  scale_color_manual(name = "Disease Aggressiveness", labels = c("Mild","Moderate/Severe"), values = meta_colors) +
  labs(x=paste0("PC1 (", round(genus_axis1.2$rel_eigen[1] * 100, digits = 2), "%)"), 
       y=paste0("PC2 (", round(genus_axis1.2$rel_eigen[2] * 100, digits = 2), "%)"),
       title = "Genus - PCoA (Bray-Curtis)") + theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = 13, color = "black"),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        plot.title = element_text(size = 13),
        legend.text = element_text(size = 13, colour = "black"),
        legend.title = element_blank(),
        strip.text = element_text(size=13, colour = "black"),
        legend.position='right') +
  annotate("text", label = "F = 4.312, p = 0.003", x = -0.4, y = -0.75, size = 4)



