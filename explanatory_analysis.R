library(ggstatsplot) ; packageVersion("ggstatsplot")


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


# FEV1 Distribution Between Clinical States-------------------------------------
processed_fev1_abund <- function(counts, meta, taxa, 
                                 taxrank = NULL,
                                 lib.size.cut.off = 1000){
  
  ps <- phyloseq(otu_table(counts,taxa_are_rows = TRUE),
                 sample_data(meta),
                 phyloseq::tax_table(as.matrix(taxa)))
  # Filter out unwanted taxa
  ps <- subset_taxa(ps,!is.na(Kingdom) & Kingdom != "" & Kingdom != "Archaea" & Kingdom != "Eukaryota" & Family != "Mitochondria" & Order != "Chloroplast" & Genus != "uncultured")
  # Aggregates taxa 
  ps.agg <- tax_glom(ps, taxrank = taxrank)
  # Filter samples with fewer than lib.size.cut.off reads
  ps.agg <- prune_samples(sample_sums(ps.agg) >= lib.size.cut.off, ps.agg)
  
  agg.meta <- data.frame(sample_data(ps.agg))
  agg.taxa <- data.frame(tax_table(ps.agg))
  agg.otu <- data.frame(t(otu_table(ps.agg)))
  agg.otu <- agg.otu[, colSums(agg.otu != 0) > 0]
  
  merge <- merge(agg.meta, agg.otu, by = 0)
  return(merge)
}

data_fev1_abund <- processed_fev1_abund(counts, meta, taxa, 
                                            taxrank = "Genus", 
                                            lib.size.cut.off = 1000)

my_comparisons <- list(
  c("Baseline", "Exacerbation"),
  c("Baseline", "Treatment"),
  c("Baseline", "Recovery"),
  c("Exacerbation", "Treatment"),
  c("Exacerbation", "Recovery"),
  c("Treatment", "Recovery")
)

ggplot(data_fev1_abund, aes(y = host_fev1, x = host_clinical_state, fill = host_clinical_state)) +
  geom_boxplot(outlier.shape = NA) +
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
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(colour = "black", size = 13),
    axis.text.x = element_text(angle = 30, vjust = 0.7, size = 13),
    strip.text = element_text(size = 13, colour = "black"),
    axis.title.y = element_text(size = 13, colour = "black")
  ) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",  
    hide.ns = FALSE, 
    tip.length = 0.01,
    vjust = 0.2
  )


# Top 5  Dominant Taxa in Clinical States --------------------------------------
process <- function(counts, meta, taxa, 
                    taxrank = NULL,
                    lib.size.cut.off = 1000){
  
  ps <- phyloseq(otu_table(counts,taxa_are_rows = TRUE),
                 sample_data(meta),
                 phyloseq::tax_table(as.matrix(taxa)))
  # Filter out unwanted taxa
  ps <- subset_taxa(ps,!is.na(Kingdom) & Kingdom != "" & Kingdom != "Archaea" & Kingdom != "Eukaryota" & Family != "Mitochondria" & Order != "Chloroplast" & Genus != "uncultured")
  # Aggregates taxa 
  ps.agg <- tax_glom(ps, taxrank = taxrank)
  # Filter samples with fewer than lib.size.cut.off reads
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

processed_for_top5_alpha <- process(counts, meta, taxa, 
                                    taxrank = "Genus", 
                                    lib.size.cut.off = 1000)

meta.tab.filt <- data.frame(sample_data(processed_for_top5_alpha))
meta.tab.filt <- meta.tab.filt %>% rownames_to_column(var = "ID")

genus.melt <- psmelt(processed_for_top5_alpha)
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
    axis.title.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    plot.title = element_blank(),
    legend.position = 'bottom',
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10.5),
    strip.text.x = element_text(size = 12, color = "black", face = "bold"),
    strip.background = element_rect(fill = NA, color = "black", linewidth = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )


# Alpha Diversities ------------------------------------------------------------
genus_alpha_div <- estimate_richness(processed_for_top5_alpha, measures = c("Observed",  "Shannon", "Simpson"))
genus_alpha_div <- round(genus_alpha_div, digits = 3)
colnames(genus_alpha_div) <- c("Observed Genera", "Shannon Index", "Simpson Index")
genus_alpha_div.meta <- genus_alpha_div %>% rownames_to_column(var = "ID") %>% 
  left_join(., meta.tab.filt %>% dplyr::select(ID, host_clinical_state), by = "ID")
genus_alpha_div.meta$host_clinical_state <- factor(genus_alpha_div.meta$host_clinical_state, levels = c("Baseline","Exacerbation", "Treatment", "Recovery"))

# Observed Genera --------------------------------------------------------------
ggbetweenstats(
  data = genus_alpha_div.meta,
  x = host_clinical_state,
  y = 'Observed Genera',
  type = "nonparametric",
  var.equal = FALSE,
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = TRUE,
  bf.message = FALSE,
  results.subtitle = FALSE,
  p.adjust.method = "BH",
  xlab = "",
  ylab = "Observed Genera",
  ggsignif.args = list(textsize = 5, tip_length = 0.01, na.rm = TRUE, vjust = 0, step_increase = 0.1),
  boxplot.args = list(width = 0.2, alpha = 0.2, na.rm = TRUE),
  violin.args = list(width = 0.4, alpha = 0.2, na.rm = TRUE, linewidth = 0.6),
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), 
                    alpha = 0.7, size = 3, stroke = 0, na.rm = TRUE), package = "ggthemes", palette = "Classic_10") +
  ylim(0,90) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        strip.text = element_text(size = 15),
        plot.title = element_text(size = 20),
        legend.position = 'none') + scale_color_manual(values=c("gold","#ff684c","#51b364","#5fa2ce"))


# Shannon Index ----------------------------------------------------------------
ggbetweenstats(
  data = genus_alpha_div.meta,
  x = host_clinical_state,
  y = 'Shannon Index',
  type = "nonparametric",
  var.equal = FALSE,
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = TRUE,
  bf.message = FALSE,
  results.subtitle = FALSE,
  p.adjust.method = "BH",
  xlab = "",
  ylab = "Shannon Diversity",
  ggsignif.args = list(textsize = 5, tip_length = 0.01, na.rm = TRUE, step_increase = 0.1),
  boxplot.args = list(width = 0.2, alpha = 0.2, na.rm = TRUE),
  violin.args = list(width = 0.4, alpha = 0.2, na.rm = TRUE, linewidth = 0.6),
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha =
                      0.7, size = 3, stroke = 0, na.rm = TRUE),
  package = "ggthemes",
  palette = "Classic_10_Medium") +
  ylim(0,5) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        strip.text = element_text(size = 15),
        plot.title = element_text(size = 20),
        legend.position = 'none') + scale_color_manual(values=c("gold","#ff684c","#51b364","#5fa2ce"))




# Simpson Index ----------------------------------------------------------------
ggbetweenstats(
  data = genus_alpha_div.meta,
  x = host_clinical_state,
  y = 'Simpson Index',
  type = "nonparametric",
  var.equal = FALSE,
  plot.type = "violin",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = TRUE,
  bf.message = FALSE,
  results.subtitle = FALSE,
  p.adjust.method = "BH",
  xlab = "",
  ylab = "Simpson Diversity",
  ggsignif.args = list(textsize = 5, tip_length = 0.01, na.rm = TRUE, step_increase = 0.1),
  boxplot.args = list(width = 0.2, alpha = 0.2, na.rm = TRUE),
  violin.args = list(width = 0.4, alpha = 0.2, na.rm = TRUE,linewidth = 0.6),
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha =
                      0.7, size = 3, stroke = 0, na.rm = TRUE),
  package = "ggthemes",
  palette = "Classic_10_Medium") +
  ylim(0.0, 1.7) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 13, colour = "black"),
        axis.text.y = element_text(size = 13, colour = "black"),
        strip.text = element_text(size = 15),
        plot.title = element_text(size = 20),
        legend.position = 'none') + scale_color_manual(values=c("gold","#ff684c","#51b364","#5fa2ce"))


# Beta Diversity ---------------------------------------------------------------
process_for_beta <- function(counts, meta, taxa, 
                                  taxrank = NULL,
                                  lib.size.cut.off = 1000){
  
  ps <- phyloseq(otu_table(counts,taxa_are_rows = TRUE),
                 sample_data(meta),
                 phyloseq::tax_table(as.matrix(taxa)))
  # Filter out unwanted taxa
  ps <- subset_taxa(ps,!is.na(Kingdom) & Kingdom != "" & Kingdom != "Archaea" & Kingdom != "Eukaryota" & Family != "Mitochondria" & Order != "Chloroplast" & Genus != "uncultured")
  # Aggregates taxa 
  ps.agg <- tax_glom(ps, taxrank = taxrank)
  # Filter samples with fewer than lib.size.cut.off reads
  ps.agg <- prune_samples(sample_sums(ps.agg) >= lib.size.cut.off, ps.agg)
  
  agg.meta <- data.frame(sample_data(ps.agg))
  agg.taxa <- data.frame(tax_table(ps.agg))
  agg.otu <- data.frame(t(otu_table(ps.agg)))
  agg.otu <- agg.otu[, colSums(agg.otu != 0) > 0]
  
  # compositional otu (total sum equals to 1)
  otu.total <- vegan::decostand(agg.otu, MARGIN = 1, method = "total")
  # Recreate the phyloseq object with the compositional OTU table
  ps <- phyloseq(otu_table(t(otu.total), taxa_are_rows = TRUE),
                     sample_data(agg.meta),
                     phyloseq::tax_table(as.matrix(agg.taxa)))
  return(ps)
}

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
        axis.title = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, colour = "black"),
        axis.text.y = element_text(size = 15, colour = "black"),
        plot.title = element_text(size = 15),
        legend.text = element_text(size = 13, colour = "black"),
        legend.title = element_text(size = 13, colour = "black"),
        strip.text = element_text(size=13, colour = "black"),
        legend.position='right') 

# Pairwise PERMANOVA -----------------------------------------------------------
dist = phyloseq::distance(processed_for_beta, method = "bray")
ordination = ordinate(processed_for_beta, method = "PCoA", distance = dist)
meta.beta <- data.frame(sample_data(processed_for_beta))
cbn <- combn(x=unique(meta.beta$host_clinical_state), m = 2)
p <- c()

set.seed(123)
for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(processed_for_beta, host_clinical_state %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = "bray") ~ host_clinical_state, 
                                data = metadata_sub, permutations = 9999)
  print(permanova_pairwise)
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)[1:6,]
p.table <- as.data.frame(p.table)

