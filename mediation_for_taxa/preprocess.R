rm(list = ls())

library(dplyr)
library(vegan)
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(nlme)
library(glmnet)
library(MASS)
library(scalreg)
library(doParallel)
library(foreach)
library(ggpubr)
library(rstatix)
################################################################################
# Data Manipulation # 
################################################################################

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
table(metadata$host_disease_aggressiveness)
metadata$host_clinical_state <- factor(metadata$host_clinical_state, levels = c("Baseline", "Exacerbation", "Treatment","Recovery"))
metadata <- metadata %>%
  mutate(
    host_disease_aggressiveness = ifelse(
      host_disease_aggressiveness == "Moderate/Severe", 
      "Moderate_Severe", 
      host_disease_aggressiveness
    )
  )
metadata$host_disease_aggressiveness <- factor(metadata$host_disease_aggressiveness, levels = c("Mild","Moderate_Severe"))
meta <- metadata

################################################################################
# Pre-processing # 
################################################################################
process_phyloseq_data <- function(counts, meta, taxa, 
                                  taxrank = NULL,
                                  lib.size.cut.off = 1000, 
                                  prevalence.cut.off = 0.1){
  
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
  # Filter genera having low abundance
  threshold <- ceiling(prevalence.cut.off * nrow(agg.otu))
  low.abundance.filt <- agg.otu[,colSums(agg.otu > 3) >= threshold]
  
  # compositional otu (total sum equals to 1)
  otu.clean.total <- vegan::decostand(low.abundance.filt, MARGIN = 1, method = "total")
  # log norm transformation
  offset <- 0.5
  lognorm_input <- otu.clean.total + offset
  lognorm <- log10(lognorm_input)
  # Recreate the phyloseq object with the compositional OTU table
  ps.lognorm <- phyloseq(otu_table(t(lognorm), taxa_are_rows = TRUE),
                         sample_data(agg.meta),
                         phyloseq::tax_table(as.matrix(agg.taxa)))
  return(ps.lognorm)
}




################################################################################
# Subsetting "host clinical state" as binary
################################################################################

subset_phyloseq_samples <- function(processed, treatment, binary) {
  if (!treatment %in% colnames(sample_data(processed))) {
    stop("The treatment column is not found in the sample data.")
  }
  
  selected_samples <- sample_data(processed)[[treatment]] %in% binary
  subset_data <- prune_samples(selected_samples, processed)
  
  otutab <- as.data.frame(t(otu_table(subset_data)))
  metatab <- data.frame(sample_data(subset_data))
  taxtab <- as.data.frame(tax_table(subset_data), stringsAsFactors = FALSE)
  
  if ("Genus" %in% colnames(taxtab)) {
    if (all(colnames(otutab) == rownames(taxtab))) {
      colnames(otutab) <- taxtab$Genus[match(colnames(otutab), rownames(taxtab))]
    } else {
      stop("Error: colnames(otutab) and rownames(taxtab) do not match.")
    }
  } else {
    warning("'Genus' column not found in tax_table. OTU IDs will be used as column names.")
  }
  
  merged_data <- cbind(otutab, metatab)
  merged_data[[treatment]] <- ifelse(merged_data[[treatment]] == binary[1], 0, 1)
  
  return(merged_data)
}

