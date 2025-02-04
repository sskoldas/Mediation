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

table(metadata$host_disease_aggressiveness)
table(metadata$host_clinical_state)
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
# Pre-processing for alpha diversity# 
################################################################################

preprocess_alpha_diversity <- function(counts, meta, taxa, 
                                      taxrank = NULL,
                                      lib.size.cut.off = 1000){
  
  ps <- phyloseq(otu_table(counts,taxa_are_rows = TRUE),
                 sample_data(meta),
                 phyloseq::tax_table(as.matrix(taxa)))
  # Filter out unwanted taxa
  ps <- subset_taxa(ps,!is.na(Kingdom) & Kingdom != "" & Kingdom != "Archaea" & Kingdom != "Eukaryota" & Family != "Mitochondria" & Order != "Chloroplast")
  # Aggregates taxa 
  ps.agg <- tax_glom(ps, taxrank = taxrank)
  # Filter samples with fewer than lib.size.cut.off reads
  ps.agg <- prune_samples(sample_sums(ps.agg) >= lib.size.cut.off, ps.agg)
  ps.agg <- prune_taxa(taxa_sums(ps.agg) > 0, ps.agg)
  return(ps.agg)
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

