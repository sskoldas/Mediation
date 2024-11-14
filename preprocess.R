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
#Import data and build phyloseq object # 
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
metadata$host_clinical_state <- factor(metadata$host_clinical_state, levels = c("Baseline", "Exacerbation", "Treatment","Recovery"))
meta <- metadata

################################################################################
# Preprocessing # 
################################################################################

process_phyloseq_data <- function(counts, meta, taxa, 
                                  lib.size.cut.off = 1000, 
                                  prevalence.cut.off = 0.1,
                                  mean.prop.cut.off = 0.0001,
                                  taxrank = NULL,
                                  relative.abundance = TRUE,
                                  clr.transform = TRUE) {
  ps <- phyloseq(
    otu_table(counts, taxa_are_rows = TRUE),
    sample_data(meta),
    phyloseq::tax_table(as.matrix(taxa))
  )
  print(table(tax_table(ps)[, "Kingdom"]))
  
  # Subset to only include Bacteria
  ps <- subset_taxa(ps,!is.na(Kingdom) & Kingdom != "" & Kingdom != "Archaea" & Kingdom != "Eukaryota" & Family != "Mitochondria" & Order != "Chloroplast")
  print(table(tax_table(ps)[, "Kingdom"]))
  
  # Filter samples with fewer than lib.size.cut.off reads
  ps <- prune_samples(sample_sums(ps) > lib.size.cut.off, ps)
  
  # Filter OTU/ASVs showing less than 10% of the samples
  presence_absence <- otu_table(ps) > 0
  samples_with_taxon <- rowSums(presence_absence)
  n_samples <- ncol(otu_table(ps))
  prevalence <- samples_with_taxon / n_samples
  taxa_to_keep_prevalence <- prevalence >= prevalence.cut.off
  ps <- prune_taxa(taxa_to_keep_prevalence, ps)
  
  # Filter OTU/ASVs based on mean proportion cut-off 
  prop.otu.tab <- apply(otu_table(ps), 2, function(x) x / sum(x))
  mean.prop <- rowMeans(prop.otu.tab)
  ps <- prune_taxa(mean.prop > mean.prop.cut.off, ps)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  # Aggregate OTUs at the specified taxonomic level, if provided
  if (!is.null(taxrank)) {
    ps <- tax_glom(ps, taxrank = taxrank)
  }
  
  if (relative.abundance) {
    # Handling zero and normalization
    otu.tab <- as.data.frame(t(otu_table(ps)))
    otu.tab[otu.tab == 0] <- 0.5
    otu.norm <- as.data.frame(decostand(otu.tab, MARGIN = 1, method = "total"))
    meta.tab <- as.data.frame(sample_data(ps))
    tax.tab <- as.data.frame(tax_table(ps))
    
    if (clr.transform){
      otu.norm <- as.data.frame(decostand(otu.norm, MARGIN = 1, method = "clr"))
    }
    
    # Recreate the phyloseq object with the compositional OTU table
    ps <- phyloseq(otu_table(t(otu.norm), taxa_are_rows = TRUE),
                   sample_data(meta.tab),
                   phyloseq::tax_table(as.matrix(tax.tab)))
  }
  
  return(ps)
}




####################################################################################################
# Subsetting "host clinical state" as binary
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

