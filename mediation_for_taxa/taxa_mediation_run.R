setwd("~/Desktop/Mediation/")
source("./mediation_for_taxa/preprocess.R")
source("./mediation_for_taxa/debiasedLasso.R")
source("./mediation_for_taxa/parameter_estimation.R")
source("./mediation_for_taxa/taxa_mediation_model.R")
source("./BCa_CI.R")


# run the preprocess
processed <- process_phyloseq_data(counts, meta, taxa, 
                                   taxrank = "Genus", 
                                   lib.size.cut.off = 1000,
                                   prevalence.cut.off = 0.05)

# effect could be "indirect_effects", "direct_effect" or "total_effect"
run_for_clinical_state <- function(a,b, effect = C("indirect_effects","direct_effect","total_effect")) {
  
  binary <- subset_phyloseq_samples(processed = processed, treatment = "host_clinical_state", binary = c(a,b))
  if (any(is.na(binary))) {
    stop("Data contains missing values. Please handle them before proceeding.")
  }
  binary$isolate <- as.character(binary$isolate)
  # Filter isolates with at least two samples
  df_filtered <- binary %>%
    group_by(isolate) %>%
    filter(n() >= 2) %>%
    ungroup()
  
  df_filtered$isolate <- droplevels(as.factor(df_filtered$isolate))
  group_sizes <- table(df_filtered$isolate)
  
  if (any(group_sizes < 2)) {
    stop("There are still groups with fewer than two observations.")
  }
  
  ID <- df_filtered$isolate
  Y <- df_filtered$host_fev1 
  X <- df_filtered$host_clinical_state
  cov <- df_filtered$Host_age
  otu <- df_filtered %>%
    dplyr::select(-c(isolate, host_fev1, Host_age, host_clinical_state)) %>%
    as.matrix()
  
  #run mediation analysis with high-dimensional, compositional data
  res <- mediation_with_multilevel(
    X=X, Y=Y, OTU = otu, COV = cov, grp = ID,
    n_boot = 3000,
    n_cores = parallel::detectCores() - 1)
  
  for (i in effect) {
    if (!is.null(res[[i]])) { 
      save_file <- res[[i]]
      file_path <- paste0(
        "./Outputs/", 
        a, "_", b, "_", i, "_table.tsv"
      )
      write.table(save_file, file_path, row.names = TRUE, col.names = NA)
    } else {
      warning(paste("Effect", i, "not found in results. Skipping."))
    }
  }
  return(res)
}


run_for_clinical_state(a="Baseline", b="Exacerbation", effect = c("indirect_effects","direct_effect","total_effect"))
run_for_clinical_state(a="Baseline", b="Treatment", effect = c("indirect_effects","direct_effect","total_effect"))
run_for_clinical_state(a="Baseline", b="Recovery", effect = c("indirect_effects","direct_effect","total_effect"))
run_for_clinical_state(a="Exacerbation", b="Treatment", effect = c("indirect_effects","direct_effect","total_effect"))
run_for_clinical_state(a="Exacerbation",b="Recovery", effect = c("indirect_effects","direct_effect","total_effect"))
run_for_clinical_state(a="Treatment", b="Recovery", effect = c("indirect_effects","direct_effect","total_effect"))

#-----------------------------------------------------------------------------------------------------------------------------------
run_for_disease_aggressiveness <- function(a,b, effect = C("indirect_effects","direct_effect","total_effect")) {
  
  binary <- subset_phyloseq_samples(processed = processed, treatment = "host_disease_aggressiveness", binary = c(a,b))
  binary <- na.omit(binary)
  binary <- binary %>% select(-host_clinical_state)
  if (any(is.na(binary))) {
    stop("Data contains missing values. Please handle them before proceeding.")
  }
  binary$isolate <- as.character(binary$isolate)
  # Filter isolates with at least two samples
  df_filtered <- binary %>%
    group_by(isolate) %>%
    filter(n() >= 2) %>%
    ungroup()
  
  df_filtered$isolate <- droplevels(as.factor(df_filtered$isolate))
  group_sizes <- table(df_filtered$isolate)
  
  if (any(group_sizes < 2)) {
    stop("There are still groups with fewer than two observations.")
  }
  
  ID <- df_filtered$isolate
  Y <- df_filtered$host_fev1 
  X <- df_filtered$host_disease_aggressiveness
  cov <- df_filtered$Host_age
  otu <- df_filtered %>%
    dplyr::select(-c(isolate, host_fev1, Host_age, host_disease_aggressiveness)) %>%
    as.matrix()
  
  #run mediation analysis with high-dimensional, compositional data
  res <- mediation_with_multilevel(
    X=X, Y=Y, OTU = otu, COV = cov, grp = ID,
    n_boot = 3000,
    n_cores = parallel::detectCores() - 1)
  
  for (i in effect) {
    if (!is.null(res[[i]])) { 
      save_file <- res[[i]]
      file_path <- paste0(
        "./Outputs/", 
        a, "_", b, "_", i, "_host_disease_aggressiveness_table.tsv"
      )
      write.table(save_file, file_path, row.names = TRUE, col.names = NA)
    } else {
      warning(paste("Effect", i, "not found in results. Skipping."))
    }
  }
  return(res)
}

run_for_disease_aggressiveness(a="Mild", b="Moderate_Severe", effect = c("indirect_effects","direct_effect","total_effect"))
