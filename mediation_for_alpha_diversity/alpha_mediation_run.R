
source("./mediation_for_alpha_diversity/preprocess_for_alpha_diversity.R")
source("./mediation_for_alpha_diversity/alpha_mediation_model.R")
source("./BCa_CI.R")

process_for_alpha <- preprocess_alpha_diversity(counts, meta, taxa, 
                                                  taxrank = "Genus", 
                                                  lib.size.cut.off = 1000)

meta.tab.filt <- data.frame(sample_data(process_for_alpha))
meta.tab.filt <- meta.tab.filt %>% rownames_to_column(var = "ID")
genus_alpha_div <- estimate_richness(process_for_alpha, measures = c("Observed",  "Shannon", "Simpson"))
genus_alpha_div <- round(genus_alpha_div, digits = 3)
colnames(genus_alpha_div) <- c("Observed_Genera", "Shannon_Index", "Simpson_Index")
genus_alpha_div.meta <- genus_alpha_div %>% rownames_to_column(var = "ID") %>% 
  left_join(., meta.tab.filt %>% dplyr::select(ID, host_clinical_state, host_disease_aggressiveness, host_fev1, Host_age, isolate), by = "ID")
genus_alpha_div.meta$host_clinical_state <- factor(genus_alpha_div.meta$host_clinical_state, levels = c("Baseline","Exacerbation", "Treatment", "Recovery"))
genus_alpha_div.meta$host_disease_aggressiveness <- factor(genus_alpha_div.meta$host_disease_aggressiveness, levels = c("Mild","Moderate_Severe"))

process <- function(data, filter_col, filter_vals) {
  data <- data %>% filter(!!sym(filter_col) %in% filter_vals)
  data[[quo_name(filter_col)]] <- ifelse(data[[quo_name(filter_col)]] == filter_vals[1], 0, 1)
  data[["isolate"]] <- as.character(data[["isolate"]])
  data <- data %>% group_by(isolate) %>% filter(n() >= 2) %>% ungroup()
  data[["isolate"]] <- droplevels(as.factor(data[["isolate"]]))
  return(data)
}

# Mild-Moderate/Severe
aggress <- process(data = genus_alpha_div.meta, filter_col = "host_disease_aggressiveness", filter_vals = c("Mild","Moderate_Severe"))
X = aggress[["host_disease_aggressiveness"]]
Y = aggress[["host_fev1"]]
COV = aggress[["Host_age"]]
grp = aggress[["isolate"]]
observed_genera = aggress[["Observed_Genera"]]
shannon = aggress[["Shannon_Index"]]
simpson = aggress[["Simpson_Index"]]

observed_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=observed_genera, n_boot = 3000, n_cores = parallel::detectCores()-1)
shannon_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=shannon, n_boot = 3000, n_cores = parallel::detectCores()-1)
simpson_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=simpson, n_boot = 3000, n_cores = parallel::detectCores()-1)

observed_me <- data.frame(Diversity = rep("Observed",3), observed_me)
shannon_me <- data.frame(Diversity = rep("Shannon",3), shannon_me)
simpson_me <- data.frame(Diversity= rep("Simpson",3),simpson_me)

aggress_alpha_diversity <- rbind(observed_me,shannon_me,simpson_me)
write.table(aggress_alpha_diversity, "./Outputs/Mild_Moderate_Severe_alpha_diversity.tsv", sep = "\t", row.names = F, col.names = T)



## baseline-exacerbation
be <- process(data = genus_alpha_div.meta, filter_col = "host_clinical_state", filter_vals = c("Baseline","Exacerbation"))
X = be[["host_clinical_state"]]
Y = be[["host_fev1"]]
COV = be[["Host_age"]]
grp = be[["isolate"]]
observed_genera = be[["Observed_Genera"]]
shannon = be[["Shannon_Index"]]
simpson = be[["Simpson_Index"]]

observed_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=observed_genera, n_boot = 3000, n_cores = parallel::detectCores()-1)
shannon_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=shannon, n_boot = 3000, n_cores = parallel::detectCores()-1)
simpson_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=simpson, n_boot = 3000, n_cores = parallel::detectCores()-1)

observed_me <- data.frame(Diversity = rep("Observed",3), observed_me)
shannon_me <- data.frame(Diversity = rep("Shannon",3), shannon_me)
simpson_me <- data.frame(Diversity= rep("Simpson",3),simpson_me)

be_alpha_diversity <- rbind(observed_me,shannon_me,simpson_me)
write.table(be_alpha_diversity, "./Outputs/baseline_exacerbation_alpha_diversity.tsv", sep = "\t", row.names = F, col.names = T)



## baseline-treatment
bt <- process(data = genus_alpha_div.meta, filter_col = "host_clinical_state", filter_vals = c("Baseline","Treatment"))
X = bt[["host_clinical_state"]]
Y = bt[["host_fev1"]]
COV = bt[["Host_age"]]
grp = bt[["isolate"]]
observed_genera = bt[["Observed_Genera"]]
shannon = bt[["Shannon_Index"]]
simpson = bt[["Simpson_Index"]]

observed_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=observed_genera, n_boot = 3000, n_cores = parallel::detectCores()-1)
shannon_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=shannon, n_boot = 3000, n_cores = parallel::detectCores()-1)
simpson_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=simpson, n_boot = 3000, n_cores = parallel::detectCores()-1)

observed_me <- data.frame(Diversity = rep("Observed",3), observed_me)
shannon_me <- data.frame(Diversity = rep("Shannon",3), shannon_me)
simpson_me <- data.frame(Diversity= rep("Simpson",3),simpson_me)

bt_alpha_diversity <- rbind(observed_me,shannon_me,simpson_me)
write.table(bt_alpha_diversity, "./Outputs/baseline_treatment_alpha_diversity.tsv", sep = "\t", row.names = F, col.names = T)



## baseline-recovery
br <- process(data = genus_alpha_div.meta, filter_col = "host_clinical_state", filter_vals = c("Baseline","Recovery"))
X = br[["host_clinical_state"]]
Y = br[["host_fev1"]]
COV = br[["Host_age"]]
grp = br[["isolate"]]
observed_genera = br[["Observed_Genera"]]
shannon = br[["Shannon_Index"]]
simpson = br[["Simpson_Index"]]

observed_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=observed_genera, n_boot = 3000, n_cores = parallel::detectCores()-1)
shannon_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=shannon, n_boot = 3000, n_cores = parallel::detectCores()-1)
simpson_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=simpson, n_boot = 3000, n_cores = parallel::detectCores()-1)

observed_me <- data.frame(Diversity = rep("Observed",3), observed_me)
shannon_me <- data.frame(Diversity = rep("Shannon",3), shannon_me)
simpson_me <- data.frame(Diversity= rep("Simpson",3),simpson_me)

br_alpha_diversity <- rbind(observed_me,shannon_me,simpson_me)
write.table(br_alpha_diversity, "./Outputs/baseline_recovery_alpha_diversity.tsv", sep = "\t", row.names = F, col.names = T)



## exacerbation-treatment
et <- process(data = genus_alpha_div.meta, filter_col = "host_clinical_state", filter_vals = c("Exacerbation","Treatment"))
X = et[["host_clinical_state"]]
Y = et[["host_fev1"]]
COV = et[["Host_age"]]
grp = et[["isolate"]]
observed_genera = et[["Observed_Genera"]]
shannon = et[["Shannon_Index"]]
simpson = et[["Simpson_Index"]]

observed_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=observed_genera, n_boot = 3000, n_cores = parallel::detectCores()-1)
shannon_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=shannon, n_boot = 3000, n_cores = parallel::detectCores()-1)
simpson_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=simpson, n_boot = 3000, n_cores = parallel::detectCores()-1)

observed_me <- data.frame(Diversity = rep("Observed",3), observed_me)
shannon_me <- data.frame(Diversity = rep("Shannon",3), shannon_me)
simpson_me <- data.frame(Diversity= rep("Simpson",3),simpson_me)

et_alpha_diversity <- rbind(observed_me,shannon_me,simpson_me)
write.table(et_alpha_diversity, "./Outputs/exacerbation_treatment_alpha_diversity.tsv", sep = "\t", row.names = F, col.names = T)



## exacerbation-recovery
er <- process(data = genus_alpha_div.meta, filter_col = "host_clinical_state", filter_vals = c("Exacerbation","Recovery"))
X = er[["host_clinical_state"]]
Y = er[["host_fev1"]]
COV = er[["Host_age"]]
grp = er[["isolate"]]
observed_genera = er[["Observed_Genera"]]
shannon = er[["Shannon_Index"]]
simpson = er[["Simpson_Index"]]

observed_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=observed_genera, n_boot = 3000, n_cores = parallel::detectCores()-1)
shannon_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=shannon, n_boot = 3000, n_cores = parallel::detectCores()-1)
simpson_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=simpson, n_boot = 3000, n_cores = parallel::detectCores()-1)

observed_me <- data.frame(Diversity = rep("Observed",3), observed_me)
shannon_me <- data.frame(Diversity = rep("Shannon",3), shannon_me)
simpson_me <- data.frame(Diversity= rep("Simpson",3),simpson_me)

er_alpha_diversity <- rbind(observed_me,shannon_me,simpson_me)
write.table(er_alpha_diversity, "./Outputs/exacerbation_recovery_alpha_diversity.tsv", sep = "\t", row.names = F, col.names = T)



## treatment-recovery
tr <- process(data = genus_alpha_div.meta, filter_col = "host_clinical_state", filter_vals = c("Treatment","Recovery"))
X = tr[["host_clinical_state"]]
Y = tr[["host_fev1"]]
COV = tr[["Host_age"]]
grp = tr[["isolate"]]
observed_genera = tr[["Observed_Genera"]]
shannon = tr[["Shannon_Index"]]
simpson = tr[["Simpson_Index"]]

observed_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=observed_genera, n_boot = 3000, n_cores = parallel::detectCores()-1)
shannon_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=shannon, n_boot = 3000, n_cores = parallel::detectCores()-1)
simpson_me <- mediation_analysis_for_alpha_diversity(X, Y, COV, grp, alpha_diversity=simpson, n_boot = 3000, n_cores = parallel::detectCores()-1)

observed_me <- data.frame(Diversity = rep("Observed",3), observed_me)
shannon_me <- data.frame(Diversity = rep("Shannon",3), shannon_me)
simpson_me <- data.frame(Diversity= rep("Simpson",3),simpson_me)

tr_alpha_diversity <- rbind(observed_me,shannon_me,simpson_me)
write.table(tr_alpha_diversity, "./Outputs/treatment_recovery_alpha_diversity.tsv", sep = "\t", row.names = F, col.names = T)

