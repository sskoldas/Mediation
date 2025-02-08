

library(tibble)
library(forestplot)
library(tidyverse)

# example
indirect <- list.files(path = "./Outputs/", pattern = "indirect")

a <- read.csv("./Outputs/Mild_Moderate_Severe_indirect_effects_host_disease_aggressiveness_table.tsv", row.names=1, sep="")
a$ci_range <- paste0("[", sprintf("%.3f", a$CI_Lower), ", ", sprintf("%.3f", a$CI_Upper), "]")
a[sapply(a, is.numeric)] <- round(a[sapply(a, is.numeric)], 3)
a <- a %>% arrange(Mediators)
a <- a %>% filter(Estimate >= 1 | Estimate <= -1 | p_value < 0.1)

output_base_data <- tibble(mean = a$Estimate,
                           lower = a$CI_Lower,
                           upper = a$CI_Upper,
                           taxa = a$Mediators,
                           estimate = as.character(a$Estimate),
                           p_value = as.character(a$p_value),
                           ci_range = a$ci_range)

output_header <- tibble(taxa = "Mediators",
                        estimate = "Estimate",
                        p_value = "P-value",
                        ci_range = "CI [Lower-Upper]",
                        summary = TRUE)

output_df <- bind_rows(output_header,
                       output_base_data)

# forest plot
the_plot <- output_df %>% 
  forestplot(labeltext = c(taxa, estimate, p_value, ci_range),
             is.summary = summary,
             xlab = "95% Confidence Interval",
             txt_gp = fpTxtGp(title = gpar(cex=1.35),
                              ticks = gpar(cex=0.75),
                              xlab = gpar(cex = 1),
                              label = gpar(fontfamily="", cex=0.9)),
             hrzl_lines = list("2" = gpar(lty = 2), "49" = gpar(lwd = 1, columns = 1:5, col = "#000044")),
             title = "Mild-Moderate/Severe",
             col = fpColors(box = "#D55E00", line="black"),
             boxsize = 0.3,
             alpha = 0.75)
the_plot



################################################################################
# create data table for direct and total effect
library(stringr)
library(tidyverse)
total <- list.files(path = "./Outputs/", pattern = "total")
direct <- list.files(path = "./Outputs/", pattern = "_direct")
total.direct <- c(total, direct)

names <- c()
read_file_list <- list()
for (i in seq_along(total.direct)){
  names[i] <- str_extract(total.direct[i], ".*effect")
  file_path <- paste0("./Outputs/", total.direct[i])
  read_file <- read.csv(file_path, row.names=1, sep ="")
  rownames(read_file) <- names[i]
  read_file_list[[i]] <- read_file
}

total_direct_merge <- do.call(rbind, read_file_list)
total_direct_merge[] <- lapply(total_direct_merge, function(x){
  if(is.numeric(x)) round(x, digits = 3) else x
})
