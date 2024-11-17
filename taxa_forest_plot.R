library(tibble)
library(forestplot)
library(tidyverse)

# example
a <- read.csv("./Outputs/Treatment_Recovery_indirect_effects_table.tsv", row.names=1, sep="")
a$ci_range <- paste0("[", sprintf("%.3f", a$CI_Lower), ", ", sprintf("%.3f", a$CI_Upper), "]")
a[sapply(a, is.numeric)] <- round(a[sapply(a, is.numeric)], 3)
a <- a %>% arrange(Mediators)

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
             hrzl_lines = list("2" = gpar(lty = 2), "39" = gpar(lwd = 1, columns = 1:5, col = "#000044")),
             title = "Treatment-Recovery",
             col = fpColors(box = "#D55E00", line="black"),
             boxsize = 0.3,
             alpha = 0.75)
the_plot
