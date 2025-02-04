setwd("~/Desktop/Mediation/")
library(grid)
library(forestplot)
library(dplyr)
library(tibble)

#example
files <- list.files(path = "./Outputs/", pattern = "alpha")

a <- read.delim("./Outputs/Mild_Moderate_Severe_alpha_diversity.tsv", row.names=NULL)
names(a)[1] <- "Alpha Diversity"
a <- a %>% mutate(p_value = ifelse(p_value == 0.000, "<0.001", p_value))

output_base_data <- tibble(mean = a$Estimate,
                           lower = a$CI_Lower,
                           upper = a$CI_Upper,
                           alpha = a$`Alpha Diversity`,
                           estimate = as.character(a$Estimate),
                           effects = a$Effects,
                           p_value = as.character(a$p_value))

output_header <- tibble(alpha = "Alpha Diversity",
                        effects = "Effects",
                        estimate = "Estimate",
                        p_value = "P-value",
                        summary = TRUE)

output_df <- bind_rows(output_header,
                       output_base_data)

# Create the forest plots for each index
observed_plot <- output_df |> filter(alpha == "Observed") |>
  forestplot(labeltext = c(effects, estimate, p_value),
             title = "Observed Genera",
             xlab = "95% Confidence Interval",
             boxsize = 0.2,
             alpha = 0.75,
             hrzl_lines = list("2" = gpar(lty = 2), "5" = gpar(lwd = 1, columns = 1:3, col = "#000044")),
             txt_gp = fpTxtGp(title = gpar(cex=1.15),
                              ticks = gpar(cex=0.75),
                              xlab = gpar(cex = 0.75),
                              label = gpar(fontfamily="", cex=0.9)),
             new_page = FALSE) |>
  fp_add_header(effects = "Effects", estimate = "Estimate", p_value = "p-value") |>
  fp_decorate_graph(box = gpar(lty=2, col="lightgray"), graph.pos =4)


shannon_plot <- output_df |> filter(alpha == "Shannon") |>
  forestplot(labeltext = c(effects, estimate, p_value),
             title = "Shannon Index",
             xlab = "95% Confidence Interval",
             boxsize = 0.2,
             alpha = 0.75,
             hrzl_lines = list("2" = gpar(lty = 2), "5" = gpar(lwd = 1, columns = 1:3, col = "#000044")),
             txt_gp = fpTxtGp(title = gpar(cex=1.15),
                              ticks = gpar(cex=0.75),
                              xlab = gpar(cex = 0.75),
                              label = gpar(fontfamily="", cex=0.9)),
             new_page = FALSE) |>
  fp_add_header(effects = "Effects", estimate = "Estimate", p_value = "p-value") |>
  fp_decorate_graph(box = gpar(lty=2, col="lightgray"), graph.pos =4)


simpson_plot <- output_df |> filter(alpha == "Simpson") |>
  forestplot(labeltext = c(effects, estimate, p_value),
             title = "Simpson Index",
             xlab = "95% Confidence Interval",
             boxsize = 0.2,
             alpha = 0.75,
             hrzl_lines = list("2" = gpar(lty = 2), "5" = gpar(lwd = 1, columns = 1:3, col = "#000044")),
             txt_gp = fpTxtGp(title = gpar(cex=1.15),
                              ticks = gpar(cex=0.75),
                              xlab = gpar(cex = 0.75),
                              label = gpar(fontfamily="", cex=0.9)),
             new_page = FALSE) |>
  fp_add_header(effects = "Effects", estimate = "Estimate", p_value = "p-value") |>
  fp_decorate_graph(box = gpar(lty=2, col="lightgray"), graph.pos =4)


# Set up the grid layout to accommodate the main title and three plots
grid.newpage()
borderWidth <- unit(4, "pt")
totalBorderWidth <- borderWidth * 2
totalPlotWidth <- unit(1, "npc") - totalBorderWidth
plotWidth <- unit(convertX(totalPlotWidth, "npc", valueOnly = TRUE) / 3, "npc")

# Adjust the grid layout to have two rows: one for the title and one for the plots
pushViewport(viewport(layout = grid.layout(
  nrow = 2,
  ncol = 5,
  heights = unit.c(unit(1.5, "lines"), unit(1, "null")),
  widths = unit.c(plotWidth, borderWidth, plotWidth, borderWidth, plotWidth)
)))

# Place the main title in the first row spanning all columns
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:5))
grid.text("Mild-Moderate/Severe", gp = gpar(fontsize = 15, fontface = "bold"), just = "center")
upViewport()

# Plot the Observed Genera in the first column (second row)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
observed_plot |>
  fp_set_style(box = "royalblue", line = "darkblue", summary = "royalblue")
upViewport()

# Add a border between plots
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.rect(gp = gpar(fill = "#dddddd", col = "#eeeeee"))
upViewport()

# Plot the Shannon Index in the third column (second row)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
shannon_plot |>
  fp_set_style(box = "royalblue", line = "darkblue", summary = "royalblue")
upViewport()

# Add another border
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 4))
grid.rect(gp = gpar(fill = "#dddddd", col = "#eeeeee"))
upViewport()

# Plot the Simpson Index in the fifth column (second row)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 5))
simpson_plot |>
  fp_set_style(box = "royalblue", line = "darkblue", summary = "royalblue")
upViewport(2)

