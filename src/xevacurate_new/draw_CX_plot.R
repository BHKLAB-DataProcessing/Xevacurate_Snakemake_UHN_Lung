library(Xeva)         # or your plotting package with plotPDX
library(gridExtra)    # for arranging plots
library(ggplot2)      # in case plotPDX returns ggplot objects

setwd("/Users/guanqiaofeng/Documents/BHK/Xeva/Xevacurate_202507all/results/explore")

# Load batch list
batch_list <- readLines("../build/all_tax_batch.txt")

# readin xevaset
x.set <- readRDS("../xevaset/UHN_TNBC_ALL_202507_DrugResponse.rds")

# Prepare output PDF
pdf("pdx_batch_plots_TAX_50d.pdf", width = 12, height = 18)  # 2 columns × 3 rows layout

# Initialize plot collector
plot_pairs <- list()

for (i in seq_along(batch_list)) {
  batch_id <- batch_list[i]
  
  # Generate the 2 plots for this batch
  p1 <- plotPDX(x.set,
                batch = batch_id,
                vol.normal = TRUE,
                control.col = "#a6611a",
                treatment.col = "#018571",
                major.line.size = 1,
                max.time = 50,
                SE.plot = "ribbon") +
    ggtitle(paste(batch_id, "- Normalized")) +
    theme(legend.position = "none")
  
  p2 <- plotPDX(x.set,
                batch = batch_id,
                vol.normal = FALSE,
                control.col = "#a6611a",
                treatment.col = "#018571",
                major.line.size = 1,
                max.time = 50,
                SE.plot = "ribbon") +
    ggtitle(paste(batch_id, "- Raw Volume")) +
    theme(legend.position = "none")
  
  # Add both plots to the collector
  plot_pairs <- c(plot_pairs, list(p1, p2))
  
  # Once 6 plots are collected (3 batches), draw and reset
  if (length(plot_pairs) == 6 || i == length(batch_list)) {
    do.call("grid.arrange", c(plot_pairs, ncol = 2, nrow = 3))
    plot_pairs <- list()  # Reset for next page
  }
}

# Close the PDF
dev.off()
