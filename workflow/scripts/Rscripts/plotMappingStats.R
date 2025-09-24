library(ggplot2)
library(cowplot)
library(plyr)
library(scales)
library(gridExtra)
library(grid)
library(ggplotify)

dir <- Sys.getenv("BASE_PROJECT_DIR")
source(paste0(dir, "/workflow/scripts/Rscripts/plotutils.R"))

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]


data <- read.table(input_file, header = TRUE, as.is = TRUE, sep = "\t")
data <- build_data(data)


plt <- ggplot(data, aes(x = sample_name, y = percentMappedReads, fill = sizeFrac)) +
    geom_bar(width = 0.75, stat = "identity", position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = {
        sizeFrac_Rpalette
    }) +
    geom_hline(aes(yintercept = 1), linetype = "dashed", alpha = 0.7, size = lineSize) +
    geom_text(aes(group = sizeFrac, y = 0.01, label = paste(sep = "", percent(percentMappedReads), "\n", "(", comma(mappedReads), ")")), angle = 90, size = geom_textSize, hjust = 0, vjust = 0.5, position = position_dodge(width = 0.9)) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    xlab("") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

plt <- facets(plt, length(unique(dat$sample_name)))
publish(plt, output_file, length(unique(dat$sample_name)))
