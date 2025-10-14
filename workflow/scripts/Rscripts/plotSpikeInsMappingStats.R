#!/usr/bin/env Rscript
library(ggplot2)
library(scales)
library(gridExtra)
library(grid)

dir <- Sys.getenv("BASE_PROJECT_DIR")
source(paste0(dir, "/workflow/scripts/Rscripts/plotutils.R"))

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]


data <- read.table(input_file, header = TRUE, as.is = TRUE, sep = "\t")
data <- build_data(data)

maxY <- max(data$percent)

plt <- ggplot(data, aes(x = sample_name, y = percent, fill = category)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(position = position_dodge(width = 0.9), size = 5, aes(y = 0, label = paste(sep = "", percent(percent), "\n", "(", comma(count), ")")), hjust = 0, vjust = 0.5, angle = 90) +
    scale_fill_manual(values = c("ERCCs" = "#e4b5ff", "SIRVs" = "#5edba9")) +
    ylab("% reads mapped on\nspike-in sequences") +
    xlab("") +
    guides(fill = guide_legend(title = "Spike-in set")) +
    scale_y_continuous(labels = percent) +
    expand_limits(y = c(0, maxY)) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

plt <- facets(plt, length(unique(data$sample_name)))
publish(plt, output_file, length(unique(data$sample_name)))
