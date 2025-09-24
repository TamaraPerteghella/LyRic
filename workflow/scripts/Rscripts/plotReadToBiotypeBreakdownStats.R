library(ggplot2)
library(cowplot)
library(plyr)
library(scales)
library(gridExtra)
library(grid)
library(ggplotify)
library(data.table)

dir <- Sys.getenv("BASE_PROJECT_DIR")
source(paste0(dir, "/workflow/scripts/Rscripts/plotutils.R"))

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

data <- read.table(input_file, header = TRUE, as.is = TRUE, sep = "\t")
data <- build_data(data)

data$biotype <- factor(data$biotype)

plt <- ggplot(data, aes(x = sample_name, y = readOverlapsPercent, fill = biotype)) +
    geom_bar(stat = "identity") +
    ylab("% read overlaps") +
    xlab("") +
    guides(fill = guide_legend(title = "Region/biotype")) +
    scale_y_continuous(labels = scales::percent) +
    coord_cartesian(ylim = c(0, 1)) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

plt <- facets(plt, length(unique(dat$sample_name)))
publish(plt, output_file, length(unique(dat$sample_name)))
