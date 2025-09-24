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

data <- read.table(input_file, header = TRUE, as.is = TRUE, sep = "\t")
data <- build_data(data)

plt <- ggplot(data = data, aes(x = sample_name, y = errorRate, fill = errorCategory)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c(deletions = "#bfbfbf", insertions = "#ffa64d", mismatches = "#1a1aff")) +
    ylab("# Errors per mapped base") +
    xlab("") +
    guides(fill = guide_legend(title = "Error class")) +
    scale_y_continuous(labels = label_scientific(digits = 1)) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

plt <- facets(plt, length(unique(dat$sample_name)))
publish(plt, args[2], length(unique(dat$sample_name)))


deletions_only <- subset(data, errorCategory == "deletions")
plt <- ggplot(data = deletions_only, aes(x = sample_name, y = errorRate, fill = errorCategory)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c(deletions = "#bfbfbf")) +
    ylab("# Errors per mapped base") +
    xlab("") +
    guides(fill = guide_legend(title = "Error class")) +
    scale_y_continuous(labels = label_scientific(digits = 1)) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

plt <- facets(plt, length(unique(dat$sample_name)))
publish(plt, args[3], length(unique(dat$sample_name)))
