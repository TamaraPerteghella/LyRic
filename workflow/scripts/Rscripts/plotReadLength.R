#!/usr/bin/env Rscript
library(ggplot2)
library(scales)
library(gridExtra)
library(grid)
library(dplyr)
library(data.table)

dir <- Sys.getenv("BASE_PROJECT_DIR")
source(paste0(dir, "/workflow/scripts/Rscripts/plotutils.R"))

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

dat <- fread(input_file, header = TRUE, sep = "\t")

datSum <- dat %>%
    group_by(sample_name) %>%
    summarise(n = n(), med = median(length))
summaryStats <- transform(datSum, LabelN = paste0("N = ", comma(n)), LabelM = paste0("Median = ", comma(med)))

plt <- ggplot(dat, aes(x = length)) +
    geom_histogram(aes(y = ..density..), binwidth = 100) +
    geom_vline(data = summaryStats, aes(xintercept = med), color = "#ff0055", linetype = "solid", size = 2) +
    geom_text(data = summaryStats, aes(label = LabelN, x = Inf, y = Inf), hjust = 1, vjust = 1, size = 5, fontface = "bold") +
    geom_text(data = summaryStats, aes(label = LabelM, x = Inf, y = Inf), hjust = 1, vjust = 2.5, size = 5, fontface = "bold", color = "#ff0055") +
    coord_cartesian(xlim = c(0, 3500)) +
    scale_y_continuous(labels = scientific) +
    scale_x_continuous(labels = comma, name = "Read length (nts)")

publish(plt, output_file, length(unique(dat$sample_name)))
