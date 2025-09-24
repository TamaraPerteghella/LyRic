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

maxY <- max(dat$percent)
plt <- ggplot(data = dat, aes(x = sample_name, y = percent, fill = category)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(position = position_dodge(width = 0.9), size = geom_textSize, aes(y = 0, label = paste(sep = "", percent(percent), "\n", "(", comma(count), ")")), hjust = 0, vjust = 0.5, angle = 90) +
    scale_fill_manual(values = c("ERCCs" = "#e4b5ff", "SIRVs" = "#5edba9")) +
    ylab("% reads mapped on\nspike-in sequences") +
    xlab("") +
    guides(fill = guide_legend(title = "Spike-in set")) +
    scale_y_continuous(labels = percent) +
    expand_limits(y = c(0, maxY)) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())


plt <- facets(plt, length(unique(dat$sample_name)))
publish(plt, output_file, length(unique(dat$sample_name)))
