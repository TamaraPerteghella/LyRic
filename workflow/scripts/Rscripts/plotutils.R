#!/usr/bin/env Rscript
library(ggplot2)

annotation <- read.table(paste0(Sys.getenv("BASE_PROJECT_DIR"), "/config/sample_annotations.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

feature1 <- colnames(annotation)[3]
feature2 <- colnames(annotation)[4]
feature3 <- colnames(annotation)[5]

build_data <- function(data) {
    if (ncol(annotation) < 3) {
        return(data)
    }

    data <- merge(data, annotation, by = "sample_name", all = TRUE)
    return(data)
}


facets <- function(plt, samples) {
    if (ncol(annotation) == 2) {
        return(plt)
    } else if (ncol(annotation) == 3) {
        facet_organisation <- paste("~", feature1)
    } else if (ncol(annotation) == 4) {
        facet_organisation <- paste(feature1, "~", feature2)
    } else {
        facet_organisation <- paste(feature1, "~", feature2, "+", feature3)
    }

    plt <- plt + facet_wrap(as.formula(facet_organisation), nrow = samples)
    return(plt)
}


publish <- function(plt, output_file, samples) {
    plt <- plt + theme_minimal(base_size = 10) +
        theme(
            axis.text = element_text(size = rel(1)),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.ticks = element_line(size = 2),
            axis.line = element_line(colour = "#595959", size = 0.5),
            axis.title = element_text(size = rel(1)),
            panel.grid.major = element_line(colour = "#d9d9d9", size = 0.3),
            panel.grid.minor = element_line(colour = "#e6e6e6", size = 0.3),
            panel.border = element_blank(), panel.background = element_blank(),
            strip.background = element_rect(colour = "#737373", fill = "white"),
            legend.key.size = unit(0.5, "line"),
            legend.title = element_text(size = rel(1)),
            legend.text = element_text(size = rel(1)),
            strip.text = element_text(size = rel(1))
        )

    ggsave(plot = plt, filename = output_file, dpi = 300, height = 5, width = 5, units = "in", limitsize = FALSE, scale = samples)
}
