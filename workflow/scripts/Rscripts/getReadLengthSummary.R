#!/usr/bin/env Rscript
library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
output_file <- args[2]

dat <- fread(input_file, header=T, sep='\t')

datSum <- dat %>%
 group_by(sample_name) %>%
 summarise(n=n(), median=median(length), mean=mean(length), max=max(length))

write.table(datSum, output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)