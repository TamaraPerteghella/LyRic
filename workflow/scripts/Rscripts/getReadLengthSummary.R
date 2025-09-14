library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
output_file <- args[2]

dat <- fread(input_file, header=T, sep='\t')

datSum <- dat %>%
 group_by(file_name) %>%
 summarise(n=n(), median=median(length), mean=mean(length), max=max(length))

write_tsv(datSum, output_file)