#!/usr/bin/env Rscript
library(ggplot2)
library(scales)
library(gridExtra)
library(grid)
library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]
output_file <- args[2]

dat <- fread(input_file, header=T, sep='\t')

datSum <- dat %>% group_by(sample_name) %>% summarise(n=n(), med=median(length))
summaryStats <- transform(datSum, LabelN = paste0('N= ', comma(n)), LabelM = paste0( 'Median= ', comma(med)))

plt <- ggplot(dat, aes(x=length)) +
geom_histogram(aes(y=..density..), binwidth=100) +
geom_vline(data = summaryStats, aes(xintercept=med), color='#ff0055', linetype='solid', size=2) +
geom_text(data = summaryStats, aes(label = LabelN, x = Inf, y = Inf), hjust=1, vjust=1, size=15, fontface = 'bold') +
geom_text(data = summaryStats, aes(label = LabelM, x = med, y = Inf), hjust=-0.1, vjust=2.5, size=15, fontface = 'bold', color='#ff0055') +
coord_cartesian(xlim=c(0, 3500)) +
scale_y_continuous(labels=scientific) +
scale_x_continuous(labels=comma, name='Read length (nts)') +
theme(axis.text= element_text(size=20), axis.ticks = element_line(size=2), 
    axis.line = element_line(colour = '#595959', size=2), axis.title=element_text(size = 30), 
    panel.grid.major = element_line(colour='#d9d9d9', size=2), panel.grid.minor = element_line(colour='#e6e6e6', size=2), 
    panel.border = element_blank(), panel.background = element_blank(), 
    strip.background = element_rect(colour='#737373',fill='white'), legend.key.size=unit(0.5,'line'), 
    legend.title=element_text(size=15), legend.text=element_text(size=15), strip.text = element_text(size = 15)) + 
theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(plot = plt, filename = output_file)
