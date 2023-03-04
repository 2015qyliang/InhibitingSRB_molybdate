# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(tidyverse)
library(microeco)
library(GUniFrac)
library(ape)
library(ggplot2)
library(ggsci)
library(doBy)
library(gridExtra)

##############################################################################

diybox = function(v.meas) {
  alphaDF = alphaDF.raw[alphaDF.raw$Measure %in% v.meas, ]
  alphaDF$Measure = factor(alphaDF$Measure, levels = v.meas)
  
  alphaDF$Group = substr(alphaDF$SampleID, 4, 4)
  alphaDF$Time = as.numeric(substr(alphaDF$SampleID, 2, 3))
  
  p = ggplot(alphaDF, aes(x = Time, y = Value, fill = Group, group = Time)) + 
    geom_boxplot(alpha = 0.7, outlier.color = 'grey50', 
                 outlier.size = 0.5, size = 0.1) + 
    facet_grid(facets = Group~Measure, rows = 1) +
    scale_fill_manual(values = c('#E69F00', '#0072B2')) + 
    scale_x_continuous(breaks = c(0, 5, 12, 21, 30)) + 
    labs(x = "", y = "") + 
    ylim(min(alphaDF$Value), 1.1*max(alphaDF$Value)) +
    theme_bw()+
    theme(panel.grid = element_blank(), 
          strip.text = element_text(size = 10, color = 'black'), 
          strip.background = element_blank(), 
          axis.text.y = element_text(size = 8, hjust = 1, vjust = 0.5, color = 'black'),
          axis.text.x = element_text(size = 8, hjust = 0.5, vjust = 0.5, angle = 0, color = 'black'),
          legend.position = "none")
  return(p)
  
}


##############################################################################

sample_info = read.table('02-samplesGroup.txt', header = T, 
                         sep = '\t', stringsAsFactors = F)
rownames(sample_info) = sample_info$SampleID
otu_table = read.table('02-OTUtable.txt', header = T, sep = '\t',
                       row.names = 1, stringsAsFactors = F)
taxonomy_table = read.table('02-comTax.txt',  header = T, sep = '\t', 
                            row.names = 1, stringsAsFactors = F)
phylo_tree = read.tree('02-FastTree.txt')
dataset = microtable$new(sample_table = sample_info, 
                         otu_table = otu_table, 
                         tax_table = taxonomy_table, 
                         phylo_tree = phylo_tree)
dataset$cal_alphadiv(PD = T)

##############################################################################

t1 = trans_alpha$new(dataset = dataset, group = "Group")
alphaDF.raw = t1$data_alpha[, c('Measure', 'Value', 'SampleID', 'Group')]

alphaDF.output = alphaDF.raw
alphaDF.output$Value = round(alphaDF.output$Value, digits = 3)
write.table(alphaDF.output, '04-AlphaDF.txt',
            append = F, quote = F, sep = '\t',
            row.names = F, col.names = T)


##############################################################################


vec.Measure = c('Observed', 'Simpson', 'Shannon')

gridList = list()
gridList[[1]] = diybox(vec.Measure[1])
gridList[[2]] = diybox(vec.Measure[2])
gridList[[3]] = diybox(vec.Measure[3])  

p = grid.arrange(gridList[[1]], gridList[[2]], gridList[[3]], nrow = 1)

ggsave("05-Alpha_main.pdf", p, units = "in", height = 2, width = 6.5)


#############
# anova TEST

for (index.alpha in vec.Measure) {
  df1 = alphaDF.output[alphaDF.output$Measure == index.alpha, ]
  posthoc = TukeyHSD(aov(Value ~ Group, data = df1), 'Group', conf.level=0.95)
  write.table(data.frame(Compares = rownames(posthoc$Group), posthoc$Group),
              paste0('06-_ANOVA_TurkeyHSDtest_', index.alpha, '.txt'),
              append = F, quote = F, sep = '\t', row.names = F, col.names = T)
}

