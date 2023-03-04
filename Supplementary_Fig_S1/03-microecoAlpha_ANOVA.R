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
  Cdf = alphaDF[alphaDF$Group %in% c('T00C', 'T05C', 'T12C', 'T21C', 'T30C'), ]
  Cdf$Group = factor(Cdf$Group,
                     levels = c('T00C', 'T05C', 'T12C', 'T21C', 'T30C'),
                     labels = c('T00C', 'T05C', 'T12C', 'T21C', 'T30C'))
  Mdf = alphaDF[alphaDF$Group %in% c('T00C', 'T05M', 'T12M', 'T21M', 'T30M'), ]
  Mdf$Group = factor(Mdf$Group,
                     levels = c('T00C', 'T05M', 'T12M', 'T21M', 'T30M'),
                     labels = c('T00C', 'T05M', 'T12M', 'T21M', 'T30M'))
  Cdf.median = summaryBy(.~Measure + Group, data = Cdf, FUN = median)
  Mdf.median = summaryBy(.~Measure + Group, data = Mdf, FUN = median)
  
  p = ggplot(alphaDF, aes(x = Group, y = Value, fill = Group)) + 
    geom_boxplot(alpha = 0.9, outlier.color = 'grey50', outlier.size = 0.5, size = 0.1) + 
    # geom_jitter(width = 0.28, alpha = 0.2, size = 1) + 
    facet_grid(facets = Measure~., rows = 3, scales = 'free_y') +
    geom_line(data = Cdf.median, aes(x = Group, y = Value.median, group = 1), 
              colour  = '#E69F00', linetype = "solid",  size = 0.6, alpha = 0.95) + 
    geom_point(data = Cdf.median, aes(x = Group, y = Value.median), 
               colour  = '#E69F00', shape = 16,  size = 1.5) + 
    geom_line(data = Mdf.median, aes(x = Group, y = Value.median, group = 1), 
              colour  = '#0072B2', linetype = "solid",  size = 0.75, alpha = 0.65) + 
    geom_point(data = Mdf.median, aes(x = Group, y = Value.median), 
               colour  = '#0072B2', shape = 17,  size = 1.5) + 
    scale_fill_npg() + labs(x = "", y = "") + 
    ylim(min(alphaDF$Value), 1.1*max(alphaDF$Value)) +
    theme_bw()+
    theme(panel.grid = element_blank(), 
          #panel.grid = element_line(size = 0.1, color = 'grey90'), 
          strip.text = element_text(size = 10), 
          strip.background = element_blank(), 
          axis.text.y = element_text(size = 8, hjust = 0.5, vjust = 0.5, angle = 90),
          axis.text.x = element_blank(),
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
##############################################################################
##############################################################################
# supplemental alpha diversity

vec.Measure = c('Chao1', 'ACE', 'InvSimpson', 'Fisher', 'PD')

gridList = list()
for (i in 1:4) {
  v.meas = vec.Measure[i]
  gridList[[i]] = diybox(v.meas)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  print(paste0(i, ' - ', v.meas))
}

gridList[[5]] = diybox(vec.Measure[5]) + 
  scale_x_discrete(labels = c("T00C" = "D00C", 
                              "T05C" = "D05C", "T05M" = "D05M", 
                              "T12C" = "D12C", "T12M" = "D12M", 
                              "T21C" = "D21C", "T21M" = "D21M", 
                              "T30C" = "D30C", "T30M" = "D30M")) 
p = grid.arrange(gridList[[1]], gridList[[2]], gridList[[3]],
                 gridList[[4]], gridList[[5]], nrow = 5)
ggsave("05-Alpha_sup.pdf", p, units = "in", height = 6, width = 4)

# # backup
# gridList[[5]] = diybox(vec.Measure[5]) + 
#   scale_x_discrete(labels = c("T00C" = "D00C", 
#                               "T05C" = "D05C", "T05M" = "D05M", 
#                               "T12C" = "D12C", "T12M" = "D12M", 
#                               "T21C" = "D21C", "T21M" = "D21M", 
#                               "T30C" = "D30C", "T30M" = "D30M")) + 
#   theme(axis.text.x = element_text(size = 8, color = 'black', angle = 90,
#                                    hjust = 0.5, vjust = 0.5))
# p = grid.arrange(gridList[[1]], gridList[[2]], gridList[[3]],
#                  gridList[[4]], gridList[[5]], nrow = 5)
# ggsave("05-Alpha_sup_backup.pdf", p, units = "in", height = 6, width = 4)


#############
# anova TEST

for (index.alpha in vec.Measure) {
  df1 = alphaDF.output[alphaDF.output$Measure == index.alpha, ]
  posthoc = TukeyHSD(aov(Value ~ Group, data = df1), 'Group', conf.level=0.95)
  write.table(data.frame(Compares = rownames(posthoc$Group), posthoc$Group), 
              paste0('06-sup_ANOVA_TurkeyHSDtest_', index.alpha, '.txt'), 
              append = F, quote = F, sep = '\t', row.names = F, col.names = T)
}

##############################################################################
# # check info -- Y axis
# # 
# gridList = list()
# for (i in 1:4) {
#   v.meas = vec.Measure[i]
#   gridList[[i]] = diybox(v.meas)+
#     theme(axis.text.y = element_text(size = 8, color = 'black', angle = 0,
#                                      hjust = 1, vjust = 0.5), 
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank())
#   print(paste0(i, ' - ', v.meas))
# }
# 
# gridList[[5]] = diybox(vec.Measure[5]) + 
#   scale_x_discrete(labels = c("T00C" = "D00C", 
#                               "T05C" = "D05C", "T05M" = "D05M", 
#                               "T12C" = "D12C", "T12M" = "D12M", 
#                               "T21C" = "D21C", "T21M" = "D21M", 
#                               "T30C" = "D30C", "T30M" = "D30M")) + 
#   theme(axis.text.y = element_text(size = 8, color = 'black', angle = 0,
#                                    hjust = 1, vjust = 0.5), 
#         axis.text.x = element_text(size = 8, color = 'black', angle = 90,
#                                    hjust = 0.5, vjust = 0.5))
# p = grid.arrange(gridList[[1]], gridList[[2]], gridList[[3]],
#                  gridList[[4]], gridList[[5]], nrow = 5)
# ggsave("08-6Alpha_sup_checkINFO.pdf", p, units = "in", height = 4, width = 4)

