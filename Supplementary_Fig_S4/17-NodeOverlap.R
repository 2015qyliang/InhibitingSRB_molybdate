# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

############################################################################

library(zetadiv)
library(tidyverse)

############################################################################

getOverlap = function(gps, raw.fns){
  nd.overlap = c()
  fn1 = raw.fns[str_detect(raw.fns, gps[1])]
  nd1.df = read.table(paste0('NetInfo/', fn1), header = T, sep = '\t', row.names = 1)
  nd.overlap = rownames(nd1.df)
  for (i in 2:length(gps)) {
    fn2 = raw.fns[str_detect(raw.fns, gps[i])]
    nd2.df = read.table(paste0('NetInfo/', fn2), header = T, sep = '\t', row.names = 1)
    nd.overlap = intersect(nd.overlap, rownames(nd2.df))
  }
  return(length(nd.overlap))
}

getComb = function(combSize, gpMarker, raw.fns) {
  m.comb = combn(gpMarker, combSize)
  overlaps = c()
  for (i in 1:ncol(m.comb)) {
    gps = m.comb[, i]
    overlaps = append(overlaps, getOverlap(gps, raw.fns))
  }
  return(data.frame(SizeComb = rep(combSize, length(overlaps)), 
                    OverlapNumber = overlaps))
}

############################################################################

nodesFns = list.files('NetInfo/', 'node')
control = c('00C', '05C', '12C', '21C', '30C')
treatment = c('00C', '05M', '12M', '21M', '30M')
comb.sizes = c(2, 3, 4, 5)

overlapDF = data.frame()
for (cs in comb.sizes) {
  overdf = getComb(cs, control, nodesFns)
  overdf$Group = 'Control'
  overlapDF = rbind(overlapDF, overdf)
  overdf = getComb(cs, treatment, nodesFns)
  overdf$Group = 'Treatment'
  overlapDF = rbind(overlapDF, overdf)
}

head(overlapDF)
overlapDF$SizeComb = factor(overlapDF$SizeComb, levels = unique(overlapDF$SizeComb))

p = ggplot(overlapDF, aes(x = SizeComb, y = OverlapNumber, color = Group)) + 
  geom_boxplot(size = 0.5, alpha = 0.75, position = 'dodge', 
               outlier.color = 'grey75', outlier.size = 1) + 
  # geom_jitter(size=0.35, width = 0.3, alpha = 0.2) + 
  scale_color_manual(values = c('#E69F00', '#0072B2'))+ 
  labs(x= NULL, y = 'The number of overlaoping nodes') +
  theme( panel.background = element_rect(fill="white",color="black"),
         panel.grid = element_blank(), 
         axis.text.y = element_text(size = 8, angle = 90, vjust =0, hjust = 0.5),
         axis.text.x = element_text(size = 10),
         axis.title.y = element_text(size = 10), 
         strip.background = element_blank(),
         legend.position = 'none')

ggsave('18-OverlapNodesOrder.pdf', p, units = 'in', width = 2, height = 2.5)

############################################################################

# Kolmogorov-Smirnov test
comb.sizes = c(2, 3, 4, 5)
cp = statisticD = Pvalue = c()
for (cs in comb.sizes) {
  c = overlapDF[overlapDF$SizeComb == cs & overlapDF$Group == 'Control', 'OverlapNumber']
  m = overlapDF[overlapDF$SizeComb == cs & overlapDF$Group == 'Treatment', 'OverlapNumber']
  kr = ks.test(c, m)
  cp = append(cp, paste0(cs, 'C_vs_T', cs))
  statisticD = append(statisticD, kr$statistic)
  Pvalue = append(Pvalue, kr$p.value)
}

write.table(data.frame(ComparePairs = cp, statisticD, Pvalue), '19-KStest.txt', 
            append = F, quote = F, sep = '\t', row.names = F, col.names = T)

# Mann-Whitney test
comb.sizes = c(2, 3, 4, 5)
cp = statisticD = Pvalue = c()
for (cs in comb.sizes) {
  c = overlapDF[overlapDF$SizeComb == cs & overlapDF$Group == 'Control', 'OverlapNumber']
  m = overlapDF[overlapDF$SizeComb == cs & overlapDF$Group == 'Treatment', 'OverlapNumber']
  wr = wilcox.test(c, m)
  cp = append(cp, paste0(cs, 'C_vs_T', cs))
  statisticD = append(statisticD, wr$statistic)
  Pvalue = append(Pvalue, wr$p.value)
}

write.table(data.frame(ComparePairs = cp, statisticD, Pvalue), '19-MWtest.txt', 
            append = F, quote = F, sep = '\t', row.names = F, col.names = T)

############################################################################
