# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

############################################################################

library(tidyverse)

############################################################################

otusDF = read.table('02-RelTable.txt', header = T, sep = '\t', stringsAsFactors = F, row.names = 1)
nodesFns = list.files('NetInfo/', 'node')

control = c('00C', '05C', '12C', '21C', '30C')
treatment = c('00C', '05M', '12M', '21M', '30M')

otus = nm.whichCol = c()
for (hfn in control) {
  fns = nodesFns[str_detect(nodesFns, hfn)]
  nodesdf = read.table(paste0('NetInfo/', fns), header = T, sep = '\t', stringsAsFactors = F)
  otus = append(otus, nodesdf$id)
  nm.whichCol = append(nm.whichCol, which(str_detect(colnames(otusDF), hfn)))
}
controlDF = otusDF[unique(otus), nm.whichCol]
controlDF = controlDF[rowSums(controlDF) != 0, ]

otus = nm.whichCol = c()
for (hfn in treatment) {
  fns = nodesFns[str_detect(nodesFns, hfn)]
  nodesdf = read.table(paste0('NetInfo/', fns), header = T, sep = '\t', stringsAsFactors = F)
  otus = append(otus, nodesdf$id)
  nm.whichCol = append(nm.whichCol, which(str_detect(colnames(otusDF), hfn)))
}
treatDF = otusDF[unique(otus), nm.whichCol]
treatDF = treatDF[rowSums(treatDF) != 0, ]

############################################################################
istart = seq(1, 100, 20)
nmcols = c('D00', 'D05', 'D12', 'D21', 'D30')
for (i in 1:5) {
  controlDF$new = apply(controlDF[, istart[i]:(istart[i]+19)], 1, mean)
  treatDF$new = apply(treatDF[, istart[i]:(istart[i]+19)], 1, mean)
  colnames(controlDF)[ncol(controlDF)] = nmcols[i]
  colnames(treatDF)[ncol(treatDF)] = nmcols[i]
}

controlDF = controlDF[, 101:ncol(controlDF)]
treatDF = treatDF[, 101:ncol(treatDF)]

############################################################################

nd.constancy = c()
for (nd in rownames(controlDF)) {
  ndmean = mean(as.numeric(controlDF[nd, ]))
  sd = sqrt(sum((as.numeric(controlDF[nd, ]) - ndmean)^2)/5)
  if (sd != 0) {nd.constancy = append(nd.constancy, ndmean/sd)}
  if (sd == 0) {nd.constancy = append(nd.constancy, 0)}
}
constancyDF = data.frame(Constancy = nd.constancy, 
                         Group = rep('Control', length(nd.constancy)))

nd.constancy = c()
for (nd in rownames(treatDF)) {
  ndmean = mean(as.numeric(treatDF[nd, ]))
  sd = sqrt(sum((as.numeric(treatDF[nd, ]) - ndmean)^2)/5)
  if (sd != 0) {nd.constancy = append(nd.constancy, ndmean/sd)}
  if (sd == 0) {nd.constancy = append(nd.constancy, 0)}
}

constancyDF = rbind(constancyDF,
                    data.frame(Constancy = nd.constancy, 
                               Group = rep('Treatment', length(nd.constancy))))

############################################################################

p = ggplot(constancyDF, aes(x = Group, y = Constancy, color = Group)) + 
  geom_boxplot(size = 0.5, alpha = 0.75, outlier.color = 'grey75', outlier.size = 1) + 
  geom_jitter(size=0.35, width = 0.3, alpha = 0.2) + 
  scale_color_manual(values = c('#E69F00', '#0072B2'))+ 
  labs(x= NULL, y = 'The constancy of nodes') +
  theme( panel.background = element_rect(fill="white",color="black"),
         panel.grid = element_blank(), 
         axis.text.y = element_text(size = 8, angle = 90, vjust =0, hjust = 0.5),
         axis.text.x = element_text(size = 10),
         axis.title.y = element_text(size = 10), 
         strip.background = element_blank(),
         legend.position = 'none')

ggsave('16-ConstancyNodes.pdf', p, units = 'in', width = 2, height = 2.5)


############################################################################

# Kolmogorov-Smirnov test
kr = ks.test(constancyDF[constancyDF$Group == 'Control', 'Constancy'], 
             constancyDF[constancyDF$Group == 'Treatment', 'Constancy'])
kr

# D = 0.50951, p-value < 2.2e-16

############################################################################

wr = wilcox.test(constancyDF[constancyDF$Group == 'Control', 'Constancy'], 
                 constancyDF[constancyDF$Group == 'Treatment', 'Constancy'])
wr
# W = 69660, p-value < 2.2e-16

