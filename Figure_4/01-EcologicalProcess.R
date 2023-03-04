# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

##################################################################

library(ggplot2)
library(ggpubr)
library(ggsci)
library(gridExtra)
library(tidyverse)
library(reshape2)
library(doBy)

##################################################################

# filter Bins by diff-taxs

binTaxDF = read.table('PD.Bin_TopTaxon.csv', 
                      header = T, sep = ',', stringsAsFactors = F)
binsDF.f = binTaxDF[, c("Bin", "TopTaxon.Family")]

# filter contribution of Bins -- HoS & DL

binContriDF = read.table('PD.BinContributeToProcess_EachGroup.csv', 
                      header = T, sep = ',', stringsAsFactors = F)
gps = unique(binContriDF$Group[1:45])
binContriDF.f = binContriDF[binContriDF$Group %in% gps, 3:ncol(binContriDF)]
colnames(binContriDF.f)[3:ncol(binContriDF.f)] = paste0('Bin', 1:(ncol(binContriDF.f)-2))
binContriDF.f = binContriDF.f[, c(1:2, which(colnames(binContriDF.f) %in% binsDF.f$Bin))]

binContriDF.f$GP = substr(binContriDF.f$Group, 4, 4)
binContriDF.f = binContriDF.f[, 2:ncol(binContriDF.f)]
binCDF.sumBy = summaryBy(.~GP+Process, data = binContriDF.f, FUN = mean)
binCDF.sumBy.sum = apply(binCDF.sumBy[, 3:ncol(binContriDF.f)], 1, sum)
contr.process.df = data.frame(binCDF.sumBy[, 1:2], binCDF.sumBy.sum)
contr.process.df$Process = factor(contr.process.df$Process, 
                                  levels = c("HeS", "HoS", "HD", "DL", "DR"))

write.table(contr.process.df, '03-RelativeImportance.txt', 
            append = F, quote = F, sep = '\t', 
            row.names = F, col.names = T)

##############
# ggpubr -- ring plot

p1 = ggdonutchart(contr.process.df[contr.process.df$GP == 'C', ], color = "white",
                  'binCDF.sumBy.sum',  label = "Process",  fill = "Process") +
  scale_fill_npg(alpha = 0.75) + 
  theme(legend.position = 'none')

p2 = ggdonutchart(contr.process.df[contr.process.df$GP == 'M', ], color = "white",
                  'binCDF.sumBy.sum',  label = "Process",  fill = "Process") +
  scale_fill_npg(alpha = 0.75) + 
  theme(legend.position = 'none')

##################################################################

icboot = readRDS('PD.Boot.rds')
term.process = c('HeS', 'HoS', 'DL', 'HD', 'DR')
Method = c()
Group = c()
values = c()
for (gp in gps) {
  Method = append(Method, rep(term.process, rep(100, 5)))
  Group = append(Group, rep(gp, 500))
  gp.m = icboot$boot.detail[gp][[1]][2:101, 1:5]
  values = append(values, as.vector(gp.m))
}
stoch.df = data.frame(Method, Group, values, 
                      GP = substr(Group, 4, 4), 
                      TP = substr(Group, 1, 3), 
                      stringsAsFactors = F)
stoch.df = stoch.df[stoch.df$TP %in% c('T05', 'T12', 'T21', 'T30'),]
stoch.df$TP = factor(stoch.df$TP, 
                     labels = c('D05', 'D12', 'D21', 'D30'),
                     levels = c('T05', 'T12', 'T21', 'T30'))
hos.df = stoch.df[stoch.df$Method == 'HoS', ]
dl.df = stoch.df[stoch.df$Method == 'DL', ]


p3 = ggplot(hos.df, aes(x = TP, y = values*100, fill = GP )) + 
  geom_boxplot(aes(group = Group), size = 0.1, outlier.color = 'grey50', 
               outlier.shape = 16, outlier.size = 1) +
  # geom_smooth(aes(group = GP), method = "lm", se = F, size = 0.5, alpha = 0.5) + 
  scale_fill_manual(values = c('#E69F00', '#0072B2')) + 
  labs(x = NULL, y = 'Relative importance (%)' ) + 
  theme_bw()+
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(),
        axis.title = element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 8, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 8), 
        plot.subtitle = element_text(size = 8, hjust = 1, vjust = 0.5),
        legend.position = 'none')

p4 = ggplot(dl.df, aes(x = TP, y = values*100, fill = GP)) + 
  geom_boxplot(aes(group = Group), size = 0.1, outlier.color = 'grey50', 
               outlier.shape = 16, outlier.size = 1) +
  # geom_smooth(aes(group = GP),method = "lm", se = F, size = 0.5, alpha = 0.5) + 
  scale_fill_manual(values = c('#E69F00', '#0072B2')) + 
  labs(x = NULL, y = 'Relative importance (%)' ) + 
  theme_bw()+
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(),
        axis.title =  element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 8, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 8), 
        plot.subtitle = element_text(size = 8, hjust = 1, vjust = 0.5),
        legend.position = 'none')

##################################################################
# # add lm smooth
# 
# hos.df$Time = as.integer(substr(hos.df$TP, 2, 3))
# dl.df$Time = as.integer(substr(dl.df$TP, 2, 3))
# 
# head(hos.df)
# 
# write(c(), '02_hos_lmSmooth.txt')
# for (gp in unique(hos.df$GP)) {
#   hos.df.f = hos.df[hos.df$GP == gp, ]
#   res = summary(lm(hos.df.f$values ~ hos.df.f$Time))
#   write('----------------------', '02_hos_lmSmooth.txt', append = T)
#   write(paste0('Group: ', gp), '02_hos_lmSmooth.txt', append = T)
#   write(paste0('r.squared: ', res$r.squared), '02_hos_lmSmooth.txt', append = T)
#   write(paste0('p_value: ', res$coefficients[,'Pr(>|t|)'][2]), '02_hos_lmSmooth.txt', append = T)
# }
# 
# write(c(), '02_dl_lmSmooth.txt')
# for (gp in unique(dl.df$GP)) {
#   dl.df.f = dl.df[dl.df$GP == gp, ]
#   res = summary(lm(dl.df.f$values ~ dl.df.f$Time))
#   write('----------------------', '02_dl_lmSmooth.txt', append = T)
#   write(paste0('Group: ', gp), '02_dl_lmSmooth.txt', append = T)
#   write(paste0('r.squared: ', res$r.squared), '02_dl_lmSmooth.txt', append = T)
#   write(paste0('p_value: ', res$coefficients[,'Pr(>|t|)'][2]), '02_dl_lmSmooth.txt', append = T)
# }

##################################################################
##################################################################
##################################################################

Method = c()
Group = c()
stochas = c()
gps = c("T00C", "T05C", "T05M", "T12C", "T12M", "T21C", "T21M", "T30C", "T30M")

###############################
icboot = readRDS('PD.Boot.rds')
for (gp in gps) {
  Method = append(Method, rep('iCAMP', 100))
  Group = append(Group, rep(gp, 100))
  gp.m = icboot$boot.detail[gp][[1]]
  print(paste0('iCAMP -- ', length(gp.m[2:101, 6])))
  stochas = append(stochas, gp.m[2:101, 6])
}

stoch.df = data.frame(Method, Group, stochas, stringsAsFactors = F)

stoch.df = stoch.df[stoch.df$Method == "iCAMP", ]

stoch.df.m = stoch.df[stoch.df$Group == 'T00C', ]
stoch.df.m$Group = 'T00M'
stoch.df = rbind(stoch.df, stoch.df.m)
stoch.df$GP = substr(stoch.df$Group, 4, 4)

p5 = ggplot(stoch.df, aes(x = Group, y = stochas*100, fill = GP)) +
  geom_boxplot(aes(group = Group),size = 0.1, outlier.color = 'grey50', alpha = 0.75,
               outlier.shape = 16, outlier.size = 1) +
  # geom_smooth(aes(group = GP), method = "lm", se = F, size = 0.5, alpha = 0.5) +
  scale_fill_manual(values = c('#E69F00', '#0072B2')) +
  labs(x = NULL, y = 'Estimated stochasticity (%)' ) +
  theme_bw()+
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(),
        axis.title =  element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 8, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        plot.subtitle = element_text(size = 8, hjust = 1, vjust = 0.5),
        legend.position = 'none')

p6 = ggplot(stoch.df, aes(x = Method, y = stochas*100, fill = GP)) + 
  geom_boxplot(size = 0.1, outlier.color = 'grey50', alpha = 0.75, 
               outlier.shape = 16, outlier.size = 1) +
  scale_fill_manual(values = c('#E69F00', '#0072B2')) + 
  labs(x = NULL, y = 'Estimated stochasticity (%)' ) + 
  theme_bw()+
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(),
        axis.title =  element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 8, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 8), 
        plot.subtitle = element_text(size = 8, hjust = 1, vjust = 0.5),
        legend.position = 'none')

# Mann-Whitney test
c = stoch.df[stoch.df$GP == 'C', 'stochas']
m = stoch.df[stoch.df$GP == 'M', 'stochas']
wr = wilcox.test(c, m)

wr$statistic
wr$p.value

###############################
# stoch.df$Time = as.integer(substr(stoch.df$Group, 2, 3))
# 
# write(c(), '02_stocha_lmSmooth.txt')
# for (gp in unique(stoch.df$GP)) {
#   hos.df.f = stoch.df[stoch.df$GP == gp, ]
#   res = summary(lm(hos.df.f$stochas ~ hos.df.f$Time))
#   write('----------------------', '02_stocha_lmSmooth.txt', append = T)
#   write(paste0('Group: ', gp), '02_stocha_lmSmooth.txt', append = T)
#   write(paste0('r.squared: ', res$r.squared), '02_stocha_lmSmooth.txt', append = T)
#   write(paste0('p_value: ', res$coefficients[,'Pr(>|t|)'][2]), '02_stocha_lmSmooth.txt', append = T)
# }

#########################################

p = grid.arrange(p1, p3, p5, p2, p4, p6, nrow = 2)
ggsave('03-RelativeImportance.pdf', p, width = 7, height = 4)

##################################################################
##################################################################
##################################################################
##################################################################

gps = c('T00C', "T05M", "T30M", "T05C", "T12C", "T12M", "T21M", "T21C", "T30C")
icboot = readRDS('PD.Boot.rds')

DFcompare = icboot$compare
DFcompare.f = DFcompare[which(DFcompare$Group1 %in% gps & DFcompare$Group2 %in% gps), ]
DFcompare.f = DFcompare.f[which(DFcompare.f$Group1 %in% 'T00C' | DFcompare.f$Group2 %in% 'T00C'), ]
DFcompare.f$Group1[which(DFcompare.f$Group1 == 'T05M')] = 'T00C'
DFcompare.f$Group2[which(DFcompare.f$Group2 == 'T00C')] = 'T05M'
DFcompare.f$Group = paste0(DFcompare.f$Group1, '_vs_', DFcompare.f$Group2)

write.table(DFcompare.f, '04-cohenD.txt', append = F, quote = F, sep = '\t', 
            row.names = F, col.names = T)
