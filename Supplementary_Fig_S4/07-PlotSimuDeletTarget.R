# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

##########################################################################

library(ggplot2)
library(ggsci)

##########################################################################

randomDF = read.table('04-simuTargetDelt.txt', header = T, sep = '\t', stringsAsFactors = F)

randomDF$TimeGP = substr(randomDF$network, 1, 3)
randomDF$Group = substr(randomDF$network, 4, 4)

##########################################################################

p = ggplot(randomDF, aes(x = network, y = remain.mean, fill = Group)) +
  geom_col() +
  geom_errorbar(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd), 
                width=.2, position=position_dodge(.9)) + 
  facet_grid(.~TimeGP, space = 'free', scales = 'free_x') +
  scale_fill_manual(values = c('#E69F00', '#0072B2'))+ 
  labs(x= NULL, y = 'Robustness') +
  theme( panel.background = element_rect(fill="white",color="black"),
         panel.grid = element_blank(), 
         axis.text.x = element_text(size = 10, angle = 90, hjust = 0.5, vjust = 0.5),
         axis.text.y = element_text(size = 8, angle = 90, hjust = 0.5, vjust = 0.5),
         axis.title.y = element_text(size = 10), 
         strip.background = element_blank(),
         strip.text = element_blank(),
         legend.position = 'none')

ggsave('11-SimulationTarget.pdf', p, units = 'in', width = 3, height = 2.5)

##########################################################################
# Kolmogorov-Smirnov test
com.pair = c("D05C", "D05M", "D12C", "D12M", "D21C", "D21M", "D30C", "D30M")
cp = statisticD = Pvalue = c()
for (i in c(1,3,5,7)) {
  set.seed(618)
  c = randomDF[randomDF$network == com.pair[i], 'remain.mean']
  m = randomDF[randomDF$network == com.pair[i+1], 'remain.mean']
  sd.c = randomDF[randomDF$network == com.pair[i], 'remain.sd']
  sd.m = randomDF[randomDF$network == com.pair[i+1], 'remain.sd']
  kr = ks.test(rnorm(100, mean = c, sd = sd.c), rnorm(100, mean = m, sd = sd.m))
  cp = append(cp, paste0(com.pair[i], '_vs_', com.pair[i+1]))
  statisticD = append(statisticD, kr$statistic)
  Pvalue = append(Pvalue, kr$p.value)
}

write.table(data.frame(ComparePairs = cp, statisticD, Pvalue), '08-KStestTarget.txt', 
            append = F, quote = F, sep = '\t', row.names = F, col.names = T)

##########################################################################
# Kolmogorov-Smirnov test
com.pair = c("D05C", "D05M", "D12C", "D12M", "D21C", "D21M", "D30C", "D30M")
cp = statisticD = Pvalue = c()
for (i in c(1,3,5,7)) {
  set.seed(618)
  c = randomDF[randomDF$network == com.pair[i], 'remain.mean']
  m = randomDF[randomDF$network == com.pair[i+1], 'remain.mean']
  sd.c = randomDF[randomDF$network == com.pair[i], 'remain.sd']
  sd.m = randomDF[randomDF$network == com.pair[i+1], 'remain.sd']
  kr = wilcox.test(rnorm(100, mean = c, sd = sd.c), rnorm(100, mean = m, sd = sd.m))
  cp = append(cp, paste0(com.pair[i], '_vs_', com.pair[i+1]))
  statisticD = append(statisticD, kr$statistic)
  Pvalue = append(Pvalue, kr$p.value)
}

write.table(data.frame(ComparePairs = cp, statisticD, Pvalue), '08-WMtestTarget.txt', 
            append = F, quote = F, sep = '\t', row.names = F, col.names = T)
