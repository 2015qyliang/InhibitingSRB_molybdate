# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(reshape2)
library(ggplot2)
library(doBy)

envDF = read.table('02-EncFac.txt', header = T, sep = '\t', stringsAsFactors = F)
envDF = na.omit(envDF)

vfaDF = read.table('02-VFAs.txt', header = T, sep = '\t', stringsAsFactors = F)
vfaDF = na.omit(vfaDF)

#################################################

DF.mean = summaryBy(.~EnvFac + EnvFacLong + TimeGroup + GP, data = envDF, FUN = mean)
DF.mean$lngroup = paste0(DF.mean$EnvFac, '_',DF.mean$GP)

##################
DF.mean = DF.mean[DF.mean$TimeGroup %in% c("D12", "D21", "D30"), ]

for (env in unique(DF.mean$EnvFac)) {
  df1 = DF.mean[DF.mean$EnvFac == env, ]
  posthoc = TukeyHSD(aov(values.mean ~ GP, data = df1), 'GP', conf.level=0.95)
  # posthoc = t.test(values.mean ~ GP, data = df1)
  # write.table(data.frame(Test.T = posthoc$statistic, Pvalue = posthoc$p.value), 
  write.table(data.frame(Compares = rownames(posthoc$GP), posthoc$GP),
              paste0('06_Env_ANOVA_TurkeyHSD_', env, '.txt'), 
              append = F, quote = F, sep = '\t', row.names = F, col.names = T)
}


#################################################

DF.mean = summaryBy(.~VFAs + VFAsLong + TimeGroup + GP, data = vfaDF, FUN = mean)
DF.mean$lngroup = paste0(DF.mean$VFAs, '_',DF.mean$GP)

##################
DF.mean = DF.mean[DF.mean$TimeGroup %in% c("D12", "D21", "D30"), ]

for (vfa in unique(DF.mean$VFAs)) {
  df1 = DF.mean[DF.mean$VFAs == vfa, ]
  posthoc = TukeyHSD(aov(values.mean ~ GP, data = df1), 'GP', conf.level=0.95)
  write.table(data.frame(Compares = rownames(posthoc$GP), posthoc$GP), 
              paste0('06_Vfa_ANOVA_TurkeyHSD_', vfa, '.txt'), 
              append = F, quote = F, sep = '\t', row.names = F, col.names = T)
}

