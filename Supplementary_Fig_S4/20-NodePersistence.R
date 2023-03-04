# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

############################################################################

library(tidyverse)

############################################################################

otusDF = read.table('02-RelTable.txt', header = T, sep = '\t', stringsAsFactors = F, row.names = 1)
ndDF.0 = read.table('NetInfo/D00C 0.80 node_attribute.txt', header = T, sep = '\t', row.names = 1)
nodesFns = list.files('NetInfo/', 'node')

control = c('05C', '12C', '21C', '30C')
treatment = c('05M', '12M', '21M', '30M')

############################################################################

persistence = c()
for (hfn in control) {
  fns = nodesFns[str_detect(nodesFns, hfn)]
  Ndf = read.table(paste0('NetInfo/', fns), header = T, sep = '\t', row.names = 1)
  totoalNDs = length(unique(c(rownames(ndDF.0), rownames(Ndf))))
  otudf.nd = otusDF[, c(rep(T, 20), str_detect(colnames(otusDF), hfn)[21:ncol(otusDF)])]
  otudf.nd[otudf.nd > 0] = 1
  persistence = append(persistence, sum(rowSums(otudf.nd) == 40)/totoalNDs)
}

for (hfn in treatment) {
  fns = nodesFns[str_detect(nodesFns, hfn)]
  Ndf = read.table(paste0('NetInfo/', fns), header = T, sep = '\t', row.names = 1)
  totoalNDs = length(unique(c(rownames(ndDF.0), rownames(Ndf))))
  otudf.nd = otusDF[, c(rep(T, 20), str_detect(colnames(otusDF), hfn)[21:ncol(otusDF)])]
  otudf.nd[otudf.nd > 0] = 1
  persistence = append(persistence, sum(rowSums(otudf.nd) == 40)/totoalNDs)
}

DFpersist = data.frame(Persistence = persistence, Group = paste0('D', c(control, treatment)))
write.table(DFpersist, '21-Persistence.txt', append = F, sep = '\t', 
            row.names = F, col.names = T)

DFpersist$TimeGP = as.integer(substr(DFpersist$Group, 2, 3))
DFpersist$GP = substr(DFpersist$Group, 4, 4)

############################################################################

p = ggplot(DFpersist, aes(x = TimeGP, y = Persistence, color = GP, group = GP)) +
  geom_point(size = 1) + 
  geom_smooth(method = "lm", se = F, size = 0.3, alpha = 0.5) + 
  scale_color_manual(values = c('#E69F00', '#0072B2'))+ 
  labs(x= NULL, y = 'Persistence') +
  theme( panel.background = element_rect(fill="white",color="black"),
         panel.grid = element_blank(), 
         axis.text.y = element_text(size = 8, angle = 90, vjust =0, hjust = 0.5),
         axis.title.y = element_text(size = 10), 
         legend.position = 'none')

ggsave('22-PersistenceNodes.pdf', p, units = 'in', width = 2, height = 2.5)

##########################################################################

DFpersist$time = as.integer(substr(DFpersist$Group, 2, 3))

DFpersist.lm = DFpersist[DFpersist$GP == 'C', ]
res = summary(lm(DFpersist.lm$Persistence~DFpersist.lm$time))
res$coefficients[,'Estimate'][1]
# (Intercept) 
# 0.16106   
res$r.squared
# [1] 0.0518919
res$adj.r.squared
# [1] -0.4221621
res$coefficients[,'Pr(>|t|)'][2]
# DFpersist.lm$time 
# 0.7722021 

print('#####################')

DFpersist.lm = DFpersist[DFpersist$GP == 'M', ]
res = summary(lm(DFpersist.lm$Persistence~DFpersist.lm$time))
res$coefficients[,'Estimate'][1]
# (Intercept) 
# 0.3564422  
res$r.squared
# [1] 0.01543
res$adj.r.squared
# [1] -0.4768607
res$coefficients[,'Pr(>|t|)'][2]
# DFpersist.lm$time 
# 0.9163921      

############################################################################
