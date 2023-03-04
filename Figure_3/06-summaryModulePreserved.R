# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(tidyverse)

##########################################################################

getNodesDF = function(fn) {
  nodes = read.table(paste0('./NetInfo/', fn,' node_attribute.txt'), 
                     header = T, sep = '\t', stringsAsFactors = F)
  nodes$No..module = paste0('M', nodes$No..module)
  rownames(nodes) = nodes$id
  tab.modules = table(nodes$No..module)
  nodes = nodes[nodes$No..module %in% names(tab.modules)[tab.modules >= 5], ]
  return(nodes)
}

# Reference:
# https://www.nature.com/articles/s41558-021-00989-9
# Climate warming enhances microbial network complexity and stability
CompareModules = function(df1, df2){
  md1 = unique(df1$No..module)
  md2 = unique(df2$No..module)
  # Modu1 = Modu2 = Both = P1A2 = P2A1 = A1A2 = Pvalues = c()
  Modu1 = Modu2 = Both = M1nodes = M2nodes = Pvalues = c()
  for (m1 in md1) {
    for (m2 in md2) {
      nds1 = df1[df1$No..module == m1, 'id']
      nds2 = df2[df2$No..module == m2, 'id']
      nds.share = intersect(nds1, nds2)
      only_1 = length(setdiff(nds1, nds.share))
      only_2 = length(setdiff(nds2, nds.share))
      overlap = length(nds.share)
      none = length(setdiff(c(df1$id, df2$id), c(nds1, nds2)))
      m.ids = matrix(c(overlap, only_1, only_2, none), nrow = 2)
      rest.p = fisher.test(m.ids, alternative = 'greater')$p.value
      Modu1 = append(Modu1, m1)
      Modu2 = append(Modu2, m2)
      Both = append(Both, overlap)
      M1nodes = append(M1nodes, length(nds1))
      M2nodes = append(M2nodes, length(nds2))
      # A1A2 = append(A1A2, none)
      Pvalues = append(Pvalues, rest.p)
    }
  }
  return(data.frame(Modu1, Modu2, Both, M1nodes, M2nodes, Pvalues, stringsAsFactors = F))
}

##########################################################################

testFisher = function(firstDN, secDNs){
  tmp = data.frame()
  for (i in secDNs) {
    df1 = df.list[[firstDN]]
    df2 = df.list[[i]]
    newdf = CompareModules(df1, df2)
    newdf$DF1vsDF2 = rep(paste0(dns[firstDN], '_', dns[i] ), nrow(newdf))
    tmp = rbind(tmp, newdf)
  }
  return(tmp)
}

##########################################################################

sampGrps = read.table('02-RMTvalue.txt', header = T,
                      sep = '\t', stringsAsFactors = F)
df.list = list()
dns = sampGrps$groups
for (i in 1:9) {
  fn = paste0(sampGrps$groups[i], ' ', sprintf('%0.2f', sampGrps$value[i]))
  df.list[[i]] = getNodesDF(fn)
}

##########################################################################

DFsum = data.frame()
DFsum = rbind(DFsum, testFisher(1, secDNs = 2:9))
DFsum = rbind(DFsum, testFisher(2, secDNs = c(3,4,6,8)))
DFsum = rbind(DFsum, testFisher(3, secDNs = c(5,7,9)))
DFsum = rbind(DFsum, testFisher(4, secDNs = c(5,6,8)))
DFsum = rbind(DFsum, testFisher(5, secDNs = c(7,9)))
DFsum = rbind(DFsum, testFisher(6, secDNs = c(7,8)))
DFsum = rbind(DFsum, testFisher(7, secDNs = c(9)))
DFsum = rbind(DFsum, testFisher(8, secDNs = c(9)))

write.table(DFsum, '07-ModulePreserve.txt', append = F, quote = F,
            sep = '\t', row.names = F, col.names = T)

##########################################################################
