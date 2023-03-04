# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(tidyverse)
library(igraph)
library(netbiov)
library(ggsci)
library(gridExtra)

NetByModule = function(fn) {
  nodes = read.table(paste0('./NetInfo/', fn,' node_attribute.txt'), 
                     header = T, sep = '\t', stringsAsFactors = F)
  nodes$No..module = paste0('M', nodes$No..module)
  isLarge = table(nodes$No..module)
  isLarge = isLarge[isLarge >= 5]
  nodesSub = nodes[nodes$No..module %in% names(isLarge), ]
  
  rownames(nodes) = nodes$id
  links = read.table(paste0('./NetInfo/', fn,' edge_attribute.txt'),
                     header = T, sep = '\t', stringsAsFactors = F)
  net.graph = graph_from_data_frame( d = links, vertices = nodes, directed = F)
  
  nmsV = names(V(net.graph))
  net.graph = induced_subgraph(net.graph, v = which(nmsV %in% nodesSub$id))
  
  # set color of vertex 
  vcol.set = colorRampPalette(pal_nejm()(8))(length(isLarge))
  names(vcol.set) = names(isLarge)
  
  v.module = data.frame(OTUID = names(V(net.graph)), 
                        module = nodes[names(V(net.graph)), 'No..module'] , 
                        stringsAsFactors = F)
  
  # module.list
  ModuList =  list()
  for (i in 1:length(isLarge)) {
    nodes.module = names(isLarge)[i]
    ModuList[[i]] = v.module[v.module$module == nodes.module, 'OTUID']
  }
  names(ModuList) = names(isLarge)
  
  set.seed(618)
  fn = function(g)layout.star(g, center=which.max(degree(g))-1)
  return(plot.abstract.module(net.graph,  
                              v.size = 1.5, sf = 0, e.sf = 1, bg = 'white',
                              lab.color = 'black', colors = 'grey90',
                              mod.lab = T,
                              mod.list = ModuList,
                              modules.color = vcol.set,
                              arrange.bydegree = T,
                              layout.function = layout.circle,
                              layout.overall = layout.kamada.kawai))
}

##########################################################################

sampGrps = read.table('02-RMTvalue.txt', header = T,
                      sep = '\t', stringsAsFactors = F)
fns = paste0(sampGrps$groups, ' ', sprintf('%0.2f', sampGrps$value))

for (i in 1:9) {
  pdf(width = 1.3, height = 1.3, file = paste0('ModulePlot/', sampGrps$groups[i],'.pdf') )
  NetByModule(fns[i])
  dev.off()
}

