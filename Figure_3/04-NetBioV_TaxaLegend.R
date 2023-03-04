# coding: utf-8
# email: qfsfdxlqy@163.com
# Github: https://github.com/2015qyliang

library(tidyverse)
library(igraph)
library(netbiov)
library(ggsci)
library(gridExtra)

NetByModule = function(fn) {
  taxDF = read.table('02-comTax.txt', header = T, sep = '\t',stringsAsFactors = F)
  colnames(taxDF)[1] = 'id'
  
  nodes = read.table(paste0('./NetInfo/', fn,' node_attribute.txt'), 
                     header = T, sep = '\t', stringsAsFactors = F)
  nodes$No..module = paste0('M', nodes$No..module)
  isLarge = table(nodes$No..module)
  isLarge = isLarge[isLarge >= 5]
  nodesSub = nodes[nodes$No..module %in% names(isLarge), ]
  nodesSub = merge(nodesSub, taxDF, by = 'id', all.x = T)
  rownames(nodesSub) = nodesSub$id
  
  rownames(nodes) = nodes$id
  links = read.table(paste0('./NetInfo/', fn,' edge_attribute.txt'),
                     header = T, sep = '\t', stringsAsFactors = F)
  net.graph = graph_from_data_frame( d = links, vertices = nodes, directed = F)

  nmsV = names(V(net.graph))
  net.graph = induced_subgraph(net.graph, v = which(nmsV %in% nodesSub$id))
  
  v.module = data.frame(OTUID = names(V(net.graph)), 
                        module = nodes[names(V(net.graph)), 'No..module'] ,
                        Phylum = nodesSub[names(V(net.graph)), 'Phylum'],
                        stringsAsFactors = F)
  
  # module.list
  ModuList =  list()
  for (i in 1:length(isLarge)) {
    nodes.module = names(isLarge)[i]
    ModuList[[i]] = v.module[v.module$module == nodes.module, 'OTUID']
  }
  names(ModuList) = names(isLarge)
  
  # node color
  # set color of vertex 
  phylums = readLines('02-Phylum.txt')
  v.module$Phylum[!(v.module$Phylum %in% phylums)] = 'other'
  vcol.set = c(colorRampPalette(pal_nejm()(8))(10), 'grey30')
  names(vcol.set) = c(phylums, 'other')
  vertexCols = vcol.set[v.module$Phylum]
  names(vertexCols) = v.module$OTUID
  
  # set.seed(618)
  return(plot.abstract.module(net.graph,  
                              v.size = 1.5, sf = 0, e.sf = 1, bg = 'white',
                              lab.color = 'black', colors = 'grey90',
                              mod.lab = T,
                              mod.list = ModuList,
                              col.grad = vertexCols[names(V(net.graph))],
                              arrange.bydegree = T,
                              layout.function = layout.circle,
                              layout.overall = layout.kamada.kawai))
}

##########################################################################

sampGrps = read.table('02-RMTvalue.txt', header = T,
                      sep = '\t', stringsAsFactors = F)
fns = paste0(sampGrps$groups, ' ', sprintf('%0.2f', sampGrps$value))

for (i in c(1:9)) {
  pdf(width = 1.3, height = 1.3, file = paste0('ModulePlot/Taxon', sampGrps$groups[i],'.pdf') )
  NetByModule(fns[i])
  dev.off()
}

##########################################################################

phylums = readLines('02-Phylum.txt')
vcol.set = c(colorRampPalette(pal_nejm()(8))(10), 'grey30')
names(vcol.set) = c(phylums, 'other')
# plot(vcol.set)

x = 1:11
y = rep(1, 11)
df = data.frame(x , y, py = c(phylums, 'other'), vcol.set)

p = ggplot(df, aes(x = py, y =y)) + 
  geom_point(color = df$vcol.set, shape = 15) + 
  theme_bw() + 
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        axis.text.y = element_text(size = 6, hjust = 0, vjust = 0.5), 
        axis.text.x = element_blank(), 
        legend.position = 'none') + 
  coord_flip()

ggsave('05-ModuleTaxonLegend.pdf', p, width = 2, height = 1.5, units = 'in')


