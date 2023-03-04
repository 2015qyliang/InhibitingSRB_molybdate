
fnsLink = list.files('NetInfo/', pattern = 'edge_attribute.txt')

fnsPears = list.files('NodesPearson/', pattern = 'Pearson Correlation.txt')

fnsOTUs = list.files('NodesPearson/', pattern = 'MV Estimated.txt')


for (i in 1:length(fnsLink)) {
  fnl = fnsLink[i]
  DFfnl = read.table(paste0('NetInfo/', fnl), header = T, sep = '\t')
  colnames(DFfnl)[3] = 'corP'
  
  fno = fnsOTUs[i]
  DFfno = read.table(paste0('NodesPearson/', fno), header = F, sep = '\t', na.strings = 'NA')
  otus = DFfno$V1
  
  fnp = fnsPears[i]
  DFfnp = read.table(paste0('NodesPearson/', fnp), header = F, sep = '\t', na.strings = 'NA')
  colnames(DFfnp) = rownames(DFfnp) = otus
  for (k in 2:nrow(DFfnp)) {
    newln  = c(rep('NA', k-1), DFfnp[k, 1:(ncol(DFfnp)-k+1)])
    DFfnp[k, ] = newln
  }
  
  for (j in 1:nrow(DFfnl)) {
    DFfnl[j, ][3] = DFfnp[as.character(DFfnl[j, ][1]), as.character(DFfnl[j, ][2])]
  }
  
  nfn = strsplit(fnl, split = ' ')[[1]][1]
  write.table(DFfnl, paste0('NodesCors/', nfn, '_linkCorPearson.txt'), 
              append = F, quote = F, sep = '\t', row.names = F, col.names = T)
  
}





