


dfFiltGet = function(fn) {
  df.r = read.table(file.path('KnockRes', fn), header = T,
                    sep = '\t', stringsAsFactors = F)
  df = df.r[df.r$consider_cascade == 'No' & df.r$weighted == 'unweighted',]
  df = df[df$Proportion.removed == 0.5, ]
  return(df[1, ])
}


fns.random = list.files('KnockRes/', pattern = '_random_deletion.txt')
fns.target = list.files('KnockRes/', pattern = '_target_deletion.txt')


randomDF = data.frame()
for (fn in fns.random) {
  randomDF = rbind(randomDF, dfFiltGet(fn))
}
write.table(randomDF, '04-simuRandomDelt.txt', append = F, quote = F, 
            sep = '\t', row.names = F, col.names = T)

####################################################

dfFiltGet = function(fn) {
  df.r = read.table(file.path('KnockRes', fn), header = T,
                    sep = '\t', stringsAsFactors = F)
  df = df.r[df.r$consider_cascade == 'No' & df.r$weighted == 'unweighted',]
  df = df[df$Number.hub.removed == 5, ]
  return(df[1, ])
}


randomDF = data.frame()
for (fn in fns.target) {
  randomDF = rbind(randomDF, dfFiltGet(fn))
}
write.table(randomDF, '04-simuTargetDelt.txt', append = F, quote = F, 
            sep = '\t', row.names = F, col.names = T)

