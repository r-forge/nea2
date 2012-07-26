linknum.list.csum <-
function(ags, fgs, netlist, gnum=TRUE, gsym)
{
  if(!is.list(fgs)) fgs = list(fgs)
  nfgs<-length(fgs)

## numbered version
  g123= 1:length(gsym)
  names(g123)= gsym

## translate ags and fgs into numbers
  agsnum = g123[ags]
    agsnum= agsnum[!is.na(agsnum)]

## output
  nlink<- matrix(0, nrow=nfgs, ncol=length(agsnum))
    rownames(nlink)<-names(fgs)
  
## sub-network
for (i in 1:nfgs){
  fgsnum= g123[fgs[[i]]]  ## numbered version of fgs
  num = unique(c(agsnum, fgsnum))
  sub = as.character(num[num %in% names(netlist)])
  subnet = netlist[sub]

## expand sub-network to a vector form
    nnet = sapply(subnet, length)
    gg1 = rep(sub, nnet)
    gg2 = unlist(subnet)
    asubnet = paste(gg1, gg2)

## count links on sub-network
    agsfgs<- matrix(outer(agsnum,fgsnum,paste), ncol=length(fgsnum))
    fgsags<- t(matrix(outer(fgsnum,agsnum,paste), ncol=length(agsnum)))
    ugsnum<-cbind(agsfgs,fgsags) ### we do not need to use "unique" function
    a = matrix(ugsnum %in% asubnet, ncol=ncol(ugsnum))
    nlink[i,]<- cumsum(rowSums(a))
}
return(nlink)
}

