linknum.list <-
function(ags, fgs, netlist, gnum=TRUE, gsym)
{
  if(!is.list(fgs)) fgs = list(fgs)
  nfgs<-length(fgs)
  nlink<-rep(0,nfgs)
    names(nlink)<-names(fgs)

## numbered version
  g123= 1:length(gsym)
  names(g123)= gsym

## translate ags and fgs into numbers
  agsnum = g123[ags]
    agsnum= agsnum[!is.na(agsnum)]
  
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
    agsfgs<-as.vector(outer(agsnum,fgsnum,paste))
    fgsags<-as.vector(outer(fgsnum,agsnum,paste)) 
    ugsnum<-c(agsfgs,fgsags) ### we do not need to use "unique" function
    nlink[i]<-sum(ugsnum %in% asubnet)
    # check using the full network:
    # net1num= list2arr(NET[[1]], 1:length(gsym))  # use numbers
    # pick= (ugsnum %in% net1num)
    # full.link= ugsnum[pick]
    ## compare with linknum(.)
    # gsym[as.numeric(unlist(strsplit(ugsnum[pick], ' ')))]
    ## ugsnum[ugsnum %in% asubnet]  # must be equal to full.link
}
  
return(nlink)
}

