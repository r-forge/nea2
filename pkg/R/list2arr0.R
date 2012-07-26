list2arr0 <-
function(network){
lvec<-unlist(lapply(network,length))

net<-matrix(0,sum(lvec),2)
net[,1]<-rep(names(network),lvec)
net[,2]<-unlist(network)

links= paste(net[,1], net[,2], sep=' ')
ulinks= unique(links)

return(ulinks)
}

