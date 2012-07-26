arr2list <-
function(anet, gsym){
  g123 = 1:length(gsym)
    names(g123) = gsym
  net = unlist(strsplit(anet,' '))
  net123 = g123[net]  ## numbered network nodes
  mnet = matrix(net123, ncol=2, byrow=TRUE)  # matrix form

  list1  = tapply(mnet[,2], mnet[,1],c)
  list2 = tapply(mnet[,1], mnet[,2],c)
  return(list(list1=list1, list2=list2))
}

