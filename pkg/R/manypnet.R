manypnet <-
function(network,nperm=50,seed=NULL){

if (is.null(seed))  {set.seed(Sys.time())}
if (!is.null(seed)) {set.seed(seed)}

ncon<-countcon(network)$ncon
NET = list()
for (i in 1:nperm){
  print(i)
  stat= 0
  while(stat==0){  # repeat if permutation fails
    netstar = netperm(ncon)
    stat= netstar$stat
  }
  NET[[i]] = netstar$pos
}
return(NET)
}

