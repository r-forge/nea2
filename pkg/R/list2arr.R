list2arr <-
function(netlist, gsym)
{
   names(gsym) = 1:length(gsym)
   g1<- names(netlist)
   deg = sapply(netlist, length)
   gg1 = gsym[rep(g1, deg)]
   gg2 = gsym[unlist(netlist)]
   arr= paste(gg1, gg2)
   return(arr)
}

