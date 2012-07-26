countcon <-
function(network){
unet = unlist(strsplit(network,' '))
ncon = sort(table(c(unet)), decr=TRUE) ## number of connections for each gene, sorted
genes = names(ncon)   #unique names
return(list(genes=genes,ncon=ncon))
}

