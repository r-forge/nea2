linknum <-
function(ags,fgs,network){

nags<-length(ags)
nfgs<-length(fgs)
nlink<-rep(0,nfgs)
names(nlink)<-names(fgs)

for (i in 1:nfgs){
    agsfgs<-as.vector(outer(ags,fgs[[i]],paste))
    fgsags<-as.vector(outer(fgs[[i]],ags,paste)) 
    ugs<-c(agsfgs,fgsags) ### we do not need to use "unique" function
    nlink[i]<-sum(ugs %in% network)
}
return(nlink)
}

