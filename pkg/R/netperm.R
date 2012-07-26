netperm <-
function(x){  
 nx <-length(x)    ## number of genes
 remain_deg <-x    ## At every step, we compute the "remaining" links
 status<- 1        ## Whether the result is a success or not.
 pos1<-list()      ### List to save the positions of 1s for each node
 
 while(max(remain_deg)>0) 
 {
    locmax = which(remain_deg== max(remain_deg))[1]  # location of maximum, take the first one
       deg<- max(remain_deg); 
       remain_deg[locmax]<-0;

    poss_pos <- which(remain_deg>0)  ### We can give 1 to positions only having positive remaining degree
    if (length(poss_pos) <deg) {status<- 0; break} ## not enough positions to sample from

    ## weighted sampling
    perm_pos <- resample(x=1:nx, size=deg, replace=FALSE, prob=remain_deg)
    remain_deg[perm_pos]<- remain_deg[perm_pos]-1
    pos1[[as.character(locmax)]]<- perm_pos
 }
return(list(status=status,pos1=pos1))
}

