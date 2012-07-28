nea <-
function(ags,fgs,fgslib=NULL,network,pnet=NULL,nperm=50,stat="F",seed=NULL){


if (is.character(fgs) && fgs!="MF" && fgs!="BP" && fgs!="CC" && fgs!="KEGG") {fgs.type=5;fgs.list<-list();fgs.list[[1]]<-fgs;fgs<-fgs.list}

if (!is.list(fgs) && length(fgs)==1){
require(fgslib, character.only = TRUE) || stop("need data package:", fgslib)
require("GOstats", character.only = TRUE) || stop("need data package:", "GOstats")
}


if (any(duplicated(ags))) print("Duplicated in AGS: This procedure makes the input(gene symbols) unique")


ags<-unique(toupper(ags))
if (is.list(fgs)) {fgs<-lapply(fgs,toupper);fgs<-lapply(fgs,unique)}
if (is.character(network)) {network<-toupper(network)}
if (is.list(network)) {network<-lapply(network,toupper)}


if (!is.null(fgslib) && is.character(fgslib)) {fgslib<-substr(fgslib,1,nchar(fgslib)-3)}

if (is.null(pnet)) {doperm<-"yes"}
if (is.list(pnet)) {doperm<-"no";nperm<-length(pnet)}

if (is.character(fgs) && fgs=="MF")   {fgs.type=1}
if (is.character(fgs) && fgs=="BP")   {fgs.type=2}
if (is.character(fgs) && fgs=="CC")   {fgs.type=3}
if (is.character(fgs) && fgs=="KEGG") {fgs.type=4;fgs.probe <- AnnotationDbi::as.list(get(paste(fgslib, "PATH2PROBE", sep = "")))}
if (is.list(fgs)) {fgs.type=5;fgs.list<-fgs}
 
if (fgs.type< 4) {fgs.probe <- AnnotationDbi::as.list(get(paste(fgslib, "GO2PROBE", sep = "")));fgsname<-names(fgs.probe);ont<-Ontology(fgsname)}
if (fgs.type==1) {fgs.probe <-  fgs.probe[which(ont=="MF")]}
if (fgs.type==2) {fgs.probe <-  fgs.probe[which(ont=="BP")]}
if (fgs.type==3) {fgs.probe <-  fgs.probe[which(ont=="CC")]}


if (fgs.type<5){
fgs.symbol<-AnnotationDbi::as.list(get(paste(fgslib, "SYMBOL", sep = "")))

nf<-length(fgs.probe)
fgs.list<-list()
for (i in 1:nf){

xnames<-unique(unlist(fgs.symbol[fgs.probe[[i]]]))
fgs.list[[i]]<-xnames

}

names(fgs.list)<-names(fgs.probe)
}


if (is.character(network)) {margin<-countcon(network);genes<-margin$genes}
if (is.list(network))      {anet<-list2arr0(network);margin<-countcon(anet);genes<-margin$genes}


if (is.character(network) && stat=="F"){netlist2<-arr2list(network,gsym=genes);nlink<-linknum.arr(ags,fgs=fgs.list,netlist2,gnum=TRUE,gsym=genes)}
if (is.list(network) && stat=="F") {nlink<-linknum.list(ags,fgs=fgs.list,network,gnum=TRUE,gsym=genes)}

if (is.character(network) && stat=="M"){netlist2<-arr2list(network,gsym=genes);nlink<-linknum.arr.csum(ags,fgs=fgs.list,netlist2,gnum=TRUE,gsym=genes)}
if (is.list(network) && stat=="M") {nlink<-linknum.list.csum(ags,fgs=fgs.list,network,gnum=TRUE,gsym=genes)}


if (is.null(seed))  {set.seed(Sys.time())}
if (!is.null(seed)) {set.seed(seed)}

if (stat=="F") {res.nlink<-matrix(0,length(fgs.list),nperm)}
if (stat=="M") {res.nlink<-vector("list",nperm)}

for (i in 1:nperm){    
   
    #if (doperm=="yes") {perm.ncon<-netperm(margin$ncon);pos1<-perm.ncon$pos1}
    #if (doperm=="no")  {pos1<-pnet[[i]]}    
    #g1<- as.numeric(names(pos1))    
    #nnet = sapply(pos1, length)
    #gg1 = genes[rep(g1, nnet)]
    #gg2 = genes[unlist(pos1)]
    #perm.network= paste(gg1, gg2)
    ### old ##res.nlink[,i]<-linknum(ags,fgs=fgs.list,network=perm.network)
    
    if (doperm=="yes" && stat=="F") {perm.ncon<-netperm(margin$ncon);pos1<-perm.ncon$pos1;res.nlink[,i]<-linknum.list(ags,fgs=fgs.list,netlist=pos1,gnum=TRUE,gsym=genes)}
    if (doperm=="no" && stat=="F")  {res.nlink[,i]<-linknum.list(ags,fgs=fgs.list,netlist=pnet[[i]],gnum=TRUE,gsym=genes)}

    if (doperm=="yes" && stat=="M") {perm.ncon<-netperm(margin$ncon);pos1<-perm.ncon$pos1;res.nlink[[i]]<-linknum.list.csum(ags,fgs=fgs.list,netlist=pos1,gnum=TRUE,gsym=genes)}
    if (doperm=="no" && stat=="M")  {res.nlink[[i]]<-linknum.list.csum(ags,fgs=fgs.list,netlist=pnet[[i]],gnum=TRUE,gsym=genes)} 
}


if (stat=="F"){

rownames(res.nlink)<-names(fgs.list)

mean<-apply(res.nlink,1,mean)
se<-apply(res.nlink,1,sd)
se<-ifelse(se==0,1,se)
zscore<-(nlink-mean)/se

#pvalue<-rep(0,length(nlink))
#for (i in 1:length(nlink)){ 
#pvalue[i]<-min(2*min(sum(nlink[i]<=res.nlink[i,])/nperm, sum(nlink[i]>=res.nlink[i,])/nperm),1)
#}

#FDR1d<-pval2FDR(pvalue)

return(list(nlink=nlink,res.nlink=res.nlink,zscore=zscore))
}


if (stat=="M"){
BSTAR<-array(0,dim=c(nrow(nlink),ncol(nlink),nperm+1)) #array(0,dim=c(length(fgs),length(ags),nperm+1))

for (i in 1:nperm){
BSTAR[,,i]<-res.nlink[[i]]
}
BSTAR[,,nperm+1]<-nlink
zscore.mat<-zstat(BSTAR)

return(list(nlink=nlink,res.nlink=res.nlink,zscore.mat=zscore.mat))
}


}

