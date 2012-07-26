gea <-
function(ags,fgs,ntotal){
if (any(duplicated(ags))) {stop("inpute gene symbols must be unique")}
if (any(duplicated(fgs))) {stop("inpute gene symbols must be unique")}
        a<-sum(ags %in% fgs)
b<-length(fgs)-a
 c<-length(ags)-a
d<-ntotal-a-b-c
        bad<-0; if (a==0||b==0||c==0||d==0) {bad<-1}
if (a==0) {a<-1/2}
if (b==0) {b<-1/2}
if (c==0) {c<-1/2}
if (d==0) {d<-1/2}
z<-log(a*d/(b*c))/(1/a+1/b+1/c+1/d)^{1/2}
hyper<-phyper(q=a,m=a+b,n=c+d,k=a+c) 
return(list(bad=bad,z=z,hyper=hyper,a=a,b=b,c=c,d=d))
}

