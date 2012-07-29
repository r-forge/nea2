\name{nea}
\alias{nea}
\alias{manypnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Network Enrichment Analysis
}
\description{
nea performs network enrichment analysis.
The number of network links between a novel altered gene set (AGS) and a known functional gene set (FGS) is computed and 
assessed by the number of network links obtained from randomly permuted network where the degree distribution is kept. 
As a measure of activation of the FGS, we provide the z-score. 

manypnet generates a specified number of randomly permuted networks. 
}
\usage{
nea(ags, fgs, fgslib = NULL, network, pnet = NULL, nperm = 50, stat="F", seed = NULL)

manypnet(network,nperm=50,seed=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{ags}{
A vector of altered genes. Gene symbols (upper case) are used as a default.
}
  \item{fgs}{
A list defined by user or a character to specify GO ontologies or KEGG pathway. 
User can provide their own functional 
gene sets as a list. 
Options to specify GO ontologies or KEGG pathway are
"CC","BP","MF" and "KEGG" (cellular component, biological process, moelcular function and KEGG pathway). 
}
  \item{fgslib}{
A character of the name of annotation data.
To use GO terms or KEGG pathways, a specific annotation data should be specified. 
}
  \item{network}{
A vector of gene pairs 
or a list representing the network link. 
For the vector case, each element has
a combined name of two gene symbols with separation ''. For example, 
"7SK ALDOA". The number of network links
is computed with respect to this network.
}
  \item{pnet}{
A list of randomly permuted networks. 
This is the output of manypnet. If you do not specify this,
the permutation procedure for the network automatically performed. 
}
  \item{nperm}{
The number of permutations
}
  \item{stat}{
Two types of network enrichment statistic: FNEA and MNEA. 
FNEA (stat="F") depends on the size of ags while MNEA (Stat="M") avoids the dependence by taking
the largest statistic across AGS. Default is FNEA.
}
  \item{seed}{
A seed number for the permutation procedure
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
 \item{nlink}{A vector (of length of fgs) of numbers (stat="F") or a matrix (length of fgs by length of ags lying in the given network) of numbers (stat="M"). The numbers are the network links between AGS and FGS under the given network.}
 \item{res.nlink}{A matrix (the length of fgs by the number of permutation) of network link numbers obtained from the randomly permuted networks (stat="F") or a list of matrix consisting of network link numbers obtained from the randomly permuted networks (stat="M"). When stat="M", each matrix has the same size as nlink and the length of the list is the same as the number of permutation.}
 \item{zscore}{A vector (of length of fgs) of numbers to assess the activation of FGS (stat="F").}
 \item{zscore.mat}{A matrix (the length of fgs by the number of permutation+1) of numbers to assess the activation of FGS (stat="M"). The first nperm columns are from the randomly permuted networks and the last column is from observed data.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Woojoo Lee and Yudi Pawitan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Example 1: artifical data.

AGS<-c("A","B","C","D","E")
FGS<-list(FGS1=c("B","E","F","G","H"),FGS2=c("F","C"))
NETWORK<-c("A B","A F","E C","Q R","S T","U V")

res<-nea(ags=AGS, fgs = FGS, fgslib = NULL, network=NETWORK, pnet = NULL, nperm = 50, stat="F", seed = 1234)
res$nlink
res$res.nlink

## Example 2: use of annotaion data.
## We use gene symbols (upper case)
## "hgu133a.db": Affymetrix Human Genome U133 Set annotation data 

AGS<-c("AIFM3","DIMT1L","ADNP","AHCYL1","EIF4H","RGL1","SEC23IP","EIF4A1","CSNK2B","NOS3")
NETWORK<-c("DNAJC6 RGL1","C1ORF156 NCBP2","AHCYL1 RTN3","PLK4 SKIV2L2","C22ORF28 MESDC2","TINP1 UTP23",    
"HEATR3 MVD","WBP11 XAB2","CSNK2B PA2G4","GCN1L1 RRM2","DIMT1L SMC1A","GPN3 THOC3",     
"DLG3 GPHN","C19ORF29 EXOSC4","AIFM3 SFXN5","HSPA1L RUVBL2","DLAT EIF4A1","ADNP XRCC5",     
"NOS2 NOS3","CIZ1 TLK2","MRPL49 RPS7","GSPT1 SLK","LUC7L2 SEC23IP","DHX8 IGF2BP3",   
"CNTROB SASS6","MRPS12 RPLP2","DHODH EIF4H","GINS3 KIF23","ANXA5 TGFBI","CDK5 PMM1")

res<-nea(ags=AGS, fgs = "KEGG", fgslib = "hgu133a.db", network=NETWORK, pnet = NULL, nperm = 5, stat="F", seed = 1234)
res$nlink
res$res.nlink
}
%## The function is currently defined as
%AA
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
