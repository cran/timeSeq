\name{timeSeq.geneset}
\alias{timeSeq.geneset}
\title{
Non-parallel Differential Expression Gene Set Detection
}
\description{
Use a negative binomial mixed-effects (NBME) model to identify both PDE and NPDE gene sets in time course RNA-seq data.
}
\usage{
timeSeq.geneset(data, group.label, genesets, gene.names, 
	geneset.names = NULL, reads = NULL, exon.length = NULL, 
	gene.length = NULL, exon.level = TRUE, gene.level = FALSE, 
	p.values = FALSE, n_cores = NULL, iterations = NULL, offset = TRUE)
}
\arguments{
\item{data}{N by p matrix of expression values, one observation per row with p time points.}
\item{group.label}{Vector indicating the experimental group of each time point.}
\item{genesets}{Genes grouping.}
\item{gene.names}{Vector including all the gene names.}
\item{geneset.names}{Vector including all the gene set names.}
\item{reads}{Library size.}
\item{exon.length}{Vector containing the length exons, only used in exon level data.}
\item{gene.length}{Vector containing the length genes, only used in gene level data.}
\item{exon.level}{Indicating if this is an exon level dataset. Default is TRUE.}
\item{gene.level}{Indicating if this is a gene level dataset. Default is FALSE.}
\item{p.values}{Indicating whether caculating p values. Default is FALSE.}
\item{n_cores}{Number of cores used to do parallel computing.}
\item{iterations}{Number of permutations for significance testing.}
\item{offset}{If effective library size can be used as offset.}

}

\value{
A list with compoents
  \item{NPDE.ratio}{Average NPDE ratios of each gene set.}
  \item{PDE.ratio}{Average PDE ratios of each gene set.}
  \item{NPDE.pvalue}{P values for detecting NPDE gene sets.}
  \item{PDE.pvalue}{P values for detecting PDE gene sets.}
  \item{geneset.names}{Vector of names for gene sets.}
}

\references{
Chong Gu. Model diagnostics for smoothing spline ANOVA models. Canadian Journal of Statistics, 32(4):347-358, 2004.

Chong Gu. Smoothing spline ANOVA models. Springer, second edition, 2013.

Chong Gu and Ping Ma. Optimal smoothing in nonparametric mixed-effect models. Annals of Statistics, pages 1357-1379, 2005a.
}


\author{Fan Gao and Xiaoxiao Sun}

\examples{
#### Simulation
set.seed(1)
geneset.data <- matrix(nrow=90, ncol=16)
mu.produce.exp <- function(x,exon.sd)
{exp(10*(1e3*(x^11*(1-x)^6)+100*(x^3*(1-x)^10)) + 5 + exon.sd)}
ti1 <- c(1:8)/10
ti2 <- c(8:1)/10
nu <- 20000
pp <- 1
exon_length <- c()
for(i in 1:10){
  for(j in 1:3){
    exon_length[pp] <- c(0.1,0.25,0.4)[j]
    exon_sd <- rnorm(1, sd=exon_length[pp])
    mu1 <- mu.produce.exp(ti1, exon_sd)
    mu2 <- mu.produce.exp(ti2, exon_sd)
    geneset.data[(i-1)*3 + j, 1:8] <- rnbinom(8,nu,mu=mu1)
    geneset.data[(i-1)*3 + j, 9:16] <- rnbinom(8,nu,mu=mu2)
    pp <- pp + 1
  }
}

for(i in 11:30){
  for(j in 1:3){
    exon_length[pp] <- c(0.1,0.25,0.4)[j]
    exon_sd <- rnorm(1, sd=exon_length[pp])
    mu1 <- mu.produce.exp(ti1, exon_sd)
    geneset.data[(i-1)*3 + j, 1:8] <- rnbinom(8,nu,mu=mu1)
    geneset.data[(i-1)*3 + j, 9:16] <- rnbinom(8,nu,mu=mu1)
    pp <- pp + 1
  }
}

### Gene names should match with the data. 
gene.names=paste("g",rep(1:30, each=3),sep="")

### create some random gene sets
gene.sets=vector("list",3)
for(i in 1:3){gene.sets[[i]]=paste("g",(1:30)[((i-1)*10+1):(i*10)],sep="")}
geneset.names=paste("set",as.character(1:3),sep="")
group.label <- c(rep(1,8), rep(2,8))

##### Model fitting (If you want to run 'timeSeq.geneset' by yourself, just remove the 
##### "## " from the head of two lines below.)
## geneset.fit <-  timeSeq.geneset(geneset.data, group.label, gene.sets, gene.names, 
## 	geneset.names, n_cores = 2, offset=FALSE)
  
###Load the object returned by timeSeq.geneset
data(object_by_timeSeq.geneset)
str(geneset.fit)  
}
