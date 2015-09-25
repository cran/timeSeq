\name{timeSeq.geneset}
\alias{timeSeq.geneset}
\title{
Differentially Expressed Gene Sets Detection
}
\description{
Use a negative binomial mixed-effects (NBME) model to identify both parallel differentially expressed (PDE) and nonparallel differentially expressed (NPDE) gene sets in time course RNA-seq data.
}
\usage{
timeSeq.geneset(data, group.label, genesets, gene.names, 
	geneset.names = NULL, reads = NULL, exon.length = NULL, 
	gene.length = NULL, exon.level = TRUE, gene.level = FALSE, 
	p.values = FALSE, n_cores = NULL, iterations = NULL, offset = TRUE)
}
\arguments{
\item{data}{a n by p matrix of expression values.}
\item{group.label}{a vector indicating the experimental conditions of each time point.}
\item{genesets}{gene set collection.}
\item{gene.names}{a vector including all the gene names.}
\item{geneset.names}{a vector including all the gene set names.}
\item{reads}{library size.}
\item{exon.length}{a vector containing the length of exons, only used in exon level data.}
\item{gene.length}{a vector containing the length of genes, only used in gene level data.}
\item{exon.level}{logical:indicating if this is an exon level dataset. Default is TRUE.}
\item{gene.level}{logical:indicating if this is a gene level dataset. Default is FALSE.}
\item{p.values}{logical: indicating whether p values should be calculated. Default is FALSE.}
\item{n_cores}{the number of cores used to do parallel computing.}
\item{iterations}{the number of permutations for significance testing.}
\item{offset}{logical: should the effective library size be used as offset. }

}

\value{
A list with compoents
  \item{NPDE.ratio}{average NPDE ratios of each gene set.}
  \item{PDE.ratio}{average PDE ratios of each gene set.}
  \item{NPDE.pvalue}{p values for detecting NPDE gene sets.}
  \item{PDE.pvalue}{p values for detecting PDE gene sets.}
  \item{geneset.names}{names of gene sets.}
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
