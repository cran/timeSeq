\name{timeSeq}
\alias{timeSeq}
\title{Non-parallel Differential Expression Genes Detection}
\description{
Accurately identifying differentially expressed (DE)
genes from time course RNA-seq data has been of tremendous
significance in creating a global picture of cellular function. DE
genes from time course RNA-seq data can be classified into two
types, parallel DE genes (PDE) and non-parallel DE (NPDE) genes.
The former are often biologically irrelevant, whereas the latter are often
biologically interesting. In this package, we propose a negative binomial
mixed-effects (NBME) model to identify both PDE and NPDE genes in
time course RNA-seq data.
}
\usage{
timeSeq(data, group.label, gene.names, reads = NULL, exon.length = NULL, 
gene.length = NULL, exon.level = TRUE, gene.level = FALSE, p.values = FALSE, 
n_cores = NULL, iterations = 10, offset = TRUE)
}
\arguments{
\item{data}{N by p matrix of expression values, one observation per row with p time points.}
\item{group.label}{Vector indicating the experimental group of each time point.}
\item{gene.names}{Vector including all the gene names.}
\item{reads}{Library size.}
\item{exon.length}{Vector containing the length exons, only used in exon level data.}
\item{gene.length}{Vector containing the length genes, only used in gene level data.}
\item{exon.level}{Indicating if this is an exon level dataset. Default is TRUE.}
\item{gene.level}{Indicating if this is a gene level dataset. Default is FALSE.}
\item{p.values}{Indicating whether caculating p values. Default is FALSE.}
\item{n_cores}{Number of cores used to do parallel computing.}
\item{iterations}{Number of permutations for significance testing. Default 10. It  must be a positive integer greater than 1.}
\item{offset}{If effective library size can be used as offset.}
}

\details{Nonparallel differential expression(NPDE) genes 
             and parallel differential expression(PDE) genes detection.}

\value{
A list with components
\item{sorted}{Object returned by NPDE.sort function. It contains sorted Kullback Leibler Ratios(KLRs) for identifying DE genes}
\item{count}{Number of exons for each gene.}
\item{NPDE.ratio}{NPDE ratios of each gene.}
\item{PDE.ratio}{PDE ratios of each gene.}
\item{genenames}{Gene names of each gene.}
\item{table}{Gene expression values for each gene.}
\item{data}{Same as the input 'data' matrix.}
\item{gene.names}{Same as the input 'gene.names' vector.}
\item{effective.lib.size}{Effective library size estimated by edgeR.}
\item{exon.length}{Lengths of each exon.}
\item{group.label}{Same as the input 'group.label' vector.}
\item{group.length}{Total number of time points.}
\item{group1.length}{Number of time points of group 1.}
\item{group2.length}{Number of time points of group 2.}
\item{exon.level}{Same as the input 'exon.level'.}
\item{p.values}{Same as the input 'p.values'.}
\item{pvalue}{P values for detecting DE genes.}
}

\references{
Chong Gu. Model diagnostics for smoothing spline ANOVA models. Canadian Journal of Statistics, 32(4):347-358, 2004.

Chong Gu. Smoothing spline ANOVA models. Springer, second edition, 2013.

Chong Gu and Ping Ma. Optimal smoothing in nonparametric mixed-effect models. Annals of Statistics, pages 1357-1379, 2005a.
}


\author{Fan Gao and Xiaoxiao Sun}

\examples{
###Data Simulation
mu.produce.exp <- function(x,exon.sd)
{exp(10*(1e3*(x^11*(1-x)^6)+100*(x^3*(1-x)^10)) + 5 + exon.sd)}
ti1 <- c(1:8)/10
ti2 <- c(8:1)/10
nu <- 1000

set.seed(1)
settingC.exon <- matrix(nrow=30, ncol=16)
pp <- 1
exon_length <- c()
for(i in 1:5){
  for(j in 1:3){
    exon_length[pp] <- sample(seq(0.1,0.4,length=3), 1)
    exon_sd <- rnorm(1, sd=exon_length[pp])
    mu1 <- mu.produce.exp(ti1, exon_sd)
    mu2 <- mu.produce.exp(ti2, exon_sd)
    p1 <- nu/(mu1+nu)
    p2 <- nu/(mu2+nu)
    settingC.exon[(i-1)*3 + j, 1:8] <- rnbinom(8,nu,p1)
    settingC.exon[(i-1)*3 + j, 9:16] <- rnbinom(8,nu,p2)
    pp <- pp + 1
  }
}

for(i in 6:10){
  for(j in 1:3){
    exon_length[pp] <- sample(seq(0.1,0.4,length=3), 1)
    exon_sd <- rnorm(1, sd=exon_length[pp])
    mu1 <- mu.produce.exp(ti1, exon_sd)
    p1 <- nu/(mu1+nu)
    settingC.exon[(i-1)*3 + j, 1:8] <- rnbinom(8,nu,p1)
    settingC.exon[(i-1)*3 + j, 9:16] <- rnbinom(8,nu,p1)
    pp <- pp + 1
  }
}

gene.names <- paste("gene", rep(1:10, each=3), sep="")
group.label <- c(rep(1,8), rep(2,8))

##### Model fitting (If you want to run 'timeSeq' by yourself, just remove the 
##### "## " from the head of two lines below.)
## model.fit <- timeSeq(settingC.exon, group.label, gene.names, exon.length=exon_length, 
## 	reads=NULL, n_cores = 2, offset=FALSE)

###Load the object returned by timeSeq
data(object_by_timeSeq)
str(model.fit)
}


