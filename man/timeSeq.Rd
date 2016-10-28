\name{timeSeq}
\alias{timeSeq}
\title{Statistical Inference for Time Course RNA-Seq Data using a Negative Binomial Mixed-Effects Model}
\description{
Accurately identifying differentially expressed (DE)
genes from time course RNA-seq data has been of tremendous
significance in creating a global picture of cellular function. DE
genes from the time course RNA-seq data can be classified into two
types, parallel DE genes (PDE) and non-parallel DE (NPDE) genes.
The former are often biologically irrelevant, whereas the latter are often
biologically interesting. In this package, we propose a negative binomial
mixed-effects (NBME) model to identify both PDE and NPDE genes in
time course RNA-seq data.
}
\usage{
timeSeq(data.count, group.label, gene.names, reads = NULL, exon.length = NULL, 
exon.level = TRUE, n.cores = NULL, offset = TRUE, iterations = 10, p.values=FALSE)
}
\arguments{
\item{data.count}{a n by p matrix of expression values.}
\item{group.label}{a vector indicating the experimental conditions of each time point.}
\item{gene.names}{a vector containing all the gene names.}
\item{reads}{library size.}
\item{exon.length}{a vector containing the length of exons, only used in exon level data.}
\item{exon.level}{logical:indicating if this is an exon level dataset. Default is TRUE.}
\item{n.cores}{the number of cores used to do parallel computing.}
\item{offset}{logical: should the library size be used as offset.}
\item{iterations}{the number of iterations times used to compute p values of NPDE genes.}
\item{p.values}{logical: should the p values be computed.}
}

\details{Nonparallel differential expression(NPDE) genes 
             and parallel differential expression(PDE) genes detection.}

\value{
A list with components
\item{sorted}{an object returned by NPDE.sort function. It contains sorted Kullback Leibler Ratios(KLRs) for identifying DE genes}
\item{count}{the number of exons for each gene.}
\item{NPDE.ratio}{the NPDE ratios.}
\item{PDE.ratio}{the PDE ratios.}
\item{genenames}{gene names.}
\item{table}{gene expression values.}
\item{data}{a n by p matrix of expression values.}
\item{gene.names}{a vector including all the gene names.}
\item{effective.lib.size}{the effective library size estimated by edgeR.}
\item{exon.length}{length of exons.}
\item{group.label}{a vector indicating the experimental conditions of each time point.}
\item{group.length}{the total number of time points.}
\item{group1.length}{the number of time points of condition one.}
\item{group2.length}{the number of time points of condition two.}
\item{exon.level}{logical:indicating if this is an exon level dataset. Default is TRUE.}
}

\references{
Sun, Xiaoxiao, David Dalpiaz, Di Wu, Jun S. Liu, Wenxuan Zhong, and Ping Ma. "Statistical inference for time course RNA-Seq data using a negative binomial mixed-effect model." BMC Bioinformatics, 17(1):324, 2016.

Chong Gu. Model diagnostics for smoothing spline ANOVA models. Canadian Journal of Statistics, 32(4):347-358, 2004.

Chong Gu. Smoothing spline ANOVA models. Springer, second edition, 2013.

Chong Gu and Ping Ma. Optimal smoothing in nonparametric mixed-effect models. Annals of Statistics, 1357-1379, 2005.
}


\author{Fan Gao and Xiaoxiao Sun}

\examples{
##Exon level data
data(pAbp)
attach(pAbp)
model.fit <- timeSeq(data.count, group.label, gene.names, reads, exon.length, n.cores = 1)
detach(pAbp)

##Gene level data (three replicates)
data(simulate.dt)
attach(simulate.dt)
model.fit <- timeSeq(data.count, group.label, gene.names, reads, exon.level=FALSE, n.cores = 1)
    
}


