\name{timeSeq.permutate}
\alias{timeSeq.permutate}
\title{Bootstrap method for Significance Testing of NPDE genes
}
\description{
The p values for identifying NPDE genes are calculated through a bootstrap procedure.
}
\usage{
timeSeq.permutate(gene.initial, gene.name, gene.names, group.length, 
	group1.length, group2.length, group.label, ratio.initial, 
	n.cores = NULL, iterations = 10, lib.size, exon.length, 
	exon.level, offset) 
}
\arguments{
  \item{gene.initial}{read counts of a gene.}
  \item{gene.name}{target gene name.}
  \item{gene.names}{a vector including all the gene names.}
  \item{group.length}{total number of time points.}
  \item{group1.length}{the number of time points for condition one.}
  \item{group2.length}{the number of time points for comdition two.}
  \item{group.label}{a vector indicating the experimental condition of each time point.}
  \item{ratio.initial}{the original Kullback-Leibler distance ratio.}
	\item{n.cores}{the number of cores used to do parallel computing.}
  \item{iterations}{the number of permutations for computing p values. Default 10. It  must be a positive integer greater than 1.}
  \item{lib.size}{library size.}
  \item{exon.length}{a vector containing the length of exons. Only used in exon level data.}
  \item{exon.level}{logical: indicating if this is an exon level dataset.}
  \item{offset}{logical: should the effective library size be used as offset.}
}
\value{
A list with components
  \item{pvalue.NPDE}{a vector of p-values for NPDE genes.}
}
\references{
Chong Gu. Model diagnostics for smoothing spline ANOVA models. Canadian Journal of Statistics, 32(4):347-358, 2004.

Chong Gu. Smoothing spline ANOVA models. Springer, second edition, 2013.

Chong Gu and Ping Ma. Optimal smoothing in nonparametric mixed-effect models. Annals of Statistics, pages 1357-1379, 2005a.
}
\author{
Fan Gao and Xiaoxiao Sun
}
