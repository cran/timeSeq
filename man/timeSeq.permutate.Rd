\name{timeSeq.permutate}
\alias{timeSeq.permutate}
\title{Permutation for Significance Testing
}
\description{
We permutate for significance testing. Shuffle the time label for each gene and recompute the test statistic -- Kullback-Leibler distance ratio -- for each permuted gene. 
}
\usage{
timeSeq.permutate(gene_initial, gene_name, gene.names, group.length, 
	group1.length, group2.length, group.label, ratio_initial, 
	n_cores = NULL, iterations = 10, effective.lib.size, exon.length, 
	exon.level, gene.length, gene.level, offset) 
}
\arguments{
  \item{gene_initial}{Expression values of one gene.}
  \item{gene_name}{Set of all the gene names.}
  \item{gene.names}{Vector including all the gene names.}
  \item{group.length}{Total number of time points.}
  \item{group1.length}{Number of time points of group 1.}
  \item{group2.length}{Number of time points of group 2.}
  \item{group.label}{Vector indicating the experimental group of each time point.}
  \item{ratio_initial}{Original Kullback-Leibler distance ratio.}
	\item{n_cores}{Number of cores used to do parallel computing.}
  \item{iterations}{Number of permutations for computing p values. Default 10. It  must be a positive integer greater than 1.}
  \item{effective.lib.size}{Effective library size.}
  \item{exon.length}{Vector containing the length exons, only used in exon level data.}
  \item{exon.level}{Indicating if this is an exon level dataset.}
  \item{gene.length}{Vector containing the length genes, only used in gene level data.}
  \item{gene.level}{Indicating if this is a gene level dataset.}
  \item{offset}{If effective library size can be used as offset.}
}
\value{
A list with components
  \item{pvalue.NPDE}{Vector of p-values for NPDE genes.}
  \item{pvalue.PDE}{Vector of p-values for PDE genes.}
}
\references{
Chong Gu. Model diagnostics for smoothing spline ANOVA models. Canadian Journal of Statistics, 32(4):347-358, 2004.

Chong Gu. Smoothing spline ANOVA models. Springer, second edition, 2013.

Chong Gu and Ping Ma. Optimal smoothing in nonparametric mixed-effect models. Annals of Statistics, pages 1357-1379, 2005a.
}
\author{
Fan Gao and Xiaoxiao Sun
}
